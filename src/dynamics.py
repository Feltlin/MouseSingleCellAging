"""Dynamical-systems models of aging trajectories (16_dynamics).

Three groups:
  * accumulation   Karin et al. 2019 saturating-removal model vs linear and exponential growth
  * compartments   thymic output -> naive -> memory T cells
  * forecasting    predict a held-out age from the others: interpolation baselines, curve fits,
                   a linear ODE, sparse nonlinear ODE discovery (SINDy) and a neural ODE

Every setting is fixed here before results are seen.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, least_squares


# ----------------------------------------------------------------------------
# 7a: accumulation (Karin et al. 2019)
# ----------------------------------------------------------------------------
def karin_mean(t, eta, beta, kappa):
    """Mean-field saturating-removal model: dX/dt = eta*t - beta*X/(kappa + X), X(0) = 0.
    Production of damaged/senescent cells rises linearly with age; removal saturates,
    so X accelerates late in life (the mechanism Karin et al. link to the Gompertz law)."""
    t = np.asarray(t, float)
    grid = np.unique(np.concatenate([[0.0], t]))
    # LSODA switches to a stiff method when removal (beta/kappa) is fast; RK45 would crawl
    sol = solve_ivp(lambda s, x: [eta * s - beta * x[0] / (kappa + max(x[0], 0.0))], (0, grid[-1]), [0.0],
                    t_eval=grid, method="LSODA", rtol=1e-6, atol=1e-9)
    return np.interp(t, sol.t, sol.y[0])


def fit_growth(t, y):
    """Fit linear, exponential and Karin curves to (age, value); return AIC and parameters.
    Each curve has a free offset a: y = a + f(t)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    n = len(y)

    def aic(resid, k):
        rss = max(float(np.sum(resid ** 2)), 1e-12)
        return n * np.log(rss / n) + 2 * k

    out = {}
    b = np.polyfit(t, y, 1)
    out["linear"] = {"aic": aic(y - np.polyval(b, t), 2), "params": b}
    try:
        p, _ = curve_fit(lambda s, a, c, r: a + c * (np.exp(r * s) - 1), t, y,
                         p0=[y.min(), 0.01, 0.1], bounds=([-np.inf, -np.inf, 0], [np.inf, np.inf, 1]), maxfev=20000)
        out["exponential"] = {"aic": aic(y - (p[0] + p[1] * (np.exp(p[2] * t) - 1)), 3), "params": p}
    except RuntimeError:
        out["exponential"] = {"aic": np.nan, "params": None}

    def resid(q):
        return y - (q[0] + karin_mean(t, q[1], q[2], q[3]))

    best = None
    for eta0, beta0, kappa0 in [(0.01, 0.1, 1.0), (0.05, 1.0, 1.0), (0.001, 0.05, 0.1)]:
        r = least_squares(resid, [y.min(), eta0, beta0, kappa0],
                          bounds=([-np.inf, 0, 0, 1e-2], [np.inf, 10, 10, 100]), max_nfev=400)
        if best is None or r.cost < best.cost:
            best = r
    out["karin"] = {"aic": aic(best.fun, 4), "params": best.x}
    return out


def karin_sde(eta, beta, kappa, eps, t_end, n=2000, dt=0.05, seed=0):
    """Euler-Maruyama simulation of dX = (eta t - beta X/(kappa+X)) dt + sqrt(2 eps) dW,
    for n individuals. Returns times and an (n x times) array; used to show that the model
    predicts the spread between individuals to grow with age."""
    rng = np.random.default_rng(seed)
    steps = int(t_end / dt)
    x = np.zeros(n)
    ts, xs = [0.0], [x.copy()]
    for k in range(1, steps + 1):
        t = k * dt
        x = np.maximum(x + (eta * t - beta * x / (kappa + x)) * dt + np.sqrt(2 * eps * dt) * rng.normal(size=n), 0)
        ts.append(t)
        xs.append(x.copy())
    return np.array(ts), np.stack(xs, axis=1)


# ----------------------------------------------------------------------------
# 7c: T-cell compartments
# ----------------------------------------------------------------------------
def tcell_fraction(t, theta, delta, k, f0):
    """Naive fraction f of a T-cell pool held at constant size (homeostasis):

        df/dt = theta * exp(-delta t) - k f,   f(0) = f0

    theta = thymic output relative to the pool, delta = rate of thymic involution,
    k = rate at which naive cells are lost (converted to memory or dying beyond
    replacement). Exact solution, so it stays finite and between the curves it mixes.
    (A first version tracked naive and memory pools separately; with only fractions
    observed it drove both pools to zero and their ratio to numerical noise.)"""
    t = np.asarray(t, float)
    if abs(k - delta) < 1e-8:
        return f0 * np.exp(-k * t) + theta * t * np.exp(-k * t)
    return f0 * np.exp(-k * t) + theta * (np.exp(-delta * t) - np.exp(-k * t)) / (k - delta)


def fit_tcell(t, f, f0, fixed_delta=None):
    """Least-squares fit of theta, delta, k (f0 fixed from the youngest mice), from
    several starting points. With fixed_delta, delta is held there (profile likelihood)."""
    t, f = np.asarray(t, float), np.asarray(f, float)

    def unpack(q):
        return (q[0], fixed_delta, q[1]) if fixed_delta is not None else tuple(q)

    best = None
    for theta0 in [0.01, 0.1, 1.0]:
        for k0 in [0.01, 0.1, 1.0]:
            x0 = [theta0, k0] if fixed_delta is not None else [theta0, 0.1, k0]
            lo = [0, 1e-4] if fixed_delta is not None else [0, 0, 1e-4]
            hi = [10, 10] if fixed_delta is not None else [10, 10, 10]
            r = least_squares(lambda q: f - tcell_fraction(t, *unpack(q), f0), x0, bounds=(lo, hi))
            if best is None or r.cost < best.cost:
                best = r
    return best


def residualize(X, covariates, train):
    """Remove covariate effects (e.g. sex, plate) from every column of X by least squares,
    fitted on the `train` rows only and applied to all rows, so held-out rows don't leak in."""
    import pandas as pd
    D = pd.get_dummies(covariates.astype(str), drop_first=True).to_numpy(float)
    D = np.hstack([np.ones((len(D), 1)), D])
    B, *_ = np.linalg.lstsq(D[train], X[train], rcond=None)
    return X - D[:, 1:] @ B[1:]


# ----------------------------------------------------------------------------
# 7e: forecasting a held-out age
# ----------------------------------------------------------------------------
def _sigmoid(t, a, b, t0, s):
    return a + b / (1 + np.exp(-(t - t0) / s))


def forecast(ages_train, means_train, t_new, method, data_t=None, data_y=None, seed=0):
    """Predict the program vector at age t_new from per-age mean vectors (ages x dims).

    nearest      the closest training age
    linear       linear interpolation (linear extrapolation outside the range)
    spline       smoothing spline per dimension (GCV-chosen smoothness)
    sigmoid      a + b / (1 + exp(-(t - t0)/s)) per dimension, fit to per-mouse data
    linear_ode   x' = A x + b from finite differences, integrated from the previous age
    sindy        sparse polynomial ODE (degree 2, STLSQ), integrated from the previous age
    neural_ode   small neural vector field fit to the whole trajectory, integrated from the previous age
    """
    T = np.asarray(ages_train, float)
    Y = np.asarray(means_train, float)
    prev = T[T < t_new].max() if (T < t_new).any() else T.min()
    y_prev = Y[T == prev][0]

    if method == "nearest":
        d = np.abs(T - t_new)
        return Y[d == d.min()].mean(axis=0)
    if method == "linear":
        if t_new > T.max():
            return Y[-1] + (Y[-1] - Y[-2]) / (T[-1] - T[-2]) * (t_new - T[-1])
        return np.array([np.interp(t_new, T, Y[:, j]) for j in range(Y.shape[1])])
    if method == "spline":
        from scipy.interpolate import make_smoothing_spline
        return np.array([make_smoothing_spline(T, Y[:, j])(t_new) for j in range(Y.shape[1])])
    if method == "sigmoid":
        out = []
        for j in range(Y.shape[1]):
            yj = data_y[:, j]
            try:
                p, _ = curve_fit(_sigmoid, data_t, yj, p0=[yj.min(), np.ptp(yj), np.median(data_t), 3.0],
                                 bounds=([-np.inf, -np.inf, 0, 0.3], [np.inf, np.inf, 40, 30]), maxfev=20000)
                out.append(_sigmoid(t_new, *p))
            except RuntimeError:
                out.append(np.polyval(np.polyfit(data_t, yj, 1), t_new))
        return np.array(out)
    if method == "linear_ode":
        dY = np.diff(Y, axis=0) / np.diff(T)[:, None]
        Xm = (Y[1:] + Y[:-1]) / 2
        D = np.hstack([Xm, np.ones((len(Xm), 1))])
        W = np.linalg.solve(D.T @ D + 1e-2 * np.eye(D.shape[1]), D.T @ dY)     # ridge
        A, b = W[:-1].T, W[-1]
        sol = solve_ivp(lambda s, x: A @ x + b, (prev, t_new), y_prev, rtol=1e-6)
        return sol.y[:, -1]
    if method == "sindy":
        import pysindy as ps
        model = ps.SINDy(feature_library=ps.PolynomialLibrary(degree=2),
                         optimizer=ps.STLSQ(threshold=0.02, alpha=0.05))
        model.fit(Y, t=T)
        try:
            out = model.simulate(y_prev, np.linspace(prev, t_new, 50))[-1]
        except ValueError:                       # the discovered ODE blew up to infinity
            return np.full(Y.shape[1], np.nan)
        # a diverging solution is a failure, reported as NaN (not replaced by a baseline)
        return out if np.all(np.isfinite(out)) and np.abs(out).max() < 1e3 else np.full(Y.shape[1], np.nan)
    if method == "neural_ode":
        import torch
        from torchdiffeq import odeint
        torch.manual_seed(seed)
        dim = Y.shape[1]
        f = torch.nn.Sequential(torch.nn.Linear(dim, 16), torch.nn.Tanh(), torch.nn.Linear(16, dim))
        tt = torch.tensor(T, dtype=torch.float32)
        yy = torch.tensor(Y, dtype=torch.float32)
        opt = torch.optim.Adam(f.parameters(), lr=1e-2, weight_decay=1e-3)
        for _ in range(400):
            pred = odeint(lambda s, x: f(x), yy[0], tt, method="rk4", options={"step_size": 0.5})
            loss = ((pred - yy) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
        with torch.no_grad():
            out = odeint(lambda s, x: f(x), torch.tensor(y_prev, dtype=torch.float32),
                         torch.tensor([prev, t_new], dtype=torch.float32), method="rk4",
                         options={"step_size": 0.5})[-1].numpy()
        return out
    raise ValueError(method)


FORECASTERS = ["nearest", "linear", "spline", "sigmoid", "linear_ode", "sindy", "neural_ode"]
