"""Age models for the model comparison (11_models), all with one interface.

A classifier here is an object with
    fit(X, y, mice)      X cells x genes (log-normalized), y 0/1 (old), mice = mouse id per cell
    predict(X) -> score  probability-like "old" score per cell
so models that need the mouse structure (pseudobulk, multiple-instance learning)
and models that don't are evaluated by exactly the same code.

Every setting is fixed here, before any result is seen; nothing is tuned on a test set.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# ----------------------------------------------------------------------------
# Classifiers (young vs. old)
# ----------------------------------------------------------------------------
class SklearnCells:
    """Wrap a scikit-learn classifier: every cell is a training example."""

    def __init__(self, make):
        self.make = make

    def fit(self, X, y, mice):
        self.model = self.make().fit(X, y)
        return self

    def predict(self, X):
        return self.model.predict_proba(X)[:, 1]


def logistic():
    """The baseline of 05a: L2 logistic regression, C = 0.1."""
    return SklearnCells(lambda: LogisticRegression(C=0.1, max_iter=2000))


def elastic_net():
    """Sparse linear model: an equal mix of L1 (drops genes) and L2 (shrinks them), at the
    baseline's C = 0.1. On Lung fibroblasts it keeps about 44 of ~2,000 genes: a short panel.
    (A first version used stochastic gradient descent, which did not converge in a smoke
    test; it was replaced by the exact saga solver before the comparison was run.)"""
    return SklearnCells(lambda: LogisticRegression(l1_ratio=0.5, C=0.1, solver="saga", max_iter=2000))


def lightgbm():
    """Gradient-boosted trees, the same settings as 05a."""
    from lightgbm import LGBMClassifier
    return SklearnCells(lambda: LGBMClassifier(n_estimators=200, learning_rate=0.05, num_leaves=15,
                                               subsample=0.8, subsample_freq=1, colsample_bytree=0.5,
                                               verbose=-1))


class MouseLevel:
    """One training example per mouse: the mean expression of its cells.
    Matches the real sample size; applied to single cells at prediction time."""

    def fit(self, X, y, mice):
        df = pd.DataFrame(X).assign(_m=mice, _y=y).groupby("_m").mean()
        self.model = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000))
        self.model.fit(df.drop(columns="_y").to_numpy(), df["_y"].round().astype(int).to_numpy())
        return self

    def predict(self, X):
        return self.model.predict_proba(X)[:, 1]


def mouse_level():
    return MouseLevel()


class AttentionMIL:
    """Attention-based multiple-instance learning (Ilse et al. 2018).

    A mouse is a bag of cells with one label. Each cell is encoded by a small network,
    an attention layer learns how much each cell should count, and the weighted
    average is classified. With ~10 training mice, each epoch draws many random
    sub-bags of `bag_size` cells per mouse, so the model sees varied bags but can
    never use more independent labels than there are mice."""

    def __init__(self, hidden=64, emb=32, bag_size=64, bags_per_mouse=8, epochs=40,
                 lr=1e-3, weight_decay=1e-3, dropout=0.3, seed=0):
        self.__dict__.update(locals())
        del self.__dict__["self"]

    def _net(self, n_genes):
        import torch.nn as nn

        class Net(nn.Module):
            def __init__(s):
                super().__init__()
                s.enc = nn.Sequential(nn.Linear(n_genes, self.hidden), nn.ReLU(), nn.Dropout(self.dropout),
                                      nn.Linear(self.hidden, self.emb), nn.ReLU())
                s.att = nn.Sequential(nn.Linear(self.emb, 16), nn.Tanh(), nn.Linear(16, 1))
                s.cls = nn.Linear(self.emb, 1)

            def forward(s, x):                       # x: bags x cells x genes
                h = s.enc(x)
                w = s.att(h).softmax(dim=1)          # attention over the cells of each bag
                return s.cls((w * h).sum(dim=1)).squeeze(-1), w.squeeze(-1)

        return Net()

    def fit(self, X, y, mice):
        import torch
        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)
        self.mu, self.sd = X.mean(axis=0), X.std(axis=0) + 1e-6
        Xs = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32)
        ids = np.unique(mice)
        idx = {m: np.where(mice == m)[0] for m in ids}
        lab = {m: float(y[mice == m][0]) for m in ids}
        pos = np.mean([lab[m] for m in ids])
        weight = {1.0: 0.5 / max(pos, 1e-6), 0.0: 0.5 / max(1 - pos, 1e-6)}   # balance old and young mice
        self.net = self._net(X.shape[1])
        opt = torch.optim.Adam(self.net.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        loss_fn = torch.nn.BCEWithLogitsLoss(reduction="none")
        self.net.train()
        for _ in range(self.epochs):
            bags, labels, weights = [], [], []
            for m in ids:
                for _ in range(self.bags_per_mouse):
                    pick = rng.choice(idx[m], size=self.bag_size, replace=len(idx[m]) < self.bag_size)
                    bags.append(Xs[pick])
                    labels.append(lab[m])
                    weights.append(weight[lab[m]])
            xb = torch.stack(bags)
            yb, wb = torch.tensor(labels), torch.tensor(weights)
            order = torch.randperm(len(xb))
            for k in range(0, len(xb), 16):
                b = order[k:k + 16]
                logit, _ = self.net(xb[b])
                loss = (loss_fn(logit, yb[b]) * wb[b]).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
        return self

    def predict(self, X):
        """Score each cell as a bag of one cell."""
        import torch
        self.net.eval()
        with torch.no_grad():
            x = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32).unsqueeze(1)
            logit, _ = self.net(x)
        return torch.sigmoid(logit).numpy()


def attention_mil():
    return AttentionMIL()


CLASSIFIERS = {
    "logistic (baseline)": logistic,
    "elastic net": elastic_net,
    "LightGBM": lightgbm,
    "mouse-level": mouse_level,
    "attention MIL": attention_mil,
}


# ----------------------------------------------------------------------------
# Clocks (age in months)
# ----------------------------------------------------------------------------
class OrdinalClock:
    """Ordinal regression by cumulative binary models (Frank & Hall 2001).

    For each adult age a_k after the youngest, a logistic model estimates P(age >= a_k);
    the predicted age is a_0 + sum_k P(age >= a_k) * (a_k - a_(k-1)). It uses only the
    order of the ages, not equal spacing, so a batch at one age can't stretch the scale."""

    def __init__(self, C=0.1):
        self.C = C

    def fit(self, X, y):
        self.ages = np.sort(np.unique(y))
        self.models = [LogisticRegression(C=self.C, max_iter=2000).fit(X, (y >= a).astype(int))
                       for a in self.ages[1:]]
        return self

    def predict(self, X):
        out = np.full(len(X), self.ages[0], dtype=float)
        for a_prev, a, m in zip(self.ages[:-1], self.ages[1:], self.models):
            out += m.predict_proba(X)[:, 1] * (a - a_prev)
        return out


def clock_models() -> dict:
    """Clock regressors with fixed settings (the 05c ridge first)."""
    from lightgbm import LGBMRegressor
    from sklearn.decomposition import PCA
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.neural_network import MLPRegressor
    return {
        "ridge (05c)": lambda: Ridge(alpha=100),
        "ordinal": lambda: OrdinalClock(C=0.1),
        "kernel ridge RBF": lambda: make_pipeline(PCA(50, random_state=0), StandardScaler(),
                                                  KernelRidge(alpha=1.0, kernel="rbf", gamma=1 / 50)),
        "LightGBM": lambda: LGBMRegressor(n_estimators=300, learning_rate=0.05, num_leaves=15, subsample=0.8,
                                          subsample_freq=1, colsample_bytree=0.5, verbose=-1),
        "MLP": lambda: make_pipeline(PCA(50, random_state=0), StandardScaler(),
                                     MLPRegressor(hidden_layer_sizes=(64, 32), alpha=1e-2, max_iter=300,
                                                  early_stopping=True, random_state=0)),
    }
