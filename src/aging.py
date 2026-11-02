"""Shared helpers for the single-cell aging analysis (Tabula Muris Senis).

Every function takes the AnnData object as its first argument, so notebooks
stay short and the same code runs on any tissue.

Expected AnnData (as written by 03_preprocess.ipynb):
    .X                         log-normalized expression
    .layers["counts"]          raw counts
    .var["highly_variable"]    from sc.pp.highly_variable_genes
    .obs                       "age", "mouse.id", "sex", "cell_ontology_class"
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
AGE_MONTHS = {"1m": 1, "3m": 3, "18m": 18, "21m": 21, "24m": 24, "30m": 30}
YOUNG = ["1m", "3m"]
OLD = ["18m", "21m", "24m", "30m"]   # ages missing from a tissue are simply ignored
ADULT = ["3m", "18m", "21m", "24m", "30m"]

# Genes likely to come from ambient RNA or blood rather than the cell itself
AMBIENT = ["S100a8", "S100a9", "Lyz2", "Ngp", "Camp", "Ltf", "Retnlg",
           "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Beta-s",   # TMS names Hbb "Beta-s"
           "Scgb1a1"]   # airway secretory gene, slips past foreign_genes (see 05b findings)

# Technical and sex genes: ribosomal proteins (incl. pseudogenes), rRNA, Malat1,
# X/Y-linked sex genes, hemoglobin
TECHNICAL_PATTERN = r"^(Rp[sl]\d|Rn45s|Malat1|Xist|Tsix|Ddx3y|Eif2s3y|Kdm5d|Uty|Beta-s|Hb[ab]-)"

# Genes induced within minutes by tissue dissociation (not aging)
DISSOCIATION = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2", "Ier3",
                "Dusp1", "Hspa1a", "Hspa1b", "Hsp90aa1", "Hspa8", "Zfp36", "Klf2", "Klf4", "Nr4a1"]

# Circadian clock genes: differences point to different times of day at collection (batch)
CIRCADIAN = ["Dbp", "Nr1d1", "Nr1d2", "Per1", "Per2", "Per3", "Arntl", "Cry1", "Tef"]

MYELOID = ["classical monocyte", "intermediate monocyte", "non-classical monocyte",
           "lung macrophage", "alveolar macrophage"]

# Mice excluded from specific cell types, with the reason (Lung findings)
EXCLUDE = {"30-M-2": {"cell_types": MYELOID, "reason": "abnormal myeloid expansion"}}


def main_kw(ct: str) -> dict:
    """Keyword arguments for the main analysis of one cell type (applies EXCLUDE)."""
    drop = [m for m, rule in EXCLUDE.items() if ct in rule["cell_types"]]
    return {"drop_mice": drop} if drop else {}


# ----------------------------------------------------------------------------
# Gene filtering
# ----------------------------------------------------------------------------
def celltype_means(adata) -> pd.DataFrame:
    """Mean log-expression of every gene in every cell type (cell types x genes)."""
    pb = sc.get.aggregate(adata, by="cell_ontology_class", func="mean")
    m = pb.layers["mean"]
    m = m.toarray() if hasattr(m, "toarray") else np.asarray(m)
    return pd.DataFrame(m, index=pb.obs["cell_ontology_class"].astype(str), columns=pb.var_names)


def foreign_genes(ct_means: pd.DataFrame, ct: str, ratio: float = 0.1, min_other: float = 1.0) -> list:
    """Genes that belong to another cell type: clearly expressed somewhere else
    (mean log-expression > min_other) but < ratio of that level in `ct`.
    In `ct` they are most likely ambient RNA (e.g. Scgb1a1 from club cells in T cells)."""
    own = ct_means.loc[ct]
    other = ct_means.drop(index=ct).max()
    return list(ct_means.columns[(other > min_other) & (own < ratio * other)])


def feature_mask(adata, level: str = "bio", ct: str | None = None,
                 ct_means: pd.DataFrame | None = None, hvg_only: bool = True) -> np.ndarray:
    """Boolean mask over adata.var_names.

    level="hvg"   highly variable genes only
    level="clean" + drop AMBIENT genes
    level="bio"   + drop technical/sex genes
    level="own"   + drop other cell types' marker genes (needs ct and ct_means)
    """
    names = adata.var_names
    mask = adata.var["highly_variable"].to_numpy().copy() if hvg_only else np.ones(len(names), bool)
    if level == "hvg":
        return mask
    mask &= ~names.isin(AMBIENT)
    if level == "clean":
        return mask
    mask &= ~names.str.match(TECHNICAL_PATTERN)
    if level == "bio":
        return mask
    if level == "own":
        if ct is None or ct_means is None:
            raise ValueError('level="own" needs ct and ct_means')
        return mask & ~names.isin(foreign_genes(ct_means, ct))
    raise ValueError(f"unknown level: {level}")


# ----------------------------------------------------------------------------
# Data prep
# ----------------------------------------------------------------------------
def prep(adata, ct, ages=None, target="binary", features=None, max_cells=300, seed=0,
         keep_mice=None, drop_mice=None):
    """Cells of one type, at most `max_cells` per mouse.

    target="binary": y = 1 for OLD ages, 0 otherwise
    target="months": y = age in months
    Returns X (dense), y, mouse ids.
    """
    ages = (YOUNG + OLD) if ages is None else ages
    obs = adata.obs
    mask = (obs.cell_ontology_class == ct) & obs.age.isin(ages)
    if keep_mice is not None:
        mask &= obs["mouse.id"].isin(keep_mice)
    if drop_mice is not None:
        mask &= ~obs["mouse.id"].isin(drop_mice)
    idx = (obs[mask].sample(frac=1, random_state=seed)
                    .groupby("mouse.id", observed=True).head(max_cells).index)
    sub = adata[idx]
    feat = adata.var["highly_variable"].to_numpy() if features is None else features
    X = sub[:, feat].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    if target == "binary":
        y = sub.obs.age.isin(OLD).to_numpy().astype(int)
    elif target == "months":
        y = sub.obs.age.astype(str).map(AGE_MONTHS).to_numpy().astype(float)
    else:
        raise ValueError(target)
    mice = sub.obs["mouse.id"].astype(str).to_numpy()
    return X, y, mice


def usable_cell_types(adata, min_cells=50, min_mice=2) -> list:
    """Cell types with >= min_cells and >= min_mice in both young and old groups."""
    obs = adata.obs[adata.obs.age.isin(YOUNG + OLD)].copy()
    obs["group"] = np.where(obs.age.isin(OLD), "old", "young")
    size = obs.groupby(["cell_ontology_class", "group"], observed=True).agg(
        cells=("mouse.id", "size"), mice=("mouse.id", "nunique")).unstack()
    ok = size[(size["cells"] >= min_cells).all(axis=1) & (size["mice"] >= min_mice).all(axis=1)]
    return list(ok.index)


# ----------------------------------------------------------------------------
# Young vs. old classifier
# ----------------------------------------------------------------------------
def logreg():
    return LogisticRegression(C=0.1, max_iter=2000)


def _shuffle_by_mouse(y, mice, seed):
    rng = np.random.default_rng(seed)
    ids = np.unique(mice)
    labels = np.array([y[mice == m][0] for m in ids])
    new = dict(zip(ids, rng.permutation(labels)))
    return np.array([new[m] for m in mice])


def age_eval(adata, ct, make_model=logreg, shuffle_seed=None, **prep_kw):
    """Mouse-held-out evaluation. Returns a dict of scores, or None if too few mice."""
    X, y, mice = prep(adata, ct, **prep_kw)
    if shuffle_seed is not None:
        y = _shuffle_by_mouse(y, mice, shuffle_seed)
    n_splits = min(5, len(np.unique(mice[y == 0])), len(np.unique(mice[y == 1])))
    if n_splits < 2:
        return None
    scores = np.zeros(len(y))
    train_acc = []
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=0)
    for tr, te in cv.split(X, y, mice):
        model = make_model().fit(X[tr], y[tr])
        scores[te] = model.predict_proba(X[te])[:, 1]
        train_acc.append(accuracy_score(y[tr], model.predict(X[tr])))
    per_mouse = pd.DataFrame({"mouse": mice, "y": y, "score": scores}).groupby("mouse").mean()
    return {
        "auroc_cell": roc_auc_score(y, scores),
        "auroc_mouse": roc_auc_score(per_mouse["y"], per_mouse["score"]),
        "bal_acc_test": balanced_accuracy_score(y, scores > 0.5),
        "acc_train": float(np.mean(train_acc)),
        "n_cells": len(y),
        "n_mice": len(per_mouse),
    }


def age_auc(adata, ct, **kw):
    r = age_eval(adata, ct, **kw)
    return np.nan if r is None else r["auroc_cell"]


def cross_age_auc(adata, ct, make_model=logreg, splits=None, **prep_kw):
    """Train on some age groups, test on different ones (batch check).
    Default splits fit Lung's ages (1m, 3m, 18m, 21m, 30m)."""
    X, y, mice = prep(adata, ct, **prep_kw)
    ages = np.array([m.split("-")[0] + "m" for m in mice])
    if splits is None:
        splits = [({"1m", "18m", "21m"}, {"3m", "30m"}),
                  ({"3m", "30m"}, {"1m", "18m", "21m"})]
    aucs = []
    for train_ages, test_ages in splits:
        tr, te = np.isin(ages, list(train_ages)), np.isin(ages, list(test_ages))
        if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
            continue
        model = make_model().fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)) if aucs else np.nan


def top_genes(adata, ct, features, n=20, **prep_kw):
    """Strongest logistic-regression weights (fit on all mice; for interpretation only)."""
    X, y, _ = prep(adata, ct, features=features, **prep_kw)
    model = logreg().fit(X, y)
    coef = pd.Series(model.coef_[0], index=adata.var_names[features]).sort_values()
    return pd.DataFrame({"up_in_old": coef.tail(n)[::-1].index, "up_in_young": coef.head(n).index})


# ----------------------------------------------------------------------------
# Age clock (regression)
# ----------------------------------------------------------------------------
def clock_eval(adata, ct, alpha=100.0, ages=None, **prep_kw):
    """Leave-one-mouse-out ridge regression on age in months (adults by default)."""
    ages = ADULT if ages is None else ages
    X, y, mice = prep(adata, ct, ages=ages, target="months", **prep_kw)
    pred = np.zeros(len(y))
    for tr, te in LeaveOneGroupOut().split(X, y, mice):
        pred[te] = Ridge(alpha=alpha).fit(X[tr], y[tr]).predict(X[te])
    cells = pd.DataFrame({"mouse": mice, "true": y, "pred": pred})
    per_mouse = cells.groupby("mouse").agg(true=("true", "first"), pred=("pred", "median"))
    return cells, per_mouse


def clock_summary(adata, ct, **kw):
    _, pm = clock_eval(adata, ct, **kw)
    naive = np.mean([abs(t - pm["true"].drop(m).mean()) for m, t in pm["true"].items()])
    return {
        "mae_months": float((pm["pred"] - pm["true"]).abs().mean()),
        "naive_mae": float(naive),
        "spearman_mice": float(spearmanr(pm["true"], pm["pred"])[0]),
        "n_mice": len(pm),
    }


def held_out_age(adata, ct, hold="21m", alpha=100.0, ages=None, **prep_kw):
    """Train without one age group, predict it. Returns median prediction per mouse."""
    ages = ADULT if ages is None else ages
    X, y, mice = prep(adata, ct, ages=ages, target="months", **prep_kw)
    te = y == AGE_MONTHS[hold]
    model = Ridge(alpha=alpha).fit(X[~te], y[~te])
    return pd.Series(model.predict(X[te]), index=mice[te]).groupby(level=0).median().round(1)


# ----------------------------------------------------------------------------
# Development vs. aging trajectories
# ----------------------------------------------------------------------------
def trajectories(adata, ct, drop_mice=None, exclude=(), min_expr=0.1):
    """Per-gene change 1m -> 3m (development) and 3m -> old (aging), from per-mouse means."""
    sub = adata[adata.obs.cell_ontology_class == ct]
    if drop_mice:
        sub = sub[~sub.obs["mouse.id"].isin(drop_mice)]
    pb = sc.get.aggregate(sub, by="mouse.id", func="mean")
    mat = pb.layers["mean"]
    mat = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
    mat = pd.DataFrame(mat, index=pb.obs["mouse.id"].astype(str), columns=pb.var_names)
    mouse_age = sub.obs.drop_duplicates("mouse.id").set_index("mouse.id")["age"]
    mouse_age.index = mouse_age.index.astype(str)
    grp = mat.groupby(mouse_age.loc[mat.index].astype(str).values).mean()
    keep = (grp.max() > min_expr) & ~grp.columns.isin(list(exclude)) \
        & ~grp.columns.str.match(TECHNICAL_PATTERN) & ~grp.columns.isin(AMBIENT)
    grp = grp.loc[:, keep]
    old_mean = grp.loc[[a for a in OLD if a in grp.index]].mean()
    return pd.DataFrame({"development": grp.loc["3m"] - grp.loc["1m"],
                         "aging": old_mean - grp.loc["3m"]})


def label_patterns(tr, big=0.5, small=0.2):
    d, a = tr["development"], tr["aging"]
    tr = tr.copy()
    tr["pattern"] = np.select(
        [(d.abs() > big) & (a.abs() < small),
         (a.abs() > big) & (d.abs() < small),
         ((d > big) & (a > big)) | ((d < -big) & (a < -big)),
         ((d > big) & (a < -big)) | ((d < -big) & (a > big))],
        ["development only", "aging only", "lifelong", "reverses"], default="")
    return tr


def pattern_top(tr, n=15) -> dict:
    """Top `n` genes of each pattern, by size of the change that defines the pattern."""
    top = {}
    for p in ["development only", "aging only", "lifelong", "reverses"]:
        sel = tr[tr.pattern == p]
        order = (sel["aging"] if p == "aging only" else sel["development"]).abs().sort_values(ascending=False)
        top[p] = list(order.index[:n])
    return top
