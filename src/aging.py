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

from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

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
           "Scgb1a1"]   # airway secretory gene; the linear foreign_genes rule now catches it too,
                        # kept here so the 05a no_ambient_genes column stays comparable

# Technical and sex genes: ribosomal proteins (incl. pseudogenes), rRNA, Malat1,
# X/Y-linked sex genes, hemoglobin
TECHNICAL_PATTERN = r"^(Rp[sl]\d|Rn45s|Malat1|Xist|Tsix|Ddx3y|Eif2s3y|Kdm5d|Uty|Beta-s|Hb[ab]-)"

# Genes induced within minutes by tissue dissociation (not aging)
DISSOCIATION = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2", "Ier3",
                "Dusp1", "Hspa1a", "Hspa1b", "Hsp90aa1", "Hspa8", "Zfp36", "Klf2", "Klf4", "Nr4a1"]

# Circadian clock genes: differences point to different times of day at collection (batch)
CIRCADIAN = ["Dbp", "Nr1d1", "Nr1d2", "Per1", "Per2", "Per3", "Arntl", "Cry1", "Tef"]

# Current symbol -> the older symbol TMS uses (TMS was annotated before these renames)
GENE_ALIASES = {"Bex3": "Ngfrap1", "Hbb-bs": "Beta-s"}


def load_gene_set(path, adata) -> tuple[list, list]:
    """Read a one-gene-per-line file, map current symbols to TMS's older ones,
    and return (genes present in adata, genes missing)."""
    genes = [g.strip() for g in Path(path).read_text().splitlines() if g.strip()]
    mapped = [GENE_ALIASES.get(g, g) for g in genes]
    present = [g for g in mapped if g in adata.var_names]
    missing = [g for g, m in zip(genes, mapped) if m not in adata.var_names]
    return present, missing

# ----------------------------------------------------------------------------
# Per-tissue settings
# ----------------------------------------------------------------------------
# A new tissue needs no entry to run; add one only after reviewing its Phase 2-4
# checks (outlier mice, pooled samples). Keys, all optional:
#   exclude            {mouse: {"cell_types": [...], "reason": str}}  dropped from those cell types
#   odd_mice           samples left out of composition trends
#   inspect_mice       mice plotted one at a time on the UMAP (03)
#   cross_age_splits   hand-picked (train ages, test ages) pairs; the generic
#                      cross_age_splits() is always reported alongside
#   composition_cell_types, example_pass, example_fail, example_clock
#                      cell types shown in 04, 05b and 05c (otherwise picked automatically)
LUNG_MYELOID = ["classical monocyte", "intermediate monocyte", "non-classical monocyte",
                "lung macrophage", "alveolar macrophage"]

TISSUE_SETTINGS = {
    "Lung": {
        "exclude": {"30-M-2": {"cell_types": LUNG_MYELOID, "reason": "abnormal myeloid expansion"}},
        "odd_mice": ["30-M-2", "3-M-5/6", "3-M-7/8"],
        "inspect_mice": ["30-M-2", "30-M-5"],
        "cross_age_splits": [({"1m", "18m", "21m"}, {"3m", "30m"}),
                             ({"3m", "30m"}, {"1m", "18m", "21m"})],
        "composition_cell_types": ["B cell", "CD8-positive, alpha-beta T cell",
                                   "CD4-positive, alpha-beta T cell", "NK cell",
                                   "bronchial smooth muscle cell", "classical monocyte",
                                   "alveolar macrophage", "fibroblast of lung"],
        "example_pass": "fibroblast of lung",
        "example_fail": "CD8-positive, alpha-beta T cell",
        "example_clock": "alveolar macrophage",
    },
}

# Compartment of each TMS cell type label, for comparing tissues (07_atlas). Ambiguous labels
# are checked against marker genes in 07_atlas before they are counted; "other" is not compared.
_COMPARTMENT_LISTS = {
    "stromal": ["adventitial cell", "bronchial smooth muscle cell", "fibroblast of cardiac tissue",
                "fibroblast of lung", "kidney mesangial cell", "mesenchymal stem cell",
                "mesenchymal stem cell of adipose", "smooth muscle cell", "stromal cell",
                "bladder cell"],                       # TMS "bladder cell" = bladder mesenchyme (checked in 07)
    "endothelial": ["endocardial cell", "endothelial cell", "endothelial cell of coronary artery",
                    "endothelial cell of hepatic sinusoid", "fenestrated cell",
                    "kidney cortex artery cell"],       # checked in 07
    "epithelial": ["basal cell", "basal cell of epidermis", "bladder urothelial cell", "brush cell",
                   "epithelial cell of proximal tubule", "hepatocyte", "keratinocyte",
                   "kidney collecting duct principal cell", "kidney distal convoluted tubule epithelial cell",
                   "kidney loop of Henle ascending limb epithelial cell",
                   "kidney loop of Henle thick ascending limb epithelial cell",
                   "kidney proximal convoluted tubule epithelial cell", "podocyte",
                   "luminal epithelial cell of mammary gland", "type II pneumocyte"],
    "immune": ["B cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "Kupffer cell",
               "Langerhans cell", "NK cell", "T cell", "alveolar macrophage", "classical monocyte",
               "intermediate monocyte", "non-classical monocyte", "granulocyte", "leukocyte", "macrophage",
               "mature NK T cell", "monocyte", "myeloid dendritic cell", "naive B cell", "naive T cell",
               "neutrophil", "plasma cell", "immature B cell"],
    "progenitor": ["erythroblast", "granulocytopoietic cell", "hematopoietic precursor cell",
                   "megakaryocyte-erythroid progenitor cell", "precursor B cell", "proerythroblast",
                   "promonocyte", "skeletal muscle satellite cell"],
}
COMPARTMENT = {ct: comp for comp, cts in _COMPARTMENT_LISTS.items() for ct in cts}


def compartment(ct: str) -> str:
    """Compartment of a cell type label ("other" if unlisted, e.g. cardiomyocyte, Schwann cell)."""
    return COMPARTMENT.get(ct, "other")


ROOT = Path(__file__).resolve().parent.parent


def settings(tissue: str) -> dict:
    """Settings for one tissue ({} if the tissue has none yet)."""
    return TISSUE_SETTINGS.get(tissue, {})


def paths(tissue: str) -> SimpleNamespace:
    """File locations for one tissue. figure() creates the tissue's figure folder.

    official   the h5ad as downloaded from figshare
    raw        the same, rewritten once by 02 (silences old-format warnings)
    processed  written by 03
    """
    data = ROOT / "data" / tissue
    fig_dir = ROOT / "figure" / tissue

    def figure(name: str) -> Path:
        fig_dir.mkdir(parents=True, exist_ok=True)
        return fig_dir / name

    return SimpleNamespace(
        data=data,
        official=data / f"tabula-muris-senis-droplet-processed-official-annotations-{tissue}.h5ad",
        raw=data / "raw.h5ad",
        processed=data / "processed.h5ad",
        figure=figure,
    )


def slug(name: str) -> str:
    """File-name-safe short form: "CD8-positive, alpha-beta T cell" -> "CD8-positive", "3-M-5/6" -> "3-M-5_6"."""
    return name.split(",")[0].replace(" ", "_").replace("/", "_")


def main_kw(ct: str, tissue: str) -> dict:
    """Keyword arguments for the main analysis of one cell type (applies the tissue's exclusions)."""
    exclude = settings(tissue).get("exclude", {})
    drop = [m for m, rule in exclude.items() if ct in rule["cell_types"]]
    return {"drop_mice": drop} if drop else {}


def ages_in(adata) -> list:
    """Ages present in the data, youngest first."""
    return sorted(adata.obs["age"].astype(str).unique(), key=AGE_MONTHS.get)


def flag_outlier_mice(adata, share: float = 0.5, min_cells: int = 100) -> pd.DataFrame:
    """Mice that supply more than `share` of a cell type's cells (cell types with >= min_cells).
    A flag is a prompt to look (as with Lung 30-M-2), not a decision to exclude."""
    tab = pd.crosstab(adata.obs["cell_ontology_class"], adata.obs["mouse.id"])
    tab = tab[tab.sum(axis=1) >= min_cells]
    frac = tab.div(tab.sum(axis=1), axis=0)
    top = frac.idxmax(axis=1)
    out = pd.DataFrame({"top_mouse": top, "share": frac.max(axis=1), "n_cells": tab.sum(axis=1),
                        "n_mice": (tab > 0).sum(axis=1)})
    return out[out["share"] > share].sort_values("share", ascending=False)


# ----------------------------------------------------------------------------
# Gene filtering
# ----------------------------------------------------------------------------
def celltype_means(adata) -> pd.DataFrame:
    """Per-cell-type gene means, and the tissue's ambient cutoff (see ambient_cutoff).

    .log     mean log-expression (cell types x genes), used to find a gene's source cell type
    .linear  mean normalized expression (counts per 10k, undoing log1p), used for the ratio
    .cutoff  linear ratio below which a gene counts as ambient in a cell type
    """
    cts = adata.obs["cell_ontology_class"].astype(str).to_numpy()
    names = sorted(set(cts))
    X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else np.asarray(adata.X)
    lin = X.copy()
    if hasattr(lin, "data"):
        lin.data = np.expm1(lin.data)
    else:
        lin = np.expm1(lin)

    def per_type(M):
        return pd.DataFrame(np.vstack([np.asarray(M[cts == c].mean(axis=0)).ravel() for c in names]),
                            index=names, columns=adata.var_names)

    log, linear = per_type(X), per_type(lin)
    hvg = adata.var["highly_variable"].to_numpy() if "highly_variable" in adata.var else None
    return SimpleNamespace(log=log, linear=linear, cutoff=ambient_cutoff(log, linear, genes=hvg))


def _foreign_ratios(ct_means, ct: str, min_other: float = 1.0) -> pd.Series:
    """For genes clearly expressed in another cell type (mean log > min_other there):
    linear level in `ct` divided by the highest linear level in any other cell type."""
    other_log = ct_means.log.drop(index=ct).max()
    cand = other_log.index[other_log > min_other]
    other_lin = ct_means.linear.drop(index=ct)[cand].max()
    return ct_means.linear.loc[ct, cand] / other_lin


def ambient_cutoff(log: pd.DataFrame, linear: pd.DataFrame, genes=None, seed: int = 0) -> float:
    """Ratio separating ambient from real expression, fit on the tissue itself.

    Over every (cell type, gene) pair where the gene is clearly expressed elsewhere,
    log10 of the linear ratio has two groups: ambient RNA at a small, near-constant
    fraction of the source (Lung: ~0.004) and real expression (Lung: ~0.12). A two-component
    Gaussian mixture places the cutoff where a pair becomes more likely real (Lung: 0.031).
    Fit per tissue because contamination differs between tissues and samples."""
    from sklearn.mixture import GaussianMixture
    ns = SimpleNamespace(log=log, linear=linear)
    if genes is not None:
        ns = SimpleNamespace(log=log.loc[:, genes], linear=linear.loc[:, genes])
    v = pd.concat([_foreign_ratios(ns, ct) for ct in ns.log.index])
    v = np.log10(v[v > 1e-6].to_numpy())[:, None]
    gm = GaussianMixture(2, random_state=seed).fit(v)
    hi = int(np.argmax(gm.means_.ravel()))
    grid = np.linspace(v.min(), v.max(), 4001)[:, None]
    return float(10 ** grid[np.argmax(gm.predict_proba(grid)[:, hi] > 0.5), 0])


def foreign_genes(ct_means, ct: str, min_other: float = 1.0) -> list:
    """Genes that belong to another cell type: clearly expressed somewhere else
    (mean log-expression > min_other) but in `ct` at a linear level below the tissue's
    ambient cutoff relative to that source. In `ct` they are most likely ambient RNA
    (e.g. Scgb1a1 from club cells in T cells, Sftpc from type II pneumocytes in neutrophils).

    The ratio is on the linear scale: the earlier log-scale ratio let abundant ambient
    genes pass, because a log compresses a 500-fold difference to about 7-fold."""
    r = _foreign_ratios(ct_means, ct, min_other)
    return list(r.index[r < ct_means.cutoff])


def feature_mask(adata, level: str = "bio", ct: str | None = None,
                 ct_means: SimpleNamespace | None = None, hvg_only: bool = True) -> np.ndarray:
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
         keep_mice=None, drop_mice=None, return_ages=False):
    """Cells of one type, at most `max_cells` per mouse.

    target="binary": y = 1 for OLD ages, 0 otherwise
    target="months": y = age in months
    Returns X (dense), y, mouse ids (and each cell's age label if return_ages).
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
    if return_ages:
        return X, y, mice, sub.obs["age"].astype(str).to_numpy()
    return X, y, mice


def usable_cell_types(adata, min_cells=50, min_mice=2, min_ages=1) -> list:
    """Cell types with >= min_cells, >= min_mice and >= min_ages distinct ages in both
    young and old groups. min_ages=2 is what the cross-age batch check needs."""
    obs = adata.obs[adata.obs.age.isin(YOUNG + OLD)].copy()
    obs["group"] = np.where(obs.age.isin(OLD), "old", "young")
    obs["age_str"] = obs["age"].astype(str)
    size = obs.groupby(["cell_ontology_class", "group"], observed=True).agg(
        cells=("mouse.id", "size"), mice=("mouse.id", "nunique"), ages=("age_str", "nunique")).unstack()
    # A group absent from the whole tissue (e.g. no young mice) must count as 0, not be skipped
    size = size.reindex(columns=pd.MultiIndex.from_product([["cells", "mice", "ages"], ["old", "young"]]),
                        fill_value=0).fillna(0)
    ok = size[(size["cells"] >= min_cells).all(axis=1) & (size["mice"] >= min_mice).all(axis=1)
              & (size["ages"] >= min_ages).all(axis=1)]
    return list(ok.index)


def cross_age_splits(ages) -> list:
    """Every (train ages, test ages) split for the batch check, from the ages present.

    One young age goes to train with a non-empty, proper subset of the old ages; the
    other young age and the remaining old ages go to test. So no age group is ever in
    both, and each side has young and old. Needs 2 young and >= 2 old ages; returns []
    otherwise (e.g. a tissue with only 3m young mice is "not testable")."""
    young = [a for a in YOUNG if a in set(ages)]
    old = [a for a in OLD if a in set(ages)]
    if len(young) < 2 or len(old) < 2:
        return []
    splits = []
    for y_train in young:
        for k in range(1, len(old)):
            for subset in combinations(old, k):
                train = {y_train, *subset}
                test = (set(young) - {y_train}) | (set(old) - set(subset))
                splits.append((train, test))
    return splits


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
    splits=None uses every split from cross_age_splits() for the ages this cell type has;
    pass a tissue's hand-picked splits (settings(tissue)["cross_age_splits"]) to reproduce them.
    NaN means not testable (no split with young and old on both sides)."""
    X, y, mice, ages = prep(adata, ct, return_ages=True, **prep_kw)
    if splits is None:
        splits = cross_age_splits(set(ages))
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


def technical_adata(adata) -> sc.AnnData:
    """Same cells, but the features are only technical: log total counts, log genes per cell,
    and the fraction of counts from ribosomal, dissociation and hemoglobin genes (z-scored).

    Run through age_auc / cross_age_auc like expression. If these alone separate young
    from old on unseen age groups, the cell type's age signal can't be told apart from
    a processing difference. Old cells can differ technically for real reasons (RNA content),
    so a high score is a warning, not proof of batch."""
    import anndata as ad
    C = adata.layers["counts"]
    C = C.tocsr() if hasattr(C, "tocsr") else np.asarray(C)
    total = np.asarray(C.sum(axis=1)).ravel()
    names = adata.var_names

    def frac(mask):
        return np.asarray(C[:, np.where(mask)[0]].sum(axis=1)).ravel() / total

    feats = pd.DataFrame({
        "log_total": np.log10(total),
        "log_genes": np.log10(np.asarray((C > 0).sum(axis=1)).ravel()),
        "ribosomal": frac(names.str.match(r"^Rp[sl]\d")),
        "dissociation": frac(names.isin(DISSOCIATION)),
        "hemoglobin": frac(names.str.match(r"^(Hb[ab]-|Beta-s)")),
    }, index=adata.obs_names)
    feats = (feats - feats.mean()) / feats.std().replace(0, 1)
    tech = ad.AnnData(feats.to_numpy(), obs=adata.obs.copy(),
                      var=pd.DataFrame(index=feats.columns))
    tech.var["highly_variable"] = True
    return tech


# ----------------------------------------------------------------------------
# Pseudobulk (mice as samples)
# ----------------------------------------------------------------------------
def pseudobulk(adata, ct, drop_mice=None, min_cells=20, ages=None):
    """Raw counts of one cell type summed per mouse (mice with >= min_cells only).

    Returns (counts: mice x genes DataFrame of ints, meta: DataFrame with age, months,
    group "old"/"young", sex, n_cells)."""
    ages = (YOUNG + OLD) if ages is None else ages
    obs = adata.obs
    mask = (obs.cell_ontology_class == ct) & obs.age.isin(ages)
    if drop_mice:
        mask &= ~obs["mouse.id"].isin(drop_mice)
    sub = adata[mask.to_numpy()]
    n = sub.obs["mouse.id"].astype(str).value_counts()
    keep = n.index[n >= min_cells]
    sub = sub[sub.obs["mouse.id"].astype(str).isin(keep).to_numpy()]
    pb = sc.get.aggregate(sub, by="mouse.id", func="sum", layer="counts")
    m = pb.layers["sum"]
    m = m.toarray() if hasattr(m, "toarray") else np.asarray(m)
    mice = pb.obs["mouse.id"].astype(str).to_numpy()
    counts = pd.DataFrame(np.round(m).astype(np.int64), index=mice, columns=pb.var_names)
    first = sub.obs.drop_duplicates("mouse.id").set_index(sub.obs.drop_duplicates("mouse.id")["mouse.id"].astype(str))
    meta = pd.DataFrame({"age": first.loc[mice, "age"].astype(str).to_numpy(),
                         "sex": first.loc[mice, "sex"].astype(str).to_numpy(),
                         "n_cells": n.loc[mice].to_numpy()}, index=mice)
    meta["months"] = meta["age"].map(AGE_MONTHS)
    meta["group"] = np.where(meta["age"].isin(OLD), "old", "young")
    return counts, meta


def deseq(counts: pd.DataFrame, meta: pd.DataFrame, design: str = "~group",
          contrast=("group", "old", "young"), min_total: int = 10) -> pd.DataFrame:
    """PyDESeq2 on a pseudobulk matrix (genes with < min_total counts dropped).
    Returns log2FoldChange, lfcSE, pvalue, padj and baseMean per gene."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    counts = counts.loc[:, counts.sum(axis=0) >= min_total]
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=list(contrast), quiet=True)
    stats.summary()
    return stats.results_df


def log_cpm_change(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.Series:
    """Fast per-gene effect: mean log2(CPM + 1) of old mice minus young mice.
    Used where DESeq2 is too slow to repeat (permutation nulls)."""
    cpm = np.log2(counts.div(counts.sum(axis=1), axis=0) * 1e6 + 1)
    old = meta["group"].eq("old").to_numpy()
    return cpm[old].mean() - cpm[~old].mean()


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


def clock_summary(adata, ct, min_mice=3, **kw):
    """Per-mouse error of the clock against the naive baseline. NaN with fewer than
    `min_mice` adult mice (leave-one-mouse-out needs others to train on)."""
    ages = kw.get("ages") or ADULT
    obs = adata.obs[(adata.obs.cell_ontology_class == ct) & adata.obs.age.isin(ages)]
    drop = kw.get("drop_mice") or []
    if obs.loc[~obs["mouse.id"].isin(drop), "mouse.id"].nunique() < min_mice:
        return {"mae_months": np.nan, "naive_mae": np.nan, "spearman_mice": np.nan, "n_mice": 0}
    _, pm = clock_eval(adata, ct, **kw)
    naive = np.mean([abs(t - pm["true"].drop(m).mean()) for m, t in pm["true"].items()])
    return {
        "mae_months": float((pm["pred"] - pm["true"]).abs().mean()),
        "naive_mae": float(naive),
        "spearman_mice": float(spearmanr(pm["true"], pm["pred"])[0]),
        "n_mice": len(pm),
    }


def held_out_age(adata, ct, hold="21m", alpha=100.0, ages=None, **prep_kw):
    """Train without one age group, predict it. Returns median prediction per mouse
    (empty if the cell type lacks that age, or fewer than two other ages are left to train on)."""
    ages = ADULT if ages is None else ages
    X, y, mice = prep(adata, ct, ages=ages, target="months", **prep_kw)
    te = y == AGE_MONTHS[hold]
    if not te.any() or len(np.unique(y[~te])) < 2:
        return pd.Series(dtype=float)
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
    old = [a for a in OLD if a in grp.index]
    nan = pd.Series(np.nan, index=grp.columns)
    # Without 1m there is no development step; without 3m, aging is measured from 1m.
    young_ref = "3m" if "3m" in grp.index else ("1m" if "1m" in grp.index else None)
    development = grp.loc["3m"] - grp.loc["1m"] if {"1m", "3m"} <= set(grp.index) else nan
    aging = grp.loc[old].mean() - grp.loc[young_ref] if (old and young_ref) else nan
    return pd.DataFrame({"development": development, "aging": aging})


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
