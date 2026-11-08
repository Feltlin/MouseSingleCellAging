"""Figure style for the README figures.

The house style of the author's other figures (gRNA candidate map, LoDEI 3'UTR rankings),
tightened: a bold flush-left title naming what is plotted ("Subject · scope"), one grey line
of facts and keys below it (no rule), bold panel titles with a grey one-line caption, value labels on the
marks, right-hand numeric columns with their own header, black axis lines with outward
ticks, and grey footnotes for method detail. Colour: brick red and slate blue (validated
for colour-vision deficiency, protan dE 16.3), warm grey for "nothing", ink for structure.

    from src import style as st
    st.lung_batch_check(res, P.figure("age_auroc.png"))
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

INK = "#111111"
INK2 = "#5f5d58"        # secondary text: subtitles, captions, numeric columns
MUTED = "#8f8c85"       # context marks, footnotes
RULE = "#d9d6ce"        # light gridlines
NEUTRAL = "#eeece6"     # "nothing": the diverging midpoint, box fill
RED = "#b8463c"
BLUE = "#2d62a3"
WHITE = "#ffffff"
RED_ARM = ["#ecc2bc", "#d4867c", RED, "#8a2f27"]          # light -> dark
BLUE_ARM = ["#a9c0de", BLUE, "#1f4677"]                    # light -> dark

SHORT = {"CD8-positive, alpha-beta T cell": "CD8 T cell", "CD4-positive, alpha-beta T cell": "CD4 T cell",
         "fibroblast of lung": "fibroblast", "fibroblast of cardiac tissue": "cardiac fibroblast",
         "kidney loop of Henle thick ascending limb epithelial cell": "loop of Henle, thick ascending",
         "kidney loop of Henle ascending limb epithelial cell": "loop of Henle, ascending",
         "kidney proximal convoluted tubule epithelial cell": "proximal convoluted tubule",
         "kidney distal convoluted tubule epithelial cell": "distal convoluted tubule",
         "kidney collecting duct principal cell": "collecting duct principal cell",
         "epithelial cell of proximal tubule": "proximal tubule epithelium",
         "endothelial cell of coronary artery": "coronary artery endothelium",
         "endothelial cell of hepatic sinusoid": "hepatic sinusoid endothelium",
         "megakaryocyte-erythroid progenitor cell": "megakaryocyte-erythroid progenitor",
         "kidney mesangial cell": "mesangial cell", "kidney cortex artery cell": "cortex artery cell",
         "skeletal muscle satellite cell": "satellite cell", "bronchial smooth muscle cell": "bronchial smooth muscle"}
TISSUE = {"Heart_and_Aorta": "Heart", "Limb_Muscle": "Muscle", "Small_Intestine": "Small intestine",
          "Brown_Fat": "Brown fat", "Gonadal_Fat": "Gonadal fat", "Mesenteric_Fat": "Mesenteric fat",
          "Subcutaneous_Fat": "Subcutaneous fat"}


def use() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": INK, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3, "ytick.major.size": 0, "xtick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
        "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE, "legend.frameon": False})


def short(name: str) -> str:
    s = SHORT.get(name, name)
    return s[0].upper() + s[1:]


class Page:
    """A figure laid out in inches from the top-left: header, body, footer."""

    def __init__(self, width: float, height: float, margin: float = 0.3):
        use()
        self.W, self.H, self.m = width, height, margin
        self.fig = plt.figure(figsize=(width, height))

    def y(self, inches_from_top: float) -> float:
        return 1 - inches_from_top / self.H

    def x(self, inches_from_left: float) -> float:
        return inches_from_left / self.W

    def text(self, left, top, s, **kw):
        kw.setdefault("ha", "left")
        kw.setdefault("va", "baseline")
        return self.fig.text(self.x(left), self.y(top), s, **kw)

    def header(self, title: str, *lines: list[str]) -> None:
        """Bold title, then one grey line per list of facts and keys (items spaced, no separator rule).
        Each extra line takes 0.18 in."""
        self.text(self.m, 0.30, title, fontsize=13, fontweight="bold")
        for i, facts in enumerate(lines):
            self.text(self.m, 0.54 + 0.18 * i, "    ".join(facts), fontsize=8.5, color=INK2)

    def footnote(self, lines: list[str]) -> None:
        for i, line in enumerate(reversed(lines)):
            self.fig.text(self.x(self.m), (0.18 + 0.15 * i) / self.H, line, ha="left", va="baseline",
                          fontsize=6.8, color=MUTED)

    def axes(self, left: float, top: float, width: float, height: float):
        return self.fig.add_axes([self.x(left), self.y(top + height), width / self.W, height / self.H])

    def save(self, path) -> None:
        self.fig.savefig(path, dpi=220)
        plt.close(self.fig)


def legend_row(page: Page, items, top: float, left: float | None = None) -> None:
    """A row of legend items: (kind, colour, label) with kind 'dot', 'ring' or 'tick'."""
    x = page.m if left is None else left
    for kind, colour, label in items:
        marker = {"dot": "o", "ring": "o", "tick": "|"}[kind]
        page.fig.add_artist(Line2D([page.x(x + 0.05)], [page.y(top)], marker=marker, ls="",
                                   ms=5.5 if kind != "tick" else 8, mfc=WHITE if kind == "ring" else colour,
                                   mec=colour, mew=1.1, transform=page.fig.transFigure))
        t = page.text(x + 0.16, top, label, va="center", fontsize=7.5)
        page.fig.canvas.draw()
        x += 0.16 + t.get_window_extent().width / page.fig.dpi + 0.28


def stepped_key(page: Page, left: float, top: float, bounds, colours, ticks, title: str, swatch: float = 0.42):
    page.text(left, top, title, fontsize=7.5, fontweight="bold")
    ax = page.axes(left, top + 0.09, swatch * len(colours), 0.10)
    for i, c in enumerate(colours):
        ax.add_patch(Rectangle((i + 0.03, 0), 0.94, 1, color=c, lw=0))
    ax.set_xlim(0, len(colours)), ax.set_ylim(0, 1)
    ax.set_xticks(range(len(bounds)))
    ax.set_xticklabels(ticks, fontsize=6.8, color=INK2)
    ax.tick_params(length=0, pad=2)
    ax.set_yticks([])
    ax.spines["bottom"].set_visible(False)


# ----------------------------------------------------------------------------
# Figure 1: ordinary test vs. batch check, one tissue
# ----------------------------------------------------------------------------
def lung_batch_check(res: pd.DataFrame, path, tissue: str = "Lung", void_above: float = 0.6) -> None:
    """Dumbbell per cell type: ring = AUROC on unseen mice (ordinary test), dot = cross-age AUROC
    (batch check; red passes 0.7, blue does not, grey = void by the shuffled control), tick =
    shuffled control. Right: the two numbers as columns."""
    d = res.sort_values("cross_age", ascending=False)
    n = len(d)
    row, lbl, plot_w, cols_w = 0.25, 1.75, 3.9, 1.3
    body_top = 1.1
    W = 0.3 + lbl + plot_w + cols_w + 0.3
    H = body_top + row * n + 0.5 + 0.45 + 0.5
    p = Page(W, H)
    p.header(f"Age classifier AUROC per cell type · {tissue}",
             [f"{n} cell types", "ring: tested on unseen mice", "dot: tested on an unseen age group",
              "0.5 = chance", "pass ≥ 0.7"])
    ax = p.axes(0.3 + lbl, body_top, plot_w, row * n)
    y = np.arange(n)[::-1]
    for v in (0.3, 0.4, 0.6, 0.8, 0.9, 1.0):
        ax.axvline(v, color=RULE, lw=0.5, zorder=0)
    ax.axvline(0.5, color=INK, lw=0.7, zorder=0)
    ax.axvline(0.7, color=INK, lw=0.7, ls=(0, (3, 2)), zorder=0)
    ax.text(0.5, n - 0.35, "chance", ha="center", va="bottom", fontsize=6.8, color=INK2)
    ax.text(0.7, n - 0.35, "pass", ha="center", va="bottom", fontsize=6.8, color=INK2)
    c1, c2 = 0.3 + lbl + plot_w + 0.7, 0.3 + lbl + plot_w + 1.3
    for cx, head in [(c1, "unseen\nmice"), (c2, "unseen\nage group")]:
        p.text(cx, body_top - 0.06, head, ha="right", va="bottom", fontsize=6.8, fontweight="bold", color=INK2,
               linespacing=1.1)
    for yi, (ct, r) in zip(y, d.iterrows()):
        void = r["shuffled"] > void_above
        passes = r["cross_age"] >= 0.7 and not void
        c = MUTED if void else (RED if passes else BLUE)
        ax.plot([r["cross_age"], r["logreg"]], [yi, yi], color=NEUTRAL, lw=3.2, zorder=1, solid_capstyle="round")
        ax.plot([r["shuffled"]] * 2, [yi - 0.25, yi + 0.25], color=MUTED, lw=1.0, zorder=2)
        ax.scatter(r["logreg"], yi, s=30, facecolor=WHITE, edgecolor=INK2, lw=1.0, zorder=3)
        ax.scatter(r["cross_age"], yi, s=36, color=c, zorder=4, lw=0)
        ty = body_top + row * (n - 1 - yi) + row / 2
        p.text(0.3, ty, short(ct), va="center", fontsize=8)
        p.text(c1, ty, f"{r['logreg']:.2f}", ha="right", va="center", fontsize=7.5, color=INK2)
        p.text(c2, ty, "void" if void else f"{r['cross_age']:.2f}", ha="right", va="center", fontsize=7.5,
               color=c, fontweight="bold" if passes else "normal")
    ax.set_xlim(0.28, 1.0)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks([])
    ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("AUROC", labelpad=4)
    legend_row(p, [("ring", INK2, "Unseen mice (ordinary test)"), ("dot", RED, "Unseen age group: passes"),
                   ("dot", BLUE, "does not pass"), ("dot", MUTED, "void"), ("tick", MUTED, "Shuffled control")],
               top=body_top + row * n + 0.62)
    p.footnote(["Batch check: trained on 1, 18 and 21 months, tested on 3 and 30 months, then the reverse; "
                "each age group was processed on its own day.",
                "Shuffled control: the whole pipeline with age shuffled between mice (mean of 5); above 0.6 the row "
                "is void. Sex-chromosome genes left out. Notebook 05a."])
    p.save(path)


# ----------------------------------------------------------------------------
# Figure 2: the atlas, cross-age AUROC per tissue x cell type
# ----------------------------------------------------------------------------
ATLAS_BOUNDS = [0.0, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0001]
ATLAS_COLOURS = [BLUE_ARM[1], BLUE_ARM[0], NEUTRAL, RED_ARM[0], RED_ARM[1], RED_ARM[3]]
COMPARTMENTS = ["stromal", "endothelial", "epithelial", "progenitor", "immune", "other"]


def atlas_heatmap(atlas: pd.DataFrame, path, value: str = "cross_age_all_splits") -> None:
    """Rows: cell types grouped by compartment, each header with its pass count; columns: tissues.
    Cell colour = cross-age AUROC in steps, value printed (truncated so it sits in its colour band);
    void rows outlined."""
    a = atlas.copy()
    tissues = sorted(a["tissue"].unique())
    groups = [(c, sorted(a.loc[a["compartment"] == c, "cell_type"].unique(), key=str.lower))
              for c in COMPARTMENTS if (a["compartment"] == c).any()]
    rows = []
    for comp, cts in groups:
        sub = a[(a["compartment"] == comp) & ~a["verdict"].isin(["void", "not testable"])]
        rows.append((comp.capitalize(), None, f"{(sub['verdict'] == 'age signal').sum()} of {len(sub)} pass"))
        rows += [(short(ct), ct, None) for ct in cts]
    n, m = len(rows), len(tissues)
    cw, ch, lbl = 0.40, 0.175, 2.15
    body_top = 1.85
    W = 0.3 + lbl + cw * m + 0.3
    H = body_top + ch * n + 0.7
    p = Page(W, H)
    p.header(f"Cross-age AUROC per cell type · {m} tissues",
             [f"{len(a)} rows", "23 mice", "mean over every split", "red: passes at 0.7",
              "blue: inverted on new age groups"])
    stepped_key(p, 0.3, 0.86, ATLAS_BOUNDS, ATLAS_COLOURS, ["0", ".4", ".5", ".7", ".8", ".9", "1"],
                "Cross-age AUROC")
    ax = p.axes(0.3 + lbl, body_top, cw * m, ch * n)
    cmap, norm = ListedColormap(ATLAS_COLOURS), BoundaryNorm(ATLAS_BOUNDS, len(ATLAS_COLOURS))
    look = a.set_index(["cell_type", "tissue"])
    for i, (label, ct, count) in enumerate(rows):
        yy = n - 1 - i
        ty = body_top + ch * i + ch / 2
        if ct is None:
            p.text(0.3, ty, label, va="center", fontsize=8, fontweight="bold")
            p.text(0.3 + lbl - 0.08, ty, count, ha="right", va="center", fontsize=6.8, color=INK2)
            if i > 0:
                p.fig.add_artist(Line2D([p.x(0.3), p.x(0.3 + lbl + cw * m)], [p.y(ty - ch / 2)] * 2, color=INK, lw=0.5))
            continue
        p.text(0.3, ty, label, va="center", fontsize=7)
        for j, t in enumerate(tissues):
            if (ct, t) not in look.index or pd.isna(look.loc[(ct, t), value]):
                continue
            r = look.loc[(ct, t)]
            v = r[value]
            if r["verdict"] == "void":
                ax.add_patch(Rectangle((j + 0.07, yy - 0.40), 0.86, 0.80, fill=False, ec=MUTED, lw=0.6))
                ax.text(j + 0.5, yy, "void", ha="center", va="center", fontsize=5.6, color=MUTED)
                continue
            ax.add_patch(Rectangle((j + 0.07, yy - 0.40), 0.86, 0.80, color=cmap(norm(v)), lw=0))
            # rounded like every number in the README; three decimals only where rounding would
            # cross a colour boundary (0.6996 would print .70 on a grey "partial" cell)
            crosses = norm(round(v, 2)) != norm(v)
            label = f"{v:.3f}" if crosses else f"{v:.2f}"
            ax.text(j + 0.5, yy, label.lstrip("0"), ha="center", va="center", fontsize=5.8 if not crosses else 5.2,
                    color=WHITE if (v >= 0.9 or v < 0.4) else INK)
    ax.set_xlim(0, m), ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks([])
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(m) + 0.5)
    ax.set_xticklabels([TISSUE.get(t, t) for t in tissues], rotation=90, fontsize=7.5)
    ax.tick_params(axis="x", length=0, pad=3)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    p.footnote(["void: shuffled control above 0.6, the row is not scored.  Blank: cell type absent or not testable there.",
                "Three decimals where rounding would cross a colour band.  Sex-chromosome genes left out. Notebook 07."])
    p.save(path)


# ----------------------------------------------------------------------------
# Figure 3: measured mature miRNA against age, per organ
# ----------------------------------------------------------------------------
MIR_BOUNDS = [-0.8, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8]
MIR_COLOURS = [BLUE_ARM[1], BLUE_ARM[0], NEUTRAL, RED_ARM[0], RED_ARM[1], RED_ARM[3]]
MIR_GROUPS = [("miR-29 family", ["miR-29a", "miR-29b", "miR-29c"]),
              ("Inflammation", ["miR-34a", "miR-146a", "miR-21", "miR-155"]),
              ("Other", ["miR-22", "miR-17"])]


def mirna_heatmap(meas: pd.DataFrame, path) -> None:
    """Rows: organs sorted by mean miR-29 rho; columns: a priori miRNAs in three groups.
    Cell colour = Spearman rho with age in steps; bottom row = organs where it rises (t > 2)."""
    piv = meas.pivot(index="organ", columns="miRNA", values="rho")
    tval = meas.pivot(index="organ", columns="miRNA", values="t")
    piv = piv.loc[piv[MIR_GROUPS[0][1]].mean(axis=1).sort_values(ascending=False).index]
    cols, xs, x, starts = [], [], 0.0, []
    for g, ms in MIR_GROUPS:
        starts.append(x)
        for mname in ms:
            cols.append(mname), xs.append(x)
            x += 1
        x += 0.5
    units = x - 0.5
    n = len(piv)
    cw, ch, lbl = 0.52, 0.215, 1.35
    head = 0.62                                   # group title + miRNA names, above the cells
    body_top = 1.32
    W = 0.3 + lbl + cw * units + 0.3
    H = body_top + head + ch * n + 0.45 + 0.55
    p = Page(W, H)
    p.header(f"Mature miRNA vs age · {n} organs",
             ["10 ages, 1 to 27 months", "small RNA, Wagner et al. 2024", "red: rises with age"])
    stepped_key(p, 0.3, 0.86, MIR_BOUNDS, MIR_COLOURS, ["−.8", "−.4", "−.2", ".2", ".4", ".6", ".8"],
                "Spearman ρ with age")
    x0 = 0.3 + lbl
    for (g, ms), s0 in zip(MIR_GROUPS, starts):
        p.text(x0 + cw * (s0 + 0.05), body_top + 0.18, g, fontsize=8, fontweight="bold")
        p.fig.add_artist(Line2D([p.x(x0 + cw * (s0 + 0.05)), p.x(x0 + cw * (s0 + len(ms) - 0.05))],
                                [p.y(body_top + 0.27)] * 2, color=INK, lw=0.8))
    for mname, xx in zip(cols, xs):
        p.text(x0 + cw * (xx + 0.5), body_top + 0.47, mname, ha="center", fontsize=7, color=INK2)
    ax = p.axes(x0, body_top + head, cw * units, ch * n)
    cmap, norm = ListedColormap(MIR_COLOURS), BoundaryNorm(MIR_BOUNDS, len(MIR_COLOURS))
    for i, organ in enumerate(piv.index):
        yy = n - 1 - i
        p.text(0.3, body_top + head + ch * i + ch / 2, TISSUE.get(organ, organ.replace("_", " ")),
               va="center", fontsize=7.5)
        for mname, xx in zip(cols, xs):
            v = piv.loc[organ, mname]
            ax.add_patch(Rectangle((xx + 0.05, yy - 0.42), 0.9, 0.84,
                                   color=cmap(norm(np.clip(v, -0.7999, 0.7999))), lw=0))
            txt = ".00" if abs(v) < 0.005 else f"{v:.2f}".replace("0.", ".").replace("-", "−")
            ax.text(xx + 0.5, yy, txt, ha="center", va="center", fontsize=6,
                    color=WHITE if (v >= 0.6 or v < -0.4) else INK)
    yb = body_top + head + ch * n + 0.25
    p.fig.add_artist(Line2D([p.x(0.3), p.x(x0 + cw * units)], [p.y(yb - 0.14)] * 2, color=INK, lw=0.5))
    p.text(0.3, yb, "Organs where it rises (t > 2)", va="center", fontsize=7, fontweight="bold", color=INK2)
    for mname, xx in zip(cols, xs):
        up = int((tval[mname] > 2).sum())
        p.text(x0 + cw * (xx + 0.5), yb, f"{up}/{len(tval)}", ha="center", va="center", fontsize=7.5,
               fontweight="bold" if up >= n / 2 else "normal", color=RED if up >= n / 2 else INK2)
    ax.set_xlim(0, units), ax.set_ylim(-0.5, n - 0.5)
    ax.set_xticks([]), ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    p.footnote(["ρ: Spearman correlation of log2(reads per million + 1) with age over the mice of one organ;  "
                "t: age slope adjusted for sex and plate.",
                "Same mice and RNA as the bulk time course. Organs sorted by mean miR-29 ρ. Notebook 13."])
    p.save(path)


# ----------------------------------------------------------------------------
# Figure 4: the replicating genes over the bulk time course
# ----------------------------------------------------------------------------
BULK_PANELS = [("Col1a1", "Collagen I"), ("H2-K1", "MHC class I"), ("C4b", "Complement C4b"), ("Cdkn1a", "p21")]
ORGAN_LABEL = {"Gonadal_Fat": "Gonadal fat", "Brown_Fat": "Brown fat", "Mesenteric_Fat": "Mesenteric fat",
               "Subcutaneous_Fat": "Subcutaneous fat", "Limb_Muscle": "Muscle", "Small_Intestine": "Small intestine",
               "White_Blood_Cells": "White blood cells"}


def bulk_trajectories(d: pd.DataFrame, tvals: pd.DataFrame, path, tvals_adult: pd.DataFrame | None = None,
                      panels=BULK_PANELS) -> None:
    """Small multiples, one per gene: each organ's mean log CPM per age, centred on that organ's own
    mean (grey), and the median organ (bold; red if at least a third of organs rise and more rise than
    fall, blue for the reverse, grey otherwise).
    Counts from `tvals` (age t, plate and sex adjusted): organs with t > 2 / t < -2.
    `d`: one row per sample with Organ, Age and one column per gene (log CPM). `tvals_adult`: the same t
    from 3 months on, so a change during growth (1-3 months) is not read as aging; the growth window is shaded."""
    k = len(panels)
    pw, ph, gap = 1.55, 1.75, 0.32
    left, top = 0.3 + 0.35, 1.47
    W = left + k * pw + (k - 1) * gap + 0.3
    H = top + ph + 0.55 + 0.6
    p = Page(W, H)
    p.header(f"Expression across the lifespan · {d['Organ'].nunique()} organs",
             [f"bulk RNA, {d['Age'].nunique()} ages, 1 to 27 months", "grey: one organ", "bold: median organ",
              "shaded: growth, 1 to 3 months"])
    ages = sorted(d["Age"].unique())
    curves = {}
    for gene, _ in panels:
        m = d.groupby(["Organ", "Age"])[gene].mean().unstack()
        curves[gene] = m.sub(m.mean(axis=1), axis=0)
    lim = max(np.nanpercentile(np.abs(c.to_numpy()), 98) for c in curves.values())
    lim = float(np.ceil(lim * 4) / 4)
    for i, (gene, name) in enumerate(panels):
        x0 = left + i * (pw + gap)
        up, down = int((tvals[gene] > 2).sum()), int((tvals[gene] < -2).sum())
        n = len(tvals)
        colour = RED if (up >= n / 3 and up > down) else (BLUE if (down >= n / 3 and down > up) else INK2)
        verdict = {RED: f"up in {up} of {n}", BLUE: f"down in {down} of {n}"}.get(colour, f"up {up}, down {down} of {n}")
        p.text(x0, 0.85, name, fontsize=9, fontweight="bold")
        p.text(x0, 1.03, f"{gene}  ·  {verdict} organs", fontsize=7.2, color=colour,
               fontweight="bold" if colour != INK2 else "normal")
        if tvals_adult is not None:
            ua, da = int((tvals_adult[gene] > 2).sum()), int((tvals_adult[gene] < -2).sum())
            p.text(x0, 1.20, f"from 3 months: up {ua}, down {da}", fontsize=7.2, color=INK2)
        ax = p.axes(x0, top, pw, ph)
        for v in np.arange(-lim, lim + 1e-9, 0.5):
            ax.axhline(v, color=RULE, lw=0.5, zorder=0)
        ax.axhline(0, color=INK, lw=0.6, zorder=0)
        ax.axvspan(0, 3, color=NEUTRAL, lw=0, zorder=0)
        if i == 0:
            ax.text(0.3, lim * 1.03, "growth", ha="left", va="bottom", fontsize=6.5, color=INK2)
        c = curves[gene]
        for organ, row in c.iterrows():
            ax.plot(row.index, row.to_numpy(), color=MUTED, lw=0.6, alpha=0.55, zorder=1)
        med = c.median(axis=0)
        ax.plot(med.index, med.to_numpy(), color=colour, lw=2.2, zorder=3, solid_capstyle="round")
        ax.set_xlim(0, 28), ax.set_ylim(-lim, lim)
        ax.set_xticks([1, 9, 18, 27])
        ax.tick_params(axis="y", length=0, labelsize=7)
        if i == 0:
            ax.spines["left"].set_visible(True)
            ax.tick_params(axis="y", length=3)
            ax.set_ylabel("change from the organ's mean\nlog(CPM + 1)", fontsize=7.2, labelpad=3)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("age, months", fontsize=7.2, labelpad=2)
    p.footnote(["Grey: one organ, mean per age, centred on its own mean. Bold: the median organ. Counts: organs with age "
                "t > 2 or < −2, plate and sex as covariates.",
                "Shaded: 1 to 3 months, still growth. Most of the collagen fall happens there. Bulk time course "
                "(Schaum et al. 2020); notebook 14."])
    p.save(path)


# ----------------------------------------------------------------------------
# Figure 5: what leaving out the sex-chromosome genes did
# ----------------------------------------------------------------------------
def sex_gene_slope(atlas: pd.DataFrame, path) -> None:
    """Slope chart, one line per testable row: cross-age AUROC with the six sex-chromosome genes in
    (left) and out (right, the default). Rows that cross 0.7 upward are red and labelled."""
    a = atlas[~atlas["verdict"].isin(["void", "not testable"])].copy()
    a["with"], a["without"] = a["with_sex_genes_cross_age"], a["cross_age_all_splits"]
    rise = a[(a["with"] < 0.7) & (a["without"] >= 0.7)].sort_values("without", ascending=False)
    fall = a[(a["with"] >= 0.7) & (a["without"] < 0.7)]
    pw, ph = 2.0, 3.6
    left, top = 0.95, 1.0
    W = left + pw + 2.9
    H = top + ph + 0.5 + 0.6
    p = Page(W, H)
    p.header("Cross-age AUROC with and without sex genes",
             [f"{len(a)} scored rows", f"{len(rise)} cross 0.7 once the genes are left out", f"{len(fall)} fall"],
             ["left: six sex-chromosome genes in (first draft)", "right: left out (the default)"])
    ax = p.axes(left, top, pw, ph)
    lo = float(np.floor(min(a["with"].min(), a["without"].min()) * 10) / 10)
    for v in np.arange(lo, 1.0001, 0.1):
        ax.axhline(v, color=RULE, lw=0.5, zorder=0)
    ax.axhline(0.7, color=INK, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for _, r in a.iterrows():
        if r.name in rise.index:
            continue
        ax.plot([0, 1], [r["with"], r["without"]], color=MUTED, lw=0.6, alpha=0.6, zorder=2)
    for _, r in rise.iterrows():
        ax.plot([0, 1], [r["with"], r["without"]], color=RED, lw=1.5, zorder=3)
        ax.scatter([0, 1], [r["with"], r["without"]], s=12, color=RED, zorder=4, lw=0)
    ax.set_xlim(-0.08, 1.08), ax.set_ylim(lo, 1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["sex genes in", "left out"], fontsize=7.5)
    ax.spines["left"].set_visible(True)
    ax.tick_params(axis="y", length=3, labelsize=7)
    ax.set_ylabel("cross-age AUROC", fontsize=7.5)
    ax.text(0.5, 0.7, "pass", ha="center", va="center", fontsize=6.8, color=INK2, zorder=5,
            bbox=dict(fc=WHITE, ec="none", pad=1.2))
    # labels for the rising rows, spread so they do not overlap
    ys = rise["without"].to_numpy().astype(float)
    step = 0.024 * (1.0 - lo) / 0.7
    placed = []
    for y in ys:
        y = min(y, placed[-1] - step) if placed else y
        placed.append(y)
    for (_, r), yl in zip(rise.iterrows(), placed):
        ax.plot([1.0, 1.1], [r["without"], yl], color=RED, lw=0.5, clip_on=False)
        fmt = lambda v: f"{v:.3f}" if (round(v, 2) >= 0.7) != (v >= 0.7) else f"{v:.2f}"
        ax.text(1.12, yl, f"{TISSUE.get(r['tissue'], r['tissue'])} · {short(r['cell_type'])}"
                          f"   {fmt(r['with'])} → {fmt(r['without'])}",
                va="center", fontsize=6.8, color=INK, clip_on=False)
    p.footnote(["Sex is tied to age group: 21-month mice all female, 24 and 30 all male. With Xist in, a split",
                "that pairs male young with female old in training, and the reverse in test, inverts the model.",
                "Void rows left out. Notebook 07, sections 11 and 12."])
    p.save(path)
