# Single-Cell Aging

Find the cell type in an aging tissue whose gene expression really tracks age,
and separate that signal from the batch it was collected in.

Give the pipeline an annotated single-cell RNA-seq atlas with donor age and
donor identity. It trains one age model per cell type, holding out whole
animals, then puts every surviving model through a second test on age groups it
has never seen. A cell type counts as having an age signal only if it passes
both. The same pipeline then estimates age in month, separates a change of
aging from a change of development, and names the gene behind each model.

Nothing in the pipeline is tissue specific. It was written for the mouse lung
of Tabula Muris Senis, and that case is carried through
[section 10](#10-result-mouse-lung) as a worked example, but no step assumes
it and none of that data is part of this repository.


## Contents

1. [Overview](#1-overview)
2. [Biological Background](#2-biological-background)
3. [Install](#3-install)
4. [Usage](#4-usage)
5. [Input](#5-input)
6. [Output](#6-output)
7. [Data Process](#7-data-process)
8. [Evaluation Design](#8-evaluation-design)
9. [File Map](#9-file-map)
10. [Result: Mouse Lung](#10-result-mouse-lung)
11. [Limitation](#11-limitation)
12. [Reference](#12-reference)


## 1. Overview

**Input.** One `.h5ad` atlas with raw count, a cell type label, and a donor age
and donor id per cell.

**Output.** A table of AUROC per cell type under every control, a per-cell-type
gene list with an artifact flag, an age clock in month, and a figure per step.

**What it does, in one pass:**

```
read the atlas
  -> rebuild raw count, QC, doublet removal, normalize, 2000 variable gene
  -> per mouse cell-type fraction vs age          (composition, not expression)
  -> per cell type, hold out whole mice:
       logistic regression young vs old  -> AUROC on unseen mice
       gate: shuffled control near 0.5
       gate: cross-age AUROC >= 0.7      (the batch check)
  -> for a model that passes, read its gene
       filter: ambient, technical, sex, other cell type's marker
       flag: circadian, dissociation
       split: change of development (1m->3m) vs change of aging (3m->old)
  -> ridge regression on age in month, leave one mouse out
```

**The headline finding on mouse lung.** Lung **fibroblast** and **adventitial
cell** carry an age signal that survives every control, including a model tested
on age groups it never trained on (AUROC 0.97 and 0.95). In the old fibroblast
the gene that rise include ***Cdkn1a*** (p21, a core senescence gene) and the
complement gene *C3* and *Cfb*, while collagen fall. Immune cell type score
well on unseen mice (AUROC up to 0.91) and then fail the batch check, which is
the finding that matters methodologically: on an atlas where each age group was
processed separately, ordinary cross-validation cannot tell aging from batch.


## 2. Biological Background

### Why One Cell Type Is Not Enough

A tissue does not age as one thing. Bulk RNA-seq grinds a whole organ together
and reports the average over million of mixed cell, so when an inflammation
gene rises in old tissue the measurement cannot say whether every cell turned it
up a little or a few immune cell moved in. Those two have opposite biological
meaning and the same bulk readout.

**Single-cell RNA-seq** tag each mRNA molecule with a barcode naming the cell it
came from, giving a cell × gene count matrix. Each cell type can then be asked
about on its own. That is what makes the question here answerable: not "does
this tissue age" but "which cell type in it ages, and how".

Two consequence of the measurement shape the whole pipeline. The matrix is
**sparse**, over 90% zero, because each cell release only a fraction of its
mRNA, so a zero usually means "not caught" rather than "not expressed". And the
cell of one animal are **not independent sample** of anything: they share that
animal's genotype, its health, the day its tissue was cut and the reagent lot it
was processed with.

### Why The Hallmark Matter Here

Aging is organized into **hallmark**, change that appear with age, accelerate
aging when worsened and slow it when fixed (López-Otín et al. 2023). This
project touches two directly.

**Cellular senescence.** A badly stressed cell can stop dividing permanently
instead of dying. It is marked by the cell-cycle blocker ***Cdkn1a*** (p21) and
***Cdkn2a*** (p16), and it secretes an inflammatory, tissue-remodeling mix
called the SASP. Senescent cell accumulate with age, and clearing them extend
healthy lifespan in mouse (Baker et al. 2011, 2016), which is why senolytic drug
are in trial. No single gene marks every senescent cell, so a gene set such as
SenMayo is scored instead (Saul et al. 2022).

**Inflammaging.** Old tissue carries chronic low-grade inflammation, partly from
the SASP and partly from immune cell accumulating in the organ. It shows up in
single-cell data as a shift in which cell are present and as immune and
interferon gene rising inside many cell type.

A per-cell-type age model is an active research direction, not a solved problem:
Buckley et al. 2023 built one per cell type in the mouse brain. This project
asks the same question of the lung, and adds a batch control that the atlas
design makes necessary.

### What Each Metric Measures

| Metric | The question it answers | Why it is here |
| --- | --- | --- |
| **AUROC on unseen mice** | Can a model rank an old cell above a young one in an animal it never saw? | The honest version of accuracy. It needs no cutoff and does not care that old cell outnumber young one |
| **Cross-age AUROC** | Does that model still work on an age group it never trained on? | Each age group was processed in its own batch, so this is the only score that a batch shortcut cannot help |
| **Shuffled control** | What does the same pipeline score when the age label is swapped between animal at random? | The floor. It says how much separation is reachable from between-animal difference alone, so a real score has to clear it |
| **Per-mouse AUROC** | Can the model rank whole animal by age? | The sample size here is the animal, not the cell, so this is the number with real degree of freedom behind it |
| **Clock error in month** | How far off is the predicted age, against a naive guess? | Turns a two-class result into a quantity, and the naive baseline keeps it honest |
| **Development vs aging change** | Did the gene change from 1m to 3m, or from 3m onward? | Growth is not aging. A gene that only moves while the animal matures is maturation, not an aging marker |

### Why They Are Kept Apart

The AUROC on unseen mice and the cross-age AUROC are deliberately **not** fused
into one rank. They disagree, and the disagreement is the entire result: a model
can separate animal perfectly and still be reading the collection batch. Collapsing
them into an average would hide exactly the case worth reporting.

A cell type is therefore reported as a row of score, with a verdict, and the
reason for the verdict is visible in the row rather than buried in a weight.


## 3. Install

```bash
conda env create -f environment.yml
conda activate scaging
```

Required: `scanpy`, `anndata`, `scikit-learn`, `pandas`, `numpy`, `seaborn`,
`matplotlib`, `igraph`, `leidenalg`. Python 3.14.

Optional:

- `lightgbm` for the non-linear comparison model. Without it, skip that one cell
  in `05a`.
- `pydeseq2` for the pseudobulk differential expression of phase 7.
- `nbdime` for a readable notebook diff.

Two git setting are worth enabling once per clone:

```bash
nbdime config-git --enable                                          # cell by cell notebook diff
git config filter.strip-output.clean "python src/strip_output.py"   # commit code only
git config filter.strip-output.required true
```

The second one matter. An output is regenerated by rerunning the notebook, and
in git it only hides the code change: a plot is a megabyte of base64 that differ
on every run, and even a printed table churns. So the committed notebook is code
and markdown only, which keeps a diff readable, while the copy on disk keeps
every output. `required true` makes a missing filter fail the commit instead of
quietly committing the outputs.


## 4. Usage

### The Run

There is no single entry point. Each phase is one notebook in `notebook/`, run
in order, and each reads the file the previous one wrote:

```bash
conda activate scaging
jupyter lab            # then run 02, 03, 04, 05a, 05b, 05c in that order
```

The order is not a style preference. `03` writes the `.h5ad` that every later
notebook load, and `05b` and `05c` both read `figure/lung_age_auroc.csv`, so
`05a` has to run before them on a fresh clone.

### Every Setting

Every shared default lives in `src/aging.py`, so a change applies to all three
phase-5 notebook at once:

| Setting | Default | Meaning |
| --- | --- | --- |
| `AGE_MONTHS` | 1, 3, 18, 21, 24, 30 | Age label to month, for the clock |
| `YOUNG` / `OLD` | 1m, 3m / 18m, 21m, 24m, 30m | The two class; an age absent from a tissue is ignored |
| `ADULT` | 3m and up | Age used by the clock, so growth is not modelled as aging |
| `max_cells` | 300 | Cap per mouse, so one animal cannot dominate a cell type |
| `AMBIENT` | *S100a8/9*, *Lyz2*, hemoglobin, *Scgb1a1* … | Gene treated as free-floating RNA rather than a cell's own |
| `TECHNICAL_PATTERN` | `Rp[sl]`, `Malat1`, `Xist`, `Ddx3y` … | Ribosomal, rRNA, and sex-linked gene |
| `DISSOCIATION` | *Fos*, *Jun*, *Hspa1a* … | Gene induced by tissue digestion, flagged not dropped |
| `CIRCADIAN` | *Dbp*, *Nr1d1*, *Per1* … | Time-of-day gene, flagged as a batch signature |
| `EXCLUDE` | 30-M-2 in myeloid | Mouse dropped from named cell type only, with a reason |
| `MYELOID` | monocyte, macrophage | The cell type that `EXCLUDE` applies to |
| `C` | 0.1 | L2 strength of the logistic regression |
| `alpha` | 10, 100, 1000 | Ridge strength of the clock; all three are reported |
| `min_cells` / `min_mice` | 50 / 2 | A cell type below this in either class is not modelled |

### Calling It Directly

Every function take the `AnnData` as its first argument, so the same code runs
on any tissue from a plain script:

```python
import scanpy as sc
from src import aging as ag

adata = sc.read_h5ad("data/lung_processed.h5ad")

for ct in ag.usable_cell_types(adata):
    kw = ag.main_kw(ct)                       # applies EXCLUDE for this cell type
    print(ct,
          round(ag.age_auc(adata, ct, **kw), 2),             # unseen mice
          round(ag.cross_age_auc(adata, ct, **kw), 2),       # unseen age group
          round(ag.age_auc(adata, ct, shuffle_seed=0, **kw), 2))   # shuffled
```


## 5. Input

### The Atlas

One `.h5ad` file. Four thing have to be there, and the pipeline checks or
rebuilds the rest:

| Needed | Where | Note |
| --- | --- | --- |
| Raw integer count | `.X` or `.raw` | Scrublet and pseudobulk both need count, not a log value |
| Cell type label | `.obs["cell_ontology_class"]` | Taken as given; this project does not re-annotate |
| Age | `.obs["age"]` | A label such as `18m`, mapped to month by `AGE_MONTHS` |
| Donor id | `.obs["mouse.id"]` | **The unit of analysis.** Without it nothing here is valid |

**The donor id is the one column that cannot be missing.** Every score in
section 6 is computed by holding out a whole animal. With age but no donor id,
a model can be trained and will report a high number, and that number means
nothing, because the cell of one animal in both train and test let the model
recognize the animal instead of the age. This is the single most common way a
single-cell age model is silently wrong.

**A log value in `.X` is handled.** The official Tabula Muris Senis file ship
with `.X` already log-normalized and the count kept in `.raw`, so `03` rebuilds
from `.raw` and keeps a copy in `layers["counts"]`. Normalizing an already
normalized matrix twice is silent and ruins every downstream number, so the
notebook assert the matrix is integer before continuing.

### What A Tissue Has To Be Checked For First

Phase 2 to 4 exist to find the thing that would otherwise be read as biology.
For a new tissue they are not optional:

- **An outlier donor.** One animal supplying most of a cell type will drive that
  cell type's result on its own. Cross-tabulate cell type against donor before
  modelling anything.
- **Sex against age.** If the young animal are male and the old female, a model
  can score well on *Xist* alone. The males-only rerun exists for this.
- **A pooled sample.** An id like `3-M-5/6` is two animal in one sample. It is
  one unit, and it cannot be split.
- **Composition versus expression.** A cell type can change in number without
  any cell changing inside, and tissue dissociation affect number. Phase 4
  answers the count question separately so phase 5 can be about expression only.
- **Which age exist.** Lung has no 24m group. The cross-age split has to be
  written for the age a tissue actually has.

### The Current Tissue

Tabula Muris Senis droplet Lung, from
[figshare](https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728):
`tabula-muris-senis-droplet-processed-official-annotations-Lung.h5ad`, placed in
`data/`. **Not committed**: it is large, not ours to redistribute, and it
changes upstream.

| | |
| --- | --- |
| Cell | 24,540, and 24,490 after doublet removal |
| Gene | 20,138 |
| Mouse | 16: six young (1m, 3m) and ten old (18m, 21m, 30m) |
| Age | 1m, 3m, 18m, 21m, 30m. **No 24m in Lung** |
| Annotation | Cell type, age, sex and mouse id, provided by the consortium |

The release carries no mitochondrial gene, so the usual mitochondrial QC cannot
be run, and the cell were already filtered by the consortium: a
`n_genes_by_counts > 300` cut remove exactly one cell. A gene name is the older
symbol of the submission, so hemoglobin beta is `Beta-s` rather than `Hbb-bs`.


## 6. Output

### The Table

`figure/lung_age_auroc.csv`, one row per cell type. Every column is an AUROC
computed with whole mice held out, so the column differ only in what was changed
about the input:

| Column | Meaning |
| --- | --- |
| `logreg` | The main number: logistic regression, unseen mice |
| `cross_age` | **The batch check.** Trained on some age group, tested on others |
| `shuffled` | The same pipeline with age randomized between animal |
| `males_only` | Rerun on male animal alone, so sex cannot stand in for age |
| `no_ambient_genes` | Rerun with the ambient gene list dropped |
| `with_excluded_mice` | Rerun keeping the excluded donor, to show its effect |
| `lgbm` | LightGBM instead of logistic regression |
| `passes_batch_check` | `cross_age >= 0.7` |

`figure/lung_age_report.csv` describes the main model instead of varying it:

| Column | Meaning |
| --- | --- |
| `auroc_cell` | AUROC over cell, same as `logreg` above |
| `auroc_mouse` | AUROC after averaging each animal's cell score, then ranking animal |
| `bal_acc_test` | Balanced accuracy on held-out cell at a 0.5 cutoff |
| `acc_train` | Training accuracy, reported to show the size of the overfit |
| `n_cells`, `n_mice` | What the row is actually based on |

The gene table (`lung_top_genes_up_in_old.csv`, `..._up_in_young.csv`,
`lung_aging_only_genes.csv`) give the strongest weight per cell type after
filtering, with a `flag` naming a circadian or dissociation gene. The clock
table (`lung_clock_summary.csv`, `lung_clock_held_out_age.csv`) give error in
month per cell type at every regularization strength tried.

Everything in `figure/` is **regenerated, not committed**, including the CSV. A
committed table would duplicate the notebook that wrote it and go stale the
moment the pipeline runs on a second tissue; the number worth keeping are quoted
in [section 10](#10-result-mouse-lung). To display one plot in this README,
un-ignore that single file in `.gitignore`.

### Reading The Number

| Column | Range | Lowest means | Highest means | Want |
| --- | --- | --- | --- | --- |
| `logreg` | 0 to 1 | **0**: every old cell ranked below every young one | **1**: perfect ranking in an animal never seen | **high**, and above `shuffled` |
| `cross_age` | 0 to 1 | **0**: what was learned is exactly inverted on a new age group | **1**: transfers perfectly | **≥ 0.7** |
| `shuffled` | 0 to 1 | **0 or 1**: both mean the pipeline can hit an extreme on noise | **0.5**: no signal left after randomizing | **near 0.5** |
| `acc_train` | 0 to 1 | — | **1**: the model memorized its training animal | expect ~1, and read the gap |
| `mae_months` | ≥ 0 | **0**: exact age | large: no better than guessing | **below `naive_mae`** |

Three caveat about the extreme.

**A `cross_age` far below 0.5 is not a weak result, it is an inverted one.** At
0.15 the model is reliably wrong on a new age group, which can only happen if
what it learned was a property of the batch that flip sign between batch. Read
it as a positive detection of a batch effect, not as a near-miss.

**`acc_train` near 1.00 everywhere is expected and is not a finding.** Any model
with 2,000 gene and a dozen animal can memorize its training set. The number is
printed so the gap against the held-out score is visible, and that gap is the
overfit.

**A per-mouse AUROC of 1.00 is weaker evidence than it looks.** With this few
animal, and age confounded with batch, perfectly ranking the animal is exactly
what a batch effect produces. It is reported next to `cross_age` for that
reason, never on its own.

### The Figure

| Figure | Shows |
| --- | --- |
| `lung_umap*.png` | The cell type, age and donor map, plus one donor at a time in colour against grey |
| `lung_qc_by_*.png` | Gene and count per cell, by age and by donor |
| `lung_composition_by_mouse.png` | One dot per animal: each cell type's fraction against age, coloured by sex |
| `lung_monocyte_state.png`, `lung_state_*.png` | Marker group per donor, for the outlier investigation |
| `lung_age_auroc.png` | AUROC against cross-age AUROC per cell type, with the shuffled control |
| `lung_dev_vs_aging_*.png` | Each gene's change of development against its change of aging |
| `lung_clock_*.png` | Predicted against true age in month, one point per cell and a median per animal |

One donor per panel, plotted in colour against everything else in grey, is worth
the extra figure: two similar default colour (30-M-2 and 30-M-5) are impossible
to tell apart on a sixteen-category map, which is how an outlier donor hides.


## 7. Data Process

Notation. A tissue has $N$ cell from $D$ donor. $x_{i,g}$ is the expression of
gene $g$ in cell $i$ after normalization, $d(i)$ the donor of cell $i$, and
$a(d)$ the age in month of donor $d$. A cell type $c$ gives the cell subset
$I_c$.

### 7.1 Step 1: Rebuild, QC And Normalize

Count are taken from `.raw` and asserted integer, because `.X` in the released
file is already log-normalized and normalizing twice is silent. A copy is kept
in `layers["counts"]` for the pseudobulk of phase 7, which needs count rather
than a log value.

Doublet are called by Scrublet **per donor** (`batch_key="mouse.id"`), on count,
before normalization. Order matter: Scrublet normalizes internally, and running
it on an already normalized matrix changed the call from 49 cell to 127 on this
tissue.

Normalization is count per 10,000 then $\log(1+x)$:

$$x_{i,g} = \log\left(1 + 10^4 \cdot \frac{n_{i,g}}{\sum_{g'} n_{i,g'}}\right),$$

which put every cell on a comparable total, so a cell that happened to be
sequenced deeper does not read as a cell expressing more. Then 2,000 highly
variable gene, PCA, a neighbour graph, Leiden cluster and UMAP. The cluster are
used as a check on the provided label, never to replace them.

### 7.2 Step 2: Composition, Before Any Expression

For donor $d$ and cell type $c$ the fraction is

$$f_{d,c} = \frac{|\lbrace i : d(i) = d,\ i \in I_c \rbrace|}{|\lbrace i : d(i) = d \rbrace|},$$

computed with a cross-tabulation so that a donor with none of a cell type keep a
0 rather than vanishing from the average. The trend against age is Spearman's
$\rho$ over donor, with Benjamini–Hochberg correction across the cell type
tested, because the count is one per animal and there are a dozen animal.

This is kept separate from everything downstream on purpose. A cell type can
change in number with no cell changing inside it, dissociation affect the
number, and a single model over all cell would learn that composition instead of
aging.

### 7.3 Step 3: Build One Matrix Per Cell Type

For cell type $c$, cell are shuffled at a fixed seed and capped at
$n_{\max} = 300$ per donor:

$$I_c^{\ast} = \bigcup_{d} \operatorname{head}_{n_{\max}}\left(\operatorname{shuffle}(\lbrace i \in I_c : d(i) = d \rbrace)\right).$$

The cap is not cosmetic. One donor on this tissue contributes a third of all
cell, and without a cap its cell alone would set the loss. The label is

$$y_i = \mathbb{1}\left[a(d(i)) \in \text{OLD}\right],$$

and a donor named in `EXCLUDE` for this cell type is dropped here, with the
reason recorded next to it in `src/aging.py` rather than in a notebook.

A cell type is modelled only with at least 50 cell and 2 donor **in both
class**, so a cell type with 400 old cell from one animal and nothing young is
never given a score.

### 7.4 Step 4: The Two Gate

**Gate 1, the shuffled control.** Age is permuted **between donor**, never
between cell, then the whole pipeline is rerun unchanged:

$$\tilde{y}_i = \mathbb{1}\left[\pi(d(i)) \in \text{OLD}\right], \qquad \pi \text{ a permutation of the donor}.$$

The reported `shuffled` column is the mean over five such permutation.

The permutation keeps every cell with its own donor and every donor on one
single label; it changes only which donor is called old. Real age is therefore
the one thing removed, so any score above 0.5 comes from between-donor
difference — sequencing depth, dissection, health — happening to line up with
the random labelling. That makes `shuffled` a **floor to clear, not a p-value**:
fibroblast at 0.98 against 0.60 has a wide margin, while neutrophil at 0.97
against 0.70 means the pipeline can manufacture 0.70 from noise on that cell
type, so the 0.97 is not a clean measurement and the row is voided.

Permuting individual cell instead would leave each donor holding a mix of label,
destroying the structure being controlled for and returning an easy 0.5 that
tests nothing. The control also says nothing about batch: age and batch are the
same grouping here, so relabelling donor breaks both at once. That is Gate 2's
job, and a cell type can clear the shuffled control and still be entirely batch
driven.

**Gate 2, the cross-age check.** Train on one set of age group and test on a
disjoint set:

$$\text{train } \lbrace 1\text{m}, 18\text{m}, 21\text{m} \rbrace \to \text{test } \lbrace 3\text{m}, 30\text{m} \rbrace, \qquad \text{and the reverse},$$

each split carrying both class, and the two AUROC averaged. The split is
written so each side holds a young and an old group; since each age group is its
own processing batch, a model leaning on batch has nothing to lean on here.

A cell type passes when $\text{cross\_age} \ge 0.7$.

### 7.5 Step 5: Score

Cross-validation is `StratifiedGroupKFold` over donor, up to five fold,
stratified so each fold carries both class. Every cell is scored by a model
that never saw its donor:

$$\widehat{p}_i = \widehat{p}^{\,(-d(i))}(x_i).$$

**AUROC over cell** is the probability that a random old cell outranks a random
young one, which needs no cutoff and is unmoved by the class being unbalanced.
**AUROC over donor** average each donor's score first,
$\bar{p}_d = \operatorname{mean}_{d(i)=d} \widehat{p}_i$, then rank donor.
**Balanced accuracy** is the mean of the two per-class rate at a 0.5 cutoff, and
**training accuracy** is kept to expose the overfit.

The main model is L2 logistic regression at $C = 0.1$, which keep the weight
small so a few hundred cell from a dozen animal cannot support 2,000 free
parameter. LightGBM is run as a comparison: a tree ensemble can express a
non-linear rule, but it needs more independent sample than a dozen animal to
find one.

### 7.6 Step 6: Filter The Gene, Then Read Them

Only a model with AUROC ≥ 0.8 **and** a passing batch check has its gene read at
all. The filter is cumulative:

| Level | Drop |
| --- | --- |
| `hvg` | Nothing: the 2,000 highly variable gene |
| `clean` | `AMBIENT`: free-floating RNA and blood gene |
| `bio` | Ribosomal, rRNA, *Malat1*, *Xist*, Y-linked, hemoglobin |
| `own` | Another cell type's marker gene |

The `own` level is the interesting one. Using the mean log expression
$\mu_{c,g}$ of gene $g$ in cell type $c$, gene $g$ is foreign to $c$ when

$$\max_{c' \neq c} \mu_{c',g} > 1 \qquad \text{and} \qquad \mu_{c,g} < 0.1 \cdot \max_{c' \neq c} \mu_{c',g},$$

that is, clearly expressed somewhere else and nearly absent here. *Scgb1a1* from
club cell showing up in a T cell is the case it was written for. **It is not
sufficient:** the comparison is on a log scale, so a very abundant ambient gene
clear the 10% bar and pass. *Sftpc* and *Mgp* still get through, and are
currently handled by name. A linear-scale ratio catch them and also drop
*Cdkn1a* in fibroblast and *Cd74* in macrophage, which is the open problem in
[section 11](#11-limitation).

A circadian or dissociation gene is **flagged, not dropped**, because the flag
is the information: *Dbp* and *Nr1d1* appearing as "young" gene say the young
and old animal were collected at a different time of day.

### 7.7 Step 7: Separate Development From Aging

Expression is averaged per donor, then per age group, and two change are taken:

$$\Delta_{\text{dev}} = \mu_{3\text{m}} - \mu_{1\text{m}}, \qquad
\Delta_{\text{age}} = \operatorname{mean}_{a \in \text{OLD}} \mu_{a} - \mu_{3\text{m}}.$$

The pair label a gene: **aging only** when
$|\Delta_{\text{age}}| > 0.5$ and $|\Delta_{\text{dev}}| < 0.2$, **development
only** when reversed, **lifelong** when both move the same way, **reverses**
when they oppose. A 1m animal is still growing, so a gene that moves only from
1m to 3m is maturation, and calling it an aging marker would be wrong.

### 7.8 Step 8: The Clock

Ridge regression on age in month, **adult only** (3m and up, so growth is not
modelled as aging), leave one donor out:

$$\widehat{a}_i = \widehat{a}^{\,(-d(i))}(x_i), \qquad
\widehat{a}_d = \operatorname{median}_{d(i) = d} \widehat{a}_i,$$

and the error is per donor, against a naive baseline that always guesses the
mean age of the other donor:

$$\text{MAE} = \frac{1}{D}\sum_d \left| \widehat{a}_d - a(d) \right|, \qquad
\text{MAE}_{\text{naive}} = \frac{1}{D}\sum_d \left| a(d) - \operatorname{mean}_{d' \neq d} a(d') \right|.$$

The baseline is the point. A clock on a few animal spread over 3m to 30m can
look accurate by predicting the middle of the range for everything, and only the
comparison exposes it. Every regularization strength tried is reported, so a
strength is not chosen after seeing the test score. A whole age group can also
be held out, which ask whether the clock orders age it never saw.


## 8. Evaluation Design

Two unit are kept deliberately separate, and confusing them is the easiest way
to report a result that does not exist.

The **cell** is the unit of measurement. There are tens of thousand, they make
every figure look convincing, and they are **not independent**: the cell of one
animal share its genotype, its health and its processing day.

The **donor** is the unit of inference. There are sixteen. Every claim in
section 10 is bounded by that number, not by the cell count.

From that, three rule hold everywhere in the code:

| Rule | Why |
| --- | --- |
| A whole donor goes into train or test, never split | Otherwise the model can identify the animal and inherit its label |
| Cell are capped per donor | Otherwise the largest donor sets the loss on its own |
| A donor is averaged before a group is compared | Otherwise a p-value counts cell as sample and comes out tiny and meaningless |

And one rule follows from this atlas specifically, which is the part most likely
to be silently wrong:

**A high AUROC on unseen donor is necessary and not sufficient.** Each age group
in Tabula Muris Senis was collected and processed at a different time, so donor
identity and batch and age all move together. A model that reads the batch
therefore score well on a held-out donor of a *trained* age group, and ordinary
cross-validation cannot see the difference. Only a test on an age group the
model never trained on can, which is why `cross_age` sits next to every AUROC in
the output and why a cell type is reported with a verdict rather than a score.


## 9. File Map

| File | Job |
| --- | --- |
| `src/aging.py` | Every shared function and every default. The one file to edit to change behaviour |
| `src/strip_output.py` | Git clean filter: commit a notebook as code and markdown, without any output |
| `notebook/01_scanpy_tutorial.ipynb` | The scanpy preprocessing and clustering tutorial, kept as a record of the learning step |
| `notebook/02_explore.ipynb` | Age, donor, sex balance and cell type per donor; the provided UMAP |
| `notebook/03_preprocess.ipynb` | Rebuild from count, QC, doublet, normalize, variable gene, PCA, Leiden, UMAP |
| `notebook/04_composition.ipynb` | Cell type fraction per donor against age, and the outlier donor investigation |
| `notebook/05a_classifier.ipynb` | Young against old per cell type, every control, and the batch check |
| `notebook/05b_genes.ipynb` | The gene behind each model, filtered and flagged; development against aging |
| `notebook/05c_clock.ipynb` | Ridge clock in month, three strength, and a held-out age group |

Inside `src/aging.py`:

| Group | Function |
| --- | --- |
| Setting | `AGE_MONTHS`, `YOUNG`, `OLD`, `ADULT`, `AMBIENT`, `TECHNICAL_PATTERN`, `DISSOCIATION`, `CIRCADIAN`, `MYELOID`, `EXCLUDE`, `main_kw` |
| Gene filter | `celltype_means`, `foreign_genes`, `feature_mask` |
| Data prep | `prep`, `usable_cell_types` |
| Classifier | `logreg`, `age_eval`, `age_auc`, `cross_age_auc`, `top_genes` |
| Clock | `clock_eval`, `clock_summary`, `held_out_age` |
| Trajectory | `trajectories`, `label_patterns`, `pattern_top` |

| Folder | Holds |
| --- | --- |
| `data/` | The atlas and the intermediate `.h5ad`. Not tracked by git |
| `figure/` | Every plot and result table. Regenerated by a run, not tracked by git |
| `notebook/` | The analysis, one notebook per phase |
| `src/` | The shared code |


## 10. Result: Mouse Lung

The pipeline is general, but it was built for one question: which cell type in
the aging mouse lung carries an age signal that is real. This section is the
worked example, on the 16 animal and 24,490 cell of
[section 5](#5-input).

### 10.1 Age Signal Per Cell Type

| Cell type | AUROC, unseen mice | Cross-age AUROC | Shuffled | Verdict |
| --- | --- | --- | --- | --- |
| Fibroblast of lung | 0.98 | **0.97** | 0.60 | Robust age signal |
| Adventitial cell | 0.97 | **0.95** | 0.53 | Robust age signal |
| Type II pneumocyte | 0.97 | **0.80** | 0.35 | Passes, on 125 cell |
| Neutrophil | 0.97 | 0.78 | 0.70 | Void: shuffled control too high |
| Alveolar macrophage | 0.95 | 0.66 | 0.40 | Partial |
| Bronchial smooth muscle | 0.92 | 0.66 | 0.47 | Partial |
| Classical monocyte | 0.91 | 0.40 | 0.34 | Batch driven |
| Intermediate monocyte | 0.88 | 0.40 | 0.50 | Batch driven |
| CD8 T cell | 0.83 | 0.31 | 0.56 | Batch driven |
| Myeloid dendritic cell | 0.79 | 0.29 | 0.51 | Batch driven |
| Non-classical monocyte | 0.76 | 0.29 | 0.38 | Batch driven |
| CD4 T cell | 0.73 | 0.24 | 0.34 | Batch driven |
| NK cell | 0.68 | 0.24 | 0.47 | Batch driven |
| B cell | 0.48 | 0.15 | 0.49 | Batch driven |

**The stromal signal holds under every control.** Fibroblast and adventitial
cell keep an AUROC of 0.96 and above when the model is restricted to male animal
and when the ambient gene are dropped (change ≤ 0.02), and logistic regression
match or beat LightGBM for every cell type in the table. Young fibroblast come
largely from the two pooled 3m sample, and the signal survives when those are
the test set.

**The immune result is a batch detection, not a near miss.** Their cross-age
AUROC sits far below 0.5, so what the model learned is reliably inverted on a new
age group. Several of them nonetheless rank the animal perfectly
(`auroc_mouse` = 1.00 for alveolar macrophage, CD8 T, intermediate monocyte and
bronchial smooth muscle), which is the clearest demonstration in the project
that perfect donor separation is not evidence of aging when age and batch move
together.

Neutrophil are voided by their own control: AUROC 0.97 with a shuffled score of
0.70 means this pipeline can reach 0.70 on that cell type from noise alone, so
0.97 cannot be read as signal.

### 10.2 Aging Gene In The Cell Type That Pass

Read only from a filtered model with AUROC ≥ 0.8 and a passing batch check:

| Cell type | Up in old | Up in young |
| --- | --- | --- |
| Fibroblast of lung | ***Cdkn1a*** (p21), *C3*, *Cfb*, *Fmo2*, *Gsta3*, *Nrp1*, *Snhg11* | *Col3a1*, *Col1a1*, *Col6a1*, *Sparc*, *Mfap5*, *Mgp* |
| Adventitial cell | *Mt1*, *Bst2*, *Cdkn1c*, *Klf9*, *Efemp1*, *Ifi27l2a* | *Col3a1*, *Col1a1*, *Col1a2*, *Col5a2*, *Sparc*, *Cd248*, *Meg3* |
| Type II pneumocyte | *Lrg1*, *Slpi*, *Lcn2*, *Hp*, ***Cdkn1a***, *Nrn1*, *B2m* | *Mgp*, *Gsn*, *Igfbp7*, *Inmt*, *Sparc*, *Cbr2* |

Three thing in that table:

- ***Cdkn1a*** (p21) rises in old fibroblast and old type II pneumocyte. It is a
  core senescence gene, and in the separate trajectory analysis it changes in
  aging ($\Delta_{\text{age}} = +0.54$) and not in development, so it is not a
  maturation effect. The complement gene *C3* and *Cfb* rising alongside it fit
  inflammaging in a stromal cell.
- **Collagen and matrix gene fall with age** in both stromal cell type
  (*Col1a1*, *Col3a1*, *Col6a1*, *Sparc*, *Mfap5*), which matches the known
  decline of extracellular matrix production. In an immune cell type the same
  gene appear and are ambient RNA from the fibroblast-rich 3m sample instead,
  which is what the `own` filter is for.
- **The acute-phase gene in type II pneumocyte** (*Lrg1*, *Slpi*, *Lcn2*, *Hp*)
  are the epithelial side of the same inflammatory shift.

Flagged rather than interpreted: ***Dbp*** and ***Nr1d1***, circadian gene,
appear as "young" gene in several cell type. They are a time-of-day signature,
so the young and old animal were very likely collected at a different hour. This
is a batch finding, and the sort of thing that would otherwise be reported as an
aging gene.

### 10.3 Cell Composition With Age

Per-donor fraction against age, over the 13 animal left after the outlier sample
are removed, Spearman:

| Cell type | $\rho$ | p | p (BH) | Read as |
| --- | --- | --- | --- | --- |
| CD8 T cell | **+0.70** | 0.008 | 0.061 | Rise with age; the strongest composition result |
| Alveolar macrophage | **−0.65** | 0.016 | 0.062 | Fall with age |
| Fibroblast of lung | −0.52 | 0.070 | 0.187 | Possible decline, not conclusive |
| B cell | +0.41 | 0.170 | 0.230 | Possible rise, not conclusive |
| Bronchial smooth muscle | +0.40 | 0.173 | 0.230 | No clear trend |
| CD4 T cell | −0.44 | 0.136 | 0.230 | No clear trend |
| NK cell | −0.30 | 0.318 | 0.363 | No clear trend |
| Classical monocyte | −0.12 | 0.707 | 0.707 | No trend once the outlier donor is out |

Both top result land at p ≈ 0.06 after correction across the eight cell type
tested, on 13 animal: **suggestive, and consistent with known lung immune
aging**, not established. The CD8 gene that change in aging only (*S100a4*,
*S100a6*, *Lgals1*, *Nkg7*, *Gzmm*, *Id2*) are effector and memory gene, which
points at the known shift from naive toward memory T cell. So part of the CD8
"age signal" is a change in which CD8 subtype is present, not every CD8 cell
aging alike.

### 10.4 The Clock

Ridge regression on adult animal, leave one donor out, at $\alpha = 100$:

| Cell type | MAE, month | Naive baseline | Spearman over donor |
| --- | --- | --- | --- |
| Alveolar macrophage | **4.4** | 9.1 | 0.83 |
| Fibroblast of lung | **5.9** | 9.1 | 0.75 |
| CD8 T cell | 5.8 | 7.6 | 0.61 |
| Intermediate monocyte | 6.4 | 7.3 | 0.44 |
| Bronchial smooth muscle | 6.7 | 9.0 | 0.61 |
| Type II pneumocyte | 7.0 | 9.7 | 0.69 |
| Classical monocyte | 7.1 | 9.1 | 0.61 |
| Adventitial cell | 7.6 | 9.0 | 0.66 |
| Neutrophil | 8.0 | 8.1 | 0.05 |

Most cell type beat the naive baseline, best in alveolar macrophage at 4.4 month
against 9.1. All three strength tried ($\alpha$ = 10, 100, 1000) are in
`05c_clock`; the ordering of cell type is stable across them, and neutrophil is
the one cell type that fails to beat the baseline at $\alpha = 1000$.

**Where it breaks, and the honest read.** Holding out a whole age group, the
clock put the 18m animal at 15 to 23 month, which is reasonable, and the 21m
animal at 10 to 16, which is younger than the 18m animal and therefore wrong.
Both 21m animal are female and form their own batch, so the clock is ordering
batch rather than resolving age inside the old range. The clock works as a
coarse young-to-old estimator on this data and no further, and the reason is the
study design rather than the model.

### 10.5 What The Audit Found First

Every claim above depends on phase 2 to 4 having found these first:

- **Donor 30-M-2 is a sample of one compartment.** It contributes 8,014 cell, a
  third of the dataset, and 97% of them are myeloid: 92% of all lung macrophage
  and 68% of all classical monocyte in the tissue come from this one animal. Its
  monocyte are in an unusual state, with low *Ccr2* and *Ly6c2*, *S100a8* and
  *S100a9* about 12× other 30m animal, some *Mki67* and *Top2a*, and **no**
  interferon or hemoglobin signal. That reads as an inflammatory, partly
  dividing, less mature myeloid expansion rather than an antiviral response or
  blood contamination, so it is excluded from the myeloid cell type, kept
  everywhere else, and reported with `with_excluded_mice` alongside. Including it
  moves alveolar macrophage from 0.95 to 0.90 and classical monocyte from 0.91 to
  0.86.
- **The two pooled 3m sample** (`3-M-5/6`, `3-M-7/8`) are about 50% fibroblast
  with almost no immune cell, the opposite imbalance, so they are left out of the
  composition trend. The NK "decline" visible on the plot came mostly from them.
- **Ambient RNA is present and only partly removable.** *S100a8* and *S100a9*
  appear in the B cell of three donor, *Scgb1a1* in many cell type, collagen in
  immune cell. SoupX and CellBender both need the unfiltered droplet data, which
  this release does not include, so the filter of
  [section 7.6](#76-step-6-filter-the-gene-then-read-them) is the available
  substitute.
- **Sex is confounded with age.** The 1m and 30m animal are all male and the 21m
  all female. Sex can hide a signal as well as fake one: B cell go from 0.48 to
  0.77 and NK cell from 0.68 to 0.89 on male animal alone.
- **Doublet are rare here** (49 cell, 0.2%), and the released data was already
  filtered by the consortium, so a `n_genes_by_counts > 300` cut removes exactly
  one cell.

### 10.6 What This Shows

- **A framework.** Mouse-held-out AUROC plus a cross-age batch check and a
  donor-permuted control separates a real age signal from a collection artifact.
  On this atlas that distinction changes the answer for 8 of 14 cell type.
- **A lung finding.** The clearest aging signature in the mouse lung is in the
  **stromal cell**: senescence (*Cdkn1a*) and complement gene rising while matrix
  gene fall, holding up on animal and age group the model never saw.
- **A limit worth naming.** The obstacle is the study design, not the model
  class: with age tied to batch and 16 animal, no amount of model capacity
  recovers what the experiment did not separate. LightGBM beating nothing here is
  the evidence for that.


## 11. Limitation

- **One tissue, sixteen donor.** Every number is bounded by that, not by the
  24,490 cell. Two donor are pooled sample treated as one unit.
- **Age is confounded with processing batch.** The cross-age check detects this
  and reports it per cell type; it cannot undo it. A cell type that fails the
  check is not shown to lack an age signal, only to lack one this design can
  demonstrate.
- **Ambient RNA cannot be removed properly.** SoupX and CellBender need the
  unfiltered droplet matrix, absent from this release. The `own` filter is a
  log-scale ratio, so an abundant ambient gene still passes: *Sftpc* and *Mgp*
  are currently handled by name, and a linear-scale ratio that catch them also
  drops *Cdkn1a* in fibroblast and *Cd74* in macrophage. Open problem.
- **No mitochondrial QC** is possible: the release carries no `mt-` gene, so the
  usual dying-cell filter is unavailable.
- **The shuffled control is five permutation, not a null distribution.** It gives
  a floor to clear, not a calibrated p-value; a few hundred permutation would be
  needed for that.
- **A gene weight is interpretation, not a test.** The weight come from a model
  fit on all animal, and the composition trend are per-cell-type correlation over
  13 animal. Neither is a differential expression test; the pseudobulk of phase 7
  is what would give a p-value with the animal as the sample.
- **The clock does not resolve age within the old range**, for the batch reason
  in [section 10.4](#104-the-clock).

### Next Step

1. **Senescence scoring** with SenMayo (Saul et al. 2022), on the cell type that
   pass the batch check, with the failing immune type shown for contrast.
2. **The ambient filter**, per the open problem above: compare the fraction of
   cell expressing rather than the mean, or test whether a candidate ambient gene
   tracks its source cell type within the same donor.
3. **More tissue.** Parameterize the tissue name, then run the pipeline on three
   to five more droplet tissue (Kidney, Liver, Thymus) and report one cell type ×
   tissue table of which signal pass the batch check in more than one organ. The
   limit to know going in: more Tabula Muris Senis tissue add cell, not
   independent animal or batch.
4. **Pseudobulk differential expression** (PyDESeq2) with the animal as the
   sample, for the cell type that pass.
5. **A design where age is not tied to batch**, which is the only real fix: a
   large human cohort with many donor pooled across run (OneK1K, CELLxGENE
   Census), where the linear baseline here can be compared against an MLP, an
   scVI representation, or a single-cell foundation model.


## 12. Reference

- Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes
  ageing tissues in the mouse. *Nature* **583**, 590–595 (2020). The dataset.
- López-Otín, C. et al. Hallmarks of aging: an expanding universe. *Cell*
  **186**, 243–278 (2023). The framing of section 2.
- Baker, D. J. et al. Clearance of p16Ink4a-positive senescent cells delays
  ageing-associated disorders. *Nature* **479**, 232–236 (2011). Why senescence
  is causal and not only correlated.
- Saul, D. et al. A new gene set identifies senescent cells and predicts
  senescence-associated pathways across tissues. *Nature Communications*
  **13**, 4827 (2022). SenMayo, used in phase 6.
- Buckley, M. T. et al. Cell-type-specific aging clocks to quantify aging and
  rejuvenation in neurogenic regions of the brain. *Nature Aging* **3**,
  121–137 (2023). The closest published form of section 10.
- Squair, J. W. et al. Confronting false discoveries in single-cell differential
  expression. *Nature Communications* **12**, 5692 (2021). Why the animal, and
  not the cell, is the sample.
- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene
  expression data analysis. *Genome Biology* **19**, 15 (2018). The toolkit.
