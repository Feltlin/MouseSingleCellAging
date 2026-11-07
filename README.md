# Single-Cell Aging

Find the cell type in an aging tissue whose gene expression really tracks age,
and separate that signal from the batch it was collected in.

Give the pipeline an annotated single-cell RNA-seq atlas with donor age and
donor identity. It trains one age model per cell type, holding out whole
animal, then puts every surviving model through a second test on age group it
has never seen. A cell type counts as having an age signal only if it passes
both. The same pipeline then estimates age in month, separates a change of
aging from a change of development, names the gene behind each model, and tests
what it finds against a second laboratory and a bulk time course.

Nothing in the pipeline is tissue specific. It was written on the mouse lung of
Tabula Muris Senis and then run on every droplet tissue of that atlas that can
be tested, nine in all. [Section 10](#10-result-mouse-aging-atlas) gives the
cross-tissue result first and the lung as the worked example, and
[section 11](#11-critical-review-and-limitation) attacks every claim the way a
thesis committee would, with the check run in response. None of the data is
part of this repository.


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
10. [Result: Mouse Aging Atlas](#10-result-mouse-aging-atlas)
11. [Critical Review And Limitation](#11-critical-review-and-limitation)
12. [Reference](#12-reference)


## 1. Overview

**The question.** Across mouse organ, which cell type carry a gene-expression
change with age that survives the batch structure of the data, which gene make
it up, and does it hold in animal processed by another laboratory?

**Why it is hard.** In the atlas used here, every age group was collected and
processed separately, so age, processing batch and (often) sex move together. An
ordinary held-out test cannot tell "this model learned aging" from "this model
learned which day the sample was processed". Most of this project is about
telling those apart, and about stating plainly when it cannot.

**What it does, in one pass:**

```
read the atlas
  -> rebuild raw count, QC, doublet removal, normalize, 2000 variable gene
  -> per mouse cell-type fraction vs age          (composition, not expression)
  -> per cell type, hold out whole mice:
       logistic regression young vs old  -> AUROC on unseen mice
       gate: shuffled control near 0.5
       gate: cross-age AUROC >= 0.7      (the batch check: unseen age group)
       control: technical feature only, sex only, male only, interval over mice
  -> for a model that passes, read its gene
       filter: ambient, technical, sex, other cell type's marker
       flag: circadian, dissociation
       split: change of development (1m->3m) vs change of aging (3m->old)
  -> ridge regression on age in month, leave one mouse out
  -> across tissue: pseudobulk DESeq2 with the mouse as sample
  -> replication: another lab (Kimmel 2019) and a 10-age bulk time course
  -> exploratory: model comparison, regulator, miRNA, dynamical model
```

**What holds up, and how strongly.**

| Finding | Evidence | Strength |
| --- | --- | --- |
| On this atlas, an ordinary mouse-held-out AUROC **cannot** be read as aging | Lung: 12 of 14 cell type reach 0.7 on unseen mice, 4 keep a signal on unseen age group; perfect per-mouse ranking with chance-level transfer | **strong**, the methodological result |
| Structural cell (stromal, endothelial, epithelial) carry a more transferable age signal than immune cell of the same organ | Within a tissue, cross-age AUROC is higher in structural row by 0.14 (within-tissue permutation p = 0.001), and by 0.11 with sex-chromosome gene removed (p = 0.002). The pass-rate difference is **not** significant within tissue (p = 0.13), and lung carries most of the gap | **moderate**, downgraded from the first draft |
| Much of the "immune cell fail the batch check" pattern outside lung is sex, not batch | Removing six sex-chromosome gene lifts 12 row above 0.7, all upward, mostly marrow; immune pass rate 0.33 → 0.59 | **new**, found in the critical review |
| MHC class I and complement rise with age | Pseudobulk in 22 of 22 robust row (*H2-D1*, *H2-K1*); another lab (*C3*, *Cfb*); bulk time course (MHC-I 6–7 organ, *C4b* 9 organ) | **strong**: three independent design |
| Collagen and matrix gene fall with age | Lung stroma per cell; bulk *Col1a1* / *Col3a1* down in 13 / 14 of 17 organ, plate balanced | **strong** |
| p21 (*Cdkn1a*) rises in old structural cell | Inside this atlas only: 14 of 15 robust structural row; **not** in another lab, **not** in bulk | **not replicated**: reported as an observation, not an aging marker |
| A weak age direction transfers to another lab | 19 of 19 matched row rank old above young (AUROC 0.51–0.76), but technical feature alone also transfer in 12 of 19; expression clearly beats them in 10 row (myeloid, epithelial, endothelial, adventitial) | **moderate** for those 10 row; **unresolved** for lymphocyte |
| Aging looks like steady drift, not acceleration | Linear beats a saturating-removal model for p16 in 10 of 11 organ, but p16 sits near detection floor | **weak**: the data cannot test acceleration |

Everything after the replication (model comparison, regulator, miRNA, dynamical
model) is labelled **exploratory**: it generates hypothesis, and none of it is
used to support the claim above.


## 2. Biological Background

This section gives the background a reader with general biology needs for the
rest of the README. Each later result section opens with a shorter
**Background** paragraph that points back here.

### 2.1 Why Aging Is Studied Per Cell Type

A tissue does not age as one thing. Bulk RNA-seq grinds a whole organ together
and reports the average over million of mixed cell, so when an inflammation
gene rises in old tissue the measurement cannot say whether every cell turned it
up a little or a few immune cell moved in. Those two have opposite biological
meaning and the same bulk readout.

Single-cell RNA-seq measures each cell separately, so each cell type can be
asked about on its own: not "does this tissue age" but "which cell type in it
ages, and how". A per-cell-type age model is an active research direction:
Buckley et al. 2023 built one per cell type in the mouse brain. The atlas's own
aging analysis (Zhang et al. 2021) tested cell rather than animal and had no
batch control, which is the gap this project is built to close.

### 2.2 How Droplet Single-Cell RNA-seq Works, And What It Gets Wrong

The tissue is cut up and digested with enzyme at 37 °C until the cell float
free. Each cell is then caught in an oil droplet with a bead carrying a unique
**cell barcode**. Inside the droplet the cell's polyadenylated mRNA is captured
and copied, every molecule tagged with the cell barcode and a random
**unique molecular identifier (UMI)**, and only the 3′ end of each transcript is
sequenced. Counting UMI per gene per barcode gives a **cell × gene count
matrix**. Four property of that matrix shape every step here:

| Property | What it is | Consequence here |
| --- | --- | --- |
| **Sparsity** | Each cell gives up only a fraction of its mRNA, so over 90% of entry are zero | A zero usually means "not caught". A rare transcript such as p16 is invisible per cell (0.8% of cell) and has to be summed per animal |
| **Depth** | Cell differ in how many molecule were caught | Expression is normalized per cell before comparison; depth itself is used as a technical control |
| **Ambient RNA** | Cell that burst during digestion release RNA into the suspension, and some lands in every droplet | A fibroblast gene can appear at low level "inside" a T cell. Filtered per cell type (section 7.6) |
| **Dissociation stress** | Digestion at 37 °C switches on immediate-early gene (*Fos*, *Jun*, *Hspa1a*) within minute | A difference in handling time looks like biology; these gene are flagged, never read as aging |

Two more practical limit: a **doublet** (two cell in one droplet) looks like a
hybrid cell and is removed; and only polyadenylated RNA is captured, so mature
**microRNA** (section 2.6) are not measured at all.

### 2.3 The Dataset And The Design Problem

**Tabula Muris Senis** (Tabula Muris Consortium 2020) profiled C57BL/6 mice at
1, 3, 18, 21, 24 and 30 month, the mouse equivalent of roughly adolescence to
late old age. This project uses the droplet (10x) part: 16 tissue, **23 mice in
total**, because one mouse gave cell from up to 12 organ.

**The design problem is confounding.** Each age group was collected and
processed at its own time. If every 3-month mouse was processed on one day and
every 18-month mouse on another, then anything that differed between those two
day (reagent lot, sequencing run, time of day, handling) differs between the two
age group as well, and a classifier can separate the age group by reading it.
On top of that, **sex is tied to age group** in every tissue: the 21-month mice
are all female, the 24- and 30-month mice all male, and most 1-month mice male.

Confounding cannot be removed by any statistical correction, because the data
contain no young and old mouse processed together to compare. It can only be
**tested around**: by asking whether a model trained on some age group works on
age group (and therefore batch) it never saw, by checking that technical
measurement or sex alone cannot do the same job, and by checking result against
data with a different design.

Two such dataset are used:

| Dataset | Design | Role |
| --- | --- | --- |
| Kimmel et al. 2019 | Another lab; lung, kidney, spleen; 3–4 young (7–8 month) and 3 old (~21 month) male mice; 10x | A model that transfers here cannot be reading TMS batch |
| TMS bulk (Schaum et al. 2020) | 17 organ, 56 mice, **10 age** from 1 to 27 month, sequencing plate balanced across age | Fills the 3–18 month gap; plate is a covariate rather than a confounder |

### 2.4 The Aging Biology Touched Here

Aging is organized into **hallmark**, process that change with age, accelerate
aging when worsened and slow it when fixed (López-Otín et al. 2023). Four of
them appear in this project.

**Cellular senescence.** A badly stressed cell can stop dividing permanently
instead of dying. Two cell-cycle brake mark it. **p16** (*Cdkn2a*) blocks the
kinase CDK4/6 and is the most specific marker of established senescence.
**p21** (*Cdkn1a*) is switched on by p53 after DNA damage and blocks CDK2; it
also rises transiently after almost any stress, so **p21 alone is not proof of
senescence** (Hernandez-Segura et al. 2018). Senescent cell secrete an
inflammatory, tissue-remodelling mix called the **SASP**. They accumulate with
age, and removing them extends healthy lifespan in mice (Baker et al. 2011). No
single gene marks every senescent cell, so a gene set such as **SenMayo** (Saul
et al. 2022) is scored instead.

**Inflammaging.** Old tissue carries chronic low-grade inflammation without
infection (Franceschi et al. 2000). Three molecular reading of it appear here:

- **NF-κB** is the master transcription factor of inflammation; **interferon
  signalling** (transcription factor IRF1, STAT1) is the antiviral alarm, which
  can also be set off by the cell's own nucleic acid.
- **MHC class I** (*H2-K1*, *H2-D1* in mouse, with *B2m*) displays fragment of a
  cell's protein on its surface so that CD8 T cell can inspect it. Interferon
  and NF-κB both raise it.
- **Complement** (*C3*, *C4b*, *Cfb*) is a cascade of blood and tissue protein
  that tags pathogen and debris for immune attack; fibroblast and liver make
  much of it.

**Extracellular matrix.** Fibroblast build the scaffold between cell, chiefly
collagen I and III (*Col1a1*, *Col3a1*). Matrix production falls with age, one
reason old tissue heals more slowly.

**Immune aging and the thymus.** New ("naive") T cell are made in the thymus,
which shrinks from early adulthood (**thymic involution**; Palmer 2013). The
T-cell pool therefore shifts with age from naive cell (marked by *Sell*, *Ccr7*,
*Lef1*, *Tcf7*) toward memory and effector cell (*Cd44*, *S100a4*, *Ccl5*).

**MicroRNA.** A microRNA (miRNA) is a ~22-nucleotide RNA that pairs with the 3′
end of a target mRNA and speeds its decay. One miRNA has hundred of target, so
when a miRNA becomes more active its target fall **together**, and that
coordinated fall can be read from mRNA data even though the miRNA itself is not
measured. The **miR-29** family represses collagen gene and is reported to rise
with age in mouse tissue (Ugalde et al. 2011).

**Organ clock and ageotype.** If the gene expression of one organ predicts the
animal's age, the error of that prediction says whether the organ looks older
or younger than expected. Whether those error agree across organ of one animal
(an "ageotype": some individual old across the body) is an open question (Ahadi
et al. 2020; Schaum et al. 2020).

### 2.5 The Statistical Idea Used Throughout

**The unit of replication.** The cell is the unit of measurement; the **animal**
is the unit of replication. Cell of one animal share its genotype, health,
processing day and sex, so 10,000 cell from 12 mice are 12 independent sample,
not 10,000. Treating cell as sample ("pseudoreplication"; Hurlbert 1984) makes
every p-value tiny and wrong, which is the most common error in single-cell
aging study (Squair et al. 2021). Every test here has the animal, or a set of
animal, as its sample.

**AUROC.** A classifier gives each cell a score for "old". The area under the
ROC curve is the probability that a randomly chosen old cell scores above a
randomly chosen young one: 0.5 is chance, 1 perfect, and below 0.5 means the
model is systematically **inverted**. It needs no cutoff and is unaffected by
one class being larger.

**Grouped cross-validation.** The data are split so that every cell of a mouse
is either in training or in testing, never both. Otherwise the model can
recognize the mouse and inherit its label.

**The cross-age check.** The model is trained on some age group and tested on
other age group. Because each age group is its own batch, a model that learned
the batch has nothing to use. This is the project's central control
([section 7.4](#74-step-4-the-two-gate)).

**Permutation null.** To ask "how large a score can arise with no age signal",
the age label is shuffled **between mice** (never between cell) and the analysis
repeated. A shuffle can also be restricted, for example only within a tissue, so
that it keeps a structure of the data that a plain test would ignore.

**Multiple testing.** Testing hundred of gene or row produces some p < 0.05 by
chance. The Benjamini–Hochberg procedure controls the expected fraction of false
discovery (FDR) among the result called significant.

**Regularization.** With 2,000 gene and a dozen mice, an unconstrained model can
fit anything. L2 regularization (logistic regression with small $C$, ridge
regression with large $\alpha$) keeps weight small, trading a little fit for
much better generalization.

**Pseudobulk.** Summing the count of one cell type per animal turns single-cell
data into one sample per animal, which standard RNA-seq statistics (DESeq2: a
negative binomial model of count; Love et al. 2014) then test correctly.

**Bootstrap interval.** Resampling the test mice with replacement many time and
recomputing a score shows how much it depends on which mice happened to be
tested; the middle 90% of the resampled value is a 90% interval.

**Confirmatory versus exploratory.** A confirmatory analysis has its rule fixed
before it sees the data and is judged by a test. An exploratory analysis
searches; it suggests hypothesis but its p-value are optimistic, because many
path were available. This README labels each result section as one or the other.

### 2.6 The Dynamical-System Idea Used In The Exploratory Part

A **dynamical system** describes how a state changes over time with a
differential equation, $\dot{x} = f(x, t)$. Here "time" is age.

**Steady state and its loss.** Where production and removal balance, $\dot{x} =
0$ and the state sits at a **fixed point**. Karin et al. 2019 model senescent
cell as produced at a rate that grows with age ($\eta t$) and removed by the
immune system at a rate that **saturates** ($\beta X/(\kappa + X)$, at most
$\beta$). While production is below $\beta$, a steady state exists; once it
passes $\beta$, no balance is possible and the burden accelerates, which they
link to the exponential rise of mortality with age (the Gompertz law). With
noise added, the model also predicts individual drifting apart.

**Critical slowing down.** Near a tipping point a system recovers slowly from
perturbation, so its fluctuation grow and become correlated (Scheffer et al.
2009). Rising variance and a dominant direction of variation are therefore
proposed **early-warning signal**. They are defined for one system observed
repeatedly, which matters for how they can be read here (section 10.11).

**Identifiability.** A model parameter is identifiable when the data pin it
down. The **profile likelihood** fixes the parameter at a grid of value,
refits everything else each time, and records the fit; a flat profile means
very different value explain the data equally well (Raue et al. 2009).

**Learning an equation from data.** **SINDy** (Brunton et al. 2016) fits
$\dot{x}$ as a sparse combination of candidate term ($x$, $x^2$, $xy$, …). A
**neural ODE** (Chen et al. 2018) lets a small neural network be $f$. Both need
many time point per parameter to be trusted.

**Optimal transport.** Given the cell of two age, optimal transport finds the
cheapest way to move one cloud of cell onto the other, and from that coupling
predicts the population in between (Waddington-OT; Schiebinger et al. 2019).
Holding out a middle age and predicting it is a direct test.

### 2.7 What Each Metric Measures

| Metric | The question it answers | Why it is here |
| --- | --- | --- |
| **AUROC on unseen mice** | Can a model rank an old cell above a young one in an animal it never saw? | The honest version of accuracy, but blind to batch |
| **Cross-age AUROC** | Does that model still work on an age group it never trained on? | Each age group was processed in its own batch, so a batch shortcut cannot help here |
| **90% interval of cross-age** | How much does that score depend on which test mice happened to be used? | 2–6 test mice per split make a single number unstable |
| **Shuffled control** | What does the same pipeline score when age is swapped between animal at random? | The floor a real score has to clear |
| **Technical-only and sex-only control** | Can depth and blood or stress content, or sex, alone pass the cross-age check? | Catches a batch or sex difference that the cross-age split does not break |
| **Per-mouse AUROC** | Can the model rank whole animal by age? | The number with real degree of freedom, but perfect ranking is exactly what batch also produces |
| **Clock error in month** | How far off is the predicted age, against a naive guess? | Turns a two-class result into a quantity; the naive baseline keeps it honest |
| **Development vs aging change** | Did the gene change from 1m to 3m, or from 3m onward? | Growth is not aging |

The AUROC on unseen mice and the cross-age AUROC are deliberately **not** fused
into one rank. They disagree, and the disagreement is the result: a model can
separate animal perfectly and still be reading the collection batch.


## 3. Install

```bash
conda env create -f environment.yml
conda activate scaging
```

Core: `scanpy`, `anndata`, `scikit-learn`, `pandas`, `numpy`, `scipy`,
`statsmodels`, `seaborn`, `matplotlib`, `igraph`, `leidenalg`. Python 3.14.

By phase:

| Package | Used in | For |
| --- | --- | --- |
| `lightgbm` | 05a, 11 | Tree-ensemble comparison model |
| `pydeseq2` | 09 | Pseudobulk differential expression |
| `decoupler`, `liana` | 12, 13 | Regulator, pathway, miRNA activity; ligand–receptor list |
| `mofapy2` | 14 | Multi-organ factor analysis |
| `pysindy`, `torch` (CPU), `torchdiffeq`, `moscot` | 16 | Equation discovery, neural ODE, optimal transport |
| `nbdime` | — | Readable notebook diff |

Two git setting are worth enabling once per clone:

```bash
nbdime config-git --enable                                          # cell by cell notebook diff
git config filter.strip-output.clean "python src/strip_output.py"   # commit code only
git config filter.strip-output.required true
```

The second one matter. An output is regenerated by rerunning the notebook, and
in git it only hides the code change: a plot is a megabyte of base64 that differ
on every run. So the committed notebook is code and markdown only, while the
copy on disk keeps every output. `required true` makes a missing filter fail the
commit instead of quietly committing the outputs. Because of this, a reader on
GitHub sees no number inside a notebook: the explanation is in its markdown and
the number are in [section 10](#10-result-mouse-aging-atlas).


## 4. Usage

### The Run

Each phase is one notebook in `notebook/`, run in order, and each reads the file
the previous one wrote. The tissue is a parameter: the first cell of `02`–`05c`
reads `TISSUE` (default `Lung`), and every path is built from it.

```bash
conda activate scaging
python src/download.py --list                  # what figshare holds, and the size
python src/download.py droplet Lung Kidney     # or: droplet all
python src/download.py external                # Kimmel 2019, bulk atlas, published gene list
python src/run_tissue.py Lung Kidney           # run 02 -> 05c for each tissue
```

`run_tissue.py` executes the notebook without opening them and saves each
executed copy, output included, to `notebook/run/{Tissue}/`; the notebook in
`notebook/` stay unchanged. `--from 05a` restarts a tissue part way. Notebook
`06` onward work across tissue and are run once every tissue is done.

The order is not a style preference. `03` writes the `.h5ad` that every later
notebook load, `05b` and `05c` both read `figure/{Tissue}/age_auroc.csv`, and
`08`–`16` read `figure/atlas/atlas.csv` from `07`.

### Every Setting

Every shared default lives in `src/aging.py`, so a change applies to every
notebook at once:

| Setting | Default | Meaning |
| --- | --- | --- |
| `AGE_MONTHS` | 1, 3, 18, 21, 24, 30 | Age label to month, for the clock |
| `YOUNG` / `OLD` | 1m, 3m / 18m, 21m, 24m, 30m | The two class; an age absent from a tissue is ignored |
| `ADULT` | 3m and up | Age used by the clock, so growth is not modelled as aging |
| `max_cells` | 300 | Cap per mouse, so one animal cannot dominate a cell type |
| `AMBIENT` | *S100a8/9*, *Lyz2*, hemoglobin, *Scgb1a1* … | Gene treated as free-floating RNA rather than a cell's own |
| `TECHNICAL_PATTERN` | `Rp[sl]`, `Malat1`, `Xist`, `Ddx3y` … | Ribosomal, rRNA, and sex-linked gene |
| `SEX_PATTERN` | *Xist*, *Tsix*, *Ddx3y*, *Eif2s3y*, *Kdm5d*, *Uty* | Sex-chromosome gene, removed in the sex check |
| `DISSOCIATION` | *Fos*, *Jun*, *Hspa1a* … | Gene induced by tissue digestion, flagged not dropped |
| `CIRCADIAN` | *Dbp*, *Nr1d1*, *Per1* … | Time-of-day gene, flagged as a batch signature |
| `TISSUE_SETTINGS` | Lung only | Every per-tissue choice, in one place (below) |
| `C` | 0.1 | L2 strength of the logistic regression |
| `alpha` | 10, 100, 1000 | Ridge strength of the clock; all three are reported |
| `min_cells` / `min_mice` | 50 / 2 | A cell type below this in either class is not modelled |
| `min_ages` | 1, or 2 for the survey | Distinct age per class; the cross-age check needs 2 |

A tissue runs with no entry in `TISSUE_SETTINGS`, on automatic choice. An entry
is added only after a person has reviewed that tissue's phase 2–4 check:

| Key | Lung | Meaning |
| --- | --- | --- |
| `exclude` | 30-M-2 in myeloid | Donor dropped from named cell type only, with a reason |
| `odd_mice` | 30-M-2, 3-M-5/6, 3-M-7/8 | Sample left out of the composition trend |
| `inspect_mice` | 30-M-2, 30-M-5 | Donor plotted one at a time on the UMAP |
| `cross_age_splits` | 1m+18m+21m vs 3m+30m | A hand-picked batch-check split, reported next to the generic one |
| `composition_cell_types`, `example_*` | 8 cell type; fibroblast, CD8 T, alveolar macrophage | Cell type shown in `04`, `05b`, `05c` |

### Calling It Directly

Every function take the `AnnData` as its first argument, so the same code runs
on any tissue from a plain script:

```python
import scanpy as sc
from src import aging as ag

TISSUE = "Lung"
adata = sc.read_h5ad(ag.paths(TISSUE).processed)    # data/Lung/processed.h5ad

for ct in ag.usable_cell_types(adata):
    kw = ag.main_kw(ct, TISSUE)                       # the tissue's reviewed exclusion
    parts = ag.cross_age_scores(adata, ct, **kw)      # held-out score per split
    print(ct,
          round(ag.age_auc(adata, ct, **kw), 2),                      # unseen mice
          round(ag.cross_age_auc(adata, ct, **kw), 2),                # unseen age group
          [round(x, 2) for x in ag.bootstrap_cross_age(parts)],       # 90% interval over test mice
          round(ag.age_auc(adata, ct, shuffle_seed=0, **kw), 2))      # shuffled
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

**The donor id is the one column that cannot be missing.** Every score is
computed by holding out a whole animal. With age but no donor id, a model can be
trained and will report a high number, and that number means nothing, because
the cell of one animal in both train and test let the model recognize the animal
instead of the age.

**A log value in `.X` is handled.** The official Tabula Muris Senis file ship
with `.X` already log-normalized and the count kept in `.raw`, so `03` rebuilds
from `.raw` and keeps a copy in `layers["counts"]`. Normalizing an already
normalized matrix twice is silent and ruins every downstream number, so the
notebook asserts the matrix is integer before continuing.

### What A Tissue Has To Be Checked For First

Phase 2 to 4 exist to find the thing that would otherwise be read as biology.
For a new tissue they are not optional:

- **An outlier donor.** One animal supplying most of a cell type will drive that
  cell type's result on its own. Cross-tabulate cell type against donor before
  modelling anything.
- **Sex against age.** If the young animal are male and the old female, a model
  can score well on *Xist* alone. Section 10.3 tests this for every row.
- **A pooled sample.** An id like `3-M-5/6` is two animal in one sample. It is
  one unit, and it cannot be split.
- **Composition versus expression.** A cell type can change in number without
  any cell changing inside, and tissue dissociation affect number. Phase 4
  answers the count question separately so phase 5 can be about expression only.
- **Which age exist.** Lung has no 24m group, and Thymus droplet has one young
  donor. `cross_age_splits` builds the split from the age a tissue actually has,
  and a tissue without two young and two old age is reported as not testable
  rather than dropped. `06_survey` counts this for every tissue before any of
  them is modelled.

### The Worked-Example Tissue

Tabula Muris Senis droplet Lung, from
[figshare](https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728):
`tabula-muris-senis-droplet-processed-official-annotations-Lung.h5ad`, placed in
`data/Lung/` by `src/download.py` (md5 checked against figshare). **Not
committed**: it is large, not ours to redistribute, and it changes upstream.

| | |
| --- | --- |
| Cell | 24,540, and 24,490 after doublet removal |
| Gene | 20,138 |
| Mouse | 16: six young (1m, 3m) and ten old (18m, 21m, 30m) |
| Age | 1m, 3m, 18m, 21m, 30m. **No 24m in Lung** |
| Annotation | Cell type, age, sex and mouse id, provided by the consortium |

The release carries no mitochondrial gene, so the usual mitochondrial QC cannot
be run, and the cell were already filtered by the consortium: a
`n_genes_by_counts > 300` cut removes exactly one cell. A gene name is the older
symbol of the submission, so hemoglobin beta is `Beta-s` rather than `Hbb-bs`.


## 6. Output

### The Table

Every output of a tissue goes to `figure/{Tissue}/`, and a table across tissue
to `figure/atlas/`. `figure/{Tissue}/age_auroc.csv` has one row per cell type.
Every column is an AUROC computed with whole mice held out, so the column differ
only in what was changed about the input:

| Column | Meaning |
| --- | --- |
| `logreg` | The main number: logistic regression, unseen mice |
| `cross_age` | **The batch check.** Trained on some age group, tested on others; the tissue's hand-picked split if it has one |
| `cross_age_all_splits` | The same check averaged over every valid split; the one compared across tissue |
| `shuffled` | The same pipeline with age randomized between animal |
| `males_only` | Rerun on male animal alone, so sex cannot stand in for age |
| `no_ambient_genes` | Rerun with the ambient gene list dropped |
| `with_excluded_mice` | Rerun keeping the excluded donor, to show its effect |
| `lgbm` | LightGBM instead of logistic regression |
| `passes_batch_check` | `cross_age >= 0.7` |

`figure/atlas/atlas.csv` adds, per tissue × cell type row: the verdict, the
compartment, the technical-only score (`tech_cross_age`, `tech_confounded`), the
dominant-donor rerun, and the stricter check of section 10.3: `ci90_low` /
`ci90_high`, `sex_only_cross_age`, `no_sex_genes_cross_age` and
`males_only_cross_age`.

`age_report.csv` describes the main model instead of varying it:

| Column | Meaning |
| --- | --- |
| `auroc_cell` | AUROC over cell, same as `logreg` above |
| `auroc_mouse` | AUROC after averaging each animal's cell score, then ranking animal |
| `bal_acc_test` | Balanced accuracy on held-out cell at a 0.5 cutoff |
| `acc_train` | Training accuracy, reported to show the size of the overfit |
| `n_cells`, `n_mice` | What the row is actually based on |

The gene table (`top_genes_up_in_old.csv`, `..._up_in_young.csv`,
`aging_only_genes.csv`) give the strongest weight per cell type after
filtering, with a `flag` naming a circadian or dissociation gene. The clock
table (`clock_summary.csv`, `clock_held_out_age.csv`) give error in month per
cell type at every regularization strength tried. `outlier_mice.csv` lists every
donor that supplies over half of a cell type, as a prompt to look.

Everything in `figure/` is **regenerated, not committed**, including the CSV.
The number worth keeping are quoted in [section 10](#10-result-mouse-aging-atlas).

### Reading The Number

| Column | Range | Lowest means | Highest means | Want |
| --- | --- | --- | --- | --- |
| `logreg` | 0 to 1 | **0**: every old cell ranked below every young one | **1**: perfect ranking in an animal never seen | **high**, and above `shuffled` |
| `cross_age` | 0 to 1 | **0**: what was learned is exactly inverted on a new age group | **1**: transfers perfectly | **≥ 0.7** |
| `ci90_low` | 0 to 1 | the score could be chance with other test mice | — | **> 0.5** at least |
| `shuffled` | 0 to 1 | **0 or 1**: both mean the pipeline can hit an extreme on noise | **0.5**: no signal left after randomizing | **near 0.5** |
| `sex_only_cross_age`, `tech_cross_age` | 0 to 1 | — | sex or technical feature alone pass | **below 0.7** |
| `acc_train` | 0 to 1 | — | **1**: the model memorized its training animal | expect ~1, and read the gap |
| `mae_months` | ≥ 0 | **0**: exact age | large: no better than guessing | **below `naive_mae`** |

Three caveat about the extreme.

**A `cross_age` far below 0.5 is not a weak result, it is an inverted one.** At
0.15 the model is reliably wrong on a new age group, which can only happen if
what it learned was a property of the batch that flips sign between batch. Read
it as a positive detection of a batch effect, not as a near-miss.

**`acc_train` near 1.00 everywhere is expected and is not a finding.** Any model
with 2,000 gene and a dozen animal can memorize its training set; the gap
against the held-out score is the overfit.

**A per-mouse AUROC of 1.00 is weaker evidence than it looks.** With this few
animal, and age confounded with batch, perfectly ranking the animal is exactly
what a batch effect produces.

### The Figure

| Figure | Shows |
| --- | --- |
| `umap*.png` | The cell type, age and donor map, plus one donor at a time in colour against grey |
| `qc_by_*.png` | Gene and count per cell, by age and by donor |
| `composition_by_mouse.png` | One dot per animal: each cell type's fraction against age, coloured by sex |
| `age_auroc.png` | AUROC against cross-age AUROC per cell type, with the shuffled control |
| `dev_vs_aging_*.png` | Each gene's change of development against its change of aging |
| `clock_*.png` | Predicted against true age in month, one point per cell and a median per animal |
| `atlas/atlas_cross_age.png` | Cross-age AUROC, cell type × tissue, grouped by compartment |
| `atlas/replication_kimmel*.png` | Transfer to another lab, per row and per mouse |
| `atlas/dynamics_*.png` | The dynamical-model view of section 10.11 |


## 7. Data Process

Notation. A tissue has $N$ cell from $D$ donor. $x_{i,g}$ is the expression of
gene $g$ in cell $i$ after normalization, $d(i)$ the donor of cell $i$, and
$a(d)$ the age in month of donor $d$. A cell type $c$ gives the cell subset
$I_c$.

### 7.1 Step 1: Rebuild, QC And Normalize

Count are taken from `.raw` and asserted integer, because `.X` in the released
file is already log-normalized and normalizing twice is silent. A copy is kept
in `layers["counts"]` for the pseudobulk, which needs count rather than a log
value.

Doublet are called by Scrublet (Wolock et al. 2019) **per donor**
(`batch_key="mouse.id"`), on count, before normalization. Order matter: Scrublet
normalizes internally, and running it on an already normalized matrix changed
the call from 49 cell to 127 on the lung.

Normalization is count per 10,000 then $\log(1+x)$:

$$x_{i,g} = \log\left(1 + 10^4 \cdot \frac{n_{i,g}}{\sum_{g'} n_{i,g'}}\right),$$

which puts every cell on a comparable total, so a cell that happened to be
sequenced deeper does not read as a cell expressing more. Then 2,000 highly
variable gene, PCA, a neighbour graph, Leiden cluster and UMAP. The cluster are
used as a check on the provided label, never to replace them.

### 7.2 Step 2: Composition, Before Any Expression

For donor $d$ and cell type $c$ the fraction is

$$f_{d,c} = \frac{|\lbrace i : d(i) = d,\ i \in I_c \rbrace|}{|\lbrace i : d(i) = d \rbrace|},$$

computed with a cross-tabulation so that a donor with none of a cell type keeps
a 0 rather than vanishing from the average. The trend against age is
Spearman's $\rho$ over donor, with Benjamini–Hochberg correction across the cell
type tested, because the count is one per animal and there are a dozen animal.

This is kept separate from everything downstream on purpose. A cell type can
change in number with no cell changing inside it, dissociation affects the
number, and a single model over all cell would learn that composition instead of
aging.

### 7.3 Step 3: Build One Matrix Per Cell Type

For cell type $c$, cell are shuffled at a fixed seed and capped at
$n_{\max} = 300$ per donor:

$$I_c^{\ast} = \bigcup_{d} \operatorname{head}_{n_{\max}}\left(\operatorname{shuffle}(\lbrace i \in I_c : d(i) = d \rbrace)\right).$$

The cap is not cosmetic. One lung donor contributes a third of all cell, and
without a cap its cell alone would set the loss. The label is

$$y_i = \mathbb{1}\left[a(d(i)) \in \text{OLD}\right],$$

and a donor named in the tissue's `exclude` setting for this cell type is
dropped here, with the reason recorded in `TISSUE_SETTINGS`. A cell type is
modelled only with at least 50 cell and 2 donor **in both class**.

### 7.4 Step 4: The Two Gate

**Gate 1, the shuffled control.** Age is permuted **between donor**, never
between cell, then the whole pipeline is rerun unchanged:

$$\tilde{y}_i = \mathbb{1}\left[\pi(d(i)) \in \text{OLD}\right], \qquad \pi \text{ a permutation of the donor}.$$

The reported `shuffled` column is the mean over five such permutation. Any
score above 0.5 here comes from between-donor difference (depth, dissection,
health) happening to line up with the random labelling, so `shuffled` is a
**floor to clear, not a p-value**: fibroblast at 0.98 against 0.60 clears it
widely, while neutrophil at 0.97 against 0.70 means the pipeline can
manufacture 0.70 from noise on that cell type, and the row is voided. The
control says nothing about batch: age and batch are the same grouping here, so
relabelling donor breaks both at once. That is Gate 2's job.

**Gate 2, the cross-age check.** Train on one set of age group and test on a
disjoint set, each side holding a young and an old group. For any tissue,
`cross_age_splits` enumerates every such split: one young age with a non-empty
proper subset $S$ of the old age in training, the other young age with the rest
in test,

$$\text{train } \lbrace y_1 \rbrace \cup S \to \text{test } \lbrace y_2 \rbrace \cup (\text{OLD} \setminus S), \qquad \text{for every } y_1 \text{ and } S,$$

twelve split for Lung. `cross_age_all_splits` averages all of them and is the
column compared across tissue; a cell type passes when it is at least 0.7. Lung
also keeps the hand-picked pair it was first tested with,
$\lbrace 1\text{m}, 18\text{m}, 21\text{m} \rbrace \leftrightarrow \lbrace
3\text{m}, 30\text{m} \rbrace$.

**What the gate cannot do.** It catches a batch that differs *between age
group*. A batch shared by every young age against every old age passes every
split, and so does a sex difference when the split keeps it aligned with age;
the technical-only and sex-only control of section 10.3 exist for those case.
It also fails a real signal that differs between age group, for example a
change specific to 1-month (still growing) mice, so failing the gate means "not
demonstrable with this design", not "no aging".

### 7.5 Step 5: Score

Cross-validation is `StratifiedGroupKFold` over donor, up to five fold,
stratified so each fold carries both class. Every cell is scored by a model
that never saw its donor:

$$\widehat{p}_i = \widehat{p}^{\,(-d(i))}(x_i).$$

**AUROC over cell** is the probability that a random old cell outranks a random
young one. **AUROC over donor** averages each donor's score first,
$\bar{p}_d = \operatorname{mean}_{d(i)=d} \widehat{p}_i$, then ranks donor.
**Balanced accuracy** is the mean of the two per-class rate at a 0.5 cutoff, and
**training accuracy** is kept to expose the overfit.

The main model is L2 logistic regression at $C = 0.1$, which keeps the weight
small so a few hundred cell from a dozen animal cannot support 2,000 free
parameter. LightGBM is run as a comparison.

**Interval over test mice.** For the cross-age score, `bootstrap_cross_age`
resamples whole test mice with replacement within each class and split (1,000
draw) and reports the 5th and 95th percentile of the mean AUROC. Training mice
are held fixed, so the interval is narrower than the full uncertainty, and a
class with a single test mouse contributes no spread.

### 7.6 Step 6: Filter The Gene, Then Read Them

Only a model with AUROC ≥ 0.8 **and** a passing batch check has its gene read at
all. The filter is cumulative:

| Level | Drop |
| --- | --- |
| `hvg` | Nothing: the 2,000 highly variable gene |
| `clean` | `AMBIENT`: free-floating RNA and blood gene |
| `bio` | Ribosomal, rRNA, *Malat1*, *Xist*, Y-linked, hemoglobin |
| `own` | Another cell type's marker gene |

With $\mu_{c,g}$ the mean log expression and $m_{c,g}$ the mean linear
expression (count per 10,000) of gene $g$ in cell type $c$, gene $g$ is foreign
to $c$ when

$$\max_{c' \neq c} \mu_{c',g} > 1 \qquad \text{and} \qquad r_{c,g} = \frac{m_{c,g}}{\max_{c' \neq c} m_{c',g}} < \tau,$$

that is, clearly expressed somewhere else and present here only at a small
fraction of that source. *Scgb1a1* from club cell showing up in a T cell is the
case it was written for.

**The ratio is linear, and $\tau$ is fit, not chosen.** Ambient RNA lands in
every droplet as a near-constant fraction of the tissue's soup, so over every
pair $(c, g)$ the value $\log_{10} r_{c,g}$ falls into two group: ambient near
0.004, real expression near 0.12. A two-component Gaussian mixture sets $\tau$
where a pair becomes more likely real, refit for each tissue (Lung:
$\tau = 0.031$, the same at every seed). An earlier version compared log mean at
10%, and a log compresses a 500-fold difference to about 7-fold, so an abundant
ambient gene passed. Two other rule were tried and rejected on the same data:
fraction of cell expressing points the wrong way (ambient *Sftpc* is detected
in 66% of neutrophil, real *Cdkn1a* in 16% of fibroblast), and tracking the
soup within a donor separates nothing (per-donor depth moves every gene
together).

The rule is checked against six case it was not fit on:

| Gene in cell type | Ratio | Want | Result |
| --- | --- | --- | --- |
| *Sftpc* in neutrophil | 0.0016 | drop | dropped |
| *Sftpc* in bronchial smooth muscle | 0.0015 | drop | dropped |
| *Mgp* in alveolar macrophage | 0.0083 | drop | dropped |
| *Scgb1a1* in CD8 T cell | 0.0044 | drop | dropped |
| *Cdkn1a* in fibroblast | 0.068 | keep | kept |
| *Cd74* in alveolar macrophage | 0.037 | keep | kept, the closest call |

A circadian or dissociation gene is **flagged, not dropped**, because the flag
is the information: *Dbp* and *Nr1d1* appearing as "young" gene say the young
and old animal were collected at a different time of day.

### 7.7 Step 7: Separate Development From Aging

Expression is averaged per donor, then per age group, and two change are taken:

$$\Delta_{\text{dev}} = \mu_{3\text{m}} - \mu_{1\text{m}}, \qquad
\Delta_{\text{age}} = \operatorname{mean}_{a \in \text{OLD}} \mu_{a} - \mu_{3\text{m}}.$$

The pair labels a gene: **aging only** when $|\Delta_{\text{age}}| > 0.5$ and
$|\Delta_{\text{dev}}| < 0.2$, **development only** when reversed, **lifelong**
when both move the same way, **reverses** when they oppose. The cutoff are a
convention, not fitted, and are used only to sort gene for reading.

### 7.8 Step 8: The Clock

Ridge regression on age in month, **adult only** (3m and up, so growth is not
modelled as aging), leave one donor out:

$$\widehat{a}_i = \widehat{a}^{\,(-d(i))}(x_i), \qquad
\widehat{a}_d = \operatorname{median}_{d(i) = d} \widehat{a}_i,$$

and the error is per donor, against a naive baseline that always guesses the
mean age of the other donor:

$$\text{MAE} = \frac{1}{D}\sum_d \left| \widehat{a}_d - a(d) \right|, \qquad
\text{MAE}_{\text{naive}} = \frac{1}{D}\sum_d \left| a(d) - \operatorname{mean}_{d' \neq d} a(d') \right|.$$

A clock on a few animal spread over 3m to 30m can look accurate by predicting
the middle of the range for everything, and only the baseline exposes it. Every
regularization strength tried is reported. A whole age group can also be held
out, which asks whether the clock orders age it never saw; leave-one-donor-out
alone cannot, because the other mice of the same age (and batch) stay in
training.

### 7.9 Step 9: Across Tissue

| Step | Method | Notebook |
| --- | --- | --- |
| Survey | Count mice, age, sex and cell per tissue from metadata only; fix the testability rule | 06 |
| Atlas | Collect every tissue's row; verdict by the fixed rule; technical-only, sex-only, male-only and dominant-donor rerun; interval over test mice; compartment comparison with a within-tissue permutation | 07 |
| Senescence | SenMayo score per cell, mean per animal; p21 and p16 per animal; Spearman with age, stable if the sign holds with any one age group left out | 08 |
| Pseudobulk | Count summed per animal and cell type; DESeq2 old vs young; agreement between row against a null that permutes the 23 mice's label once for every tissue alike | 09 |
| Replication | Train on TMS, score Kimmel unchanged; same for the technical-only model; key gene per Kimmel mouse | 10 |
| Bulk | Per gene and organ, age effect with plate and sex as covariate; MOFA+; per-organ clock and its residual | 14 |

### 7.10 Every Learning Method Used

Every model in the project, grouped by what it learns. Each setting is fixed in
code (`src/aging.py`, `src/models.py`, `src/dynamics.py`) before the result it
produces is seen; where one was changed after a smoke test, the notebook says
so.

**Unsupervised: finding structure without label**

| Method | Where | Role |
| --- | --- | --- |
| Principal component analysis (PCA) | 03, 14, 16 | Compress 2,000 gene to a few axis: neighbour graph, MOFA input, the 3 "program" of a forecast |
| k-nearest-neighbour graph, Leiden clustering, UMAP | 03 | Cluster and draw the cell, as a check on the provided label |
| Scrublet | 03 | Doublet detection per donor |
| Two-component Gaussian mixture | 05b (`ambient_cutoff`) | Set the ambient cutoff per tissue |
| MOFA+ (multi-omics factor analysis, Bayesian; Argelaguet et al. 2020) | 14 | Latent factor shared across 17 organ of the same animal |

**Supervised classification: young against old**

| Method | Where | Role |
| --- | --- | --- |
| L2 logistic regression ($C = 0.1$) | 05a, 05b, 07, 10, 11 | The main model; its weight are the gene list of 05b |
| LightGBM (gradient-boosted tree; Ke et al. 2017) | 05a, 11 | Non-linear comparison |
| Elastic net logistic regression (L1 + L2; Zou and Hastie 2005) | 11 | Sparse model, about 40 gene per cell type |
| Animal-level (pseudobulk) logistic regression | 11 | One training example per animal, the true sample size |
| Attention-based multiple-instance learning (Ilse et al. 2018) | 11 | Neural network: an animal is a bag of cell with one label |
| Logistic regression on technical feature or sex only | 07, 10 | Control: can depth, blood or stress fraction, or sex, alone do the job? |

**Regression: age in month**

| Method | Where | Role |
| --- | --- | --- |
| Ridge regression ($\alpha$ = 10, 100, 1000) | 05c, 10, 14 | The clock; per-organ clock in bulk |
| Ordinal regression (Frank and Hall 2001) | 11 | Uses only the order of the age, not their spacing |
| Kernel ridge (RBF), LightGBM regressor, multilayer perceptron | 11 | Non-linear clock |

**Statistical model with the animal as sample**

| Method | Where | Role |
| --- | --- | --- |
| Pseudobulk DESeq2 (PyDESeq2; Muzellec et al. 2023) | 09 | Differential expression, old against young, per cell type |
| Ordinary least squares with covariate | 07, 09, 14 | Tissue fixed effect; power check; plate- and sex-adjusted age effect in bulk |
| Univariate and multivariate linear model (`decoupler`; Badia-i-Mompel et al. 2022) | 12, 13 | Transcription factor, pathway, hallmark and miRNA activity from target gene |
| Gene-set scoring against expression-matched control gene | 08, 16 | SenMayo; naive and memory T-cell score |
| Spearman, Fisher, Mann–Whitney, Wilcoxon, sign and binomial test; Benjamini–Hochberg | throughout | Trend and comparison, over animal or row |
| Permutation null: between animal, within tissue, or across tissue consistently | 05a, 07, 09, 14 | A floor, and a null that keeps the shared-animal structure |
| Bootstrap over test mice | 07 | 90% interval of the cross-age AUROC |

**Dynamical model: age as time (exploratory)**

| Method | Where | Role |
| --- | --- | --- |
| Saturating-removal ODE and SDE (Karin et al. 2019) | 16 | Does a senescence proxy accelerate; does the spread between animal grow |
| Homeostatic naive-fraction ODE with profile likelihood | 16 | Can the rate of thymic involution be identified |
| Smoothing spline, sigmoid curve fit, linear ODE | 16 | Forecast a held-out age |
| SINDy (STLSQ) | 16 | Forecast only; its equation are not interpreted (too few time point) |
| Neural ODE (`torchdiffeq`) | 16 | Forecast; 5 seed where a claim rests on it |
| Optimal transport (`moscot`; Klein et al. 2025) | 16 | Rebuild a held-out age's cell population |

A model is never ranked on ordinary held-out accuracy alone, because in this
atlas batch inflates it ([section 8](#8-evaluation-design)).


## 8. Evaluation Design

Two unit are kept deliberately separate. The **cell** is the unit of
measurement: there are tens of thousand, they make every figure look
convincing, and they are not independent. The **donor** is the unit of
inference: the lung has sixteen and the whole droplet atlas twenty-three,
because one animal gave up to twelve organ. Every claim is bounded by that
number, and a pattern seen in several tissue is not several independent
observation.

| Rule | Why |
| --- | --- |
| A whole donor goes into train or test, never split | Otherwise the model can identify the animal and inherit its label |
| Cell are capped per donor | Otherwise the largest donor sets the loss on its own |
| A donor is averaged before a group is compared | Otherwise a p-value counts cell as sample and comes out tiny and meaningless |
| Row of one tissue are compared within that tissue | They share mice; a test over row would treat them as independent |

**A high AUROC on unseen donor is necessary and not sufficient.** Each age
group was processed at a different time, so donor identity, batch and age move
together, and a model that reads the batch scores well on a held-out donor of a
*trained* age group. Only a test on an age group the model never trained on can
see the difference.

**Every threshold is fixed before a new tissue is looked at.** The lung was the
development tissue: the 0.7 cross-age cut and the 0.6 shuffled cut were set
while building the pipeline on it. They were then written down in `06_survey`
before any other tissue was modelled:

| Stage | Rule |
| --- | --- |
| Classifiable | ≥ 2 young and ≥ 2 old animal, ≥ 50 cell per class |
| Batch-testable | also ≥ 2 distinct young age and ≥ 2 distinct old age |
| Age signal | all-split cross-age ≥ 0.7 and shuffled ≤ 0.6 |
| Void | shuffled > 0.6 |
| Partial / batch driven | cross-age 0.5–0.7 / below 0.5 |
| Technically confounded | the technical-only model also reaches cross-age ≥ 0.7; reported, not deleted |
| Gene read | AUROC ≥ 0.8, a passing batch check, and the `own` filter |

**Which check were planned and which were added.** The survey rule, the verdict
rule, the pseudobulk, the replication and the bulk test were planned before
their data were seen. The technical-only control was added after bladder looked
uniformly perfect; the dominant-donor rerun after the outlier flag; and the
within-tissue permutation, the interval over test mice, the sex control and the
technical transfer test of section 11 were added during the critical review. An
added check can only remove or downgrade a claim here, never create one.


## 9. File Map

| File | Job |
| --- | --- |
| `src/aging.py` | Every shared function, every default and every per-tissue setting. The one file to edit to change behaviour |
| `src/models.py` | The classifier and clock of 11, all with one `fit` / `predict` interface |
| `src/dynamics.py` | The ODE, SDE and forecasting model of 16 |
| `src/download.py` | Fetch the per-tissue file and the external table from figshare through its API, md5 checked; `--list` shows the size first |
| `src/run_tissue.py` | Run notebook `02`–`05c` top to bottom for one or more tissue, saving each executed copy to `notebook/run/{Tissue}/` |
| `src/strip_output.py` | Git clean filter: commit a notebook as code and markdown, without any output |
| `notebook/01_scanpy_tutorial.ipynb` | The scanpy preprocessing and clustering tutorial, kept as a record of the learning step |
| `notebook/02_explore.ipynb` | Age, donor, sex balance and cell type per donor; the provided UMAP |
| `notebook/03_preprocess.ipynb` | Rebuild from count, QC, doublet, normalize, variable gene, PCA, Leiden, UMAP |
| `notebook/04_composition.ipynb` | Donor that dominate a cell type, then cell type fraction per donor against age |
| `notebook/04b_lung_outlier.ipynb` | Lung only: is donor 30-M-2 sick? The investigation behind its exclusion |
| `notebook/05a_classifier.ipynb` | Young against old per cell type, every control, and the batch check |
| `notebook/05b_genes.ipynb` | The gene behind each model, filtered and flagged; development against aging |
| `notebook/05c_clock.ipynb` | Ridge clock in month, three strength, and a held-out age group |
| `notebook/06_survey.ipynb` | Across every downloaded tissue: donor per age, sex balance, which cell type are classifiable and batch-testable |
| `notebook/07_atlas.ipynb` | Every tissue's result in one table: verdict, compartment, technical-only, sex-only, male-only and dominant-donor check, interval over test mice, structural against immune within tissue |
| `notebook/08_senescence.ipynb` | SenMayo, p21 and p16 per animal in every row, trend with age and its stability |
| `notebook/09_pseudobulk.ipynb` | DESeq2 per row with the animal as sample, agreement between row against a shared-mouse null, power check, replication of the published global aging gene |
| `notebook/10_replication.ipynb` | Train on this atlas, test on Kimmel et al. 2019; technical-only transfer; key gene per mouse; is the p21 rise handling stress |
| `notebook/11_models.ipynb` | Five classifier and five clock compared on unseen age group and on Kimmel transfer only |
| `notebook/12_mechanism.ipynb` | Transcription factor, pathway and hallmark activity, with a circularity and an ambient check; ligand–receptor change |
| `notebook/13_mirna.ipynb` | miRNA activity from target gene; check against primary miRNA in bulk |
| `notebook/14_multiomics.ipynb` | Bulk time course (10 age, 17 organ): the single-cell finding re-tested, p21 timing, MOFA+ factor, organ age residual |
| `notebook/16_dynamics.ipynb` | Accumulation ODE, heterogeneity, T-cell ODE and its identifiability, forecasting a held-out age (ODE, SINDy, neural ODE, optimal transport) |

Notebook `02`–`05c` run per tissue; `06` onward work across tissue. (Number 15
was a planned transposable-element and RNA-editing module, dropped.)

Inside `src/aging.py`:

| Group | Function |
| --- | --- |
| Setting | `AGE_MONTHS`, `YOUNG`, `OLD`, `ADULT`, `AMBIENT`, `TECHNICAL_PATTERN`, `SEX_PATTERN`, `DISSOCIATION`, `CIRCADIAN` |
| Per tissue | `TISSUE_SETTINGS`, `settings`, `paths`, `main_kw`, `ages_in`, `flag_outlier_mice`, `slug`, `COMPARTMENT`, `compartment` |
| Gene set | `GENE_ALIASES`, `load_gene_set` |
| Gene filter | `celltype_means`, `ambient_cutoff`, `foreign_genes`, `feature_mask` |
| Data prep | `prep`, `usable_cell_types` |
| Classifier | `logreg`, `age_eval`, `age_auc`, `cross_age_splits`, `cross_age_scores`, `cross_age_auc`, `bootstrap_cross_age`, `top_genes` |
| Control | `technical_adata`, `sex_adata` |
| Pseudobulk | `pseudobulk`, `deseq`, `log_cpm_change` |
| Clock | `clock_eval`, `clock_summary`, `held_out_age` |
| Trajectory | `trajectories`, `label_patterns`, `pattern_top` |

| Folder | Holds |
| --- | --- |
| `data/{Tissue}/` | The official file, `raw.h5ad` and `processed.h5ad` of one tissue. Not tracked by git |
| `data/external/` | Table and dataset from other study, for comparison and replication. Not tracked by git |
| `figure/{Tissue}/`, `figure/atlas/` | Every plot and result table, per tissue and across tissue. Regenerated by a run, not tracked by git |
| `notebook/` | The analysis, one notebook per phase |
| `notebook/run/` | Executed copy written by `run_tissue.py`. Not tracked by git |
| `src/` | The shared code |


## 10. Result: Mouse Aging Atlas

The question across the atlas: in which cell type, in which organ, does gene
expression carry an age signal that survives the batch check, and what part of
it holds outside this atlas?

| Section | Status | One-line answer |
| --- | --- | --- |
| [10.1](#101-which-tissue-can-be-tested) Survey | confirmatory | 9 of 16 tissue testable; 23 mice in all |
| [10.2](#102-age-signal-across-tissue) Atlas | confirmatory | structural row score higher within tissue; the pass-rate gap is not significant within tissue |
| [10.3](#103-check-behind-the-atlas) Check | confirmatory | technical batch in bladder; sex never fakes a pass but hides immune signal; every pass has an interval above 0.5 |
| [10.4](#104-senescence-across-tissue) Senescence | confirmatory | SenMayo flat; p21 rises inside this atlas only |
| [10.5](#105-gene-with-the-animal-as-sample) Pseudobulk | confirmatory | a small shared core: MHC class I, inflammation |
| [10.6](#106-replication-in-another-lab) Another lab | confirmatory | weak direction transfers; complement replicates, p21 does not |
| [10.7](#107-a-bulk-time-course) Bulk time course | confirmatory | collagen down, MHC-I and complement up; p21 flat |
| [10.8](#108-which-model-generalizes) Model | exploratory | nothing beats logistic regression; ordinal clock helps |
| [10.9](#109-which-regulator-move-with-the-shared-signal) Regulator | exploratory | co-varies with interferon / NF-κB activity; not circular, not ambient |
| [10.10](#1010-microrna) miRNA | exploratory | miR-29 target set falls in old stroma |
| [10.11](#1011-aging-as-a-dynamical-system) Dynamics | exploratory | steady drift; most model question not answerable with this data |
| [10.12](#1012-lung-age-signal-per-cell-type)–[10.16](#1016-lung-what-the-audit-found-first) Lung | worked example | every step on one tissue |

### 10.1 Which Tissue Can Be Tested

**Background.** The cross-age check needs at least two young age and two old
age in a tissue, each with enough cell, so whether a tissue can be tested at
all is a property of the study design. `06_survey` decided this from metadata
alone, before any tissue was modelled.

| Tissue | Young / old mice | Batch-testable cell type | Read as |
| --- | --- | --- | --- |
| Lung | 6 / 10 | 14 | tested |
| Marrow | 4 / 13 | 14 | tested |
| Kidney | 5 / 11 | 12 | tested |
| Limb muscle | 4 / 12 | 8 | tested |
| Heart and aorta | 2 / 9 | 6 | tested, thin: one 1m and one 3m animal |
| Spleen | 3 / 10 | 6 | tested |
| Bladder | 4 / 4 | 4 | tested, technically confounded (10.3) |
| Liver | 5 / 7 | 4 | tested |
| Tongue | 5 / 4 | 3 | tested |
| Mammary gland | 2 / 4 | 0 | not testable: every young animal is 3m |
| Thymus | 1 / 5 | 0 | not testable: one young animal |
| Trachea | 2 / 0 | 0 | not testable: no old animal |
| Fat, pancreas, skin, large intestine | 0 / 3–8 | 0 | not testable: no young animal |

**The whole droplet atlas is 23 animal.** The median animal gave 7 tissue and
one gave 12, so tissue are not independent replicate: a quirk of one animal can
reach every organ it gave.

### 10.2 Age Signal Across Tissue

**Background.** If structural cell (fibroblast, endothelium, epithelium) and
immune cell aged equally visibly, they should pass the batch check equally
often. The lung suggested they do not. Here the same test is applied to 70
tissue × cell type row, with the threshold of section 8.

Voided and not-testable row are left out; the right half also drops the 8 row
that technical feature alone can classify (10.3):

| Compartment | Row | Age signal | Rate | Without technical row | Median cross-age |
| --- | --- | --- | --- | --- | --- |
| Stromal | 8 | 6 | 0.75 | 4 of 6, 0.67 | 0.84 |
| Epithelial | 12 | 9 | 0.75 | 8 of 10, 0.80 | 0.77 |
| Endothelial | 6 | 3 | 0.50 | 3 of 4, 0.75 | 0.68 |
| Immune | 28 | 9 | 0.32 | 9 of 27, 0.33 | 0.63 |
| Progenitor | 8 | 2 | 0.25 | 2 of 7, 0.29 | 0.54 |

The first version of this comparison used Fisher and Mann–Whitney test over
row (structural 0.75 vs immune 0.33, p = 0.008 and 0.0009). Those test treat
every row as independent, but row of one tissue share mice, and marrow and
spleen contribute **only** immune row: a tissue where everything scores low for
a tissue-wide reason then lowers the immune average without saying anything
about immune cell. The stricter version compares row **within** a tissue
(`07` section 10):

| Test (without technical row) | Structural − immune | p | Read as |
| --- | --- | --- | --- |
| Row as independent sample (Fisher; Mann–Whitney) | pass rate +0.42; cross-age higher | 0.008; 0.0009 | **superseded**: overstates independence |
| Within-tissue permutation, pass rate | +0.42 | 0.13 | **not significant** |
| Within-tissue permutation, mean cross-age AUROC | +0.14 | 0.001 | **holds** |
| Tissue fixed effect, adjusted for cell and mouse count | +0.14 | 0.0009 | holds; cell count (p 0.70) and mouse count (p 0.36) do not explain it |
| Tissue with both kind where structural is higher | 5 of 7 | sign test 0.23 | consistent in direction, not decisive |
| Same comparison, sex-chromosome gene removed (10.3) | pass rate 0.80 vs 0.59; cross-age +0.11 | 0.23; 0.002 | the continuous difference holds; the pass-rate gap nearly closes |

**Read as:** inside an organ, structural cell keep a more transferable age
signal than the immune cell of the same mice, by about 0.1 AUROC. The binary
"pass rate" contrast, the headline of the first draft, is not supported once
tissue is accounted for: the lung (structural 0.85 against immune 0.54) carries
most of it, and outside the lung much of the immune shortfall was sex (10.3).

- **Structural row that pass:** lung fibroblast (0.95) and adventitial cell
  (0.94), limb muscle mesenchymal stem cell (0.90), six of eight kidney
  epithelial segment (0.73–0.86), kidney mesangial cell (0.77), kidney and
  muscle endothelium (0.88, 0.79).
- **Immune aging is organ dependent, not batch-only.** Lung and marrow immune
  cell mostly sit at 0.38–0.62, but limb muscle B (0.83), T (0.86) and
  macrophage (0.76) pass, as do spleen T cell. Without sex-chromosome gene,
  marrow plasma cell, monocyte, macrophage and naive T cell also reach 0.71–0.77
  (10.3); the lung immune row stay where they are.
- **Caveat.** Within-tissue cross-age falls as the number of mice rises (ρ
  −0.36), because the large tissue (lung, marrow) are the ones with many immune
  row; the fixed-effect model shows mouse count itself does not explain the
  compartment gap.

### 10.3 Check Behind The Atlas

**Background.** The cross-age split only breaks a batch that differs *between*
age group. Three thing could still pass it: a processing difference shared by
all young against all old mice, sex (tied to age group in every tissue), and a
single animal dominating a cell type. Each is tested per row; the uncertainty of
each score is shown as a 90% interval over test mice.

**Technical batch.** Every row rerun with only five technical feature per cell
(log count, log gene, ribosomal, dissociation and hemoglobin fraction):

| Row | Cross-age, expression | Cross-age, technical only | Read as |
| --- | --- | --- | --- |
| Bladder urothelial cell | 0.96 | 0.86 | technically confounded |
| Bladder cell (mesenchyme) | 0.96 | 0.80 | technically confounded |
| Heart fibroblast | 0.77 | 0.74 | technically confounded |
| Lung fibroblast | 0.95 | 0.60 | clean |
| Limb muscle mesenchymal stem cell | 0.90 | 0.68 | clean, but close |

26 of the 29 age-signal row survive it. An old cell can differ technically for a
real reason, such as RNA content, so a confounded row is a warning rather than a
disproof; it is reported, not deleted.

**Sex.** The main model sees every variable gene, sex-chromosome gene
included, while sex is tied to age group. Three test per row (`07` §11–12):

| Test | Result | Read as |
| --- | --- | --- |
| Sex as the only feature, cross-age | never ≥ 0.7 in any row (highest 0.68, lung fibroblast) | sex alone cannot pass the check |
| Main check without the six sex-chromosome gene | all 29 age-signal row stay ≥ 0.7 (median change 0.03) | no pass depends on *Xist* or a Y gene |
| Male mice only, where two young and two old age remain | 17 of 20 testable age-signal row stay ≥ 0.7; kidney brush cell 0.68, collecting duct 0.59, spleen T 0.69 drop | most pass survive without female mice |
| Verdict that change without sex gene | **12 row rise** to ≥ 0.7, none fall: 11 immune or progenitor row (marrow 7, spleen 2, liver Kupffer cell, heart leukocyte; one spleen row is void by its shuffled control) and lung bronchial smooth muscle | sex gene **hid** signal: in some split the training age pair male young with female old and the test age the reverse, inverting a model that leaned on *Xist* |

The last row is the most consequential finding of the review. Without the sex
gene the immune pass rate rises from 0.33 to 0.59 (structural 0.75 → 0.80), so
outside the lung a large part of the "immune fails the batch check" pattern was
sex misalignment rather than batch. The lung immune row do not change. The
verdict table and the downstream module (08–16) keep the original verdict,
fixed before this check; removing sex-chromosome gene by default is the first
item of section 11.4.

**Uncertainty.** All 29 age-signal row have a 90% interval (test mice resampled)
entirely above 0.5, and 23 entirely above 0.7. The six that reach below 0.7 sit
near the cutoff: lung alveolar macrophage (0.74, interval 0.69–0.78), kidney
macrophage, brush cell, collecting duct, loop of Henle ascending limb, and marrow
erythroblast. The median interval width is 0.05; training mice are held fixed,
so the true uncertainty is larger.

**Animal that dominate a cell type.** A cell type that takes over half its cell
from one animal was rerun without that animal. Two verdict move, neither of
them a main-table pass: liver sinusoidal endothelium goes from partial (0.57) to
an age signal (0.82) without 1-M-63, and heart coronary-artery endothelium from
partial (0.53) to batch driven (0.49) without 30-M-2. **Donor 30-M-2**, excluded
from the lung myeloid cell type, also dominates cell type in heart and marrow,
which fits a body-wide condition in that animal.

### 10.4 Senescence Across Tissue

**Background.** If senescent cell accumulate with age (section 2.4), a
senescence signature should rise per animal in cell type that carry a real age
signal. Each cell is scored with SenMayo (115 of its 117 mouse gene present),
averaged per animal, and correlated with age in month (`08_senescence`). A
trend is `stable` when it keeps its sign with any one age group left out.

| Measure | Robust structural row, ρ > 0 | BH < 0.1 | Median ρ | Read as |
| --- | --- | --- | --- | --- |
| SenMayo score | 8 of 15 | 0 | 0.00 | no age trend |
| SenMayo-high fraction | 9 of 15 | 0 | 0.20 | weak, not significant |
| ***Cdkn1a*** (p21) | **14 of 15**, all 5 tissue | 5 | **0.45** | rises inside this atlas; does not replicate (10.6, 10.7) |
| *Cdkn2a* (p16) | — | — | — | detected in 0.8% of cell: unusable per cell |

The composite score does not move: most SenMayo gene are secreted factor,
weakly captured per cell, and averaging 115 of them dilutes the few that do.
p21 alone rises in structural cell of every tissue that has one: bladder
endothelium ρ 0.95 (BH 0.03), lung fibroblast 0.78 (0.04), lung adventitial 0.77
(0.09), kidney collecting duct 0.73 (0.05), tongue basal epithelium 0.72,
muscle endothelium 0.66. The kidney gives seven row from the same animal, so the
honest count is five tissue of five, not fourteen result.

**Caveat.** p21 also rises after any acute stress (section 2.4), so even inside
the atlas a p21 rise is not by itself senescence.

### 10.5 Gene With The Animal As Sample

**Background.** The classifier's weight say which gene it used, not which gene
changed significantly. Pseudobulk DESeq2 (section 2.5) tests each gene with the
animal as sample, on 55 row (`09_pseudobulk`), and asks two further question:
do different cell type share their change, and does the atlas's own published
list of "global aging gene" (Zhang et al. 2021, tested on cell) hold with the
animal as sample?

| Question | Result | Read as |
| --- | --- | --- |
| Gene changed per row (padj < 0.1, own gene) | robust median 175, other 58; same median mouse count (13); ×4.9 after adjusting for mice (p = 0.0006) | not a power artifact; but "robust" was defined by an age difference, so this is consistency, not new evidence |
| Row pair agreeing beyond a shared-mouse null | 116 of 1,485 (7.8%; 5% by chance) | **little shared signal is detectable** at 5–17 mice per row |
| Correlation of effect, between tissue | median 0.06 (null 95th percentile 0.22) | no general program visible |
| Correlation of effect, within one tissue | median 0.28 (null 0.39) | what shared animal produce on their own |
| Published "global aging gene", up half | 60 of 93 replicate (p = 0.007) | **holds**: MHC class I and II lead |
| Published "global aging gene", down half | 68 of 190 replicate (p = 0.0001, below chance) | **does not hold**; ribosomal gene reverse |

The null keeps the shared-mouse structure: the 23 animal's old/young label are
permuted once and applied to every tissue alike, so two row can only beat it by
agreeing about age, not about which animal they came from.

The shared core across four or five tissue is small and coherent: MHC class I
(*H2-D1*, *H2-K1*, *B2m*) and inflammation (*Tnf*, *Casp1*, *Casp4*, *Gbp3*) up.
*H2-D1* and *H2-K1* rise in all 22 robust row, and they are also the top
replicating gene of the published list. Of that list, 47 of 330 gene have no
direction even in its own data and are left out. *Tsc22d3* falls in five tissue
but is glucocorticoid responsive, so it is flagged as possible time-of-day
signal. *Cdkn2a* (p16), invisible per cell, rises in 6 row across 4 tissue once
count are summed per animal (log2FC 1.8–5.2, low count).

**Caveat.** A young/old pseudobulk contrast in TMS is still a contrast between
batch. A gene list from here counts as aging only for the gene that also hold in
10.6 or 10.7.

### 10.6 Replication In Another Lab

**Background.** Kimmel et al. 2019 profiled lung, kidney and spleen of young
(7–8 month) and old (about 21 month) male mice in another lab. A model trained
on TMS cannot know Kimmel's batch, so ranking Kimmel's old cell above its young
one is evidence of a shared age signal. Each TMS row was matched to a Kimmel cell
type (19 row, label checked by marker: TMS lung fibroblast =
*Npnt*/*Inmt*/*Tcf21* stromal cell, adventitial = *Dcn*/*Pi16*), trained on TMS
and tested on Kimmel unchanged (`10_replication`).

| Question | Result | Read as |
| --- | --- | --- |
| Does the direction of aging transfer? | 19 of 19 row score old above young; cell AUROC 0.51–0.76; per-mouse 1.00 in 13 | a weak age signal is shared across lab |
| Could that transfer be technical? | a model on five technical feature only (depth, gene count, ribosomal, dissociation, hemoglobin fraction) also ranks old above young in 12 of 19 row (median 0.54). Expression beats it by > 0.05 in 10: alveolar macrophage (0.72 vs 0.44), kidney endothelium (0.64 vs 0.38), collecting duct (0.63 vs 0.37), non-classical and classical monocyte, proximal tubule, kidney macrophage, type II pneumocyte, spleen B, adventitial cell. In the other 9 (lung B, NK, CD4 and CD8 T; spleen NK, T, macrophage; kidney mesangial; lung fibroblast) technical feature do as well or better | **partly**: myeloid and structural transfer is not explained by technical feature; lymphocyte transfer cannot be told from a shared RNA-content or handling difference |
| Does the TMS clock transfer? | old placed older in 19 of 19, compressed (0.7–5 month vs ~14) | yes, in order, not in scale |
| Does the TMS batch check predict transfer? | ρ = 0.17 (p = 0.48); robust 0.64 vs other 0.59 | no detectable relation (19 row, low power) |
| Lung fibroblast, the strongest TMS row | 0.51 (0.52 without the pooled 3m sample) | **does not transfer** |
| Complement (*C3*, *Cfb*) | *C3* +0.61 in every old adventitial cell; *Cfb* +0.43 in every old mesangial cell | **replicates** |
| p21 (*Cdkn1a*) in structural cell | lung fibroblast −0.03, adventitial −0.46, AT2 −0.21 | **does not replicate** |

**p21 is therefore not called an aging marker.** Inside TMS its rise survives
removing each animal's dissociation-stress level (lung adventitial ρ 0.86,
p = 0.003), so it is not ordinary dissociation stress. In Kimmel, p21 instead
rises in old lung immune cell (classical monocyte +0.34, CD4 T +0.30), the
reverse of the TMS pattern.

**What the batch check is, after this.** It removes TMS-specific artifact (the
technical and dominant-animal case of 10.3 are real), but it does not rank which
signal travel: myeloid cell it fails (lung monocyte, alveolar macrophage, kidney
macrophage) carry a weak age signal that another lab reproduces and that
technical feature do not explain. It is a filter, not a measure of portability.

**Caveat.** The 19 row are not 19 independent replication: they come from 3
tissue of the same 6–7 Kimmel mice, and the Kimmel paper does not say whether
its young and old mice shared processing day. A technical difference can itself
be biological (an old cell may hold less RNA), so a row whose transfer technical
feature match is unresolved, not disproved.

### 10.7 A Bulk Time Course

**Background.** The same consortium sequenced bulk RNA from 17 organ of 56 mice
at **10 age** (1 to 27 month), with sequencing plate balanced across age
(Schaum et al. 2020). Bulk cannot resolve cell type, but it has what the
single-cell atlas lacks: age not tied to plate, and the 3–18 month window
filled (`14_multiomics`).

| Gene (age effect, plate and sex adjusted) | Organ up | Organ down | Read as |
| --- | --- | --- | --- |
| Collagen *Col1a1* / *Col3a1* | 0 | **13 / 14 of 17** | **matrix loss is the most universal aging change** |
| MHC class I (*H2-K1*, *H2-D1*, *B2m*) | 6–7 | 2 | replicates in most organ |
| Complement *C4b* / *C3* | 9 / 5 | 1 / 2 | replicates |
| p21 (*Cdkn1a*) | 1 | 5 | **does not rise**, neither from 1 to 9 month (median ρ −0.08) nor from 9 to 27 (−0.13) |
| p16 (*Cdkn2a*) | 3 | 0 | weak |

**p21 does not replicate in bulk either.** A bulk sample is a whole organ, so a
rise confined to fibroblast could be diluted (lung bulk p21 falls slightly, t =
−1.9). Together with Kimmel's sorted cell type the reading is: **not
replicated, cause unknown**, rather than a proven artifact.

**Collagen.** Bulk cannot say whether fibroblast make less collagen or there are
fewer fibroblast; the single-cell lung stroma shows the per-cell fall (10.13).

**A body-wide aging axis.** MOFA+ over all 17 organ finds one factor that tracks
age at ρ = 0.95 over mice and explains variation in 15 of them, most in fat
depot (16–18%), muscle (14%) and lung (9%). **Organ of one mouse look old
together:** per-organ clock reach 2.1–4.7 month of error, and after removing
regression to the mean, plate and sex, the residual of 93% of organ pair
correlate positively across mice (median 0.23; shuffled-mouse null 95th
percentile 0.03), in male (0.22) and female (0.33) alike. This is the
"ageotype" pattern, with one caveat that cannot be resolved here: every organ of
a mouse was dissected together, so a per-mouse handling effect would look the
same.

### 10.8 Which Model Generalizes

**Background.** A natural objection to section 10.2: a linear model may simply
be too weak, and a non-linear or neural model might find the aging signal the
linear one misses. Five classifier and five clock (`11_models`) were compared,
scored **only** on unseen age group and on the other lab, over the 19 row
matched to Kimmel. Exploratory: the comparison was not planned before 10.2 was
seen.

| Classifier | Unseen age group | Row below logistic | Other lab | Row below logistic | Read as |
| --- | --- | --- | --- | --- | --- |
| **L2 logistic regression** | **0.747** | — | 0.618 | — | the baseline holds |
| Elastic net (~40 gene) | 0.720 | 15 of 19 | 0.616 | 8 of 19 | sparsity does not help |
| LightGBM | 0.705 | 16 of 19 | 0.619 | 8 of 19 | non-linearity does not help |
| Attention MIL (neural network) | 0.661 | 14 of 19 | 0.597 | 11 of 19 | the most flexible model, among the worst |
| Animal-level (pseudobulk) | 0.657 | 16 of 19 | 0.568 | 15 of 19 | too few training example |

| Clock (52 row, every tissue) | Error, month | Lower error than ridge | Held-out age placed between its neighbour |
| --- | --- | --- | --- |
| Ridge | 6.0 | — | 0.51 |
| **Ordinal regression** | **5.3** | **45 of 52 row; majority in 8 of 9 tissue** | 0.49 |
| LightGBM | 5.8 | 26 of 52 | 0.54 |
| Kernel ridge (RBF), multilayer perceptron | 6.7, 6.8 | 16, 15 of 52 | 0.53 |

The row count are the result; a Wilcoxon test over row would treat row of one
tissue as independent. With 10 to 16 animal per row and age tied to batch,
flexibility buys nothing that generalizes. Ordinal regression, which uses only
the **order** of the age, is the one improvement to the clock (it loses in
both liver row), but no model places a held-out age group between its neighbour
more often than a coin: the inversion of [section 10.15](#1015-lung-the-clock)
is the data, not the model.

### 10.9 Which Regulator Move With The Shared Signal

**Background.** A transcription factor's activity can be estimated from whether
its known target gene move together (section 2.4; CollecTRI network,
Müller-Dott et al. 2023), and likewise for signalling pathway (PROGENy; Schubert
et al. 2018) and hallmark gene set (Liberzon et al. 2015). Run on the
animal-level statistic of 10.5 (`12_mechanism`). Exploratory, and correlational:
a regulator whose target rise is **associated** with the change, not shown to
cause it.

| Layer | Up in robust row of 6–7 tissue | Down | Read as |
| --- | --- | --- | --- |
| Transcription factor | NF-κB (*Nfkb1* 18 row, *Rela* 17), interferon (*Irf1* 18, *Stat1* 16), *Stat3*, *Ciita* | nothing comparable (max 8–9 row) | an inflammatory program |
| Pathway | JAK-STAT (t 4.0), NF-κB, hypoxia | TGF-β (8 row) | |
| Hallmark | interferon-α, IL-6/JAK/STAT3, inflammatory response, complement, p53, apoptosis, mTORC1 (t 3.0–4.1) | epithelial–mesenchymal transition, oxidative phosphorylation | inflammaging up, matrix and mitochondrial gene down |
| MHC class I statistic vs regulator activity, across row | *Irf1* ρ 0.61, *Stat1* 0.58, *Rela* 0.54 | | co-varies with interferon and NF-κB activity |
| Same, MHC gene removed from every target set | 0.60, 0.56, 0.52 (*Ciita*, 14 MHC target: 0.50 → 0.36) | | **not circular** |
| Regulator recomputed without each row's foreign (ambient) gene | *Irf1* 18 → 17 row, *Stat1* 16 → 17, *Rela* 17 → 15, *Nfkb1* 18 → 13, all still 6 tissue; *Ciita* 17 → 9 | | the interferon / NF-κB rise is the cell's own; the MHC class II driver is half ambient immune RNA |

The glucocorticoid receptor does not change, so the fall of *Tsc22d3* is not a
stress-hormone effect through it. AP-1 (*Jun*, *Fos*, *Egr1*) does rise, and
AP-1 is also the dissociation-stress response, so part of the NF-κB/AP-1 signal
may be how old tissue reacts to being dissociated. A ligand–receptor scan
(LIANA resource; Dimitrov et al. 2022) adds little beyond 10.5: its score is
built from the same ligand change (*C3*, *H2-K1*, *B2m* rising in stroma), so it
restates that old stroma display more immune-recognition molecule rather than
showing new communication.

### 10.10 MicroRNA

**Background.** Mature miRNA are not captured by droplet sequencing (section
2.2), so activity is inferred from their target falling together (section 2.4;
`13_mirna`; 11,333 experimentally supported pair from human, mapped to mouse
ortholog). Exploratory, and an inference: the miRNA itself is never measured.

| miRNA | Robust row, target falling together | Without matrix gene, stromal row | Read as |
| --- | --- | --- | --- |
| **miR-29a/b/c** | 10–13 row, 5–6 tissue; stromal activity 6.3–7.7 | **2.8–4.5** | target-set signal **consistent with** more miR-29 activity in aging stroma; independent of NF-κB (ρ ≈ 0) |
| miR-34a, miR-146a, miR-155 | slightly less active | unchanged | not readable: their target are inflammatory gene that NF-κB drives up |

miR-29 represses collagen and is reported to rise with age (Ugalde et al. 2011),
which fits the matrix loss of 10.7, and the signal survives removing the matrix
gene. The check against measured precursor fails for a technical reason: the
bulk library captures polyadenylated RNA, and *Mir29a* has any read in only 27%
of sample. A measured small-RNA dataset of the same cohort (Wagner et al. 2024)
is the open check.

### 10.11 Aging As A Dynamical System

**Background.** Section 2.6 introduces the model. Here age is treated as time,
and model are fitted to the bulk time course (10 age) and to single cell (4–6
age) in `16_dynamics`. Exploratory: every model has few time point per
parameter, so each row is reported as what the data **can and cannot**
distinguish.

| Question | Model | Result | Read as |
| --- | --- | --- | --- |
| Does a senescence proxy accelerate? | Karin saturating removal, $\dot{X} = \eta t - \beta X/(\kappa + X)$, vs linear and exponential | a straight line wins for p16 in 10 of 11 organ; Karin never (2.4–6.2 AIC behind; 14.7 for the MOFA aging factor) | **not testable here**: bulk p16 sits near its detection floor |
| Do animal drift apart? | spread between mice at each age | grows in 8 of 12 series (median ρ 0.16) | consistent with the stochastic model, but almost any model predicts it |
| Do cell of one type become more heterogeneous? | per-mouse distance to centroid, depth-matched | rises in 16 of 26 robust row (10 stable to dropping any age) | aged population are more varied; batch can contribute |
| Early-warning signal of a tipping point? | leading eigenvalue share, gene–gene correlation | no direction (6–8 stable rise, 6–8 stable fall) | no sign; this cross-sectional design is a weak test of it |
| Can thymic involution be measured? | homeostatic naive-fraction ODE, $\dot{f} = \theta e^{-\delta t} - k f$, profile likelihood | involution rate 100-fold apart fit almost equally (profile varies 1–13% over δ = 0.003–3 per month) | **not identifiable**; no rate is reported |
| Forecast an age inside the range (8 held-out age × 17 organ) | nearest age, linear, spline, sigmoid, linear ODE, SINDy, neural ODE | no model beats the nearest age: sigmoid wins 74 of 136 (sign p 0.35), spline 61; neural ODE (57), linear ODE (33) and SINDy (20 of 86 that ran, 50 diverged) are worse (BH Wilcoxon ≤ 0.001) | the spread between mice of one age exceeds any model difference |
| Forecast beyond the range (27 from 1–24 month) | same; neural ODE from 5 seed | **no model beats carrying 24 month forward**: neural ODE wins 8 of 17 (sign p 1.0; median skill −0.27, −0.37 to 0.29 by seed), sigmoid 5, spline 3, linear 1; SINDy diverges in 10 of 17 | the earlier single-seed neural ODE "win" (0.21) is withdrawn |
| Rebuild a held-out age's cell population | optimal transport (`moscot`) | for muscle endothelium and mesenchymal stem cell the interpolation beats random and lands below the distance between two mice of that age; for 5 immune and satellite row random does as well | 2 of 7 row; sex differs between the age involved, so part of every distance is sex |

Figure (`figure/atlas/dynamics_*.png`): the rate balance of the Karin model and
the loss of its steady state (`karin_phase`); p16 with the three fit
(`p16_fits`); a phase portrait in which the gene and axis are chosen on one half
of the mice and only the other half is drawn, so that its ordering by age is not
built in; the drawn half still orders by age (|ρ| of the first program: kidney
0.92, limb muscle 0.83, mesenteric fat 0.79, lung 0.58) (`phase_portrait`); the T-cell fit beside three involution rate that fit
equally well (`tcell`); the forecast of one held-out age (`forecast_example`);
and the optimal-transport distance (`ot`). An earlier figure printing the
equation SINDy "discovered" for kidney was removed: 10–30 coefficient fitted to
about nine time point are not interpretable.

### 10.12 Lung: Age Signal Per Cell Type

The lung, on the 16 animal and 24,490 cell of [section 5](#5-input), as a
worked example of every step. Two cross-age column: the hand-picked pair of
split the lung was first tested with, and the all-split rule the atlas uses.

| Cell type | AUROC, unseen mice | Cross-age, 2 split | Cross-age, all split | Shuffled | Verdict (all split) |
| --- | --- | --- | --- | --- | --- |
| Fibroblast of lung | 0.98 | **0.97** | **0.95** | 0.60 | Age signal |
| Adventitial cell | 0.97 | **0.95** | **0.94** | 0.53 | Age signal |
| Type II pneumocyte | 0.97 | **0.80** | **0.83** | 0.35 | Age signal, on 125 cell |
| Alveolar macrophage | 0.95 | 0.66 | **0.74** | 0.40 | Age signal under all split only |
| Neutrophil | 0.97 | 0.78 | 0.69 | 0.70 | Void: shuffled control too high |
| Bronchial smooth muscle | 0.92 | 0.66 | 0.68 | 0.47 | Partial |
| Classical monocyte | 0.91 | 0.40 | 0.56 | 0.34 | Partial |
| Intermediate monocyte | 0.88 | 0.40 | 0.55 | 0.50 | Partial |
| CD8 T cell | 0.83 | 0.31 | 0.54 | 0.56 | Partial |
| Myeloid dendritic cell | 0.79 | 0.29 | 0.52 | 0.51 | Partial |
| Non-classical monocyte | 0.76 | 0.29 | 0.51 | 0.38 | Partial |
| CD4 T cell | 0.73 | 0.24 | 0.52 | 0.34 | Partial, technically confounded |
| NK cell | 0.68 | 0.24 | 0.46 | 0.47 | Batch driven |
| B cell | 0.48 | 0.15 | 0.40 | 0.49 | Batch driven |

**The stromal signal holds under every control.** Fibroblast and adventitial
cell keep an AUROC of 0.96 and above when the model is restricted to male animal
and when the ambient gene are dropped (change ≤ 0.02), and logistic regression
matches or beats LightGBM for every cell type in the table. Young fibroblast
come largely from the two pooled 3m sample, and the signal survives when those
are the test set.

**The immune result is a batch detection on one split and chance on the rest.**
On the two hand-picked split the immune cross-age AUROC sits far below 0.5
(0.15–0.40): what the model learned is inverted on that particular new age
group, a positive detection of a batch difference. Averaged over all twelve
split it sits at 0.40–0.56, which is chance: whatever separates the animal does
not carry to an unseen batch. Alveolar macrophage is the exception, passing
under the all-split rule only (0.74), and stays open. Several immune type still
rank the animal perfectly (`auroc_mouse` = 1.00 for alveolar macrophage, CD8 T,
intermediate monocyte), the clearest demonstration in the project that perfect
donor separation is not evidence of aging when age and batch move together.

Neutrophil are voided by their own control: AUROC 0.97 with a shuffled score of
0.70 means this pipeline can reach 0.70 on that cell type from noise alone.

### 10.13 Lung: Aging Gene In The Cell Type That Pass

Read only from a filtered model with AUROC ≥ 0.8 and a passing batch check.
The weight are interpretation, not a test; section 10.5 tests them.

| Cell type | Up in old | Up in young |
| --- | --- | --- |
| Fibroblast of lung | ***Cdkn1a*** (p21), *C3*, *Cfb*, *Fmo2*, *Gsta3*, *Nrp1*, *Snhg11* | *Col3a1*, *Col1a1*, *Col6a1*, *Sparc*, *Mfap5*, *Mgp* |
| Adventitial cell | *Mt1*, *Bst2*, *Cdkn1c*, *Klf9*, *Efemp1*, *Ifi27l2a* | *Col3a1*, *Col1a1*, *Col1a2*, *Col5a2*, *Sparc*, *Cd248*, *Meg3* |
| Type II pneumocyte | *Nrn1*, *Lrg1*, *Hp*, ***Cdkn1a***, *Nupr1*, *Cyp2b10* | *Igfbp7*, *Sparc*, *Cbr2*, *Serpinh1*, *Egfl6* |

- ***Cdkn1a*** (p21) rises in old fibroblast and old type II pneumocyte, and in
  the trajectory analysis it changes between 3m and old
  ($\Delta_{\text{age}} = +0.54$) and not between 1m and 3m. It does not rise in
  another lab ([10.6](#106-replication-in-another-lab)) nor at any age in bulk
  ([10.7](#107-a-bulk-time-course)), so it is reported as an observation of this
  atlas that did not replicate. The complement gene *C3* and *Cfb* rising
  alongside it **do** replicate.
- **Collagen and matrix gene fall with age** in both stromal cell type
  (*Col1a1*, *Col3a1*, *Col6a1*, *Sparc*, *Mfap5*), per cell, which matches the
  bulk result. In an immune cell type the same gene appear and are ambient RNA
  from the fibroblast-rich 3m sample instead, which is what the `own` filter is
  for.
- **The acute-phase gene in type II pneumocyte** (*Lrg1*, *Hp*) are the
  epithelial side of the same inflammatory shift.

Flagged rather than interpreted: ***Dbp*** and ***Nr1d1***, circadian gene,
appear as "young" gene in several cell type. They are a time-of-day signature,
so the young and old animal were very likely collected at a different hour.

### 10.14 Lung: Cell Composition With Age

Per-donor fraction against age, over the 13 animal left after the outlier sample
are removed, Spearman:

| Cell type | $\rho$ | p | p (BH) | Read as |
| --- | --- | --- | --- | --- |
| CD8 T cell | **+0.70** | 0.008 | 0.061 | Rise with age; suggestive |
| Alveolar macrophage | **−0.65** | 0.016 | 0.062 | Fall with age; suggestive |
| Fibroblast of lung | −0.52 | 0.070 | 0.187 | Possible decline, not conclusive |
| B cell | +0.41 | 0.170 | 0.230 | Possible rise, not conclusive |
| Bronchial smooth muscle | +0.40 | 0.173 | 0.230 | No clear trend |
| CD4 T cell | −0.44 | 0.136 | 0.230 | No clear trend |
| NK cell | −0.30 | 0.318 | 0.363 | No clear trend |
| Classical monocyte | −0.12 | 0.707 | 0.707 | No trend once the outlier donor is out |

Both top result land at p ≈ 0.06 after correction, on 13 animal:
**suggestive, and consistent with known lung immune aging**, not established.
The CD8 gene that change in aging only (*S100a4*, *S100a6*, *Lgals1*, *Nkg7*,
*Gzmm*, *Id2*) are effector and memory gene, the naive-to-memory shift of
section 2.4. So part of the CD8 "age signal" is a change in which CD8 subtype is
present. Dissociation can also change which cell survive, so composition is the
least batch-proof measure here.

### 10.15 Lung: The Clock

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
`05c_clock`; the ordering of cell type is stable across them.

**Where it breaks.** Holding out a whole age group, the clock put the 18m animal
at 15 to 23 month, which is reasonable, and the 21m animal at 10 to 16, younger
than the 18m animal and therefore wrong. Both 21m animal are female and form
their own batch, so the clock is ordering batch rather than resolving age inside
the old range. Leave-one-donor-out error is optimistic for the same reason: the
held-out donor's batch-mates stay in training. The clock works as a coarse
young-to-old estimator on this data and no further.

### 10.16 Lung: What The Audit Found First

Every claim above depends on phase 2 to 4 having found these first:

- **Donor 30-M-2 is a sample of one compartment.** It contributes 8,014 cell, a
  third of the dataset, and 97% of them are myeloid: 92% of all lung macrophage
  and 68% of all classical monocyte in the tissue come from this one animal. Its
  monocyte are in an unusual state (low *Ccr2* and *Ly6c2*, *S100a8* and
  *S100a9* about 12× other 30m animal, some *Mki67* and *Top2a*, **no**
  interferon or hemoglobin signal): an inflammatory, partly dividing, less
  mature myeloid expansion rather than an antiviral response or blood
  contamination. It is excluded from the myeloid cell type, kept everywhere
  else, and reported with `with_excluded_mice` alongside (alveolar macrophage
  0.95 → 0.90 with it).
- **The two pooled 3m sample** (`3-M-5/6`, `3-M-7/8`) are about 50% fibroblast
  with almost no immune cell, so they are left out of the composition trend.
- **Ambient RNA is present and only partly removable.** *S100a8* and *S100a9*
  appear in the B cell of three donor, *Scgb1a1* in many cell type, collagen in
  immune cell. SoupX and CellBender need the unfiltered droplet data, which this
  release does not include, so the filter of
  [section 7.6](#76-step-6-filter-the-gene-then-read-them) is the available
  substitute.
- **Sex is confounded with age.** The 1m and 30m animal are all male and the 21m
  all female. Sex can hide a signal as well as fake one: B cell go from 0.48 to
  0.77 and NK cell from 0.68 to 0.89 on male animal alone.
- **Doublet are rare here** (49 cell, 0.2%).

### 10.17 What This Shows

- **A method result.** On an atlas where age and batch move together,
  mouse-held-out AUROC is not evidence of aging: in the lung 12 of 14 cell type
  reach 0.7 on unseen mice and 4 keep a signal on unseen age group. The
  cross-age check, a donor-permuted control, a technical-only and a sex-only
  control together separate what survives.
- **An atlas result, smaller than first claimed.** Within an organ, structural
  cell keep a more transferable age signal than immune cell of the same mice;
  whether they *pass* more often is not established once tissue is accounted
  for, and outside the lung part of the immune shortfall was sex-chromosome
  gene inverting the model between age group, not batch. Immune aging is
  visible in some organ (limb muscle, spleen, and marrow once sex gene are
  removed) and not in the lung.
- **What replicates.** Across a second lab and a bulk time course with ten age:
  **complement and MHC class I rise**, **collagen falls**. The p21 rise and the
  falling half of the published global aging gene do not hold.
- **What the exploratory part adds.** Ordinal regression is a better clock than
  ridge; the shared signal co-varies with interferon and NF-κB activity; the
  miR-29 target set falls in old stroma. The dynamical model mostly show the
  limit of the data: acceleration, tipping point and thymic rate cannot be
  tested with this many age and animal.
- **The limit.** The obstacle is the study design, not the model class: with
  age tied to batch and sex, and 23 animal across the whole atlas, no amount of
  model capacity recovers what the experiment did not separate.


## 11. Critical Review And Limitation

This section reviews the project the way a thesis committee would: every
objection that a careful reader could raise, what was checked in response, and
what changed. A claim that did not survive is listed, not quietly removed.

### 11.1 Objection And Response

| # | Objection | Check run | Outcome |
| --- | --- | --- | --- |
| 1 | "Mouse-held-out AUROC in TMS is just batch." | Cross-age check on every row; shuffled control | Agreed; it is the project's premise. Only the cross-age score is used for any verdict |
| 2 | "Your structural-vs-immune p-value treats 70 row as independent, but they share 23 mice." | Within-tissue permutation; tissue fixed effect; direction per tissue (07 §10) | **Pass-rate claim withdrawn** (p 0.13). Continuous difference holds (+0.14, p 0.001) |
| 3 | "Structural row may just be bigger, so easier to classify." | Fixed-effect model with log cell count and mouse count | Neither explains the gap (p 0.70, 0.36) |
| 4 | "Sex is tied to age group, and the model sees *Xist*." | Sex-only cross-age; sex-chromosome gene removed; male-only cross-age (07 §11–12) | Sex never fakes a pass (sex alone ≤ 0.68; all 29 pass without sex gene; 17 of 20 in male only). But it **hid** signal: 12 row rise without sex gene, immune pass rate 0.33 → 0.59 |
| 5 | "With 2–6 test mice, a cross-age score of 0.72 could be 0.55 with other mice." | Bootstrap over test mice (07 §11) | Every pass has its interval above 0.5; 6 near-cutoff row reach below 0.7 |
| 6 | "A batch shared by all young against all old passes every split." | Technical-only model | Bladder (2 row) and heart fibroblast flagged; 26 of 29 age-signal row clean. A batch that only shifts biology-like gene would still pass (limitation 11.3) |
| 7 | "Failing the cross-age check could be real biology that differs between age group, e.g. growth at 1 month." | — (design) | Agreed: failing means "not demonstrable", stated in 7.4 |
| 8 | "Robust row change more gene only because they have more mice." | Regression on mouse count (09 §9) | No: same median mice, ×4.9 after adjustment (p = 0.0006) |
| 9 | "Low overlap between cell type does not prove aging is private." | — (wording) | Agreed: reworded to "little shared signal detectable at this sample size" |
| 10 | "The transfer to another lab could be a shared technical difference." | Technical-only model trained on TMS, scored on Kimmel (10 §5b) | **Partly yes**: technical feature transfer in 12 of 19 row; expression clearly beats them in 10 (myeloid, epithelial, endothelial, adventitial). Lymphocyte transfer is unresolved |
| 11 | "Absence in bulk does not prove the p21 rise is an artifact: bulk dilutes fibroblast." | — (wording) | Agreed: "not replicated, cause unknown" replaces "artifact" |
| 12 | "Your model comparison p-values treat row as independent." | Row count per tissue | p-value removed; "lower in 14–16 of 19 row" reported. Ordinal clock wins in 8 of 9 tissue |
| 13 | "MHC class I 'follows' IRF1 activity partly because MHC gene are IRF1 target." | MHC gene removed from every target set (12 §3) | Not circular: *Irf1* 0.61 → 0.60, *Rela* 0.54 → 0.52; only *Ciita* drops (0.50 → 0.36). Wording changed from "driven by" to "co-varies with" |
| 14 | "The inflammatory signal in stromal cell may be ambient immune RNA, since old tissue hold more immune cell." | Activity recomputed without each row's foreign gene (12 §3b) | Interferon / NF-κB stay up in the same 6 tissue; the MHC class II driver *Ciita* halves (17 → 9 row), so class II in non-immune row is partly ambient |
| 15 | "SINDy equation from nine time point are overfit." | Count of coefficient against point | Agreed: the equation figure was removed; SINDy is kept as a forecaster that diverges in about a third of case |
| 16 | "The T-cell model reports a thymic half-life next to 'not identifiable'." | Profile likelihood; three fit at δ 100-fold apart | Agreed: no rate is reported; the figure now shows the non-identifiability |
| 17 | "The neural ODE 'win' is one seed, one held-out age, best of six model." | 5 seed; sign and Wilcoxon test with BH over model | **Claim withdrawn**: wins 8 of 17 (sign p 1.0); median skill swings from −0.37 to 0.29 with the seed alone |
| 18 | "Choosing gene by their age correlation and then showing they order by age is circular." | Gene chosen on half of the mice, the other half drawn | Phase portrait redrawn; it is labelled an illustration |
| 19 | "Forecast gene were chosen with the held-out mice included." | Variable gene chosen on training mice only | Fixed before the forecast was rerun |
| 20 | "Early-warning signal need repeated measurement of one system." | — (design) | Agreed: reported as cross-sectional heterogeneity, not a test of critical slowing down |
| 21 | "Optimal transport: the held-out 21-month mice are all female, their neighbour mostly male." | — (design) | Agreed: stated; the 2-of-7 result is reported as a count |
| 22 | "Thresholds were tuned on the lung." | — (disclosure) | Yes: lung was the development tissue; the rule was frozen before the other eight (section 8) |

### 11.2 Claim Withdrawn Or Downgraded

| First draft | Now | Why |
| --- | --- | --- |
| Structural pass rate 0.75 vs immune 0.33, Fisher p = 0.008 | Continuous cross-age higher within tissue; pass-rate difference not significant | Row are not independent (objection 2) |
| "Immune cell in marrow and elsewhere are batch driven" | Partly sex: 12 row pass once sex-chromosome gene are removed | Sex misaligned between age group (objection 4) |
| "Aging is mostly private to each cell type" | Little shared signal is detectable at this sample size | Low power produces low overlap (9) |
| "This settles p21: a dataset artifact" | Not replicated, cause unknown | Bulk dilutes a fibroblast change (11) |
| "Two regulator drive the shared signal" | Interferon and NF-κB activity co-vary with it | Activity inference is correlational (13) |
| "Old stroma becomes more visible to the immune system" (ligand–receptor) | Restates the MHC-I and *C3* rise | The score is built from the same gene change |
| "Only the neural ODE beats the nearest age" | No forecaster beats the nearest age, inside or beyond the range | One seed, best of six (17) |
| "19 of 19 row transfer to another lab" | 10 row transfer beyond what technical feature explain | Technical-only transfer (10) |
| Immune cell "carry a weak but real age signal another lab reproduces" | Myeloid cell do; lymphocyte unresolved | Technical-only transfer (10) |
| Kidney SINDy equation | Removed | Overfit (15) |
| Lung thymic half-life ~14 month | Removed | Not identifiable (16) |
| Spread grows "as the Karin model predicts" | Consistent with, not evidence for | Almost any model predicts it |
| Optimal transport "the structural-versus-immune split again" | 2 of 7 row, sex-confounded | Small count, sex (21) |

### 11.3 Limitation That Remain

- **Twenty-three donor for the whole atlas.** Every number is bounded by that,
  not by the cell count. Two lung donor are pooled sample treated as one unit.
- **Age is confounded with processing batch and with sex.** The cross-age, sex
  and technical control detect this and report it per row; none can undo it. A
  cell type that fails the check is not shown to lack an age signal, only to
  lack one this design can demonstrate.
- **A batch that shifts only biology-like gene is invisible.** The
  technical-only control sees depth, gene count, ribosomal, dissociation and
  blood fraction; a processing difference that moves, say, a stress-response
  program in a gene-specific way would still pass every check.
- **Seven of sixteen droplet tissue cannot be tested.**
- **Ambient RNA cannot be removed, only filtered.** SoupX and CellBender need
  the unfiltered droplet matrix, absent from this release. *Cd74* in alveolar
  macrophage sits close to the fitted cutoff (0.037 against 0.031).
- **No mitochondrial QC** is possible: the release carries no `mt-` gene.
- **The shuffled control is five permutation**, a floor rather than a calibrated
  p-value; the bootstrap interval holds training mice fixed and is therefore
  too narrow.
- **A gene weight is interpretation, not a test**, and the pseudobulk test has
  5 to 17 animal per row, so it detects only large change.
- **The single-cell atlas has no animal between 3 and 18 month.** The bulk time
  course fills it for whole organ only.
- **The replication is small**: 3–4 young and 3 old male mice per tissue, about
  10,000 gene (no *Cdkn2a*), and no statement on processing day.
- **The dynamical model are fitted to few time point** and are exploratory
  throughout.
- **Many analysis were run.** Exploratory p-value are not corrected across the
  whole project, which is why no exploratory result is used as support for a
  confirmatory claim.

### 11.4 Next Step

1. **Remove sex-chromosome gene by default and rerun 07–16.** Section 10.3
   shows they hide immune signal; the downstream module still use verdict set
   with them in. This is a code change of one default (`features` in
   `age_eval` and `cross_age_auc`) and a rerun, no new data.
2. **A design where age is not tied to batch**, the only real fix: a large
   human cohort with many donor pooled across run (OneK1K: about 1,000 donor aged
   19 to 97, blood; Yazar et al. 2022), where the immune signal TMS could not
   judge can be tested and the clock and forecasting model of 10.8 and 10.11 have
   enough independent donor to be told apart.
3. **Measured miRNA.** Check the miR-29 inference of 10.10 against mature miRNA
   measured in the same cohort (Wagner et al. 2024).
4. **A human check of the replicating gene.** MHC class I, complement and
   collagen loss against age in a human tissue atlas (GTEx).
5. **A direct senescence measure.** Acceleration (Karin) could not be tested on
   bulk p16; a dataset with a stronger senescence readout is needed.


## 12. Reference

**Data**

- Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes
  ageing tissues in the mouse. *Nature* **583**, 590–595 (2020).
- Schaum, N. et al. Ageing hallmarks exhibit organ-specific temporal signatures.
  *Nature* **583**, 596–602 (2020). The bulk time course.
- Kimmel, J. C. et al. Murine single-cell RNA-seq reveals cell-identity- and
  tissue-specific trajectories of aging. *Genome Research* **29**, 2088–2103
  (2019). The independent replication.
- Zhang, M. J., Pisco, A. O., Darmanis, S. & Zou, J. Mouse aging cell atlas
  analysis reveals global and cell type-specific aging signatures. *eLife*
  **10**, e62293 (2021). The published global aging gene.
- Wagner, V. et al. Characterizing expression changes in noncoding RNAs during
  aging and heterochronic parabiosis across mouse tissues. *Nature
  Biotechnology* **42**, 109–118 (2024).
- Yazar, S. et al. Single-cell eQTL mapping identifies cell type–specific genetic
  control of autoimmune disease. *Science* **376**, eabf3041 (2022). OneK1K.

**Aging biology**

- López-Otín, C. et al. Hallmarks of aging: an expanding universe. *Cell*
  **186**, 243–278 (2023).
- Baker, D. J. et al. Clearance of p16Ink4a-positive senescent cells delays
  ageing-associated disorders. *Nature* **479**, 232–236 (2011).
- Hernandez-Segura, A., Nehme, J. & Demaria, M. Hallmarks of cellular
  senescence. *Trends in Cell Biology* **28**, 436–453 (2018).
- Saul, D. et al. A new gene set identifies senescent cells and predicts
  senescence-associated pathways across tissues. *Nature Communications*
  **13**, 4827 (2022). SenMayo.
- Franceschi, C. et al. Inflamm-aging: an evolutionary perspective on
  immunosenescence. *Annals of the New York Academy of Sciences* **908**,
  244–254 (2000).
- Palmer, D. B. The effect of age on thymic function. *Frontiers in Immunology*
  **4**, 316 (2013).
- Ugalde, A. P. et al. Aging and chronic DNA damage response activate a
  regulatory pathway involving miR-29 and p53. *EMBO Journal* **30**,
  2219–2232 (2011).
- Ahadi, S. et al. Personal aging markers and ageotypes revealed by deep
  longitudinal profiling. *Nature Medicine* **26**, 83–90 (2020).
- Buckley, M. T. et al. Cell-type-specific aging clocks to quantify aging and
  rejuvenation in neurogenic regions of the brain. *Nature Aging* **3**,
  121–137 (2023).

**Single-cell method and statistics**

- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene
  expression data analysis. *Genome Biology* **19**, 15 (2018).
- Wolock, S. L., Lopez, R. & Klein, A. M. Scrublet: computational identification
  of cell doublets in single-cell transcriptomic data. *Cell Systems* **8**,
  281–291 (2019).
- Young, M. D. & Behjati, S. SoupX removes ambient RNA contamination from
  droplet-based single-cell RNA sequencing data. *GigaScience* **9**, giaa151
  (2020).
- van den Brink, S. C. et al. Single-cell sequencing reveals
  dissociation-induced gene expression in tissue subpopulations. *Nature
  Methods* **14**, 935–936 (2017).
- Hurlbert, S. H. Pseudoreplication and the design of ecological field
  experiments. *Ecological Monographs* **54**, 187–211 (1984).
- Squair, J. W. et al. Confronting false discoveries in single-cell differential
  expression. *Nature Communications* **12**, 5692 (2021).
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550 (2014).
- Muzellec, B. et al. PyDESeq2: a python package for bulk RNA-seq differential
  expression analysis. *Bioinformatics* **39**, btad547 (2023).
- Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate.
  *Journal of the Royal Statistical Society B* **57**, 289–300 (1995).

**Machine learning and prior knowledge**

- Zou, H. & Hastie, T. Regularization and variable selection via the elastic
  net. *Journal of the Royal Statistical Society B* **67**, 301–320 (2005).
- Ke, G. et al. LightGBM: a highly efficient gradient boosting decision tree.
  *NeurIPS* **30** (2017).
- Ilse, M., Tomczak, J. & Welling, M. Attention-based deep multiple instance
  learning. *ICML* (2018).
- Frank, E. & Hall, M. A simple approach to ordinal classification. *ECML*,
  145–156 (2001).
- Argelaguet, R. et al. MOFA+: a statistical framework for comprehensive
  integration of multi-modal single-cell data. *Genome Biology* **21**, 111
  (2020).
- Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to
  infer biological activities from omics data. *Bioinformatics Advances* **2**,
  vbac016 (2022).
- Müller-Dott, S. et al. Expanding the coverage of regulons from high-confidence
  prior knowledge for accurate estimation of transcription factor activities.
  *Nucleic Acids Research* **51**, 10934–10949 (2023). CollecTRI.
- Schubert, M. et al. Perturbation-response genes reveal signaling footprints in
  cancer gene expression. *Nature Communications* **9**, 20 (2018). PROGENy.
- Liberzon, A. et al. The Molecular Signatures Database hallmark gene set
  collection. *Cell Systems* **1**, 417–425 (2015).
- Dimitrov, D. et al. Comparison of methods and resources for cell-cell
  communication inference from single-cell RNA-Seq data. *Nature
  Communications* **13**, 3224 (2022). LIANA.

**Dynamical system**

- Karin, O., Agrawal, A., Porat, Z., Krizhanovsky, V. & Alon, U. Senescent cell
  turnover slows with age providing an explanation for the Gompertz law.
  *Nature Communications* **10**, 5495 (2019).
- Scheffer, M. et al. Early-warning signals for critical transitions. *Nature*
  **461**, 53–59 (2009).
- Raue, A. et al. Structural and practical identifiability analysis of partially
  observed dynamical models by exploiting the profile likelihood.
  *Bioinformatics* **25**, 1923–1929 (2009).
- Brunton, S. L., Proctor, J. L. & Kutz, J. N. Discovering governing equations
  from data by sparse identification of nonlinear dynamical systems. *PNAS*
  **113**, 3932–3937 (2016). SINDy.
- Chen, R. T. Q., Rubanova, Y., Bettencourt, J. & Duvenaud, D. Neural ordinary
  differential equations. *NeurIPS* **31** (2018).
- Schiebinger, G. et al. Optimal-transport analysis of single-cell gene
  expression identifies developmental trajectories in reprogramming. *Cell*
  **176**, 928–943 (2019). Waddington-OT.
- Klein, D. et al. Mapping cells through time and space with moscot. *Nature*
  **638**, 1065–1075 (2025).
