# Activity-By-Captured-Contact (ABCC)

ABCC (also known as CHiC-ABC) is a computational model for predicting the effects of enhancers on gene expression based on chromatin activity information and chromosomal contact data from high-resolution Promoter Capture Hi-C (PCHi-C). It is a modification of the Activity-By-Contact (ABC) model developed by Jesse Engreitz's lab [1][2] that uses PCHi-C data — rather than conventional Hi-C — as a source of chromosomal contact information.

Promoter Capture Hi-C (PCHi-C) is a method of enriching Hi-C libraries for contacts involving (at least on one end) gene promoters using hybridisation probes. Please see Freire-Pritchett et al. [3] for information about PCHi-C and its analysis tools.

ABCC leverages the statistical modelling of PCHi-C data by CHiCAGO software [3][4] for normalisation and imputation, and uses the results as a source of contact information, as described in more detail below.

The ABC code used for this implementation was downloaded from https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction on 22 September 2022 (latest commit: a74aa73).

---

## Requirements

For each cell type, the inputs to ABCC are:

* BAM file for DNase-seq (indexed and sorted) **or** tagAlign file for ATAC-seq
* BAM file for H3K27ac ChIP-seq (indexed and sorted) *(optional)*
* CHiCAGO-processed PCHi-C data — a CHiCAGO RDS object produced at Step 2 of the CHiCAGO pipeline (see the [PCHi-C imputation section](#step-2-running-pchic-imputation) below)
* A measure of gene expression (gene-level TPM) *(required only if using the original ABC algorithm for enhancer/gene definition)*

The following genome annotation files are also required (not cell-type-specific):

* BED file containing gene annotations (may vary across cell types if using cell-type-specific TSSs)
* BED file containing chromosome sizes
* A CHiCAGO design folder matching the experimental design of the PCHi-C experiment (contains `.baitmap` and `.rmap` files)

---

## The Principle of the ABC(C) Model

The Activity-by-Contact (ABC) model [1][2] represents a mechanistic model in which enhancers activate gene transcription upon enhancer–promoter contact. Quantitatively, the effect of an enhancer depends on how frequently it contacts a promoter multiplied by the intrinsic strength of the enhancer (i.e., its ability to activate transcription upon contact). The contribution of any specific enhancer to a gene's expression depends on the surrounding context — the activity and contact frequencies of all other nearby enhancers.

The ABC score is defined as:

> **ABC score for effect of element E on gene G** = Activity(E) × Contact(E, G) / Σ [Activity(e) × Contact(e, G)] over all candidate elements within 5 Mb

**Activity (A)** is the geometric mean of DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq read counts at element E. **Contact (C)** is the contact frequency between E and the promoter of gene G. In the original ABC model, KR-normalised Hi-C contact frequency is used. In ABCC, CHiCAGO-normalised and imputed contact frequencies are used instead (see below).

Candidate elements are ~500 bp regions centred on DHS/ATAC peaks.

---

## Running the PCHiC-ABC Model: Step-by-Step

Running ABCC consists of three sequential steps:

1. Define candidate enhancer and gene regions *(same as in the original ABC model)*
2. Run PCHi-C imputation *(new step introduced in ABCC)*
3. Compute ABC scores *(using a modified `predict.py`)*

---

### Step 1. Defining Candidate Enhancers and Genes

This step quantifies chromatin accessibility and H3K27ac activity at candidate enhancer regions, and computes analogous measures at gene promoters and bodies. It produces two files — `EnhancerList.txt` and `GeneList.txt` — that feed into the imputation step.

Either the original ABC scripts in this repository (`src_mod/makeCandidateRegions.py` and `src_mod/run.neighborhoods.py`) or the newer ABC Snakemake pipeline can be used, as described at https://abc-enhancer-gene-prediction.readthedocs.io/. A worked example for the original scripts applied to GM12878 is available at https://hoellin.github.io/eg/notes_ABC/generic_notebooks/turnkey_notebook_to_run_ABC_with_example_over_GM12878.html.

#### Substep 1a: Call peaks with MACS2

MACS2 must be run with `--call-summits` to produce summit coordinates, which `makeCandidateRegions.py` uses to centre candidate regions. A lenient significance threshold is recommended at this stage, because the next substep selects the strongest peaks by absolute read count rather than by p-value.

```bash
macs2 callpeak \
  -t /path/to/dnase_or_atac.bam \
  -f BAM \
  -n CELLTYPE \
  --call-summits \
  --nolambda \
  --outdir /path/to/macs2_output/
```

#### Substep 1b: Define candidate enhancer regions

`makeCandidateRegions.py` takes the MACS2 peaks, ranks them by absolute read count in the accessibility BAM, retains the top N (default 175,000), and extends each summit by `--peakExtendFromSummit` bp in each direction (default 250 bp, giving 500 bp regions). Promoter regions should be force-included via `--regions_includelist` even if they do not rank among the strongest peaks, ensuring all genes of interest have a promoter element in the candidate set. Block-listed regions should be excluded via `--regions_blocklist` (see https://sites.google.com/site/anshulkundaje/projects/blacklists).

```bash
python src_mod/makeCandidateRegions.py \
  --narrowPeak /path/to/macs2_output/CELLTYPE_peaks.narrowPeak \
  --bam /path/to/dnase_or_atac.bam \
  --chrom_sizes /path/to/chrom.sizes \
  --outDir /path/to/candidate_regions/ \
  --nStrongestPeaks 175000 \
  --regions_includelist /path/to/gene_promoters.bed \
  --regions_blocklist /path/to/blocklist.bed
```

Key output: `candidate_regions/CandidateRegions.bed` — a BED file of ~500 bp candidate enhancer regions.

#### Substep 1c: Count reads over candidate regions and gene bodies

`run.neighborhoods.py` counts DNase/ATAC-seq and H3K27ac ChIP-seq reads over the candidate regions from the previous substep, and over gene bodies and promoters. It produces the two files consumed by the imputation step.

```bash
python src_mod/run.neighborhoods.py \
  --candidate_enhancer_regions /path/to/candidate_regions/CandidateRegions.bed \
  --outdir /path/to/neighborhoods/ \
  --genes /path/to/gene_annotations.bed \
  --DHS /path/to/dnase.bam \
  --H3K27ac /path/to/h3k27ac.bam \
  --chrom_sizes /path/to/chrom.sizes \
  --cellType CELLTYPE
```

Use `--ATAC` in place of `--DHS` for ATAC-seq input. `--H3K27ac` may be omitted if only one accessibility assay is available (activity will then be computed from accessibility alone). The `--genes` file must be in BED-6 format (chr, start, end, name, score, strand); the name field is used to match genes to the PCHi-C baitmap in Step 2.

#### Gene expression data (optional)

The `--expression_table` flag accepts a path to a gene expression file, or a comma-separated list of paths if you have expression data from multiple samples. Each file must be a **headerless, two-column, tab-separated** file:

```
GENE_SYMBOL<TAB>TPM_VALUE
BRCA1           12.4
TP53            8.7
...
```

The gene identifiers in column 1 must match the `--primary_gene_identifier` (default: `symbol`). If a gene appears more than once (e.g. multiple transcript isoforms), the maximum expression value is retained. If multiple files are provided, the per-gene expression values are averaged across files after the per-file maximum has been taken.

**If `--expression_table` is omitted**, `Expression` is set to `NaN` for all genes. In this case, Step 3 (`predict.py`) falls back to using promoter chromatin activity as a proxy for expression: genes whose `PromoterActivityQuantile` is at or above 0.4 (the top 60% by promoter accessibility/H3K27ac activity) are treated as expressed and included in predictions. This fallback means expression data is not strictly required, but providing it gives more accurate control over which genes are modelled.

Key outputs:

* `neighborhoods/EnhancerList.txt` — candidate enhancer regions with read counts, RPKM, and quantile-normalised activity scores
* `neighborhoods/GeneList.txt` — gene regions with TSS coordinates, activity measures, and expression values (or `NaN` if no expression file was provided)

These two files are required by the imputation step.

---

### Step 2. Running PCHi-C Imputation

#### How imputation works

PCHi-C data is inherently sparse: it only measures contacts involving captured (baited) promoters, and even within the captured set, not every bait–other-end pair will have sequencing reads. ABCC needs a contact frequency estimate for *every* candidate enhancer–gene pair within 5 Mb, so the imputation step fills in the gaps using CHiCAGO's statistical model.

The imputation procedure applies a three-tier hierarchy to each candidate enhancer–gene pair:

**Tier 1 — Observed contact frequency.** If the pair appears in the CHiCAGO object with read counts exceeding the expected Brownian noise level (N > Bmean), the CHiCAGO-normalised observed count is used directly:

> N\_imp = N / (s\_i × s\_j)

where N is the observed read count and s\_i, s\_j are the bait- and other-end-specific scaling factors estimated by CHiCAGO.

**Tier 2 — Brownian noise floor.** If the pair appears in the CHiCAGO object but the observed count is at or below the Brownian noise level (Bmean ≥ N), the expected count based on the Brownian collision model is used instead:

> N\_imp = Bmean / (s\_i × s\_j)

This prevents low observed counts from underestimating the contact frequency expected by random diffusion alone.

**Tier 3 — Distance-function imputation.** If the pair does not appear in the CHiCAGO object at all — because the promoter was not captured in the experimental design, reads at that bait were filtered during QC, or the specific bait–other-end combination has no recorded reads — the expected contact frequency is estimated purely from the interaction distance using CHiCAGO's fitted distance decay function f(d):

> N\_imp = f(|distance from enhancer midpoint to TSS|)

This function is estimated from the full CHiCAGO object (by `chic_split.r`) and saved to the `distout_rds` file for use during imputation.

**Important caveat on Tier 3 (whole-bait imputation).** When an entire promoter bait is absent from the CHiCAGO object — most commonly because the promoter was not included in the PCHi-C capture design — *all* of its enhancer contacts are estimated from the genome-wide distance function with no locus-specific information. This function represents a reasonable average expectation across all baits but cannot capture locus-specific variation in contact frequency. Predictions for such genes should be interpreted with greater caution than those backed by observed PCHi-C data. The imputation script emits diagnostic messages partitioning the distance-function-imputed contacts into three sub-cases — promoter not in baitmap at all; promoter in baitmap but absent from the CHiCAGO object (likely QC-filtered); or the specific bait–other-end contact absent from the CHiCAGO object despite both fragments being present — which can help diagnose the source and scale of the imputation.

**Contact capping.** For enhancers very close to the TSS (within ~1,500 bp, approximately one median restriction fragment), contact frequencies are capped at the value predicted by the distance function at 1,500 bp. This prevents artefactually high contact scores for elements in the immediate vicinity of the promoter.

After imputation, the final contact value assigned to each enhancer–gene pair is the maximum across all applicable tier estimates. Where a candidate enhancer overlaps multiple restriction fragments, the fragment yielding the highest contact is selected.

#### Running the imputation script

`imputation_script_whole_genome.r` orchestrates the full imputation process. It automatically calls two helper scripts internally:

* `chic_split.r` — splits the CHiCAGO object by chromosome and extracts the distance function parameters
* `chr_split.r` — splits `EnhancerList.txt` and `GeneList.txt` by chromosome

**You do not need to run these helper scripts manually.**

```bash
pchic_rds=/path/to/K562_fres_merged_reweighted.Rds        # CHiCAGO Step 2 RDS object
design=/path/to/Design/Human_DpnII_binsize1500_maxL75K/    # CHiCAGO design folder (must end with /)
enhancerdir=/path/to/neighborhoods/                         # directory with EnhancerList.txt and GeneList.txt (must end with /)
split_pchic=/path/to/neighborhoods/K562_frag_dpnii_         # full path prefix for intermediate per-chromosome CHiCAGO RDS files
distout_rds=/path/to/neighborhoods/K562_frag_dpnii_dist.rds # output path for distance parameters RDS
imputed_pchic_prefix=/path/to/neighborhoods/imputed_contact/ # output directory (must end with /, must not already exist)

Rscript chic_imputation/imputation_script_whole_genome.r \
  ${pchic_rds} \
  ${design} \
  ${enhancerdir} \
  ${split_pchic} \
  ${distout_rds} \
  ${imputed_pchic_prefix}
```

**Arguments in order:**

| # | Variable | Description |
|---|----------|-------------|
| 1 | `pchic_rds` | CHiCAGO Step 2 RDS object for your cell type |
| 2 | `design` | Path to CHiCAGO design directory; must end with `/` |
| 3 | `enhancerdir` | Directory containing `EnhancerList.txt` and `GeneList.txt` from Step 1; must end with `/` |
| 4 | `split_pchic` | Full path **prefix** (including filename prefix) for intermediate per-chromosome CHiCAGO RDS files; chromosome name and `.rds` are appended automatically |
| 5 | `distout_rds` | Output path for the distance function parameters RDS file |
| 6 | `imputed_pchic_prefix` | Output directory for imputed contact files; must end with `/` and **must not already exist** |

#### Intermediate files

During processing, per-chromosome splits of the CHiCAGO object (`{split_pchic}{N}.rds`) and of the enhancer/gene lists (`EnhancerList_chr{N}.txt`, `GeneList_chr{N}.txt`) are created and then automatically deleted once each chromosome has been processed, to keep disk usage manageable. The distance function parameters file (`{distout_rds}`) is retained; it can be reused if you need to re-run imputation with different candidate regions against the same PCHi-C dataset.

#### Output files

For each chromosome N, two output files are written to `imputed_contact/chr{N}/`:

**`chr{N}.bedpe.gz`** — the primary contact file used by Step 3. Eight tab-separated columns, no header:

| Column | Name | Description |
|--------|------|-------------|
| 1 | `bait_chr` | Chromosome of the gene (with `chr` prefix) |
| 2 | `bait_start` | TSS position (bp) |
| 3 | `bait_end` | TSS position (same as `bait_start`; the bait is represented as a point at the TSS) |
| 4 | `otherEnd_chr` | Chromosome of the enhancer (with `chr` prefix) |
| 5 | `otherEnd_start` | Enhancer midpoint (integer bp) |
| 6 | `otherEnd_end` | Enhancer midpoint + 1 |
| 7 | `otherEnd_name` | Name of the restriction fragment overlapping the enhancer (from the `.rmap` file) |
| 8 | `contact` | Imputed contact frequency on the CHiCAGO-normalised scale |

**`CG_score_chr{N}.bedpe`** — the same 8 columns plus three additional columns retained for diagnostic purposes (not used by Step 3):

| Column | Name | Description |
|--------|------|-------------|
| 9 | `CG_score` | BED score field from the gene annotation file (column 5 of the input gene BED) |
| 10 | `otherEndID` | Numeric restriction fragment ID of the other end (from `.rmap`) |
| 11 | `baitID` | Numeric restriction fragment ID of the bait (from `.baitmap`) |

All diagnostic plots are written to `imputed_contact/plots/`.

#### Diagnostic plots

The script produces 14 PDF plots in total: a set of per-chromosome plots produced inside the chromosome loop, and a set of genome-wide summary plots produced once after all chromosomes have been processed.

**Per-chromosome plots** (filename suffix `_chr{N}.pdf`, saved in `imputed_contact/plots/`)

| Filename prefix | What it shows |
|-----------------|---------------|
| `observed.only.distr_N_v_Bmean_for_chr` | *Pre-imputation sanity check.* Boxplot of normalised observed contact frequency (asinh scale) vs distance bins up to 200 kb, coloured by whether N > Bmean or Bmean ≥ N. Useful for verifying that CHiCAGO's Brownian noise estimate is behaving as expected in the raw data before imputation. |
| `distr_imputed.data_N_v_Bmean_for_chr` | *Post-imputation contact distribution, up to 200 kb.* Same layout but using the final imputed `contact` values, coloured by imputation tier. Contacts inferred from observed data should generally sit above distance-function-imputed contacts at the same distance. |
| `distr_imputed.data_v.close.range_N_v_Bmean_for_chr` | Same as above, zoomed to ±25 kb. Useful for examining imputation behaviour at short range where the distance function may diverge most from the observed data. |
| `distr_imputed.data_v.v.close.range_N_v_Bmean_for_chr` | Same, zoomed further to ±5 kb. |
| `distance_distr_imputed_contacts_for_chr` | Two-panel histogram (count and density) of enhancer–gene pair distances, split by tier: observed vs distance-function-imputed. A large proportion of distance-function-imputed contacts at all ranges suggests limited PCHi-C capture coverage. |
| `count_distr_imputed_contacts_for_chr` | Boxplot of final contact frequency (asinh scale) vs distance bins up to 5 Mb, coloured by tier. |
| `count_distr_close_range_imputed_contacts_for_chr` | Same, zoomed to 1 Mb. |
| `count_distr_v_close_range_imputed_contacts_for_chr` | Same, zoomed to 100 kb. |

**Genome-wide summary plots** (saved directly in `imputed_contact/plots/`)

| Filename | What it shows |
|----------|---------------|
| `genomewide_contact_distribution.pdf` | Boxplot of imputed contact frequency vs distance (up to 1 Mb), genome-wide, coloured by tier. The top-level overview of imputation quality across the full dataset. |
| `genomewide_contact_distribution_100kb.pdf` | Same, zoomed to 100 kb. |
| `genomewide_contact_distribution_25kb.pdf` | Same, zoomed to 25 kb. |
| `genomewide_contact_distribution_5kb.pdf` | Same, zoomed to 5 kb. |
| `genomewide_distance_histogram_count.pdf` | Genome-wide histogram (raw count) of pair distances split by tier: observed vs distance-function-imputed. |
| `genomewide_distance_histogram_density.pdf` | Same, density-normalised. |

**Note on compute time.** Imputation is computationally intensive. The script processes chromosomes sequentially; for large datasets, consider running on a compute cluster.

---

### Step 3. Computing the ABC Score

Combine Activity (from Step 1) with the imputed PCHi-C contacts using the modified `predict.py`. You must use `src_mod/predict.py` from this repository — it has been modified to read BEDPE-format imputed contact files. The `--hic_type bedpe` flag is required; without it, `predict.py` will attempt to read Juicebox-format Hi-C matrices instead.

```bash
python src_mod/predict.py \
  --enhancers /path/to/neighborhoods/EnhancerList.txt \
  --genes /path/to/neighborhoods/GeneList.txt \
  --HiCdir /path/to/neighborhoods/imputed_contact/ \
  --chrom_sizes reference/chrom.sizes.og \
  --threshold .02 \
  --cellType K562 \
  --hic_type bedpe \
  --outdir ABC_scores \
  --make_all_putative
```

**Key flags:**

| Flag | Description |
|------|-------------|
| `--enhancers` | `EnhancerList.txt` from Step 1 |
| `--genes` | `GeneList.txt` from Step 1 |
| `--HiCdir` | The `imputed_pchic_prefix` output directory from Step 2; the script looks for `chr{N}/chr{N}.bedpe.gz` within this directory |
| `--hic_type bedpe` | **Required** — instructs `predict.py` to read imputed BEDPE files rather than Hi-C matrices |
| `--threshold` | Minimum ABC score to include in the main output (default 0.02; see Tips below) |
| `--cellType` | Label for the cell type, used in output file naming |
| `--outdir` | Directory for output ABC score files |
| `--make_all_putative` | Output scores for all enhancer–gene pairs, not just those above the threshold |
| `--expression_cutoff` | Genes with `Expression` below this value are excluded from predictions (default 1; applies only if expression data was provided in Step 1) |
| `--promoter_activity_quantile_cutoff` | When `Expression` is `NaN` (no expression data provided), genes with `PromoterActivityQuantile` at or above this quantile are treated as expressed (default 0.4, i.e. top 60%) |
| `--run_all_genes` | Skip expression filtering entirely and make predictions for all genes regardless of expression or promoter activity |

---

## Human Genome Assembly Conversion (hg19 ↔ hg38)

If your PCHi-C data and ABC candidate regions are on different genome assemblies, use `chic_imputation/enhancerliftoverwrapper.R` to liftover `EnhancerList.txt` and `GeneList.txt` between hg19 and hg38. The `.chain` file must be obtained separately (e.g. from the UCSC genome browser) and its path edited into the script before running.

```bash
Rscript chic_imputation/enhancerliftoverwrapper.R EnhancerList.txt GeneList.txt
```

This produces `EnhancerListhg38.txt` and `GeneListhg38.txt` alongside the originals.

---

## Tips

* **Choosing an ABCC score threshold.** The optimal threshold is highly data-dependent. One approach: plot the correlation between ABCC numerators (unnormalised ABCC scores) per gene and gene expression levels across a range of cutoffs, and choose the lowest cutoff at which an appreciable correlation is achieved.
* **Checking imputation coverage.** After imputation, examine the diagnostic messages printed to stdout — these report the percentage of enhancer–gene pairs falling into each imputation tier, including sub-categories of Tier 3. A high fraction of distance-function-imputed contacts may indicate limited PCHi-C capture coverage, and predictions for the affected genes should be treated with appropriate caution.
* **Output directory must not already exist.** The imputation script will stop with an error if `imputed_pchic_prefix` already exists, to prevent accidental overwriting of results.

---

## Contributors

* ABCC development was led by Pavel Artemov supervised by Mikhail Spivakov, with valuable contributions by Joseph Ellaway, Valeriya Malysheva and Helen Ray-Jones.

## Acknowledgements

* We thank the Engreitz lab for developing the original ABC model and making the codebase accessible.
* We thank Jonathan Cairns for leading the development of CHiCAGO in Mikhail Spivakov's lab and for helpful discussions about CHiCAGO-based imputation.
* We thank Stephen Rong for testing the ABCC software and for helpful suggestions on improving its usability.

## Contact

Please submit a GitHub issue with any questions or if you encounter any bugs.

---

## Dependencies

**Imputation step (R):** CHiCAGO, data.table, dplyr, stringr, ggplot2, ggpubr

**Modified `predict.py` (Python):** pyranges, numpy, pandas, scipy, pysam

**Liftover script (R):** GenomicRanges, rtracklayer, data.table

**If using the original ABC code for Step 1:** samtools (0.1.19), bedtools (2.26.0), Tabix (0.2.5), MACS2 (2.1.1.20160309), Java (1.8), Juicer Tools (1.7.5), pyBigWig (0.3.2)

Conda environments are provided in `r43.yml` (full environment) and `r43_nohistory.yml` (without history). The conda environment for the Snakemake ABC pipeline is described in the ABC docs at https://abc-enhancer-gene-prediction.readthedocs.io/.

---

## References

[1] Fulco CP, Nasser J, Jones TR, et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nat. Genet.* 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

[2] Nasser J, Bergman DT, Fulco CP, et al. Genome-wide enhancer maps link risk variants to disease genes. *Nature* 593, 238–243 (2021). https://doi.org/10.1038/s41586-021-03446-x

[3] Freire-Pritchett P, Ray-Jones H, Della Rosa M, et al. Detecting chromosomal interactions in Capture Hi-C data with CHiCAGO and companion tools. *Nat. Protoc.* 16, 4144–4176 (2021). https://doi.org/10.1038/s41596-021-00567-5

[4] Cairns J, Freire-Pritchett P, Wingett SW, et al. CHiCAGO: robust detection of DNA looping interactions in Capture Hi-C data. *Genome Biol.* 17, 127 (2016). https://doi.org/10.1186/s13059-016-0992-2

## Citation
Malysheva, V., Ray-Jones, H., Lakes, N. et al. High-resolution promoter interaction analysis implicates genes involved in activation of type 3 innate lymphoid cells in immune disease risk. Nat Genet 58, 1953–1966 (2026). https://doi.org/10.1038/s41588-026-02681-0
