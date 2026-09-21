# Functional GAP Pipelines

Functional GAP is a SLURM-based workflow for processing paired-end shotgun metagenomic reads from raw FASTQ files to quality-controlled reads, metagenome assemblies, a non-redundant gene catalog, gene abundance estimates, functional and taxonomic annotation, contig-level ancient DNA damage assessment, and final integrated analysis tables.

The workflow is designed for HPC environments and is particularly suitable for environmental and sedimentary ancient DNA (sedaDNA) metagenomic datasets.

## Contents

1. [Workflow overview](#workflow-overview)
2. [Before you start](#before-you-start)
3. [Get the pipeline](#get-the-pipeline)
4. [Project setup](#project-setup)
5. [Input reads and sample naming](#input-reads-and-sample-naming)
6. [R environment](#r-environment)
7. [Databases and cluster-specific paths](#databases-and-cluster-specific-paths)
8. [Running the pipeline](#running-the-pipeline)
9. [Metadata](#metadata)
10. [Main outputs](#main-outputs)
11. [Output directory structure](#output-directory-structure)
12. [Important interpretation notes](#important-interpretation-notes)
13. [Quick QC checks](#quick-qc-checks)
14. [Troubleshooting](#troubleshooting)
15. [Recommended run order](#recommended-run-order)
16. [Notes for shared HPC installations](#notes-for-shared-hpc-installations)

---

## Workflow overview

The pipeline contains two complementary branches after assembly.

```text
Paired-end raw FASTQ
        |
        v
Step 1  FastQC -> optional deduplication -> fastp
        |
        v
Step 2  Tadpole error correction -> MEGAHIT assembly
        |
        +------------------------------------------------------+
        |                                                      |
        | Gene-centric branch                                  | Contig/aDNA branch
        |                                                      |
        v                                                      v
Step 3  Prodigal                                      Step 8.1 Bowtie2
        |                                                      |
        v                                                      v
Step 4  CD-HIT non-redundant gene catalog             Step 8.2 pyDamage 1.x
        |                                                      |
        +-----------------------+                              v
        |                       |                      Step 8.3 Kraken2
        v                       v                              |
Step 5  Salmon abundance   Step 6 eggNOG                       v
        |                       |                      Step 8.4 pyDamage
        |                       |                      + Kraken2
        |                       |                      + NCBI lineage
        |                       |                      + optional metadata
        |                       |
        |                       v
        |                  Step 7 MMseqs2 taxonomy
        |                       |
        +-----------+-----------+
                    |
                    v
             Step 9 final gene table
             abundance + function
             + taxonomy
             + optional metadata
```

Step 8.4 and Step 9 produce different final products:

- **Step 8.4** produces a **contig-level** table combining pyDamage, Kraken2, NCBI taxonomy and optional sample metadata.
- **Step 9** produces a **gene-level** table combining Salmon abundance, eggNOG functional annotation, MMseqs2 taxonomy and optional sample metadata.

Step 8 does not feed into Step 9. It is a separate contig-level authentication/taxonomy branch.

---

## Before you start

The workflow assumes:

- paired-end shotgun sequencing reads,
- a Linux HPC environment,
- SLURM for job submission,
- access to the required software modules,
- sufficient storage for assemblies, gene catalogs and databases,
- and local access to eggNOG, MMseqs2/NR and Kraken2 databases.

All pipeline scripts use the directory from which the job is submitted as the working directory:

```bash
WORK=${PWD}
```

**Always submit jobs from the project/pipeline root directory.**

For example:

```bash
cd /path/to/my_project/Functional_GAP_pipelines
sbatch ...
```

Do not mix unrelated datasets in the same working directory. Several steps automatically discover all matching sample files and combine them into a single catalog.

---

## Get the pipeline

### From GitHub

```bash
git clone https://github.com/PolarTerrestrialEnvironmentalSystems/Functional_GAP_pipelines.git
cd Functional_GAP_pipelines
```

### From a shared server installation

If the pipeline is installed in a shared directory, it is recommended to copy the full pipeline into your own project directory before running it. This allows you to change database paths, SLURM settings and metadata without modifying the shared copy.

For example:

```bash
mkdir -p /path/to/my_project

rsync -a \
    /path/to/shared/Functional_GAP_pipelines/ \
    /path/to/my_project/Functional_GAP_pipelines/

cd /path/to/my_project/Functional_GAP_pipelines
```

Keep the directory structure intact because the SLURM scripts call files from `R_script/`, `Python_script/`, `environment/` and `examples/`.

---

## Project setup

A typical starting directory can look like this:

```text
Functional_GAP_pipelines/
├── raw_reads/
│   ├── sample01_R1_001.fastq.gz
│   ├── sample01_R2_001.fastq.gz
│   ├── sample02_R1_001.fastq.gz
│   ├── sample02_R2_001.fastq.gz
│   └── ...
├── Python_script/
├── R_script/
├── environment/
├── examples/
├── step1_preprocessing.sl
├── step2_de_novo_assembly.sl
├── step3_prodigal_gene_pred.sl
├── step4_cd_hit.sl
├── step5.1_salmon_index.sl
├── step5.2_abundance_salmon.sl
├── step5.3_salmon_merge.sl
├── step6_eggnog_annotation.sl
├── step7_mmseqs2_taxonomy.sl
├── step8.1_bowtie_mapping.sl
├── step8.2_pydamage.sl
├── step8.3_kraken2.sl
├── step8.4_pydamage_annotation_merge.sl
└── step9_annotation_merge.sl
```

The raw-read directory name does not have to be `raw_reads/`; set the correct path in `step1_preprocessing.sl`.

The pipeline creates an `output/` directory and its subdirectories as the workflow progresses.

---

## Input reads and sample naming

Each sample must have one R1 and one R2 file.

For example:

```text
sample01_R1_001.fastq.gz
sample01_R2_001.fastq.gz
sample02_R1_001.fastq.gz
sample02_R2_001.fastq.gz
```

or:

```text
sample01_R1.fastq.gz
sample01_R2.fastq.gz
```

Before running Step 1, open `step1_preprocessing.sl` and check:

```bash
INDIR=...
R1_ENDING=...
R2_ENDING=...
```

The endings must match your input filenames exactly.

The sample identifier is derived from the R1 filename after removing `R1_ENDING`, and the same identifier is propagated through the pipeline.

Avoid spaces in sample names.

Recommended examples:

```text
Lake01_001
KL77_sample03
Sediment_A12
```

---

## R environment

The final merge steps use R.

A reproducible Conda environment is provided:

```text
environment/r_environment.yml
```

Create it once:

```bash
conda env create -f environment/r_environment.yml
```

The environment name is:

```text
functional_gap_r
```

Test it with:

```bash
conda activate functional_gap_r

Rscript -e \
"library(data.table); library(dplyr); library(tidyr); library(stringr); cat('R environment OK\n')"

conda deactivate
```

`step8.4_pydamage_annotation_merge.sl` and `step9_annotation_merge.sl` activate this environment automatically.

Check the following variable in those scripts:

```bash
CONDA_SH="${HOME}/miniforge3/etc/profile.d/conda.sh"
```

If Conda is installed elsewhere on your system, change this path.

---

## Databases and cluster-specific paths

Before the first run, check the database paths in the relevant scripts.

### eggNOG-mapper

In:

```text
step6_eggnog_annotation.sl
```

set the local eggNOG database directory to the correct path for your system.

### MMseqs2 taxonomy

In:

```text
step7_mmseqs2_taxonomy.sl
```

set the pre-built MMseqs2 NR taxonomy database:

```bash
NR_DB="/path/to/mmseqs2/NR_ncbi"
```

### Kraken2

In:

```text
step8.3_kraken2.sl
```

set:

```bash
KRAKEN_DB="/path/to/kraken2/database"
```

The standard workflow uses:

```text
confidence = 0
```

unless you deliberately want to change the classification behaviour.

### NCBI taxonomy for Step 8.4

The user does **not** need to manually prepare an `ncbi_df.csv`.

Step 8.4 automatically downloads the NCBI `new_taxdump` archive on first use and extracts:

```text
rankedlineage.dmp
merged.dmp
```

into:

```text
resources/ncbi_taxonomy/
```

These files are reused on later runs.

`merged.dmp` is used to resolve deprecated Kraken TaxIds to current NCBI TaxIds before lineage annotation.

If compute nodes cannot access the internet, download/extract these two files once beforehand and place them in:

```text
resources/ncbi_taxonomy/
```

Step 8.4 will detect and reuse them.

### SLURM settings

The `#SBATCH` headers are templates. Check at least:

```text
--account
--partition
--qos
--time
--cpus-per-task
--mem
```

before running on another cluster.

The module names and versions used by the workflow are defined directly in the individual SLURM scripts.

---

## SLURM arrays

Sample-level steps are submitted as SLURM arrays.

You can supply the array at submission time without editing the file:

```bash
sbatch --array=1-38%6 step1_preprocessing.sl
```

Here:

- `1-38` = 38 samples,
- `%6` = run at most 6 samples simultaneously.

Change these values for your dataset and available resources.

The following steps are **sample-level array jobs**:

```text
Step 1
Step 2
Step 3
Step 5.2
Step 8.1
Step 8.2
Step 8.3
```

The following are **single/global jobs**:

```text
Step 4
Step 5.1
Step 5.3
Step 6
Step 7
Step 8.4
Step 9
```

Do not start a downstream step until its required upstream step has completed successfully.

---

# Running the pipeline

## Step 1 — Read preprocessing

Script:

```text
step1_preprocessing.sl
```

Step 1 performs:

1. FastQC on raw reads,
2. optional read deduplication with Clumpify,
3. fastp filtering,
4. paired-read merging with fastp,
5. FastQC on processed reads.

Before submission, check the user-editable variables in the script, particularly:

```bash
INDIR=...
R1_ENDING=...
R2_ENDING=...
DEDUP=...
FILTER=...
```

If deduplication is not wanted:

```bash
DEDUP="FALSE"
```

Submit one array task per sample:

```bash
sbatch --array=1-N%M step1_preprocessing.sl
```

Replace `N` with the number of samples and `M` with the maximum number of concurrent jobs.

Main outputs:

```text
output/out.fastqc/pre/
output/out.fastqc/post/
output/out.dedup/
output/out.fastp/
```

For each sample, fastp produces files including:

```text
<sample>_fastp_R1.fq.gz
<sample>_fastp_R2.fq.gz
<sample>_fastp_merged_R2.fq.gz
```

The merged-read filename is used by later scripts to discover sample IDs.

---

## Step 2 — Error correction and assembly

Script:

```text
step2_de_novo_assembly.sl
```

Step 2:

1. performs read error correction with Tadpole,
2. assembles corrected paired and merged reads with MEGAHIT.

Submit:

```bash
sbatch --array=1-N%M step2_de_novo_assembly.sl
```

Main directories:

```text
output/out.tadpole/
output/out.megahit/
```

The main assembly for each sample is:

```text
output/out.megahit/<sample>/final.contigs.fa
```

This file is reused by Prodigal, Bowtie2 and Kraken2.

---

## Step 3 — Gene prediction

Script:

```text
step3_prodigal_gene_pred.sl
```

Prodigal is run in metagenomic mode on each sample assembly.

Submit:

```bash
sbatch --array=1-N%M step3_prodigal_gene_pred.sl
```

Main outputs:

```text
output/out.prodigal/<sample>.faa
output/out.prodigal/<sample>.fna
output/out.prodigal/<sample>.gff
```

Where:

- `.faa` = predicted proteins,
- `.fna` = predicted coding sequences,
- `.gff` = gene coordinates/annotations.

---

## Step 4 — Non-redundant gene catalog

Script:

```text
step4_cd_hit.sl
```

Run Step 4 only after **all Step 3 samples are complete**.

The script:

1. combines predicted proteins from all samples,
2. combines predicted nucleotide coding sequences,
3. clusters proteins with CD-HIT,
4. extracts the nucleotide sequences corresponding to the non-redundant protein catalog.

Submit once:

```bash
sbatch step4_cd_hit.sl
```

Important outputs:

```text
output/out.prodigal/redundant_prot_catalog.faa
output/out.prodigal/redundant_pCDS_catalog.fna
output/out.prodigal/non_redundant_prot_catalog.faa
output/out.prodigal/non_redundant_prot_catalog.faa.clstr
output/out.prodigal/non_redundant_prot_catalog.id
output/out.prodigal/non_redundant_pCDS_catalog.fna
```

The two key catalogs used downstream are:

```text
non_redundant_prot_catalog.faa
non_redundant_pCDS_catalog.fna
```

---

## Step 5 — Gene abundance with Salmon

Step 5 contains three parts.

### Step 5.1 — Build the Salmon index

Script:

```text
step5.1_salmon_index.sl
```

The index is built from:

```text
output/out.prodigal/non_redundant_pCDS_catalog.fna
```

Submit once:

```bash
sbatch step5.1_salmon_index.sl
```

The generic index is stored under:

```text
output/out.salmon/non_redundant_pCDS_catalog.index/
```

The index name is independent of site, core or dataset name.

### Step 5.2 — Quantify each sample

Script:

```text
step5.2_abundance_salmon.sl
```

Salmon quantifies the corrected unmerged paired reads and the corrected merged reads against the same non-redundant pCDS catalog.

Submit one array task per sample:

```bash
sbatch --array=1-N%M step5.2_abundance_salmon.sl
```

Per-sample outputs are stored under:

```text
output/out.salmon/
```

with paired and merged quantification directories.

### Step 5.3 — Combine Salmon results

Script:

```text
step5.3_salmon_merge.sl
```

The helper script:

```text
Python_script/sum_up_qc_merged_paired.py
```

combines paired and merged quantification results for each sample.

Submit once:

```bash
sbatch step5.3_salmon_merge.sl
```

Final merged abundance files are written to:

```text
output/out.salmon_merged_paired/
```

Important files:

```text
output_all_tpm_quant.sf
output_all_raw_quant.sf
output_gene_quant.raw.count.len
```

Interpretation:

- `output_all_tpm_quant.sf` contains Salmon TPM values.
- `output_all_raw_quant.sf` contains Salmon `NumReads` values and is the default input for the final annotation merge and count-based downstream analyses.
- `output_gene_quant.raw.count.len` contains target lengths.

For this gene-catalog workflow, TPM is best interpreted as a **gene-length- and library-size-normalized relative abundance measure**, not as transcript expression.

For count-based statistical methods, use the raw `NumReads` table.

---

## Step 6 — Functional annotation with eggNOG-mapper

Script:

```text
step6_eggnog_annotation.sl
```

Input:

```text
output/out.prodigal/non_redundant_prot_catalog.faa
```

Submit once:

```bash
sbatch step6_eggnog_annotation.sl
```

Main output directory:

```text
output/out.eggnog/
```

The main annotation file used by Step 9 is:

```text
output/out.eggnog/non_redundant_prot_catalog_eggNOG.emapper.annotations
```

---

## Step 7 — Gene taxonomy with MMseqs2

Script:

```text
step7_mmseqs2_taxonomy.sl
```

MMseqs2 assigns taxonomy to the non-redundant protein catalog against a locally prepared NR taxonomy database.

Submit once:

```bash
sbatch step7_mmseqs2_taxonomy.sl
```

Main output directory:

```text
output/out.mmseqs2/
```

The taxonomy table used in Step 9 is:

```text
output/out.mmseqs2/non_redundant_prot_catalog.faa.result.tsv
```

The workflow retains taxonomic lineage information including:

```text
species
genus
family
order
class
phylum
kingdom
superkingdom
```

when available.

---

## Step 8 — Ancient DNA and contig taxonomy

Step 8 is a contig-level branch for mapping, damage assessment and taxonomy.

It can be run after Step 2 because it uses the sample assemblies and processed reads directly.

### Step 8.1 — Bowtie2 mapping

Script:

```text
step8.1_bowtie_mapping.sl
```

Processed reads are mapped back to each sample's MEGAHIT assembly:

```text
output/out.megahit/<sample>/final.contigs.fa
```

Merged and paired mappings are combined into a sorted BAM.

Submit:

```bash
sbatch --array=1-N%M step8.1_bowtie_mapping.sl
```

Main output:

```text
output/out.bowtie/<sample>.merge_paired.sorted.bam
```

and its BAM index.

### Step 8.2 — pyDamage

Script:

```text
step8.2_pydamage.sl
```

pyDamage 1.x is run on the mapped BAM for each sample.

Submit:

```bash
sbatch --array=1-N%M step8.2_pydamage.sl
```

Main output:

```text
output/out.pydamage/<sample>/pydamage_results.csv
```

The native pyDamage output is retained.

The merge workflow recognizes pyDamage 1.x model fields including:

```text
predicted_accuracy
null_model_p0
null_model_p0_stdev
damage_model_p
damage_model_p_stdev
damage_model_pmin
damage_model_pmin_stdev
damage_model_pmax
damage_model_pmax_stdev
pvalue
qvalue
RMSE
nb_reads_aligned
coverage
reflen
```

Position-specific damage columns are detected dynamically, including:

```text
CtoT-0, CtoT-1, ...
GtoA-0, GtoA-1, ...
```

The number of positions is not hard-coded.

Step 8.2 itself does not apply downstream biological thresholds such as a minimum predicted accuracy or minimum contig length.

### Step 8.3 — Kraken2 contig classification

Script:

```text
step8.3_kraken2.sl
```

Kraken2 classifies the MEGAHIT contigs.

Submit:

```bash
sbatch --array=1-N%M step8.3_kraken2.sl
```

Outputs are stored separately from the assembly:

```text
output/out.kraken/
```

For each sample, the workflow writes the raw Kraken2 output, a Kraken2 report and a normalized TSV used by Step 8.4.

The normalized table contains fields such as:

```text
sample_id
kraken_status
contig_id
kraken_taxid
sequence_length
kraken_hitlist
```

`kraken_status` indicates whether Kraken2 classified the sequence (`C`) or left it unclassified (`U`).

### Step 8.4 — Merge pyDamage, Kraken2, NCBI lineage and metadata

Script:

```text
step8.4_pydamage_annotation_merge.sl
```

R script:

```text
R_script/pydamage_annotation_merge.R
```

Submit once after all Step 8.2 and Step 8.3 array tasks are complete:

```bash
sbatch step8.4_pydamage_annotation_merge.sl
```

This step:

1. reads all `pydamage_results.csv` files,
2. reads all normalized Kraken2 tables,
3. downloads or reuses NCBI ranked-lineage resources,
4. resolves deprecated Kraken TaxIds using `merged.dmp`,
5. adds NCBI taxonomic lineage,
6. joins pyDamage and taxonomy using `sample_id + contig_id`,
7. optionally adds sample metadata,
8. writes a final contig-level table and QC reports.

pyDamage is the primary table. A pyDamage contig is retained even when Kraken2 did not assign taxonomy.

Possible TaxId states include:

```text
current
merged
unclassified
unresolved
no_kraken_match
```

Main outputs:

```text
output/out.pydamage_annotated/pydamage_annotation_full.tsv
output/out.pydamage_annotated/pydamage_annotation_results.rds
```

QC:

```text
output/out.pydamage_annotated/qc/contig_overlap_report.tsv
output/out.pydamage_annotated/qc/unmatched_pydamage_contigs.tsv
output/out.pydamage_annotated/qc/unmatched_kraken_contigs.tsv
output/out.pydamage_annotated/qc/unresolved_taxids.tsv
```

---

## Step 9 — Final gene annotation merge

Script:

```text
step9_annotation_merge.sl
```

R script:

```text
R_script/annotation_merge.R
```

Step 9 integrates:

```text
Salmon raw abundance
        +
eggNOG functional annotation
        +
MMseqs2 taxonomy
        +
optional metadata
```

Submit once:

```bash
sbatch step9_annotation_merge.sl
```

Default inputs are:

```text
output/out.salmon_merged_paired/output_all_raw_quant.sf
output/out.eggnog/non_redundant_prot_catalog_eggNOG.emapper.annotations
output/out.mmseqs2/non_redundant_prot_catalog.faa.result.tsv
```

Step 9 performs a **strict inner merge** on gene IDs.

A gene must therefore be present in:

1. the abundance table,
2. eggNOG annotation,
3. MMseqs2 taxonomy

to appear in the final integrated annotation table.

Main outputs:

```text
output/out.annotation_merged/full_annotation_raw.tsv
output/out.annotation_merged/full_annotation_relative.tsv
output/out.annotation_merged/functional_gap_results.rds
```

QC outputs:

```text
output/out.annotation_merged/qc/ID_overlap_report.tsv
output/out.annotation_merged/qc/ID_mismatch_report.tsv
```

`functional_gap_results.rds` contains the integrated R objects for convenient downstream use.

---

## Metadata

Metadata is optional for both Step 8.4 and Step 9.

An example is provided:

```text
examples/metadata_example.tsv
```

The only required column name is:

```text
sample_id
```

All additional columns are user-defined.

Example:

```text
sample_id\tsite\tage_yr_bp\tperiod\tenvironment
sample01\tSite_A\t8500\tHolocene\tterrestrial
sample02\tSite_A\t25000\tGlacial\tterrestrial
sample03\tSite_B\t145000\tMIS5\tmarine
```

Important rules:

- `sample_id` values must be unique.
- Empty sample IDs are not allowed.
- IDs must match the sample names produced by the pipeline.
- Do not rename samples to ages or other metadata values.
- Additional columns such as `site`, `core`, `age`, `depth`, `period`, `environment`, `location` or treatment information can be added freely.

By default, the merge scripts use:

```bash
METADATA_FILE="NONE"
```

To use metadata, change the relevant SLURM script to:

```bash
METADATA_FILE="${WORK}/metadata.tsv"
```

or another valid path.

For Step 9, metadata sample IDs must exactly match the sample columns in the Salmon abundance table.

For Step 8.4, metadata sample IDs must exactly match the pyDamage sample set.

---

## Main outputs

For most downstream work, the most important files are:

### Gene abundance

```text
output/out.salmon_merged_paired/output_all_raw_quant.sf
output/out.salmon_merged_paired/output_all_tpm_quant.sf
```

Use raw `NumReads` for count-based statistical analyses and TPM for descriptive length-normalized gene abundance.

### Functional annotation

```text
output/out.eggnog/non_redundant_prot_catalog_eggNOG.emapper.annotations
```

### Gene taxonomy

```text
output/out.mmseqs2/non_redundant_prot_catalog.faa.result.tsv
```

### Final integrated gene table

```text
output/out.annotation_merged/full_annotation_raw.tsv
output/out.annotation_merged/full_annotation_relative.tsv
```

### Final pyDamage/taxonomy contig table

```text
output/out.pydamage_annotated/pydamage_annotation_full.tsv
```

---

## Output directory structure

A completed run will approximately look like:

```text
output/
├── out.fastqc/
│   ├── pre/
│   └── post/
├── out.dedup/
├── out.fastp/
├── out.tadpole/
├── out.megahit/
│   └── <sample>/
│       └── final.contigs.fa
├── out.prodigal/
│   ├── <sample>.faa
│   ├── <sample>.fna
│   ├── <sample>.gff
│   ├── non_redundant_prot_catalog.faa
│   └── non_redundant_pCDS_catalog.fna
├── out.salmon/
│   ├── non_redundant_pCDS_catalog.index/
│   ├── <sample>_paired/
│   └── <sample>_merged/
├── out.salmon_merged_paired/
│   ├── <sample>/
│   │   └── quant.sf
│   ├── output_all_tpm_quant.sf
│   ├── output_all_raw_quant.sf
│   └── output_gene_quant.raw.count.len
├── out.eggnog/
├── out.mmseqs2/
├── out.bowtie/
├── out.pydamage/
├── out.kraken/
├── out.pydamage_annotated/
│   ├── pydamage_annotation_full.tsv
│   ├── pydamage_annotation_results.rds
│   └── qc/
└── out.annotation_merged/
    ├── full_annotation_raw.tsv
    ├── full_annotation_relative.tsv
    ├── functional_gap_results.rds
    └── qc/
```

NCBI taxonomy files used by Step 8.4 are stored outside `output/`:

```text
resources/ncbi_taxonomy/
├── rankedlineage.dmp
└── merged.dmp
```

---

## Important interpretation notes

### Salmon raw counts and TPM are not interchangeable

The workflow retains both Salmon `NumReads` and TPM.

```text
output_all_raw_quant.sf
```

contains raw Salmon estimated read counts and is the default abundance input for Step 9.

```text
output_all_tpm_quant.sf
```

contains length- and library-size-normalized abundance.

Although the targets are genes/pCDS rather than transcripts, TPM is used here as a normalized gene-catalog abundance measure because target length differs among genes.

Do not call these TPM values CPM.

For DESeq2 or other count-based statistical models, use raw `NumReads`.

### Step 9 relative abundance is different from Salmon TPM

`full_annotation_relative.tsv` is calculated from the raw gene abundance table.

For each sample, abundance is first converted to percentage across **all quantified genes**. Functional and taxonomic annotation filtering is applied afterwards.

Therefore, after the strict annotation merge, the sample columns in:

```text
full_annotation_relative.tsv
```

can sum to **less than 100%**.

This is expected and reflects genes that were quantified but did not survive the strict annotation intersection.

### Step 9 uses a strict gene intersection

The final integrated gene table contains only genes shared by:

```text
Salmon abundance
∩ eggNOG annotation
∩ MMseqs2 taxonomy
```

Inspect:

```text
qc/ID_overlap_report.tsv
qc/ID_mismatch_report.tsv
```

to understand how many genes were retained or lost at each integration stage.

### Step 8.4 does not define an authentication threshold

The full pyDamage table is retained.

The workflow does not automatically enforce thresholds such as:

```text
predicted_accuracy >= 0.5
reflen >= 600
```

or a specific taxonomic group.

These are downstream analytical decisions and should be selected according to the study design.

---

## Quick QC checks

### Check the number of raw read pairs

For example:

```bash
find /path/to/raw_reads \
    -maxdepth 1 \
    -name '*_R1_001.fastq.gz' \
    | wc -l
```

The number should match the SLURM array range.

### Check compressed FASTQ integrity

```bash
gzip -t /path/to/raw_reads/*.fastq.gz
```

No output generally indicates that the gzip files passed the integrity check.

### Check that all assemblies exist

```bash
find output/out.megahit \
    -name 'final.contigs.fa' \
    -type f \
    | sort
```

### Check the main non-redundant catalog

```bash
ls -lh \
    output/out.prodigal/non_redundant_prot_catalog.faa \
    output/out.prodigal/non_redundant_pCDS_catalog.fna
```

### Check the Salmon index

```bash
ls output/out.salmon/non_redundant_pCDS_catalog.index/versionInfo.json
```

### Check final abundance files

```bash
ls -lh \
    output/out.salmon_merged_paired/output_all_raw_quant.sf \
    output/out.salmon_merged_paired/output_all_tpm_quant.sf
```

### Check pyDamage/NCBI integration QC

```bash
column -t output/out.pydamage_annotated/qc/contig_overlap_report.tsv
```

Also inspect:

```text
unmatched_pydamage_contigs.tsv
unmatched_kraken_contigs.tsv
unresolved_taxids.tsv
```

### Check final gene integration QC

```bash
column -t output/out.annotation_merged/qc/ID_overlap_report.tsv
```

---

## Troubleshooting

### `SLURM_ARRAY_TASK_ID is not set`

The script is a sample-level array job but was submitted without an array.

Use:

```bash
sbatch --array=1-N%M script_name.sl
```

### Array task is outside the sample range

The array contains more task IDs than matching samples. Count the input samples and submit the correct range.

### A database cannot be found

Check the cluster-specific database variables in:

```text
step6_eggnog_annotation.sl
step7_mmseqs2_taxonomy.sl
step8.3_kraken2.sl
```

These database paths are system-dependent.

### NCBI taxonomy download fails in Step 8.4

The compute node may not have outbound internet access.

Manually place:

```text
rankedlineage.dmp
merged.dmp
```

under:

```text
resources/ncbi_taxonomy/
```

and rerun Step 8.4.

### Metadata/sample mismatch

Both merge workflows deliberately stop when metadata sample IDs do not match the pipeline samples.

The required column header is:

```text
sample_id
```

Check the first metadata column with:

```bash
cut -f1 metadata.tsv
```

### Step 9 contains fewer genes than the abundance table

This is normally caused by the strict inner merge.

Inspect:

```text
output/out.annotation_merged/qc/ID_overlap_report.tsv
output/out.annotation_merged/qc/ID_mismatch_report.tsv
```

Genes without both eggNOG and MMseqs2 annotation are not included in the final integrated table.

### pyDamage and Kraken2 contigs do not all overlap

This does not necessarily indicate an error.

Inspect:

```text
output/out.pydamage_annotated/qc/contig_overlap_report.tsv
output/out.pydamage_annotated/qc/unmatched_pydamage_contigs.tsv
output/out.pydamage_annotated/qc/unmatched_kraken_contigs.tsv
```

The final Step 8.4 table keeps all pyDamage contigs and adds taxonomy where available.

---

## Recommended run order

For a complete run:

```text
1.  step1_preprocessing.sl
2.  step2_de_novo_assembly.sl
3.  step3_prodigal_gene_pred.sl
4.  step4_cd_hit.sl
5.  step5.1_salmon_index.sl
6.  step5.2_abundance_salmon.sl
7.  step5.3_salmon_merge.sl
8.  step6_eggnog_annotation.sl
9.  step7_mmseqs2_taxonomy.sl
10. step8.1_bowtie_mapping.sl
11. step8.2_pydamage.sl
12. step8.3_kraken2.sl
13. step8.4_pydamage_annotation_merge.sl
14. step9_annotation_merge.sl
```

Step 8 can also be run independently after Step 2 while the gene-catalog branch is running.

### Minimal complete example

Assume there are 38 samples and at most 6 array jobs should run at the same time.

```bash
cd /path/to/my_project/Functional_GAP_pipelines

sbatch --array=1-38%6 step1_preprocessing.sl
# wait for successful completion

sbatch --array=1-38%6 step2_de_novo_assembly.sl
# wait for successful completion

sbatch --array=1-38%6 step3_prodigal_gene_pred.sl
# wait for successful completion

sbatch step4_cd_hit.sl
# wait for successful completion

sbatch step5.1_salmon_index.sl
# wait for successful completion

sbatch --array=1-38%6 step5.2_abundance_salmon.sl
# wait for successful completion

sbatch step5.3_salmon_merge.sl
sbatch step6_eggnog_annotation.sl
sbatch step7_mmseqs2_taxonomy.sl
# wait for all three required gene-level products before Step 9

sbatch --array=1-38%6 step8.1_bowtie_mapping.sl
# wait for successful completion

sbatch --array=1-38%6 step8.2_pydamage.sl
sbatch --array=1-38%6 step8.3_kraken2.sl
# wait for both arrays before Step 8.4

sbatch step8.4_pydamage_annotation_merge.sl
sbatch step9_annotation_merge.sl
```

**Do not submit dependent steps all at once unless you explicitly add SLURM dependencies.**

---

## Repository helper files

```text
Python_script/sum_up_qc_merged_paired.py
```

Combines paired and merged Salmon quantification results before cross-sample merging.

```text
R_script/pydamage_annotation_merge.R
```

Creates the final contig-level pyDamage + Kraken2 + NCBI taxonomy dataset.

```text
R_script/annotation_merge.R
```

Creates the final gene-level abundance + functional annotation + taxonomy dataset.

```text
environment/r_environment.yml
```

Defines the R environment needed for the merge scripts.

```text
examples/metadata_example.tsv
```

Provides a template for optional metadata.

---

## Notes for shared HPC installations

If this repository is placed in a common server directory for multiple users:

1. Keep the shared copy clean and preferably read-only.
2. Let each user copy the repository into their own project/work directory.
3. Run jobs from the user's project copy, not from the shared installation directory.
4. Keep large runtime resources such as NCBI taxdump files outside version control.
5. Check database paths and SLURM account/partition settings before the first run.
6. Use one working directory per dataset to avoid accidentally combining samples from different projects.

This keeps the shared installation reproducible while allowing each project to maintain its own outputs, metadata and cluster-specific settings.
