# Changelog

All notable changes to Functional GAP Pipelines are documented in this file.

The workflow was originally introduced in:

Çabuk U, Herzschuh U, Harms L, von Hippel B, Stoof-Leichsenring KR (2025).
*Functional annotation of eukaryotic genes from sedimentary ancient DNA.*
Frontiers in Ecology and Evolution, 13:1459690.
https://doi.org/10.3389/fevo.2025.1459690

---

## [2.0.0] - Unreleased

### Added

- Added a generic end-to-end workflow from raw paired-end shotgun metagenomic reads to final functional, taxonomic and ancient-DNA-aware output tables.
- Added a dedicated Bowtie2 mapping step for mapping processed reads back to sample-specific MEGAHIT assemblies.
- Added pyDamage 1.x support while preserving the native pyDamage output schema.
- Added dynamic handling of pyDamage C-to-T and G-to-A position-specific damage columns.
- Added support for current pyDamage model fields, including:
  - `damage_model_p`
  - `damage_model_pmin`
  - `damage_model_pmax`
  - associated standard-deviation fields
  - `RMSE`
  - `reflen`
- Added Kraken2 contig classification as a separate workflow step.
- Added normalized Kraken2 output tables for downstream integration.
- Added automated retrieval of the NCBI `new_taxdump` taxonomy archive.
- Added NCBI ranked lineage annotation using `rankedlineage.dmp`.
- Added resolution of deprecated NCBI TaxIds using `merged.dmp`.
- Added a generic contig-level integration workflow combining:
  - pyDamage
  - Kraken2
  - NCBI taxonomy
  - optional sample metadata
- Added QC reports for pyDamage/Kraken2 contig overlap and unresolved TaxIds.
- Added a generic gene-level annotation merge workflow combining:
  - Salmon abundance
  - eggNOG functional annotation
  - MMseqs2 taxonomy
  - optional sample metadata
- Added gene-ID overlap and mismatch QC reports for the final annotation merge.
- Added optional metadata support using a standardized `sample_id` column.
- Added `examples/metadata_example.tsv`.
- Added a reproducible Conda R environment in `environment/r_environment.yml`.
- Added `CITATION.cff` for software citation metadata.
- Added a complete user-oriented README covering the workflow from raw FASTQ input to final output tables.

### Changed

- Standardized the main assembly workflow around Tadpole error correction followed by MEGAHIT.
- Standardized gene prediction using Prodigal in metagenomic mode.
- Generalized sample discovery and removed dataset-specific naming assumptions from the main pipeline.
- Standardized the MEGAHIT assembly filename as:

  `final.contigs.fa`

- Standardized the non-redundant nucleotide gene catalog as:

  `non_redundant_pCDS_catalog.fna`

- Standardized the Salmon index name as:

  `non_redundant_pCDS_catalog.index`

- Updated Salmon processing so paired and merged read quantifications are combined consistently across samples.
- Renamed the normalized Salmon cross-sample abundance output from CPM to TPM to reflect the actual Salmon `quantmerge` output.
- Raw Salmon `NumReads` are retained separately for count-based downstream analyses.
- Generalized eggNOG-mapper annotation to operate directly on the non-redundant protein catalog.
- Generalized MMseqs2 taxonomy assignment to operate directly on the non-redundant protein catalog.
- Separated ancient-DNA processing into independent steps:
  - Bowtie2 mapping
  - pyDamage
  - Kraken2
  - final pyDamage/taxonomy integration
- pyDamage results are now read directly from native `pydamage_results.csv` files instead of reconstructing a fixed output schema.
- pyDamage and Kraken2 results are merged directly using `sample_id + contig_id`.
- Metadata are added only once during the final merge rather than independently to multiple intermediate tables.
- Kraken2 outputs are now stored in a dedicated:

  `output/out.kraken/`

  directory rather than alongside assembly outputs.
- The final pyDamage integration retains all pyDamage contigs even when Kraken2 taxonomy is unavailable.
- Final gene integration now uses a strict intersection of Salmon abundance, eggNOG annotation and MMseqs2 taxonomy.
- Relative gene abundance is calculated before annotation filtering so that quantified but unannotated genes remain represented in the normalization denominator.
- SLURM headers and comments were standardized across workflow scripts.
- Sample-level and global workflow steps are now clearly separated.
- Pipeline documentation was rewritten for use both from GitHub and from shared HPC installations.

### Removed

- Removed Tiara-based sequence classification from the core workflow.
- Removed MetaEuk-based gene prediction from the core workflow.
- Removed the metaSPAdes assembly branch.
- Removed BWA mapping from the ancient-DNA workflow.
- Removed dataset-specific hard-coded sample, site and core names.
- Removed dataset-specific Salmon index names.
- Removed the requirement for a manually prepared `ncbi_df.csv`.
- Removed the fixed pyDamage 0.x column definition used by the previous downstream R workflow.
- Removed project-specific downstream plotting and exploratory analysis scripts from the core repository.
- Removed project-specific CSV, RData and metadata input files.
- Removed the legacy `pydamage_workflow_v1.R`.
- Removed the legacy `ngc_calculation.py` helper.

### Fixed

- Fixed inconsistent MEGAHIT contig filename usage across Prodigal, Bowtie2 and Kraken2 steps.
- Fixed dataset-specific Salmon index path mismatches.
- Fixed the normalized Salmon abundance output being incorrectly labelled as CPM.
- Fixed pyDamage module version consistency in the SLURM workflow.
- Fixed previous assumptions about a fixed number of pyDamage C-to-T positions.
- Added support for G-to-A damage-position columns produced by pyDamage 1.x.
- Fixed duplicate intermediate catalog inclusion during repeated CD-HIT catalog construction.
- Improved handling of duplicated eggNOG and taxonomy gene IDs during final integration.
- Added explicit validation of metadata sample identifiers.
- Added explicit validation of pyDamage and Kraken2 sample/contig identifiers.
- Added QC reporting for genes and contigs that do not match across workflow components.
- Added handling of deprecated Kraken/NCBI TaxIds before lineage assignment.

---

## [1.0.0]

Initial archived version of Functional GAP Pipelines.

Zenodo:
https://doi.org/10.5281/zenodo.14810762

Associated publication:
https://doi.org/10.3389/fevo.2025.1459690
