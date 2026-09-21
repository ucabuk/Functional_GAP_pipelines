#=============================================================================#
# R-script
# author: Ugur Cabuk
# Description:
#   Read eggNOG, MMseqs2 taxonomy and Salmon gene abundance
#   Merge the full gene-level dataset
#   Save raw and relative abundance outputs
#=============================================================================#

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})


#=============================================================================#
# Safe gene ID cleaning
#=============================================================================#

clean_gene_id_safe <- function(x) {

  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace("^>", "")
}


#=============================================================================#
# Read eggNOG
#=============================================================================#

read_eggnog_clean <- function(eggnog_file) {

  cat("\n============================================================\n")
  cat("Reading eggNOG file:\n")
  cat(eggnog_file, "\n")
  cat("============================================================\n")

  eggnog_input <- fread(
    eggnog_file,
    skip = "#query",
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    data.table = TRUE
  )

  colnames(eggnog_input)[1] <- "gene_id"

  eggnog_input[, gene_id := clean_gene_id_safe(gene_id)]

  eggnog_input <- eggnog_input[
    !is.na(gene_id) &
      gene_id != "" &
      !grepl("^#", gene_id)
  ]

  cat("\neggNOG dimensions before unique:\n")
  print(dim(eggnog_input))

  eggnog_dup <- eggnog_input[, .N, by = gene_id][N > 1][order(-N)]

  cat("\nDuplicated eggNOG gene IDs:\n")
  print(nrow(eggnog_dup))

  if (nrow(eggnog_dup) > 0) {
    warning(
      nrow(eggnog_dup),
      " duplicated eggNOG gene IDs found. Keeping the first row per gene_id."
    )
  }

  eggnog_input <- unique(eggnog_input, by = "gene_id")

  cat("\neggNOG dimensions after unique:\n")
  print(dim(eggnog_input))

  return(
    list(
      data = eggnog_input,
      duplicated_ids = eggnog_dup
    )
  )
}


#=============================================================================#
# Clean taxonomy table
#=============================================================================#

clean_taxonomy <- function(taxonomy_file) {

  cat("\n============================================================\n")
  cat("Reading taxonomy file:\n")
  cat(taxonomy_file, "\n")
  cat("============================================================\n")

  taxon_NCBI <- fread(
    taxonomy_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    data.table = FALSE
  )

  if (ncol(taxon_NCBI) < 5) {
    stop(
      "Taxonomy file has fewer than 5 columns. ",
      "The MMseqs2 taxonomy output with taxonomic lineage is required."
    )
  }

  taxon_NCBI_2 <- separate(
    data = taxon_NCBI,
    col = V5,
    into = c(
      "species",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "kingdom",
      "superkingdom"
    ),
    sep = ";",
    fill = "right",
    extra = "merge"
  )

  taxon_NCBI_3 <- taxon_NCBI_2 %>%
    select(
      V1,
      V4,
      species,
      genus,
      family,
      order,
      class,
      phylum,
      kingdom,
      superkingdom
    )

  colnames(taxon_NCBI_3)[1] <- "gene_id"
  colnames(taxon_NCBI_3)[2] <- "taxonomy_assignment"

  taxon_cols <- c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom"
  )

  taxon_NCBI_3[taxon_cols] <- lapply(
    taxon_NCBI_3[taxon_cols],
    as.character
  )

  taxon_NCBI_3[taxon_NCBI_3 == ""] <- "Unclassified"

  for (tc in taxon_cols) {

    taxon_NCBI_3[[tc]] <- sub(
      "^uc_.*|unknown",
      "Unclassified",
      taxon_NCBI_3[[tc]]
    )
  }

  taxon_NCBI_4 <- replace(
    taxon_NCBI_3,
    is.na(taxon_NCBI_3),
    "Unclassified"
  )

  taxon_NCBI_4 <- as.data.table(taxon_NCBI_4)

  taxon_NCBI_4[, gene_id := clean_gene_id_safe(gene_id)]

  cat("\nTaxonomy dimensions before unique:\n")
  print(dim(taxon_NCBI_4))

  taxon_dup <- taxon_NCBI_4[, .N, by = gene_id][N > 1][order(-N)]

  cat("\nDuplicated taxonomy gene IDs:\n")
  print(nrow(taxon_dup))

  if (nrow(taxon_dup) > 0) {
    warning(
      nrow(taxon_dup),
      " duplicated taxonomy gene IDs found. Keeping the first row per gene_id."
    )
  }

  taxon_NCBI_4 <- unique(taxon_NCBI_4, by = "gene_id")

  cat("\nTaxonomy dimensions after unique:\n")
  print(dim(taxon_NCBI_4))

  return(
    list(
      data = taxon_NCBI_4,
      duplicated_ids = taxon_dup
    )
  )
}


#=============================================================================#
# Read and check metadata
#=============================================================================#

read_metadata_clean <- function(metadata_file, sample_cols) {

  if (toupper(metadata_file) == "NONE") {

    cat("\nNo metadata file supplied.\n")

    return(NULL)
  }

  cat("\n============================================================\n")
  cat("Reading metadata file:\n")
  cat(metadata_file, "\n")
  cat("============================================================\n")

  metadata_table <- fread(
    metadata_file,
    data.table = TRUE
  )

  if (!"sample_id" %in% colnames(metadata_table)) {
    stop("Metadata must contain a column named 'sample_id'.")
  }

  metadata_table[, sample_id := trimws(as.character(sample_id))]

  if (any(is.na(metadata_table$sample_id) | metadata_table$sample_id == "")) {
    stop("Metadata contains empty sample_id values.")
  }

  duplicated_samples <- metadata_table[
    duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE),
    unique(sample_id)
  ]

  if (length(duplicated_samples) > 0) {
    stop(
      "Duplicated sample_id values found in metadata: ",
      paste(head(duplicated_samples, 10), collapse = ", "),
      if (length(duplicated_samples) > 10) " ..." else ""
    )
  }

  missing_in_metadata <- setdiff(sample_cols, metadata_table$sample_id)
  missing_in_abundance <- setdiff(metadata_table$sample_id, sample_cols)

  if (length(missing_in_metadata) > 0 || length(missing_in_abundance) > 0) {

    cat("\nMetadata/sample mismatch:\n")

    if (length(missing_in_metadata) > 0) {
      cat(
        "Samples in abundance but missing from metadata:",
        length(missing_in_metadata),
        "\n"
      )
      print(head(missing_in_metadata, 20))
    }

    if (length(missing_in_abundance) > 0) {
      cat(
        "Samples in metadata but missing from abundance:",
        length(missing_in_abundance),
        "\n"
      )
      print(head(missing_in_abundance, 20))
    }

    stop(
      "Metadata sample_id values must exactly match the Salmon abundance sample columns."
    )
  }

  metadata_table <- metadata_table[
    match(sample_cols, sample_id)
  ]

  cat("\nMetadata rows:\n")
  print(nrow(metadata_table))

  return(metadata_table)
}


#=============================================================================#
# ID overlap and mismatch reports
#=============================================================================#

write_id_reports <- function(
    eggnog_input,
    quant_merge,
    taxon_clean,
    eggnog_dup,
    taxon_dup,
    qc_dir
) {

  ids_eggnog <- unique(eggnog_input$gene_id)
  ids_quant  <- unique(quant_merge$gene_id)
  ids_tax    <- unique(taxon_clean$gene_id)

  overlap_report <- data.table(
    comparison = c(
      "unique_eggNOG_IDs",
      "unique_abundance_IDs",
      "unique_taxonomy_IDs",
      "duplicated_eggNOG_IDs_before_unique",
      "duplicated_taxonomy_IDs_before_unique",
      "abundance_intersect_eggNOG",
      "abundance_intersect_taxonomy",
      "eggNOG_intersect_taxonomy",
      "abundance_intersect_eggNOG_intersect_taxonomy",
      "abundance_missing_eggNOG",
      "abundance_missing_taxonomy",
      "eggNOG_missing_abundance",
      "taxonomy_missing_abundance"
    ),
    n = c(
      length(ids_eggnog),
      length(ids_quant),
      length(ids_tax),
      nrow(eggnog_dup),
      nrow(taxon_dup),
      length(intersect(ids_quant, ids_eggnog)),
      length(intersect(ids_quant, ids_tax)),
      length(intersect(ids_eggnog, ids_tax)),
      length(Reduce(intersect, list(ids_quant, ids_eggnog, ids_tax))),
      length(setdiff(ids_quant, ids_eggnog)),
      length(setdiff(ids_quant, ids_tax)),
      length(setdiff(ids_eggnog, ids_quant)),
      length(setdiff(ids_tax, ids_quant))
    )
  )

  cat("\n================ ID OVERLAP REPORT ================\n")
  print(overlap_report)
  cat("===================================================\n\n")

  fwrite(
    overlap_report,
    file.path(qc_dir, "ID_overlap_report.tsv"),
    sep = "\t"
  )

  mismatch_ids <- unique(
    c(
      setdiff(ids_quant, ids_eggnog),
      setdiff(ids_quant, ids_tax),
      setdiff(ids_eggnog, ids_quant),
      setdiff(ids_tax, ids_quant)
    )
  )

  mismatch_report <- data.table(
    gene_id = mismatch_ids,
    in_abundance = mismatch_ids %chin% ids_quant,
    in_eggnog = mismatch_ids %chin% ids_eggnog,
    in_taxonomy = mismatch_ids %chin% ids_tax
  )

  if (nrow(mismatch_report) > 0) {
    setorder(mismatch_report, gene_id)
  }

  fwrite(
    mismatch_report,
    file.path(qc_dir, "ID_mismatch_report.tsv"),
    sep = "\t"
  )

  return(
    list(
      overlap_report = overlap_report,
      mismatch_report = mismatch_report
    )
  )
}


#=============================================================================#
# Main function
#=============================================================================#

analyze_full_annotation <- function(
    eggnog_file,
    gene_abundance_file,
    taxonomy_file,
    metadata_file = "NONE",
    output_dir
) {

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  qc_dir <- file.path(output_dir, "qc")

  dir.create(
    qc_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  cat("\n============================================================\n")
  cat("Functional GAP annotation merge\n")
  cat("============================================================\n")


  #------------------------------------------------------------#
  # 1. Read eggNOG
  #------------------------------------------------------------#

  eggnog_result <- read_eggnog_clean(eggnog_file)

  eggnog_input <- eggnog_result$data
  eggnog_dup <- eggnog_result$duplicated_ids


  #------------------------------------------------------------#
  # 2. Read gene abundance
  #------------------------------------------------------------#

  quant_merge <- fread(
    gene_abundance_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    data.table = TRUE
  )

  colnames(quant_merge)[1] <- "gene_id"

  quant_merge[, gene_id := clean_gene_id_safe(gene_id)]

  abundance_dup <- quant_merge[, .N, by = gene_id][N > 1]

  if (nrow(abundance_dup) > 0) {
    stop(
      "Duplicated gene IDs found in the abundance table. ",
      "The Salmon merged abundance table must contain one row per gene_id."
    )
  }

  sample_cols <- setdiff(
    colnames(quant_merge),
    "gene_id"
  )

  if (length(sample_cols) == 0) {
    stop("No sample columns found in the abundance file.")
  }

  quant_merge[
    ,
    (sample_cols) := lapply(.SD, as.numeric),
    .SDcols = sample_cols
  ]

  cat("\nAbundance dimensions:\n")
  print(dim(quant_merge))

  cat("\nNumber of samples:\n")
  print(length(sample_cols))


  #------------------------------------------------------------#
  # 3. Read metadata
  #------------------------------------------------------------#

  metadata_table <- read_metadata_clean(
    metadata_file = metadata_file,
    sample_cols = sample_cols
  )


  #------------------------------------------------------------#
  # 4. Read taxonomy
  #------------------------------------------------------------#

  taxonomy_result <- clean_taxonomy(taxonomy_file)

  taxon_clean <- taxonomy_result$data
  taxon_dup <- taxonomy_result$duplicated_ids


  #------------------------------------------------------------#
  # 5. ID reports
  #------------------------------------------------------------#

  id_reports <- write_id_reports(
    eggnog_input = eggnog_input,
    quant_merge = quant_merge,
    taxon_clean = taxon_clean,
    eggnog_dup = eggnog_dup,
    taxon_dup = taxon_dup,
    qc_dir = qc_dir
  )

  overlap_report <- id_reports$overlap_report


  #------------------------------------------------------------#
  # 6. RAW abundance
  #------------------------------------------------------------#

  quant_raw <- copy(quant_merge)


  #------------------------------------------------------------#
  # 7. Relative abundance
  #    Calculated against all quantified genes before annotation filtering.
  #------------------------------------------------------------#

  quant_rel <- copy(quant_merge)

  quant_rel[
    ,
    (sample_cols) := lapply(
      .SD,
      function(x) {

        total <- sum(x, na.rm = TRUE)

        if (is.na(total) || total == 0) {
          rep(0, length(x))
        } else {
          (x / total) * 100
        }
      }
    ),
    .SDcols = sample_cols
  ]

  cat("\nRelative abundance sums before annotation filtering:\n")
  print(colSums(quant_rel[, ..sample_cols], na.rm = TRUE))


  #------------------------------------------------------------#
  # 8. Merge RAW abundance + eggNOG + taxonomy
  #------------------------------------------------------------#

  cat("\n================ STRICT INNER MERGE ================\n")

  cat("Abundance rows before merge:\n")
  print(nrow(quant_raw))

  cat("eggNOG rows after unique:\n")
  print(nrow(eggnog_input))

  cat("Taxonomy rows after unique:\n")
  print(nrow(taxon_clean))

  full_raw_step1 <- merge(
    quant_raw,
    eggnog_input,
    by = "gene_id",
    all = FALSE
  )

  cat("After abundance + eggNOG:\n")
  print(nrow(full_raw_step1))

  full_raw <- merge(
    full_raw_step1,
    taxon_clean,
    by = "gene_id",
    all = FALSE
  )

  cat("After abundance + eggNOG + taxonomy:\n")
  print(nrow(full_raw))

  cat("====================================================\n")


  #------------------------------------------------------------#
  # 9. Merge RELATIVE abundance + eggNOG + taxonomy
  #------------------------------------------------------------#

  full_relative_step1 <- merge(
    quant_rel,
    eggnog_input,
    by = "gene_id",
    all = FALSE
  )

  full_relative <- merge(
    full_relative_step1,
    taxon_clean,
    by = "gene_id",
    all = FALSE
  )


  #------------------------------------------------------------#
  # 10. Replace missing taxonomy after merge
  #------------------------------------------------------------#

  taxon_cols <- c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom"
  )

  for (tc in taxon_cols) {

    if (tc %in% colnames(full_raw)) {
      full_raw[[tc]][
        is.na(full_raw[[tc]]) |
          full_raw[[tc]] == ""
      ] <- "Unclassified"
    }

    if (tc %in% colnames(full_relative)) {
      full_relative[[tc]][
        is.na(full_relative[[tc]]) |
          full_relative[[tc]] == ""
      ] <- "Unclassified"
    }
  }


  #------------------------------------------------------------#
  # 11. Save outputs
  #------------------------------------------------------------#

  fwrite(
    full_raw,
    file.path(output_dir, "full_annotation_raw.tsv"),
    sep = "\t"
  )

  fwrite(
    full_relative,
    file.path(output_dir, "full_annotation_relative.tsv"),
    sep = "\t"
  )

  results <- list(
    full_raw = full_raw,
    full_relative = full_relative,
    metadata = metadata_table,
    overlap_report = overlap_report
  )

  saveRDS(
    results,
    file.path(output_dir, "functional_gap_results.rds")
  )


  #------------------------------------------------------------#
  # 12. Final report
  #------------------------------------------------------------#

  cat("\n============================================================\n")
  cat("Functional GAP annotation merge finished\n")
  cat("============================================================\n")

  cat("Genes in abundance :", nrow(quant_merge), "\n")
  cat("Genes in eggNOG    :", nrow(eggnog_input), "\n")
  cat("Genes in taxonomy  :", nrow(taxon_clean), "\n")
  cat("Common genes       :", nrow(full_raw), "\n")
  cat("Samples            :", length(sample_cols), "\n")

  cat("\nMain outputs:\n")
  cat("  full_annotation_raw.tsv\n")
  cat("  full_annotation_relative.tsv\n")
  cat("  functional_gap_results.rds\n")

  cat("\nQC:\n")
  cat("  qc/ID_overlap_report.tsv\n")
  cat("  qc/ID_mismatch_report.tsv\n")

  cat("\nOutput directory:\n")
  cat(output_dir, "\n")

  cat("============================================================\n")

  return(results)
}


#=============================================================================#
# Command line
#=============================================================================#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {

  stop(
    paste0(
      "\nUsage:\n",
      "Rscript annotation_merge.R ",
      "<eggnog_file> ",
      "<gene_abundance_file> ",
      "<taxonomy_file> ",
      "<metadata_file|NONE> ",
      "<output_dir>\n"
    )
  )
}


eggnog_file <- args[1]
gene_abundance_file <- args[2]
taxonomy_file <- args[3]
metadata_file <- args[4]
output_dir <- args[5]


required_files <- c(
  eggnog_file,
  gene_abundance_file,
  taxonomy_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (
  toupper(metadata_file) != "NONE" &&
    !file.exists(metadata_file)
) {
  stop("Metadata file not found: ", metadata_file)
}


result <- analyze_full_annotation(
  eggnog_file = eggnog_file,
  gene_abundance_file = gene_abundance_file,
  taxonomy_file = taxonomy_file,
  metadata_file = metadata_file,
  output_dir = output_dir
)
