rm(list = ls())

suppressPackageStartupMessages({
    library(data.table)
})

clean_text <- function(x) {
    trimws(as.character(x))
}

stop_if_duplicate_keys <- function(dt, key_cols, label) {
    dup <- dt[, .N, by = key_cols][N > 1]
    if (nrow(dup) > 0) {
        stop(
            label, " contains duplicated ",
            paste(key_cols, collapse = " + "),
            " combinations. Example:\n",
            paste(capture.output(print(head(dup, 10))), collapse = "\n")
        )
    }
}

read_pydamage_directory <- function(pydamage_dir) {

    files <- list.files(
        pydamage_dir,
        pattern = "^pydamage_results\\.csv$",
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(files) == 0) {
        stop("No pydamage_results.csv files were found in: ", pydamage_dir)
    }

    cat("\n============================================================\n")
    cat("Reading pyDamage results\n")
    cat("============================================================\n")
    cat("Files found:", length(files), "\n")

    tables <- lapply(files, function(f) {
        x <- fread(
            f,
            header = TRUE,
            sep = ",",
            quote = "\"",
            fill = TRUE,
            data.table = TRUE
        )

        sample_id <- basename(dirname(f))

        if ("reference" %in% names(x) && !"contig_id" %in% names(x)) {
            setnames(x, "reference", "contig_id")
        }

        if (!"contig_id" %in% names(x)) {
            stop(
                "pyDamage file does not contain a 'reference' or 'contig_id' column:\n",
                f
            )
        }

        x[, sample_id := sample_id]
        x[, contig_id := clean_text(contig_id)]
        setcolorder(x, c("sample_id", "contig_id"))
        x
    })

    result <- rbindlist(tables, use.names = TRUE, fill = TRUE)

    stop_if_duplicate_keys(
        result,
        c("sample_id", "contig_id"),
        "Combined pyDamage table"
    )

    ctot_cols <- grep("^CtoT-[0-9]+$", names(result), value = TRUE)
    gtoa_cols <- grep("^GtoA-[0-9]+$", names(result), value = TRUE)

    cat("Rows:", nrow(result), "\n")
    cat("Samples:", uniqueN(result$sample_id), "\n")
    cat("C->T position columns:", length(ctot_cols), "\n")
    cat("G->A position columns:", length(gtoa_cols), "\n")

    expected_model_cols <- c(
        "predicted_accuracy",
        "damage_model_p",
        "damage_model_pmin",
        "damage_model_pmax",
        "pvalue",
        "qvalue",
        "nb_reads_aligned",
        "coverage"
    )

    missing_model_cols <- setdiff(expected_model_cols, names(result))
    if (length(missing_model_cols) > 0) {
        warning(
            "Some expected pyDamage 1.x columns were not found: ",
            paste(missing_model_cols, collapse = ", ")
        )
    }

    if (length(ctot_cols) == 0) {
        warning("No CtoT-N columns were found in the pyDamage output.")
    }

    if (length(gtoa_cols) == 0) {
        warning(
            "No GtoA-N columns were found. ",
            "This can occur with older pyDamage output or when G->A was disabled."
        )
    }

    list(
        data = result,
        files = files,
        ctot_cols = ctot_cols,
        gtoa_cols = gtoa_cols
    )
}

read_kraken_directory <- function(kraken_dir) {

    files <- list.files(
        kraken_dir,
        pattern = "\\.kraken\\.tsv$",
        recursive = FALSE,
        full.names = TRUE
    )

    if (length(files) == 0) {
        stop(
            "No normalized *.kraken.tsv files were found in: ",
            kraken_dir,
            "\nRun step8.3_kraken2.sl before this step."
        )
    }

    cat("\n============================================================\n")
    cat("Reading Kraken2 results\n")
    cat("============================================================\n")
    cat("Files found:", length(files), "\n")

    tables <- lapply(files, function(f) {
        fread(
            f,
            header = TRUE,
            sep = "\t",
            quote = "",
            fill = TRUE,
            data.table = TRUE
        )
    })

    result <- rbindlist(tables, use.names = TRUE, fill = TRUE)

    required <- c(
        "sample_id",
        "kraken_status",
        "contig_id",
        "kraken_taxid"
    )

    missing <- setdiff(required, names(result))
    if (length(missing) > 0) {
        stop(
            "Kraken2 table is missing required column(s): ",
            paste(missing, collapse = ", ")
        )
    }

    result[, sample_id := clean_text(sample_id)]
    result[, contig_id := clean_text(contig_id)]
    result[, kraken_taxid := clean_text(kraken_taxid)]

    stop_if_duplicate_keys(
        result,
        c("sample_id", "contig_id"),
        "Combined Kraken2 table"
    )

    cat("Rows:", nrow(result), "\n")
    cat("Samples:", uniqueN(result$sample_id), "\n")

    result
}

read_ranked_lineage <- function(rankedlineage_file) {

    cat("\n============================================================\n")
    cat("Reading NCBI rankedlineage.dmp\n")
    cat("============================================================\n")

    x <- fread(
        rankedlineage_file,
        sep = "|",
        header = FALSE,
        quote = "",
        fill = TRUE,
        strip.white = TRUE,
        data.table = TRUE,
        showProgress = TRUE
    )

    if (ncol(x) < 10) {
        stop(
            "rankedlineage.dmp has fewer than 10 fields. ",
            "Unexpected NCBI taxonomy format."
        )
    }

    x <- x[, 1:10]

    setnames(
        x,
        c(
            "taxid",
            "tax_name",
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
            "superkingdom"
        )
    )

    for (col in names(x)) {
        set(x, j = col, value = clean_text(x[[col]]))
    }

    x <- x[!is.na(taxid) & taxid != ""]

    if (anyDuplicated(x$taxid)) {
        warning(
            "Duplicated taxids were detected in rankedlineage.dmp. ",
            "Keeping the first row per taxid."
        )
        x <- unique(x, by = "taxid")
    }

    cat("NCBI taxids:", nrow(x), "\n")
    x
}

read_merged_taxids <- function(merged_file) {

    cat("\n============================================================\n")
    cat("Reading NCBI merged.dmp\n")
    cat("============================================================\n")

    x <- fread(
        merged_file,
        sep = "|",
        header = FALSE,
        quote = "",
        fill = TRUE,
        strip.white = TRUE,
        data.table = TRUE,
        showProgress = FALSE
    )

    if (ncol(x) < 2) {
        stop("merged.dmp has fewer than 2 fields.")
    }

    x <- x[, 1:2]
    setnames(x, c("old_taxid", "new_taxid"))

    x[, old_taxid := clean_text(old_taxid)]
    x[, new_taxid := clean_text(new_taxid)]

    x <- x[
        !is.na(old_taxid) &
        old_taxid != "" &
        !is.na(new_taxid) &
        new_taxid != ""
    ]

    if (anyDuplicated(x$old_taxid)) {
        x <- unique(x, by = "old_taxid")
    }

    cat("Merged taxid mappings:", nrow(x), "\n")
    x
}

resolve_taxids <- function(kraken, ranked, merged) {

    current_ids <- ranked$taxid
    merged_map <- setNames(merged$new_taxid, merged$old_taxid)

    original <- kraken$kraken_taxid
    resolved <- original

    is_unclassified <- is.na(original) | original == "" | original == "0"
    is_current <- !is_unclassified & original %chin% current_ids

    status <- rep("unresolved", length(original))
    status[is_unclassified] <- "unclassified"
    status[is_current] <- "current"

    needs_resolution <- !is_unclassified & !is_current

    for (iteration in seq_len(10)) {
        idx <- which(needs_resolution & !(resolved %chin% current_ids))
        if (length(idx) == 0) {
            break
        }

        mapped <- unname(merged_map[resolved[idx]])
        can_map <- !is.na(mapped) & mapped != ""

        if (!any(can_map)) {
            break
        }

        resolved[idx[can_map]] <- mapped[can_map]
    }

    became_current <- needs_resolution & resolved %chin% current_ids
    status[became_current] <- "merged"

    resolved[is_unclassified] <- NA_character_

    kraken[, kraken_taxid_original := original]
    kraken[, resolved_taxid := resolved]
    kraken[, taxid_status := status]

    kraken
}

read_metadata <- function(metadata_file, sample_ids) {

    if (toupper(metadata_file) == "NONE") {
        cat("\nNo metadata file supplied.\n")
        return(NULL)
    }

    cat("\n============================================================\n")
    cat("Reading metadata\n")
    cat("============================================================\n")

    metadata <- fread(metadata_file, data.table = TRUE)

    if (!"sample_id" %in% names(metadata)) {
        stop("Metadata must contain a column named 'sample_id'.")
    }

    metadata[, sample_id := clean_text(sample_id)]

    if (any(is.na(metadata$sample_id) | metadata$sample_id == "")) {
        stop("Metadata contains empty sample_id values.")
    }

    if (anyDuplicated(metadata$sample_id)) {
        stop("Metadata contains duplicated sample_id values.")
    }

    missing_in_metadata <- setdiff(sample_ids, metadata$sample_id)
    extra_in_metadata <- setdiff(metadata$sample_id, sample_ids)

    if (length(missing_in_metadata) > 0 || length(extra_in_metadata) > 0) {

        if (length(missing_in_metadata) > 0) {
            cat(
                "Samples in pyDamage but missing from metadata:",
                length(missing_in_metadata),
                "\n"
            )
            print(head(missing_in_metadata, 20))
        }

        if (length(extra_in_metadata) > 0) {
            cat(
                "Samples in metadata but missing from pyDamage:",
                length(extra_in_metadata),
                "\n"
            )
            print(head(extra_in_metadata, 20))
        }

        stop(
            "Metadata sample_id values must exactly match the pyDamage sample set."
        )
    }

    metadata[match(sample_ids, sample_id)]
}

build_overlap_reports <- function(pydamage, kraken, qc_dir) {

    pydamage_keys <- unique(
        pydamage[, .(sample_id, contig_id)]
    )

    kraken_keys <- unique(
        kraken[, .(sample_id, contig_id)]
    )

    pydamage_keys[, key := paste(sample_id, contig_id, sep = "\t")]
    kraken_keys[, key := paste(sample_id, contig_id, sep = "\t")]

    common_keys <- intersect(pydamage_keys$key, kraken_keys$key)

    unmatched_pydamage <- pydamage_keys[!key %chin% kraken_keys$key]
    unmatched_kraken <- kraken_keys[!key %chin% pydamage_keys$key]

    overlap_report <- data.table(
        comparison = c(
            "pydamage_contigs",
            "kraken_contigs",
            "common_contigs",
            "pydamage_without_kraken",
            "kraken_without_pydamage"
        ),
        n = c(
            nrow(pydamage_keys),
            nrow(kraken_keys),
            length(common_keys),
            nrow(unmatched_pydamage),
            nrow(unmatched_kraken)
        )
    )

    fwrite(
        overlap_report,
        file.path(qc_dir, "contig_overlap_report.tsv"),
        sep = "\t"
    )

    fwrite(
        unmatched_pydamage[, .(sample_id, contig_id)],
        file.path(qc_dir, "unmatched_pydamage_contigs.tsv"),
        sep = "\t"
    )

    fwrite(
        unmatched_kraken[, .(sample_id, contig_id)],
        file.path(qc_dir, "unmatched_kraken_contigs.tsv"),
        sep = "\t"
    )

    overlap_report
}

main <- function(
    pydamage_dir,
    kraken_dir,
    rankedlineage_file,
    merged_file,
    metadata_file,
    output_dir
) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    qc_dir <- file.path(output_dir, "qc")
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

    pydamage_result <- read_pydamage_directory(pydamage_dir)
    pydamage <- pydamage_result$data

    kraken <- read_kraken_directory(kraken_dir)

    overlap_report <- build_overlap_reports(
        pydamage = pydamage,
        kraken = kraken,
        qc_dir = qc_dir
    )

    ranked <- read_ranked_lineage(rankedlineage_file)
    merged <- read_merged_taxids(merged_file)

    kraken <- resolve_taxids(
        kraken = kraken,
        ranked = ranked,
        merged = merged
    )

    ranked_join <- copy(ranked)
    setnames(ranked_join, "taxid", "resolved_taxid")

    kraken_annotated <- merge(
        kraken,
        ranked_join,
        by = "resolved_taxid",
        all.x = TRUE,
        sort = FALSE
    )

    full <- merge(
        pydamage,
        kraken_annotated,
        by = c("sample_id", "contig_id"),
        all.x = TRUE,
        sort = FALSE
    )

    full[
        is.na(taxid_status) | taxid_status == "",
        taxid_status := "no_kraken_match"
    ]

    taxonomy_cols <- c(
        "tax_name",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "superkingdom"
    )

    for (col in taxonomy_cols) {
        if (col %in% names(full)) {
            full[
                is.na(get(col)) | get(col) == "",
                (col) := "Unclassified"
            ]
        }
    }

    sample_ids <- unique(pydamage$sample_id)

    metadata <- read_metadata(
        metadata_file = metadata_file,
        sample_ids = sample_ids
    )

    if (!is.null(metadata)) {

        overlapping_metadata_cols <- intersect(
            setdiff(names(metadata), "sample_id"),
            names(full)
        )

        if (length(overlapping_metadata_cols) > 0) {
            stop(
                "Metadata contains column name(s) already present in the merged table: ",
                paste(overlapping_metadata_cols, collapse = ", ")
            )
        }

        full <- merge(
            full,
            metadata,
            by = "sample_id",
            all.x = TRUE,
            sort = FALSE
        )
    }

    unresolved_taxids <- kraken[
        taxid_status == "unresolved",
        .(n_contigs = .N),
        by = .(kraken_taxid_original)
    ][order(-n_contigs)]

    fwrite(
        unresolved_taxids,
        file.path(qc_dir, "unresolved_taxids.tsv"),
        sep = "\t"
    )

    fwrite(
        full,
        file.path(output_dir, "pydamage_annotation_full.tsv"),
        sep = "\t"
    )

    results <- list(
        full = full,
        metadata = metadata,
        contig_overlap = overlap_report,
        unresolved_taxids = unresolved_taxids,
        pydamage_files = pydamage_result$files,
        ctot_columns = pydamage_result$ctot_cols,
        gtoa_columns = pydamage_result$gtoa_cols
    )

    saveRDS(
        results,
        file.path(output_dir, "pydamage_annotation_results.rds")
    )

    cat("\n============================================================\n")
    cat("pyDamage annotation merge finished\n")
    cat("============================================================\n")
    cat("pyDamage rows        :", nrow(pydamage), "\n")
    cat("Kraken2 rows         :", nrow(kraken), "\n")
    cat("Final rows           :", nrow(full), "\n")
    cat("Samples              :", uniqueN(full$sample_id), "\n")
    cat("C->T columns         :", length(pydamage_result$ctot_cols), "\n")
    cat("G->A columns         :", length(pydamage_result$gtoa_cols), "\n")
    cat("Current TaxIds       :", kraken[taxid_status == "current", .N], "\n")
    cat("Merged TaxIds        :", kraken[taxid_status == "merged", .N], "\n")
    cat("Unclassified TaxIds  :", kraken[taxid_status == "unclassified", .N], "\n")
    cat("Unresolved TaxIds    :", kraken[taxid_status == "unresolved", .N], "\n")
    cat("\nOutputs:\n")
    cat("  pydamage_annotation_full.tsv\n")
    cat("  pydamage_annotation_results.rds\n")
    cat("  qc/contig_overlap_report.tsv\n")
    cat("  qc/unmatched_pydamage_contigs.tsv\n")
    cat("  qc/unmatched_kraken_contigs.tsv\n")
    cat("  qc/unresolved_taxids.tsv\n")
    cat("============================================================\n")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
    stop(
        paste0(
            "\nUsage:\n",
            "Rscript pydamage_annotation_merge.R ",
            "<pydamage_dir> ",
            "<kraken_dir> ",
            "<rankedlineage.dmp> ",
            "<merged.dmp> ",
            "<metadata_file|NONE> ",
            "<output_dir>\n"
        )
    )
}

required_paths <- args[1:4]
missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0) {
    stop(
        "Required input path(s) not found:\n",
        paste(missing_paths, collapse = "\n")
    )
}

if (toupper(args[5]) != "NONE" && !file.exists(args[5])) {
    stop("Metadata file not found: ", args[5])
}

main(
    pydamage_dir = args[1],
    kraken_dir = args[2],
    rankedlineage_file = args[3],
    merged_file = args[4],
    metadata_file = args[5],
    output_dir = args[6]
)
