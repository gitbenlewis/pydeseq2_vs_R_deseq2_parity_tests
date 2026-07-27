#!/usr/bin/env Rscript

usage <- function() {
    cat(paste(
        "Usage: Rscript scripts/run_r_reference.R",
        "--dataset srp254919|pasilla|pickrell",
        "--mode tximport|matrix",
        "--cache-dir DIR",
        "--design FORMULA",
        "--factor NAME",
        "--reference-level LEVEL",
        "--test-level LEVEL",
        "--contrast-factor NAME",
        "--contrast-numerator LEVEL",
        "--contrast-denominator LEVEL",
        "--alpha NUMBER",
        "--fit-type parametric",
        "--size-factor-fit-type ratio",
        "--refit-cooks true|false",
        "--n-cpus 1",
        "--output-dir DIR",
        "--expected-r-version VERSION",
        "--expected-deseq2-version VERSION",
        "--expected-pasilla-version VERSION",
        "--expected-tweedeseqcountdata-version VERSION",
        "[--counts PATH --samples PATH --sample-column NAME",
        "--samples-sep auto|csv|tsv --lengths PATH --prepare-only true|false]\n"
    ))
}

parse_args <- function(args) {
    defaults <- list(
        sample_column = "sample",
        samples_sep = "auto",
        prepare_only = "false"
    )
    allowed <- c(
        "dataset", "mode", "cache_dir", "counts", "samples", "sample_column",
        "samples_sep", "lengths", "design", "factor", "reference_level",
        "test_level", "contrast_factor", "contrast_numerator",
        "contrast_denominator", "alpha", "fit_type", "size_factor_fit_type",
        "refit_cooks", "n_cpus", "output_dir", "expected_r_version",
        "expected_deseq2_version", "expected_pasilla_version",
        "expected_tweedeseqcountdata_version", "prepare_only"
    )
    opts <- defaults
    seen <- character()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (key %in% c("--help", "-h")) {
            usage()
            quit(status = 0L)
        }
        if (!startsWith(key, "--") || i == length(args)) {
            stop("Arguments must be provided as --name VALUE pairs; invalid argument: ", key)
        }
        name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
        if (!name %in% allowed) {
            stop("Unknown argument: ", key)
        }
        if (name %in% seen) {
            stop("Argument provided more than once: ", key)
        }
        opts[[name]] <- args[[i + 1L]]
        seen <- c(seen, name)
        i <- i + 2L
    }

    required <- c(
        "dataset", "mode", "cache_dir", "design", "factor", "reference_level",
        "test_level", "contrast_factor", "contrast_numerator",
        "contrast_denominator", "alpha", "fit_type", "size_factor_fit_type",
        "refit_cooks", "n_cpus", "output_dir", "expected_r_version",
        "expected_deseq2_version", "expected_pasilla_version",
        "expected_tweedeseqcountdata_version"
    )
    missing <- required[!vapply(required, function(name) {
        !is.null(opts[[name]]) && nzchar(opts[[name]])
    }, logical(1))]
    if (length(missing) > 0L) {
        stop("Missing required arguments: --", paste(gsub("_", "-", missing), collapse = ", --"))
    }
    opts
}

parse_boolean <- function(value, name) {
    normalized <- tolower(trimws(value))
    if (!normalized %in% c("true", "false")) {
        stop("--", name, " must be true or false")
    }
    identical(normalized, "true")
}

require_installed <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Missing required R package: ", package, call. = FALSE)
    }
}

format_tsv_column <- function(column) {
    if (is.factor(column)) {
        column <- as.character(column)
    }
    if (is.integer(column)) {
        output <- as.character(column)
    } else if (is.numeric(column)) {
        output <- sprintf("%.17g", column)
    } else {
        output <- as.character(column)
    }
    output[is.na(column)] <- "NA"
    output
}

write_tsv <- function(value, path) {
    table <- as.data.frame(value, check.names = FALSE, stringsAsFactors = FALSE)
    serialized <- as.data.frame(
        lapply(table, format_tsv_column),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    utils::write.table(
        serialized,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "NA"
    )
}

read_samples <- function(path, sample_column, separator = "auto") {
    if (!file.exists(path)) {
        stop("Sample metadata file does not exist: ", path)
    }
    if (separator == "auto") {
        separator <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "csv" else "tsv"
    }
    if (!separator %in% c("csv", "tsv")) {
        stop("--samples-sep must be one of: auto, csv, tsv")
    }
    samples <- if (separator == "csv") {
        utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
        utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    }
    if (!sample_column %in% colnames(samples)) {
        stop("Sample metadata is missing the sample ID column '", sample_column, "'")
    }
    sample_ids <- as.character(samples[[sample_column]])
    if (
        length(sample_ids) == 0L || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
            anyDuplicated(sample_ids)
    ) {
        stop("Sample IDs must be non-empty, non-missing, and unique")
    }
    rownames(samples) <- sample_ids
    samples
}

read_gene_matrix <- function(path, sample_ids, label) {
    if (!file.exists(path)) {
        stop(label, " file does not exist: ", path)
    }
    table <- utils::read.delim(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if (!"gene_id" %in% colnames(table)) {
        stop(label, " table is missing the 'gene_id' column")
    }
    gene_ids <- as.character(table[["gene_id"]])
    if (
        length(gene_ids) == 0L || anyNA(gene_ids) || any(!nzchar(gene_ids)) ||
            anyDuplicated(gene_ids)
    ) {
        stop(label, " gene IDs must be non-empty, non-missing, and unique")
    }
    missing_samples <- setdiff(sample_ids, colnames(table))
    if (length(missing_samples) > 0L) {
        stop(
            label, " table is missing sample columns: ",
            paste(missing_samples, collapse = ", ")
        )
    }
    values <- table[, sample_ids, drop = FALSE]
    if (!all(vapply(values, is.numeric, logical(1)))) {
        stop(label, " sample columns must all be numeric")
    }
    matrix <- as.matrix(values)
    rownames(matrix) <- gene_ids
    if (any(!is.finite(matrix)) || any(matrix < 0)) {
        stop(label, " values must be finite and non-negative")
    }
    matrix
}

normalise_gender <- function(values) {
    gender <- tolower(trimws(as.character(values)))
    gender[gender %in% c("m", "man", "male")] <- "male"
    gender[gender %in% c("f", "woman", "female")] <- "female"
    if (any(!gender %in% c("female", "male"))) {
        stop(
            "Pickrell gender metadata contains unsupported labels: ",
            paste(sort(unique(gender)), collapse = ", ")
        )
    }
    factor(gender, levels = c("female", "male"))
}

prepare_pasilla <- function(cache_dir) {
    counts_path <- system.file(
        "extdata", "pasilla_gene_counts.tsv.gz",
        package = "DESeq2", mustWork = TRUE
    )
    samples_path <- system.file(
        "extdata", "pasilla_sample_annotation.csv",
        package = "DESeq2", mustWork = TRUE
    )
    counts <- as.matrix(utils::read.csv(
        counts_path,
        sep = "\t",
        row.names = "gene_id",
        check.names = FALSE
    ))
    samples <- utils::read.csv(
        samples_path,
        row.names = 1L,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    samples <- samples[, c("condition", "type"), drop = FALSE]
    rownames(samples) <- sub("fb$", "", rownames(samples))
    if (!all(rownames(samples) %in% colnames(counts))) {
        stop("Pasilla sample labels do not match the count matrix")
    }
    counts <- counts[, rownames(samples), drop = FALSE]
    counts <- counts[rowSums(counts >= 10) >= 3, , drop = FALSE]
    if (!identical(dim(counts), c(8148L, 7L))) {
        stop(
            "Pinned Pasilla preparation produced ", nrow(counts), " genes and ",
            ncol(counts), " samples; expected 8148 genes and 7 samples"
        )
    }
    samples_out <- data.frame(
        sample = rownames(samples),
        condition = as.character(samples$condition),
        type = as.character(samples$type),
        check.names = FALSE
    )
    write_tsv(
        data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
        file.path(cache_dir, "prepared_counts.tsv")
    )
    write_tsv(samples_out, file.path(cache_dir, "prepared_samples.tsv"))
    list(counts = counts, samples = samples_out)
}

prepare_pickrell <- function(cache_dir) {
    for (package in c("Biobase", "tweeDEseqCountData")) {
        require_installed(package)
    }
    data("pickrell", package = "tweeDEseqCountData", envir = environment())
    if (!exists("pickrell.eset", inherits = FALSE)) {
        stop("The tweeDEseqCountData package did not provide pickrell.eset")
    }

    counts_all <- as.matrix(Biobase::exprs(pickrell.eset))
    phenotype <- Biobase::pData(pickrell.eset)
    gender_column <- intersect(c("gender", "Gender", "sex", "Sex"), colnames(phenotype))
    if (length(gender_column) != 1L) {
        stop("Expected exactly one gender or sex column in Pickrell metadata")
    }
    if (!identical(colnames(counts_all), rownames(phenotype))) {
        stop("Pickrell count columns and metadata rows are not identically ordered")
    }
    if (
        any(!is.finite(counts_all)) || any(counts_all < 0) ||
            any(counts_all != round(counts_all))
    ) {
        stop("pickrell.eset must contain finite non-negative integer raw counts")
    }

    counts <- counts_all[rowSums(counts_all) > 0, , drop = FALSE]
    if (!identical(dim(counts), c(12531L, 69L))) {
        stop(
            "Pinned Pickrell preparation produced ", nrow(counts), " genes and ",
            ncol(counts), " samples; expected 12531 genes and 69 samples"
        )
    }
    samples_out <- data.frame(
        sample = colnames(counts),
        gender = as.character(normalise_gender(phenotype[[gender_column]])),
        check.names = FALSE
    )
    write_tsv(
        data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
        file.path(cache_dir, "prepared_counts.tsv")
    )
    write_tsv(samples_out, file.path(cache_dir, "prepared_samples.tsv"))
    list(counts = counts, samples = samples_out)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
options(digits = 17, scipen = 0)
invisible(Sys.setlocale("LC_NUMERIC", "C"))

if (!opts$dataset %in% c("srp254919", "pasilla", "pickrell")) {
    stop("--dataset must be one of: srp254919, pasilla, pickrell")
}
expected_mode <- if (opts$dataset == "srp254919") "tximport" else "matrix"
if (!identical(opts$mode, expected_mode)) {
    stop("Dataset ", opts$dataset, " requires --mode ", expected_mode)
}
alpha <- suppressWarnings(as.numeric(opts$alpha))
if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("--alpha must be a finite number strictly between zero and one")
}
prepare_only <- parse_boolean(opts$prepare_only, "prepare-only")
refit_cooks <- parse_boolean(opts$refit_cooks, "refit-cooks")
n_cpus <- suppressWarnings(as.integer(opts$n_cpus))
if (length(n_cpus) != 1L || is.na(n_cpus) || n_cpus != 1L) {
    stop("--n-cpus must be 1 for this parity suite")
}
if (!identical(opts$fit_type, "parametric")) {
    stop("--fit-type must be parametric for this parity suite")
}
if (!identical(opts$size_factor_fit_type, "ratio")) {
    stop("--size-factor-fit-type must be ratio for this parity suite")
}
min_replicates_for_replace <- if (refit_cooks) 7L else Inf

thread_variables <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS"
)
if (any(Sys.getenv(thread_variables, unset = "") != "1")) {
    stop(
        "Set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, and ",
        "VECLIB_MAXIMUM_THREADS to 1 before launching Rscript"
    )
}

for (package in c(
    "DESeq2", "Biobase", "BiocParallel", "pasilla", "tweeDEseqCountData"
)) {
    require_installed(package)
}
actual_r_version <- as.character(getRversion())
actual_deseq2_version <- as.character(utils::packageVersion("DESeq2"))
actual_biobase_version <- as.character(utils::packageVersion("Biobase"))
actual_biocparallel_version <- as.character(utils::packageVersion("BiocParallel"))
actual_pasilla_version <- as.character(utils::packageVersion("pasilla"))
actual_tweedeseqcountdata_version <- as.character(
    utils::packageVersion("tweeDEseqCountData")
)
if (!identical(actual_r_version, opts$expected_r_version)) {
    stop(
        "R version mismatch: expected ", opts$expected_r_version,
        ", found ", actual_r_version
    )
}
if (!identical(actual_deseq2_version, opts$expected_deseq2_version)) {
    stop(
        "DESeq2 version mismatch: expected ", opts$expected_deseq2_version,
        ", found ", actual_deseq2_version
    )
}
if (!identical(actual_pasilla_version, opts$expected_pasilla_version)) {
    stop(
        "pasilla version mismatch: expected ", opts$expected_pasilla_version,
        ", found ", actual_pasilla_version
    )
}
if (!identical(
    actual_tweedeseqcountdata_version,
    opts$expected_tweedeseqcountdata_version
)) {
    stop(
        "tweeDEseqCountData version mismatch: expected ",
        opts$expected_tweedeseqcountdata_version,
        ", found ", actual_tweedeseqcountdata_version
    )
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opts$cache_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(opts$output_dir, mustWork = TRUE)
dataset_cache_dir <- normalizePath(opts$cache_dir, mustWork = TRUE)

prepared_counts_path <- ""
prepared_samples_path <- ""
input_lengths_path <- ""
if (opts$dataset == "srp254919") {
    for (name in c("counts", "samples", "lengths")) {
        if (is.null(opts[[name]]) || !nzchar(opts[[name]])) {
            stop("Dataset srp254919 requires --", name)
        }
    }
    samples <- read_samples(opts$samples, opts$sample_column, opts$samples_sep)
    sample_ids <- rownames(samples)
    counts <- read_gene_matrix(opts$counts, sample_ids, "Count")
    lengths <- read_gene_matrix(opts$lengths, sample_ids, "Transcript-length")
    if (!identical(rownames(counts), rownames(lengths))) {
        stop("Count and transcript-length gene labels are not identically ordered")
    }
    input_counts_path <- normalizePath(opts$counts, mustWork = TRUE)
    input_samples_path <- normalizePath(opts$samples, mustWork = TRUE)
    input_lengths_path <- normalizePath(opts$lengths, mustWork = TRUE)
} else {
    has_counts <- !is.null(opts[["counts"]]) && nzchar(opts[["counts"]])
    has_samples <- !is.null(opts[["samples"]]) && nzchar(opts[["samples"]])
    if (xor(has_counts, has_samples)) {
        stop("Matrix datasets require both --counts and --samples, or neither")
    }
    if (has_counts) {
        samples <- read_samples(
            opts$samples,
            opts$sample_column,
            opts$samples_sep
        )
        sample_ids <- rownames(samples)
        counts <- read_gene_matrix(opts$counts, sample_ids, "Count")
        prepared_counts_path <- normalizePath(opts$counts, mustWork = TRUE)
        prepared_samples_path <- normalizePath(opts$samples, mustWork = TRUE)
    } else {
        prepared <- if (opts$dataset == "pasilla") {
            prepare_pasilla(dataset_cache_dir)
        } else {
            prepare_pickrell(dataset_cache_dir)
        }
        counts <- prepared$counts
        samples <- prepared$samples
        rownames(samples) <- as.character(samples$sample)
        sample_ids <- rownames(samples)
        prepared_counts_path <- file.path(dataset_cache_dir, "prepared_counts.tsv")
        prepared_samples_path <- file.path(dataset_cache_dir, "prepared_samples.tsv")
    }
    input_counts_path <- prepared_counts_path
    input_samples_path <- prepared_samples_path
}

if (prepare_only) {
    if (opts$dataset == "srp254919") {
        stop("--prepare-only is only valid for generated matrix datasets")
    }
    utils::capture.output(
        utils::sessionInfo(),
        file = file.path(output_dir, "r_preparation_session_info.txt")
    )
    quit(save = "no", status = 0L)
}

if (!identical(colnames(counts), sample_ids)) {
    stop("Count columns and sample metadata rows are not identically ordered")
}
if (!opts$factor %in% colnames(samples)) {
    stop("Sample metadata is missing factor column '", opts$factor, "'")
}
if (!identical(opts$factor, opts$contrast_factor)) {
    stop("--contrast-factor must match --factor for this single-factor parity suite")
}
observed_levels <- as.character(samples[[opts$factor]])
allowed_levels <- c(opts$reference_level, opts$test_level)
if (
    identical(opts$reference_level, opts$test_level) ||
        anyNA(observed_levels) || any(!observed_levels %in% allowed_levels) ||
        !all(allowed_levels %in% observed_levels)
) {
    stop(
        "Factor '", opts$factor, "' must contain both and only the configured ",
        "reference and test levels"
    )
}
if (
    !identical(opts$contrast_numerator, opts$test_level) ||
        !identical(opts$contrast_denominator, opts$reference_level)
) {
    stop(
        "Contrast direction must be configured test level versus reference level"
    )
}

samples[[opts$factor]] <- factor(
    observed_levels,
    levels = c(opts$reference_level, opts$test_level)
)
design <- tryCatch(
    stats::as.formula(opts$design),
    error = function(error) {
        stop("Invalid design formula: ", conditionMessage(error))
    }
)
design_variables <- all.vars(design)
if (!opts$factor %in% design_variables) {
    stop("Design formula does not contain configured factor '", opts$factor, "'")
}
missing_design_variables <- setdiff(design_variables, colnames(samples))
if (length(missing_design_variables) > 0L) {
    stop(
        "Design formula references missing sample columns: ",
        paste(missing_design_variables, collapse = ", ")
    )
}

coldata <- samples
rownames(coldata) <- sample_ids
source_counts <- counts
construction_start <- proc.time()[["elapsed"]]
if (opts$mode == "tximport") {
    txi <- list(
        counts = counts,
        length = lengths,
        countsFromAbundance = "no"
    )
    dds <- DESeq2::DESeqDataSetFromTximport(
        txi,
        colData = coldata,
        design = design
    )
} else {
    if (any(counts != round(counts)) || max(counts) > .Machine$integer.max) {
        stop("Matrix-mode counts must be non-negative integers representable by R")
    }
    storage.mode(counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = design
    )
}
construction_seconds <- proc.time()[["elapsed"]] - construction_start

fit_start <- proc.time()[["elapsed"]]
dds <- DESeq2::DESeq(
    dds,
    test = "Wald",
    fitType = opts$fit_type,
    sfType = opts$size_factor_fit_type,
    minReplicatesForReplace = min_replicates_for_replace,
    parallel = FALSE,
    quiet = TRUE
)
fit_seconds <- proc.time()[["elapsed"]] - fit_start
observed_dispersion_fit_type <- as.character(
    attr(DESeq2::dispersionFunction(dds), "fitType")
)
if (!identical(observed_dispersion_fit_type, opts$fit_type)) {
    stop(
        "R DESeq2 used dispersion fit ", observed_dispersion_fit_type,
        "; requested ", opts$fit_type
    )
}

results_start <- proc.time()[["elapsed"]]
result <- DESeq2::results(
    dds,
    contrast = c(
        opts$contrast_factor,
        opts$contrast_numerator,
        opts$contrast_denominator
    ),
    alpha = alpha,
    independentFiltering = TRUE,
    cooksCutoff = TRUE,
    parallel = FALSE
)
results_seconds <- proc.time()[["elapsed"]] - results_start

rounded_counts <- DESeq2::counts(dds, normalized = FALSE)
normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
if (
    !identical(rownames(rounded_counts), rownames(source_counts)) ||
        !identical(colnames(rounded_counts), colnames(source_counts))
) {
    stop("DESeq2 changed count matrix labels or order")
}

write_tsv(
    data.frame(gene_id = rownames(rounded_counts), rounded_counts, check.names = FALSE),
    file.path(output_dir, "r_rounded_counts.tsv")
)
write_tsv(
    data.frame(gene_id = rownames(normalized_counts), normalized_counts, check.names = FALSE),
    file.path(output_dir, "r_normalized_counts.tsv")
)
write_tsv(
    data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE),
    file.path(output_dir, "r_results.tsv")
)

normalization_deviation <- NA_real_
if (opts$mode == "tximport") {
    normalization_factors <- DESeq2::normalizationFactors(dds)
    normalization_deviation <- max(
        abs(exp(rowMeans(log(normalization_factors))) - 1)
    )
    write_tsv(
        data.frame(
            gene_id = rownames(normalization_factors),
            normalization_factors,
            check.names = FALSE
        ),
        file.path(output_dir, "r_normalization_factors.tsv")
    )
} else {
    size_factors <- DESeq2::sizeFactors(dds)
    write_tsv(
        data.frame(
            sample = names(size_factors),
            size_factor = unname(size_factors),
            check.names = FALSE
        ),
        file.path(output_dir, "r_size_factors.tsv")
    )
}

rounded_difference <- abs(source_counts - rounded_counts)
metadata <- data.frame(
    key = c(
        "engine",
        "r_version",
        "deseq2_version",
        "biobase_version",
        "biocparallel_version",
        "pasilla_version",
        "tweedeseqcountdata_version",
        "dataset",
        "mode",
        "genes",
        "samples",
        "design",
        "factor",
        "reference_level",
        "test_level",
        "contrast",
        "alpha",
        "dispersion_fit_type",
        "size_factor_type",
        "refit_cooks",
        "min_replicates_for_replace",
        "n_cpus",
        "parallel",
        "input_counts_path",
        "input_samples_path",
        "input_lengths_path",
        "prepared_counts_path",
        "prepared_samples_path",
        "rounded_count_cells",
        "rounded_count_max_abs_change",
        "normalization_factor_gene_geomean_max_abs_deviation",
        "construction_seconds",
        "fit_seconds",
        "results_seconds"
    ),
    value = c(
        "R_DESeq2",
        actual_r_version,
        actual_deseq2_version,
        actual_biobase_version,
        actual_biocparallel_version,
        actual_pasilla_version,
        actual_tweedeseqcountdata_version,
        opts$dataset,
        opts$mode,
        nrow(dds),
        ncol(dds),
        opts$design,
        opts$factor,
        opts$reference_level,
        opts$test_level,
        paste(
            opts$contrast_factor,
            opts$contrast_numerator,
            "vs",
            opts$contrast_denominator,
            sep = "_"
        ),
        alpha,
        observed_dispersion_fit_type,
        opts$size_factor_fit_type,
        refit_cooks,
        min_replicates_for_replace,
        n_cpus,
        FALSE,
        input_counts_path,
        input_samples_path,
        input_lengths_path,
        prepared_counts_path,
        prepared_samples_path,
        length(rounded_difference),
        max(rounded_difference),
        normalization_deviation,
        construction_seconds,
        fit_seconds,
        results_seconds
    ),
    check.names = FALSE
)
write_tsv(metadata, file.path(output_dir, "r_metadata.tsv"))
utils::capture.output(
    utils::sessionInfo(),
    file = file.path(output_dir, "r_session_info.txt")
)
