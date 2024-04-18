## helper function to add more info to a simulated SCE ----
add_sce_data <- function(sce, batch) {
  # add some coldata columns
  colData(sce)[["sum"]] <- runif(ncol(sce))
  colData(sce)[["detected"]] <- runif(ncol(sce))
  colData(sce)[["total"]] <- runif(ncol(sce))

  # add some rowdata columns
  rowData(sce)[["gene_names"]] <- rownames(sce)
  rowData(sce)[["other_column"]] <- runif(nrow(sce))

  # add some metadata
  library_id <- glue::glue("library-{batch}")
  sample_id <- glue::glue("sample-{batch}")

  metadata(sce)$library_id <- library_id
  metadata(sce)$sample_id <- sample_id
  metadata(sce)$total_reads <- ceiling(runif(1, 100, 1000))
  metadata(sce)$sample_metadata <- data.frame(
    sample_id = sample_id,
    library_id = library_id
  )

  # Copy counts -> logcounts just to make sure the assay is retained
  logcounts(sce) <- counts(sce)
  return(sce)
}

# Generate some shared data for testing `prepare_sce_for_merge()` -----
set.seed(1665)
total_cells <- 24 # divisible by 3
total_genes <- 12 # number of months intentionally.
sce <- add_sce_data(
  sim_sce(n_cells = total_cells, n_genes = total_genes, n_empty = 0),
  batch = "1"
)
sce_name <- "sce_object"
batch_column <- "batch" # not the default
cell_id_column <- "cell_id"
shared_features <- rownames(sce)[1:10]
retain_coldata_cols <- c("sum", "detected")
preserve_rowdata_cols <- "gene_names"
expected_coldata_cols <- c("sum", "detected", batch_column, cell_id_column)

sce1 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce2 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce3 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)

# Generate some shared data for testing `merge_sce_list()` ------
sce_list <- list(
  "sce1" = sce1,
  "sce2" = sce2,
  "sce3" = sce3
) |>
  purrr::imap(
    add_sce_data
  )


# Tests without altexps ----------------------------------------------

test_that("`prepare_sce_for_merge` works as expected when all columns are present, no altexps", {
  result_sce <- prepare_sce_for_merge(
    sce,
    sce_name,
    batch_column,
    cell_id_column,
    shared_features,
    retain_coldata_cols,
    preserve_rowdata_cols
  )

  expect_equal(ncol(result_sce), total_cells) # cells
  expect_equal(nrow(result_sce), length(shared_features)) # genes

  # colData names and contents:
  expect_setequal(
    names(colData(result_sce)),
    expected_coldata_cols
  )
  expect_equal(unique(result_sce[[batch_column]]), sce_name)

  expect_equal(
    rownames(colData(result_sce)),
    c(glue::glue("{sce_name}-{rownames(colData(sce))}"))
  )

  expect_equal(
    colnames(sce),
    colData(result_sce)$cell_id
  )



  # rowData names and contents:
  expect_setequal(
    names(rowData(result_sce)),
    c("gene_names", paste(sce_name, "other_column", sep = "-"))
  )
})



test_that("`prepare_sce_for_merge` works as expected when an expected column is missing, no altexps", {
  # REMOVE "detected" column first
  # It should get re-added in as all NAs
  colData(sce)$detected <- NULL

  result_sce <- prepare_sce_for_merge(
    sce,
    sce_name,
    batch_column,
    cell_id_column,
    shared_features,
    retain_coldata_cols, # will add detected back in as NA
    preserve_rowdata_cols
  )

  # `detected` column should exist as all NA:
  expect_true(all(is.na(colData(result_sce)$detected)))
})


test_that("merging SCEs with matching genes works as expected, no altexps", {
  # First a warning for this early-return scenario
  expect_warning(merge_sce_list(sce_list = list("sce1" = sce1)))

  # Works as expected:
  merged_sce <- merge_sce_list(
    sce_list,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )


  # correct number of genes and cells:
  expect_equal(nrow(merged_sce), total_genes)
  expect_equal(ncol(merged_sce), total_cells)

  # colData names and contents:
  expect_setequal(
    names(colData(merged_sce)),
    expected_coldata_cols
  )
  expect_equal(
    merged_sce[[batch_column]],
    c(
      rep("sce1", total_cells / 3),
      rep("sce2", total_cells / 3),
      rep("sce3", total_cells / 3)
    )
  )
  expect_setequal(
    rownames(colData(merged_sce)),
    c(
      glue::glue("sce1-{rownames(colData(sce1))}"),
      glue::glue("sce2-{rownames(colData(sce2))}"),
      glue::glue("sce3-{rownames(colData(sce3))}")
    )
  )

  # rowData names and contents:
  expect_setequal(
    names(rowData(merged_sce)),
    c("gene_names", "sce1-other_column", "sce2-other_column", "sce3-other_column")
  )
  expect_equal(
    merged_sce[[batch_column]],
    c(
      rep("sce1", total_cells / 3),
      rep("sce2", total_cells / 3),
      rep("sce3", total_cells / 3)
    )
  )

  # assays
  expect_setequal(
    assayNames(merged_sce),
    c("counts", "logcounts")
  )

  # metadata names check
  expect_contains(
    names(metadata(merged_sce)),
    c("library_id", "sample_id", "library_metadata", "sample_metadata")
  )

  # check that sample metadata is a data frame
  expect_s3_class(metadata(merged_sce)$sample_metadata, "data.frame")

  # library metadata should contain a list of library metadata with all three libraries
  expect_length(
    metadata(merged_sce)$library_metadata,
    3
  )

  # check that contents of library id and sample id are correct
  expect_setequal(
    metadata(merged_sce)$library_id,
    c("library-sce1", "library-sce2", "library-sce3")
  )

  expect_setequal(
    metadata(merged_sce)$sample_id,
    c("sample-sce1", "sample-sce2", "sample-sce3")
  )

  # check contents of sample metadata
  expect_setequal(
    metadata(merged_sce)$library_id,
    metadata(merged_sce)$sample_metadata$library_id
  )

  expect_setequal(
    metadata(merged_sce)$sample_id,
    metadata(merged_sce)$sample_metadata$sample_id
  )
})

test_that("merging SCEs with multiple sample ids per library to mirror cellhash works as expected", {
  # add multiple sample ids per library into the metadata to mirror cellhash data
  sce_list <- sce_list |>
    purrr::map2(
      c("sample4", "sample5", "sample6"),
      \(sce, other_sample)      {
        metadata(sce)$sample_id <- c(metadata(sce)$sample_id, other_sample)
        metadata(sce)$sample_metadata <- data.frame(
          sample_id = paste0(metadata(sce)$sample_id, collapse = ","),
          library_id = metadata(sce)$library_id
        )
        return(sce)
      }
    )


  # merge
  merged_sce <- merge_sce_list(
    sce_list,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names"),
    # explicit false here -
    include_altexp = FALSE
  )

  # Check that all sample ids are present
  expect_setequal(
    metadata(merged_sce)$sample_id,
    sce_list |> purrr::map(\(x) metadata(x)$sample_id) |> purrr::reduce(c)
  )
})


test_that("merging SCEs with different genes among input SCEs works as expected, no altexps", {
  # rename sce2 and sce3 genes so that only 1-6 are overlapping
  # hence, we started with 12 genes.
  rownames(sce_list[[2]]) <- c(
    rownames(sce_list[[2]])[1:6],
    month.name[1:6]
  )
  rownames(sce_list[[3]]) <- c(
    rownames(sce_list[[3]])[1:6],
    month.name[7:12]
  )

  # Works as expected:
  merged_sce <- merge_sce_list(sce_list)


  # correct number of genes and cells:
  expect_equal(nrow(merged_sce), 6)
  expect_equal(ncol(merged_sce), total_cells)
})




test_that("merging SCEs with no matching genes fails as expected, no altexps, no altexps", {
  # ensure different gene names entirely
  rownames(sce_list[[1]]) <- rownames(sce_list[[2]])
  rownames(sce_list[[1]]) <- paste0(rownames(sce_list[[1]]), "-new")

  expect_error(
    merge_sce_list(
      sce_list = list(
        "sce1" = sce_list[[1]],
        "sce2" = sce_list[[2]]
      )
    )
  )
})



test_that("merging SCEs without names works as expected, no altexps", {
  # First make sure it generates a warning -
  expect_warning(
    merged_sce <- merge_sce_list(
      unname(sce_list),
      batch_column = batch_column
    )
  )

  # The resulting batch names should be 1,2,3
  expect_equal(
    colData(merged_sce)[, batch_column],
    c(
      rep("1", total_cells / 3),
      rep("2", total_cells / 3),
      rep("3", total_cells / 3)
    )
  )
})

test_that("merging SCEs with library metadata fails as expected, no altexps", {
  # add library metadata to one of the objects in the list
  metadata(sce_list$sce1)$library_metadata <- "library_metadata"

  expect_error(
    merge_sce_list(
      sce_list
    )
  )
})


# Tests with altexps -------------------------------------------------------

## helper function to add an altExp to a simulated SCE ----
add_sce_altexp <- function(
    sce,
    batch,
    altexp_name,
    num_altexp_features,
    n_cells) {
  sce_alt <- sim_sce(
    n_genes = num_altexp_features,
    n_cells = n_cells,
    n_empty = 0
  )

  # ensure matching barcodes
  colnames(sce_alt) <- colnames(sce)

  # change feature names
  rownames(sce_alt) <- stringr::str_replace(
    rownames(sce_alt),
    "^GENE",
    toupper(altexp_name)
  )

  # add some rowdata columns
  rowData(sce_alt)[["target_type"]] <- "target" # should be retained
  rowData(sce_alt)[["feature_column"]] <- rownames(sce_alt)
  rowData(sce_alt)[["other_column"]] <- runif(nrow(sce_alt))

  # add a coldata column
  colData(sce_alt)[["coldata_column"]] <- runif(ncol(sce_alt))

  # add logcounts
  logcounts(sce_alt) <- counts(sce_alt)

  # add metadata
  library_id <- glue::glue("library-{batch}")
  sample_id <- glue::glue("sample-{batch}")

  metadata(sce_alt)$library_id <- library_id
  metadata(sce_alt)$sample_id <- sample_id
  metadata(sce_alt)$mapped_reads <- 100
  metadata(sce_alt)$ambient_profile <- runif(num_altexp_features)

  # Add some more columns to retain to the SCE based on this altExp
  colData(sce)[[glue::glue("altexps_{altexp_name}_sum")]] <- runif(ncol(sce))
  colData(sce)[[glue::glue("altexps_{altexp_name}_percent")]] <- runif(ncol(sce))
  colData(sce)[[glue::glue("altexps_{altexp_name}_detected")]] <- runif(ncol(sce))

  # add sce_alt as sce's altExp
  altExp(sce, altexp_name) <- sce_alt

  return(sce)
}

num_altexp_features <- 5
altexp_name <- "adt"
sce_list_with_altexp <- sce_list |>
  purrr::imap(
    add_sce_altexp,
    altexp_name,
    num_altexp_features,
    total_cells / 3
  )
# vector of all expected names
full_altexp_features <- rownames(altExp(sce_list_with_altexp[[1]]))



test_that("prepare_sce_for_merge() works as expected with is_altexp=TRUE", {
  test_altexp <- altExp(sce_list_with_altexp[[1]])
  prepared_altexp <- prepare_sce_for_merge(
    test_altexp,
    "test",
    batch_column = "batch",
    cell_id_column = "cell",
    full_altexp_features,
    retain_coldata_cols = NULL,
    preserve_rowdata_cols = "target_type",
    is_altexp = TRUE
  )

  expect_equal(
    colnames(rowData(prepared_altexp)),
    c("target_type", "test-feature_column", "test-other_column")
  )

  # column names should be unchanged
  expect_equal(
    colnames(prepared_altexp), colnames(test_altexp)
  )

  expect_equal(
    colnames(colData(prepared_altexp)), c("batch", "cell")
  )
})


test_that("merging SCEs with altExps works as expected when include_altexps = FALSE", {
  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names"),
    include_altexp = FALSE
  )

  # there should not be any altExps
  expect_length(altExpNames(merged_sce), 0)
})


test_that("merging SCEs with altExps has correct altExp colData names when retaining altExps columns", {
  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names"),
    retain_altexp_coldata_cols = list("not_present" = c("coldata_column")),
    preserve_altexp_rowdata_cols = list("adt" = c("target_type")),
    include_altexp = TRUE
  )

  # test correct altExp colData names
  expected_cols <- c(batch_column, cell_id_column)
  observed_cols <- altExp(merged_sce) |>
    colData() |>
    names()
  expect_setequal(
    expected_cols,
    observed_cols
  )

  # test correct altExp rowData names
  expected_cols <- c(
    "target_type",
    "sce1-feature_column", "sce1-other_column",
    "sce2-feature_column", "sce2-other_column",
    "sce3-feature_column", "sce3-other_column"
  )
  observed_cols <- altExp(merged_sce) |>
    rowData() |>
    names()
  expect_setequal(
    expected_cols,
    observed_cols
  )
})


test_that("merging SCEs with altExps has correct altExp colData names when retaining _irrelevant_ altExps columns", {
  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names"),
    retain_altexp_coldata_cols = list("adt" = c("coldata_column")),
    include_altexp = TRUE
  )

  # test correct altExp colData names
  expected_cols <- c(batch_column, cell_id_column, "coldata_column")
  observed_cols <- altExp(merged_sce) |>
    colData() |>
    names()
  expect_setequal(
    expected_cols,
    observed_cols
  )
})


test_that("merging SCEs with 1 altexp and same features works as expected, with altexps", {
  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = c(
      retain_coldata_cols,
      "altexps_adt_sum",
      "altexps_adt_detected",
      "altexps_adt_percent"
    ),
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )

  merged_altexp <- altExp(merged_sce)


  expect_true(altExpNames(merged_sce) == altexp_name)
  expect_equal(dim(merged_altexp), c(num_altexp_features, total_cells))
  expect_equal(rownames(merged_altexp), full_altexp_features)

  expected_colnames <- sce_list_with_altexp |>
    purrr::imap(
      \(sce, sce_name) glue::glue("{sce_name}-{colnames(sce)}")
    ) |>
    unlist() |>
    unname()
  expect_equal(colnames(merged_altexp), expected_colnames)

  # Check colData columns
  expected_coldata <- c(
    "sum",
    "detected",
    "altexps_adt_sum",
    "altexps_adt_detected",
    "altexps_adt_percent",
    batch_column,
    "cell_id"
  )

  expect_setequal(
    colnames(colData(merged_sce)),
    expected_coldata
  )
})




test_that("merging SCEs with 1 altexp but different features fails as expected, with altexps", {
  # keep only the first 3 features from the first SCE
  altExp(sce_list_with_altexp[[1]]) <- altExp(sce_list_with_altexp[[1]])[1:3, ]


  expect_error(
    merge_sce_list(
      sce_list_with_altexp,
      batch_column = batch_column,
      # "total" should get removed
      retain_coldata_cols = retain_coldata_cols,
      # this row name should not be modified:
      preserve_rowdata_cols = c("gene_names")
    )
  )
})




test_that("merging SCEs where 1 altExp is missing works as expected, with altexps", {
  sce_list_with_altexp$sce4 <- sce_list[[1]]

  # update the metdata list with sce4 name
  metadata(sce_list_with_altexp$sce4) <- list(
    library_id = "library-sce4",
    sample_id = "sample-sce4",
    total_reads = 100,
    sample_metadata = data.frame(
      sample_id = "sample-sce4",
      library_id = "library-sce4"
    )
  )

  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )

  expect_equal(altExpNames(merged_sce), "adt")

  merged_altexp <- altExp(merged_sce, "adt")

  # check merged_altexp metadata
  altexp_metadata <- metadata(merged_altexp)
  expect_setequal(
    names(altexp_metadata),
    c("library_id", "sample_id", "library_metadata")
  )
  expect_setequal(
    altexp_metadata$library_id,
    glue::glue("library-{names(sce_list_with_altexp)}")
  )
  expect_setequal(
    altexp_metadata$sample_id,
    glue::glue("sample-{names(sce_list_with_altexp)}")
  )


  expect_setequal(
    # all but sce4 contain all metadata components
    altexp_metadata$library_metadata[-4] |>
      purrr::map(names) |>
      purrr::reduce(intersect),
    c("library_id", "sample_id", "mapped_reads", "ambient_profile")
  )

  expect_setequal(
    # sce4 has only library id and sample id, as it was missing the altExp
    names(altexp_metadata$library_metadata$sce4),
    c("library_id", "sample_id")
  )


  expect_true(
    is.null(altexp_metadata$library_metadata$sce4$ambient_profile) &
      is.null(altexp_metadata$library_metadata$sce4$mapped_reads)
  )
  expect_true(
    all(
      is.numeric(altexp_metadata$library_metadata$sce1$ambient_profile),
      is.numeric(altexp_metadata$library_metadata$sce2$ambient_profile),
      is.numeric(altexp_metadata$library_metadata$sce3$ambient_profile)
    )
  )
  expect_true(
    all(
      altexp_metadata$library_metadata$sce1$library_id == "library-sce1",
      altexp_metadata$library_metadata$sce2$library_id == "library-sce2",
      altexp_metadata$library_metadata$sce3$library_id == "library-sce3",
      altexp_metadata$library_metadata$sce4$library_id == "library-sce4"
    )
  )

  # check dimensions
  expect_equal(
    dim(merged_sce),
    c(total_genes, total_cells * 4 / 3)
  )
  expect_equal(
    dim(merged_altexp),
    c(num_altexp_features, total_cells * 4 / 3)
  )

  # check colData names are as expected
  expected_coldata <- c(
    "sum",
    "detected",
    batch_column,
    cell_id_column
  )
  expect_setequal(
    colnames(colData(merged_sce)),
    expected_coldata
  )

  # check that the NAs are as expected
  counts_mat <- counts(merged_altexp)
  sce4_counts <- counts_mat[, merged_sce[[batch_column]] == "sce4"]
  expect_true(
    all(is.na(sce4_counts))
  )
  numeric_counts <- counts_mat[, merged_sce[[batch_column]] != "sce4"]
  expect_true(
    all(is.finite(numeric_counts))
  )
})


test_that("merging SCEs with different altExps works as expected; each SCE has 1 altExp of a different name", {
  sce1 <- sce_list_with_altexp[[1]]
  sce2 <- sce_list_with_altexp[[2]]
  other_altexp_name <- "other"
  altExpNames(sce2) <- other_altexp_name
  rownames(altExp(sce2)) <- c("ADT-A", "ADT-B", "ADT-C", "ADT-D", "ADT-E")

  sce_list <- list(
    "sce1" = sce1,
    "sce2" = sce2
  )

  # Merge:
  merged_sce <- merge_sce_list(
    sce_list,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )

  # Correct names
  expect_setequal(
    altExpNames(merged_sce),
    c(altexp_name, other_altexp_name)
  )

  # Check the "adt" (`altexp_name`) altexp
  adt_merged <- altExp(merged_sce, altexp_name)
  expect_equal(
    dim(adt_merged),
    c(num_altexp_features, ncol(merged_sce))
  )
  expect_equal(
    rownames(adt_merged),
    rownames(altExp(sce1, altexp_name))
  )
  expect_equal(
    colnames(adt_merged),
    colnames(merged_sce)
  )
  expect_true(
    all(is.na(counts(adt_merged)[, merged_sce[[batch_column]] == "sce2"]))
  )

  # Check the "other"  altexp
  other_merged <- altExp(merged_sce, other_altexp_name)
  expect_equal(
    dim(other_merged),
    c(num_altexp_features, ncol(merged_sce))
  )
  expect_equal(
    rownames(other_merged),
    rownames(altExp(sce2))
  )
  expect_equal(
    colnames(adt_merged),
    colnames(merged_sce)
  )
  expect_true(
    all(is.na(counts(other_merged)[, merged_sce[[batch_column]] == "sce1"]))
  )
})



test_that("merging SCEs with different altExps works as expected; each SCE has 2 different altExps", {
  other_altexp_name <- "other"
  other_n_features <- 3
  sce1 <- sce_list_with_altexp[[1]]
  sce2 <- sce_list_with_altexp[[2]]

  sce_list <- list(
    "sce1" = sce_list_with_altexp[[1]],
    "sce2" = sce_list_with_altexp[[2]]
  ) |>
    # add "other" altExp to each sce
    purrr::imap(
      add_sce_altexp,
      other_altexp_name,
      other_n_features,
      ncol(sce1)
    )

  # Merge
  merged_sce <- merge_sce_list(
    sce_list,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )

  # Correct names
  expect_setequal(
    altExpNames(merged_sce),
    c(altexp_name, other_altexp_name)
  )

  # Check the "adt" (`altexp_name`) altexp
  adt_merged <- altExp(merged_sce, altexp_name)
  expect_equal(
    dim(adt_merged),
    c(num_altexp_features, ncol(merged_sce))
  )
  expect_equal(
    rownames(adt_merged),
    rownames(altExp(sce1, altexp_name))
  )
  expect_equal(
    colnames(adt_merged),
    colnames(merged_sce)
  )

  # next two tests check matrix values
  expect_equal(
    counts(adt_merged)[, merged_sce[[batch_column]] == "sce1"] |>
      as.numeric(),
    counts(altExp(sce_list[[1]], altexp_name)) |>
      as.numeric()
  )
  expect_equal(
    counts(adt_merged)[, merged_sce[[batch_column]] == "sce2"] |>
      as.numeric(),
    counts(altExp(sce_list[[2]], altexp_name)) |>
      as.numeric()
  )

  # Check the "other"  altexp
  other_merged <- altExp(merged_sce, other_altexp_name)
  expect_equal(
    dim(other_merged),
    c(other_n_features, ncol(merged_sce))
  )
  expect_equal(
    rownames(other_merged),
    rownames(altExp(sce_list[[1]], other_altexp_name))
  )
  expect_equal(
    colnames(other_merged),
    colnames(merged_sce)
  )
  # next two tests check matrix values
  expect_equal(
    counts(other_merged)[, merged_sce[[batch_column]] == "sce1"] |>
      as.numeric(),
    counts(altExp(sce_list[[1]], other_altexp_name)) |>
      as.numeric()
  )
  expect_equal(
    counts(other_merged)[, merged_sce[[batch_column]] == "sce2"] |>
      as.numeric(),
    counts(altExp(sce_list[[2]], other_altexp_name)) |>
      as.numeric()
  )
})




## Other tests ------------------

test_that("get_altexp_attributes passes when it should pass", {
  attribute_list <- get_altexp_attributes(sce_list_with_altexp)
  expect_equal(
    attribute_list[["adt"]][["assays"]], c("counts", "logcounts")
  )
  expect_equal(
    attribute_list[["adt"]][["features"]], full_altexp_features
  )
})


test_that("get_altexp_attributes throws an error as expected when features do not match", {
  altExp(sce_list_with_altexp[[1]]) <- altExp(sce_list_with_altexp[[1]])[1:3, ]
  expect_error(get_altexp_attributes(sce_list_with_altexp))
})


test_that("check_metadata throws an error when a field is missing", {
  metadata(sce_list[[1]])$library_id <- NULL
  expect_error(check_metadata(sce_list))
  expect_no_error(check_metadata(sce_list, expected_fields = "sample_id"))
})

test_that("get_altexp_metadata returns the correct values", {
  expected_list <- list(
    "library_id" = "library-sce1",
    "sample_id"  = "sample-sce1"
  )
  observed_list <- get_altexp_metadata(sce_list[[1]], "MISSING_ALTEXP_NAME")
  expect_equal(observed_list, expected_list)
})

test_that("prepare_merged_metadata works as expected, with sample_metadata present", {
  observed_metadata <- sce_list |>
    purrr::map(metadata) |>
    prepare_merged_metadata()

  expect_setequal(
    names(observed_metadata),
    c("library_id", "sample_id", "library_metadata", "sample_metadata")
  )

  expected_library_ids <- glue::glue("library-{names(sce_list)}")
  expected_sample_ids <- glue::glue("sample-{names(sce_list)}")
  expect_setequal(
    observed_metadata$library_id,
    expected_library_ids
  )
  expect_setequal(
    observed_metadata$sample_id,
    expected_sample_ids
  )
  expect_setequal(
    names(observed_metadata$library_metadata),
    names(sce_list)
  )
  expect_setequal(
    observed_metadata$library_metadata |>
      purrr::map(names) |>
      purrr::reduce(intersect),
    c("library_id", "sample_id", "total_reads")
  )
  expect_equal(
    data.frame(
      sample_id = expected_sample_ids,
      library_id = expected_library_ids
    ),
    observed_metadata$sample_metadata
  )
})


test_that("prepare_merged_metadata works as expected from altExps metadata (aka, sample_metadata NOT present)", {
  observed_metadata <- sce_list_with_altexp |>
    purrr::map(altExp) |>
    purrr::map(metadata) |>
    prepare_merged_metadata()

  expect_setequal(
    names(observed_metadata),
    c("library_id", "sample_id", "library_metadata")
  )

  expect_setequal(
    observed_metadata$library_id,
    glue::glue("library-{names(sce_list)}")
  )
  expect_setequal(
    observed_metadata$sample_id,
    glue::glue("sample-{names(sce_list)}")
  )
  expect_setequal(
    observed_metadata$library_metadata |>
      purrr::map(names) |>
      purrr::reduce(intersect),
    c("library_id", "sample_id", "mapped_reads", "ambient_profile")
  )
})
