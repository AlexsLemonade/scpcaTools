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
expected_coldata_cols <- sort(c("sum", "detected", batch_column, cell_id_column))

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

test_that("`update_sce_metadata()` returns the expected list", {

  metadata_list <- metadata(sce_list[[1]])

  new_metadata <- update_sce_metadata(metadata_list)

  expect_equal(
    names(new_metadata),
    c("library_id", "sample_id", "library_metadata", "sample_metadata")
  )

  expect_equal(
    names(new_metadata$library_metadata),
    c("library_id",  "sample_id", "total_reads")
  )
})


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
  expect_equal(
    sort(names(colData(result_sce))),
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

  # metadata names check
  expect_contains(
    names(metadata(result_sce)),
    c("library_id", "sample_id", "library_metadata", "sample_metadata")
  )

  # check that sample metadata is a data frame
  expect_s3_class(metadata(result_sce)$sample_metadata, "data.frame")

  # check that contents of library id and sample id are correct
  expect_equal(
    metadata(result_sce)$library_id,
    "library-1"
  )

  expect_equal(
    metadata(result_sce)$sample_id,
    "sample-1"
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
    n_empty = 0)

  # ensure matching barcodes
  colnames(sce_alt) <- colnames(sce)

  # add some rowdata columns
  rowData(sce_alt)[["feature_column"]] <- rownames(sce_alt)
  rowData(sce_alt)[["other_column"]] <- runif(nrow(sce_alt))

  # add a coldata columns
  colData(sce_alt)[["coldata_column"]] <- runif(ncol(sce_alt))

  # add logcounts
  logcounts(sce_alt) <- counts(sce_alt)

  # add metadata
  library_id <- glue::glue("library-{batch}")
  sample_id <- glue::glue("sample-{batch}")

  metadata(sce_alt)$library_id <- library_id
  metadata(sce_alt)$sample_id <- sample_id
  metadata(sce_alt)$mapped_reads <- 100

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

test_that("merging SCEs with 1 altexp and same features works as expected, with altexps", {

  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
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

})




test_that("merging SCEs with 1 altexp but different features fails as expected, with altexps", {

  # second list where 1 is missing all
  altExp(sce_list_with_altexp[[1]]) <- altExp(sce_list_with_altexp[[1]])[1:3,]

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

  sce_list_with_altexp[["sce4"]] <- removeAltExps(sce_list_with_altexp[[1]])

  merged_sce <- merge_sce_list(
    sce_list_with_altexp,
    batch_column = batch_column,
    # "total" should get removed
    retain_coldata_cols = retain_coldata_cols,
    # this row name should not be modified:
    preserve_rowdata_cols = c("gene_names")
  )

  expect_true(altExpNames(merged_sce) == "adt")

})
