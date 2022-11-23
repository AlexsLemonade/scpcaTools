# helper function to add more info to a simulated SCE ----
add_sce_data <- function(sce) {

  # add some coldata columns
  colData(sce)[["sum"]] <- runif(nrow(colData(sce)))
  colData(sce)[["detected"]] <- runif(nrow(colData(sce)))
  colData(sce)[["total"]] <- runif(nrow(colData(sce)))

  # add some rowdata columns
  rowData(sce)[["gene_names"]] <- rownames(rowData(sce))
  rowData(sce)[["other_column"]] <- runif(nrow(rowData(sce)))

  # Copy counts -> logcounts just to make sure the assay is retained
  logcounts(sce) <- counts(sce)
  return(sce)
}


# Generate some shared data for testing `prepare_sce_for_merge()` ---
set.seed(1665)
total_cells <- 24 # divisible by 3
total_genes <- 12 # number of months intentionally.
sce <- add_sce_data(
  sim_sce(n_cells = total_cells, n_genes = total_genes, n_empty = 0)
)
sce_name <- "sce_object"
batch_column <- "batch" # not the default
barcode_column <- "barcode_col" # not the default
shared_features <- rownames(sce)[1:10]
retain_coldata_cols <- c("sum", "detected")
preserve_rowdata_cols <- "gene_names"
expected_coldata_names <- c("sum", "detected", "total")

# Generate some shared data for testing `merge_sce_list()` ------
sce1 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
sce2 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
sce3 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)








test_that("`prepare_sce_for_merge` works as expected when all columns are present", {

  # Barcode is NOT explicitly included in `retain_coldata_cols` or `expected_coldata_names`
  result_sce <- prepare_sce_for_merge(sce,
                                      sce_name,
                                      batch_column,
                                      barcode_column,
                                      shared_features,
                                      retain_coldata_cols,
                                      preserve_rowdata_cols,
                                      expected_coldata_names)

  expect_equal(ncol(result_sce), total_cells) # cells
  expect_equal(nrow(result_sce), length(shared_features)) # genes

  # colData names and contents:
  expect_equal(
    sort(names(colData(result_sce))),
    sort(c(batch_column, barcode_column, retain_coldata_cols))
  )
  expect_equal(unique(result_sce[[batch_column]]), sce_name)
  expect_equal(
    sort(rownames(colData(result_sce))),
    sort(c(glue::glue("{rownames(colData(sce))}-{sce_name}")))
  )

  # rowData names and contents:
  expect_equal(
    sort(names(rowData(result_sce))),
    c("gene_names", paste("other_column", sce_name, sep = "-"))
  )

})



test_that("`prepare_sce_for_merge` works as expected when all an expected column is missing", {

  # REMOVE "detected" column first -
  colData(sce)$detected <- NULL

  result_sce <- prepare_sce_for_merge(sce,
                                      sce_name,
                                      batch_column,
                                      barcode_column,
                                      shared_features,
                                      retain_coldata_cols,
                                      preserve_rowdata_cols,
                                      expected_coldata_names)

  # `detected` column should exist as all NA:
  expect_true(all(is.na(colData(result_sce)$detected)))

})



test_that("`prepare_sce_for_merge` works as expected when barcode column already exists", {

  # create barcode
  barcodes <- rownames(colData(sce))
  colData(sce)$barcode_col <- barcodes

  result_sce <- prepare_sce_for_merge(sce,
                                      sce_name,
                                      batch_column,
                                      barcode_column,
                                      shared_features,
                                      retain_coldata_cols,
                                      preserve_rowdata_cols,
                                      expected_coldata_names)

  # `barcode_col` column should have been added
  expect_equal(colData(result_sce)$barcode_col,
               barcodes)

})


####################################################################
####################################################################


test_that("merging SCEs with matching genes works as expected", {

  # Set up list for input to function
  sce_list <- purrr::map(
    list(
      "sce1" = sce1,
      "sce2" = sce2,
      "sce3" = sce3
    ),
    add_sce_data
  )

  # First expect some errors or warnings:
  expect_error(merge_sce_list(sce_list = unname(sce_list))) # List must be named
  expect_warning(merge_sce_list(sce_list = list("sce1" = sce1))) # this is an early return situation- List must have >=2 SCEs for no warning

  # Works as expected:
  merged_sce <- merge_sce_list(sce_list,
                               # use a non-default batch and barcode column names
                               batch_column = batch_column,
                               barcode_column = barcode_column,
                               # "detected" should get removed:
                               retain_coldata_cols = c("sum", "total"),
                               # this row name should not be modified:
                               preserve_rowdata_cols = c("gene_names"))


  # correct number of genes and cells:
  expect_equal(nrow(merged_sce), total_genes)
  expect_equal(ncol(merged_sce), total_cells)

  # colData names and contents:
  expect_equal(
    sort(names(colData(merged_sce))),
    sort(c(batch_column, barcode_column, "sum", "total"))
  )
  expect_equal(merged_sce[[batch_column]],
               c(rep("sce1", total_cells/3),
                 rep("sce2", total_cells/3),
                 rep("sce3", total_cells/3))
  )
  expect_equal(
    sort(rownames(colData(merged_sce))),
    sort(
      c(
      glue::glue("{rownames(colData(sce1))}-sce1"),
      glue::glue("{rownames(colData(sce2))}-sce2"),
      glue::glue("{rownames(colData(sce3))}-sce3")
    ))
)

  # rowData names and contents:
  expect_equal(
    sort(names(rowData(merged_sce))),
    c("gene_names", "other_column-sce1", "other_column-sce2", "other_column-sce3")
  )
  expect_equal(merged_sce[[batch_column]],
               c(rep("sce1", total_cells/3),
                 rep("sce2", total_cells/3),
                 rep("sce3", total_cells/3))
  )

  # assays
  expect_equal(
    sort(assayNames(merged_sce)),
    c("counts", "logcounts")
  )

})



test_that("merging SCEs with different genes among input SCEs works as expected", {

  # rename sce2 and sce3 genes so that only 1-6 are overlapping
  # hence, we started with 12 genes.
  rownames(sce2) <- c(rownames(sce2)[1:6],
                      month.name[1:6])
  rownames(sce3) <- c(rownames(sce3)[1:6],
                      month.name[7:12])

  # Set up list for input to function
  sce_list <- purrr::map(
    list(
      "sce1" =  sce1,
      "sce2" = sce2,
      "sce3" = sce3
    ),
    add_sce_data
  )

  # Works as expected:
  merged_sce <- merge_sce_list(sce_list,
                               # "detected" should get removed:
                               retain_coldata_cols = c("sum", "total"),
                               # this row name should not be modified:
                               preserve_rowdata_cols = c("gene_names"))


  # correct number of genes and cells:
  expect_equal(nrow(merged_sce), 6)
  expect_equal(ncol(merged_sce), total_cells)
})




test_that("merging SCEs with no matching genes fails as expected", {

  rownames(sce1) <- month.name # ensure different gene names entirely

  # Set up list for input to function
  sce_list <- purrr::map(
    list(
      "sce1" = sce1,
      "sce2" = sce2),
    add_sce_data
  )

  expect_error(merge_sce_list(sce_list = sce_list,
                              # use this arg to prevent other errors
                              retain_coldata_cols = "sum"))
})
