# helper function to add more info to a simulated SCE ----
add_sce_data <- function(sce) {
  # add some coldata columns
  colData(sce)[["sum"]] <- runif(ncol(sce))
  colData(sce)[["detected"]] <- runif(ncol(sce))
  colData(sce)[["total"]] <- runif(ncol(sce))

  # add some rowdata columns
  rowData(sce)[["gene_names"]] <- rownames(sce)
  rowData(sce)[["other_column"]] <- runif(nrow(sce))

  # add some metadata
  metadata(sce)$library_id <- "library1"
  metadata(sce)$sample_id <- "sample1"
  metadata(sce)$sample_metadata <- data.frame(
    diagnosis = "diagnosis",
    ontology_id = "ontology"
  )

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
cell_id_column <- "cell_id"
shared_features <- rownames(sce)[1:10]
retain_coldata_cols <- c("sum", "detected")
preserve_rowdata_cols <- "gene_names"
expected_coldata_cols <- sort(c("sum", "detected", batch_column, cell_id_column))

# Generate some shared data for testing `merge_sce_list()` ------
sce1 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce2 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce3 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce_list <- purrr::map(
  list(
    "sce1" = sce1,
    "sce2" = sce2,
    "sce3" = sce3
  ),
  add_sce_data
)



test_that("`prepare_sce_for_merge` works as expected when all columns are present", {
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
  expect_equal(
    sort(names(rowData(result_sce))),
    c("gene_names", paste(sce_name, "other_column", sep = "-"))
  )

  # metadata names check
  expect_true(
    all(c("library_id", "sample_id", "library_metadata", "sample_metadata") %in% names(metadata(result_sce)))
  )

  # check that sample metadata is a data frame
  expect_true(
    is.data.frame(metadata(result_sce)$sample_metadata)
  )

})



test_that("`prepare_sce_for_merge` works as expected when an expected column is missing", {
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


####################################################################
####################################################################


test_that("merging SCEs with matching genes works as expected", {
  # First a warning for this early-return scenario
  expect_warning(merge_sce_list(sce_list = list("sce1" = sce1)))

  # Works as expected:
  merged_sce <- merge_sce_list(sce_list,
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
  expect_equal(
    sort(names(colData(merged_sce))),
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
  expect_equal(
    sort(rownames(colData(merged_sce))),
    sort(
      c(
        glue::glue("sce1-{rownames(colData(sce1))}"),
        glue::glue("sce2-{rownames(colData(sce2))}"),
        glue::glue("sce3-{rownames(colData(sce3))}")
      )
    )
  )

  # rowData names and contents:
  expect_equal(
    sort(names(rowData(merged_sce))),
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
  expect_equal(
    sort(assayNames(merged_sce)),
    c("counts", "logcounts")
  )
})



test_that("merging SCEs with different genes among input SCEs works as expected", {
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




test_that("merging SCEs with no matching genes fails as expected", {
  # ensure different gene names entirely
  rownames(sce_list[[1]]) <- rownames(sce_list[[2]])
  rownames(sce_list[[1]]) <- paste0(rownames(sce_list[[1]]), "-new")

  expect_error(merge_sce_list(sce_list = list(
    "sce1" = sce_list[[1]],
    "sce2" = sce_list[[2]]
  )))
})



test_that("merging SCEs without names works as expected", {
  # First make sure it generates a warning -
  expect_warning(
    merged_sce <- merge_sce_list(unname(sce_list),
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
