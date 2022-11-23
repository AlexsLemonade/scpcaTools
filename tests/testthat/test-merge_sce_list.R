# helper function for tests:
# Add more info to the SCEs
add_sce_data <- function(sce) {
  colData(sce) <- colData(sce) |>
    as.data.frame() |>
    dplyr::mutate(
      sum = runif(nrow(colData(sce))),
      detected = runif(nrow(colData(sce))),
      total = runif(nrow(colData(sce)))
    ) |>
    S4Vectors::DataFrame()

  rowData(sce) <- rowData(sce) |>
    as.data.frame() |>
    dplyr::mutate(
      gene_names = rownames(sce),
      other_column = runif(nrow(rowData(sce)))
    ) |>
    S4Vectors::DataFrame()

  # Copy counts -> logcounts just to make sure the assay is retained
  logcounts(sce) <- counts(sce)
  return(sce)
}


test_that("merging SCEs with matching genes works as expected", {

  # create three SCEs with the same number of genes
  set.seed(1665)
  total_cells <- 75
  total_genes <- 50
  sce1 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
  sce2 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
  sce3 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)

  # Set up list for input to function
  sce_list <- purrr::map(
    list(
      "sce1" =  sce1,
      "sce2" = sce2,
      "sce3" = sce3
    ),
    add_sce_data
  )

  # First expect some errors or warnings:
  expect_error(merge_sce_list(sce_list = unname(sce_list))) # List must be named
  expect_warning(merge_sce_list(sce_list = list("sce1" = sce1))) # List must have >=2 SCEs for no warning
  expect_error(merge_sce_list(sce_list = sce_list)) # Default columns aren't in these sces


  # Works as expected:
  batch_column <- "nondefault_batch"
  merged_sce <- merge_sce_list(sce_list,
                               # use a non-default batch column
                               batch_column = batch_column,
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
    c(batch_column, "sum", "total")
  )
  expect_equal(merged_sce[[batch_column]],
               c(rep("sce1", total_cells/3),
                 rep("sce2", total_cells/3),
                 rep("sce3", total_cells/3))
  )
  expect_equal(
    sort(rownames(colData(merged_sce))),
    sort(
      c(glue::glue("{rownames(colData(sce1))}-sce1"),
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

  # create three SCEs with the same number of genes
  set.seed(1665)
  total_cells <- 75
  total_genes <- 12 # number of months, wait for it!
  sce1 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
  sce2 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)
  sce3 <- sim_sce(n_cells = total_cells/3, n_genes = total_genes, n_empty = 0)

  # rename sce2 and sce3 genes so that only 1-6 are overlapping
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
