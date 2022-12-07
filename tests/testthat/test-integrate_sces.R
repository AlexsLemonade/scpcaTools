# generate testing data
set.seed(1665)
merged_sce <- sim_sce(n_cells = 300, n_genes = 100, n_empty = 0)
batches <- c(rep("a", 100),
             rep("b", 100),
             rep("c", 100))
barcodes <- rownames(colData(merged_sce))

# set up colData rownames as batch-barcode
new_rownames <- tibble::tibble(
  batches = batches,
  barcodes = barcodes,
  ids = glue::glue("{batches}-{barcodes}")) %>%
  dplyr::pull(ids)
rownames(colData(merged_sce)) <- new_rownames

# Add some info to colData:
# Add "sample" to colData
colData(merged_sce)$sample <- batches
# Add in a "covariate" column for testing
colData(merged_sce)$covariate <- sample(c("patient-a", "patient-b", "patient-c", "patient-d"),
                                          size = 300,
                                          replace = TRUE)

# add a logcounts assay for testing (numbers don't matter)
logcounts(merged_sce) <- counts(merged_sce)

# add PCs for testing (again numbers don't matter)
# make a 300x100 matrix
reducedDim(merged_sce, "PCA") <- matrix(runif(300*100,min=0,max=100),nrow = 300)

# Other variables that will be used in tests:
batch_column <- "sample"

################################################################################



test_that("`integrate_fastmnn` works as expected", {
  suppressWarnings(
    # warnings are supressed here b/c simulated data plays poorly enough with
    # algorithms to trigger warnings like:
    ### Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
    ### You're computing too large a percentage of total singular values, use a standard svd instead.
    integrated_sce <- integrate_fastmnn(merged_sce,
                                        batch_column)
  )
  # Result should have a `reconstructed` assay
  expect_true("reconstructed" %in% assayNames(integrated_sce))

  # Result should have a `corrected` reducedDim
  expect_true("corrected" %in% reducedDimNames(integrated_sce))

})


test_that("`integrate_fastmnn` fails as expected when logcounts is missing", {
  logcounts(merged_sce) <- NULL
  expect_error(integrate_fastmnn(merged_sce,
                                 batch_column))
})


test_that("`integrate_harmony` works as expected", {
  harmony_result <- integrate_harmony(merged_sce,
                                      batch_column)
  # check type
  expect_true(all(class(harmony_result) == c("matrix", "array")))

  # check dimensions
  expect_true(nrow(harmony_result) == 300)

  # rownames should equal the cell_id values
  expect_true(all(rownames(harmony_result) == new_rownames))
})



test_that("`integrate_harmony` fails when PCs are missing", {
  reducedDim(merged_sce, "PCA") <- NULL
  expect_error(
    integrate_harmony(merged_sce, batch_column)
  )
})

test_that("`integrate_harmony` fails when covariate columns are missing", {
  expect_error(
    integrate_harmony(merged_sce,
                      batch_column,
                      harmony_covariate_cols = "not_a_column"
    )
  )
})


################################################################################
################################################################################


test_that("`integrate_sces` fail as expected", {


  # bad sce
  expect_error(
    integrate_sces(batch_column, batch_column)
  )

  # bad integration method
  expect_error(
    integrate_sces(merged_sce, "not_a_real_method")
  )

  # missing batch column
  expect_error(
    integrate_sces(merged_sce,
                   "fastMNN",
                   batch_column = "not_a_column")
  )

  # insufficient batches
  merged_sce$sample <- "a"
  expect_error(
    integrate_sces(merged_sce,
                   "fastMNN")
  )

})


test_that("`integrate_sces` works as expected for fastmnn defaults", {

  suppressWarnings({
    # simulated-data related numerical warnings
    integrated_sce <- integrate_sces(merged_sce,
                                     "fastMNN")
  })

  expect_equal(assayNames(integrated_sce),
               c("counts", "logcounts", "fastMNN_corrected")
  )

  expect_true(
    "fastMNN_PCA" %in% reducedDimNames(integrated_sce)
  )

})


test_that("`integrate_sces` works as expected for return_corrected_expression=FALSE", {

  suppressWarnings({
    # simulated-data related numerical warnings
    integrated_sce <- integrate_sces(merged_sce,
                                     "fastMNN",
                                     return_corrected_expression = FALSE)
  })

  expect_equal(assayNames(integrated_sce),
               c("counts", "logcounts")
  )

})





test_that("`integrate_sces` works as expected with fastmnn extra arguments", {

  expect_no_error(
    suppressWarnings({
      # simulated-data related numerical warnings
      integrated_sce <- integrate_sces(merged_sce,
                                       "fastMNN",
                                       cos.norm = FALSE)
    })
  )

})



test_that("`integrate_sces` works as expected for harmony defaults", {

  integrated_sce <- integrate_sces(merged_sce,
                                   "harmony")
  expect_true(
    "harmony_PCA" %in% reducedDimNames(integrated_sce)
  )
})


test_that("`integrate_sces` works as expected with harmony extra arguments", {

  expect_no_error(
    integrated_sce <- integrate_sces(merged_sce,
                                     "harmony",
                                     batch_column,
                                     lambda = 2)
  )

})



test_that("`integrate_sces` works as expected with a harmony covariate", {
  expect_no_error(
    integrate_sces(merged_sce,
                   "harmony",
                   batch_column,
                   covariate_cols = "covariate")
  )
})

