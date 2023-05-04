test_that("vireo cell assignment works as expected", {
  # create an sce with cellhash data
  set.seed(1665)
  sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
  metadata(sce)$sample_id <- paste0("sample", 1:4)
  # construct cell assignments for half of cells, order arbitrary
  vireo_cells <- sample(colnames(sce), size = 50)
  donor_ids <- sample(c(metadata(sce)$sample_id, "doublet", "unassigned"),
    size = length(vireo_cells), replace = TRUE
  )

  vireo_df <- data.frame(cell = vireo_cells, donor_id = donor_ids)

  vireo_sce <- add_demux_vireo(sce, vireo_df)

  # dimensions should not change
  expect_equal(dim(vireo_sce), dim(sce))

  # check that columns were added to colData
  vireo_cols <- c(
    "vireo_sampleid",
    "vireo_donor_id"
  )
  expect_true(all(vireo_cols %in% colnames(colData(vireo_sce))))

  # test with different sample_ids
  sce2 <- sce
  metadata(sce2)$sample_id <- paste0("other", 1:4)
  expect_warning(add_demux_vireo(sce2, vireo_df))
})
