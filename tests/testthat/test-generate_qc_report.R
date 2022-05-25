set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 10)
filt_sce <- filter_counts(sce)

test_that("generating a qc report works", {
  qc_file <- generate_qc_report(library_id = "TEST",
                                unfiltered_sce = sce,
                                filtered_sce = filt_sce, )
  expect_true(file.exists(qc_file))
  # clean up
  file.remove(qc_file)
})
