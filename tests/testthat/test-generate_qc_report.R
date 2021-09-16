set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
filt_sce <- sce[,1:50]

test_that("generating a qc report works", {
  qc_file <- generate_qc_report(sample_name = "TEST",
                                unfiltered_sce = sce,
                                filtered_sce = filt_sce, )
  expect_true(file.exists(qc_file))
})
