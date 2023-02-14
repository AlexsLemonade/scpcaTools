set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 10)
filt_sce <- filter_counts(sce)

test_that("generating a qc report works", {
  qc_file <- generate_qc_report(library_id = "TEST",
                                unfiltered_sce = sce,
                                filtered_sce = filt_sce)
  expect_true(file.exists(qc_file))
  # clean up
  file.remove(qc_file)
})

test_that("generating a qc report fails with a missing template", {
  expect_error(generate_qc_report(library_id = "TEST",
                                unfiltered_sce = sce,
                                filtered_sce = filt_sce,
                                rmd_file = "missing.rmd"))
})


test_that("generating a qc report works with a warning", {
  expect_warning(generate_qc_report(library_id = "TEST",
                                    unfiltered_sce = sce,
                                    filtered_sce = filt_sce,
                                    extra_params = list(param = "test_param")))
})
