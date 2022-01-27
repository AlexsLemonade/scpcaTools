test_that("cellhash functions work", {
  # create an sce with cellhash data
  set.seed(1665)
  sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
  # create barcode and sample ids
  hashsample_table <- data.frame(barcode_id = paste0("tag", 1:4),
                              sample_id = paste0("sample", 1:4))
  # create hash reads matrix (random order for barcodes)
  hash_reads <- matrix(ncol = 100, nrow = 4,
                       dimnames = list(sample(hashsample_table$barcode_id), colnames(sce)),
                       c(rep(c(100,0,0,0), 25),
                         rep(c(0,100,0,0), 25),
                         rep(c(0,0,100,0), 25),
                         rep(c(0,0,0,100), 25)))
  # add as altExp
  hash_sce <- SingleCellExperiment(list(counts = hash_reads))
  altExp(sce, "cellhash") <- hash_sce

  # add barcode table to sce
  sce_barcode <- add_hashsample_table(sce, hashsample_table, altexp_id = "cellhash")
  # pull out barcodes and sort to match original
  extracted_barcodes <- rowData(altExp(sce_barcode, "cellhash"))|>
    as.data.frame() |>
    tibble::remove_rownames() |>
    dplyr::arrange(barcode_id)
  expect_equal(hashsample_table, extracted_barcodes)

  # test when not all barcodes are present
  expect_warning(add_hashsample_table(sce, hashsample_table[1:2,]))
  # test with the wrong altexp_id
  expect_error(add_hashsample_table(sce, hashsample_table, altexp_id = "foo"))
  # test with bad barcode tables
  expect_error(add_hashsample_table(sce, data.frame(a = 1:4, b = 1:4)))
  expect_error(add_hashsample_table(sce, hashsample_table[,1]))
})


