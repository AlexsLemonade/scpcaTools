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
                       rpois(400, 10)) # random background
  hash_reads[cbind(1:4, 1:100)] <- 1000 # one barcode higher for each column
  # add as altExp
  hash_sce <- SingleCellExperiment(list(counts = hash_reads))
  altExp(sce, "cellhash") <- hash_sce

  # add barcode table to sce
  sce_hashtable <- add_cellhash_ids(sce, hashsample_table, altexp_id = "cellhash")
  # pull out barcodes and sort to match original
  extracted_barcodes <- rowData(altExp(sce_hashtable, "cellhash"))|>
    as.data.frame() |>
    tibble::remove_rownames() |>
    dplyr::arrange(barcode_id)
  expect_equal(hashsample_table, extracted_barcodes)

  # test when not all barcodes are present
  expect_warning(add_cellhash_ids(sce, hashsample_table[1:2,], remove_unlabeled = FALSE))
  reduced_sce <- add_cellhash_ids(sce, hashsample_table[1:2,], remove_unlabeled = TRUE)
  expect_equal(sort(rownames(altExp(reduced_sce))), sort(hashsample_table$barcode_id[1:2]))
  # test rowname replacements
  reduced_sce <- add_cellhash_ids(sce, hashsample_table[1:2,], remove_unlabeled = TRUE, replace_rownames = TRUE)
  expect_equal(sort(rownames(altExp(reduced_sce))), sort(hashsample_table$sample_id[1:2]))
  # test with the wrong altexp_id
  expect_error(add_cellhash_ids(sce, hashsample_table, altexp_id = "foo"))
  # test with bad barcode tables
  expect_error(add_cellhash_ids(sce, data.frame(a = 1:4, b = 1:4)))
  expect_error(add_cellhash_ids(sce, rbind(hashsample_table, c("tag1", "sample5"))))
  expect_error(add_cellhash_ids(sce, hashsample_table[,1]))


  hash_sce <- add_demux_hashedDrops(sce_hashtable)
  hash_cols <- c("hashedDrops_sampleid",
                 "hashedDrops_bestsample",
                 "hashedDrops_LogFC",
                 "hashedDrops_Confident") # a subset of the expected columns for the altExp rowData
  expect_true(all(hash_cols %in% colnames(colData(altExp(hash_sce)))))
  expect_false(is.null(hash_sce$hashedDrops_sampleid))
  # check the results are by sample_id
  expect_true(all(hash_sce$hashedDrops_sampleid[!is.na(hash_sce$hashedDrops_sampleid)] %in% hashsample_table$sample_id))
  # results with no sample table present
  expect_warning(hash_sce_nosample <- add_demux_hashedDrops(sce))
  expect_true(all(hash_sce_nosample$hashedDrops_sampleid[!is.na(hash_sce_nosample$hashedDrops_sampleid)] %in% rownames(altExp(hash_sce_nosample))))

  # test seurat functions
  hto_sce <- add_demux_seurat(sce_hashtable)
  hto_cols <- c("HTODemux_maxID",
                "HTODemux_secondID",
                "HTODemux_margin",
                "HTODemux_classification") # a subset of the expected columns for the altExp rowData
  expect_true(all(hto_cols %in% colnames(colData(altExp(hto_sce)))))
  expect_false(is.null(hto_sce$HTODemux_sampleid))
  # check the results are by sample_id
  expect_true(all(hto_sce$HTODemux_sampleid[!is.na(hto_sce$HTODemux_sampleid)] %in% hashsample_table$sample_id))
  # results with no sample table present
  expect_warning(hto_sce_nosample <- add_demux_seurat(sce))
  expect_true(all(hto_sce_nosample$hashedDrops_sampleid[!is.na(hto_sce_nosample$hashedDrops_sampleid)] %in% rownames(altExp(hto_sce_nosample))))

})


