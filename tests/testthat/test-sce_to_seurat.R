set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
colData(sce) <- data.frame("test_column" = sample(0:10, 100, rep= TRUE),
                          "barcodes" = colnames(sce)) %>%
  tibble::column_to_rownames("barcodes") %>%
  DataFrame()
rowData(sce) <- data.frame("test_row" = sample(0:10, 200, rep = TRUE),
                           "gene_id" = rownames(sce)) %>%
  tibble::column_to_rownames("gene_id") %>%
  DataFrame()
metadata(sce)$test <- "test"
alt_sce <- sim_sce(n_cells = 100, n_genes = 10, n_empty = 0)
# rename the alt_exp so it has the same column names as the base sce
colnames(alt_sce) <- colnames(sce)
# merge alt into main
alt_name = "alt"
sce <- merge_altexp(sce, alt_sce, alt_name)
rowData(altExp(sce, alt_name)) <- data.frame("alt_test_row" = sample(0:5, 10, rep = TRUE),
                                             "alt_id" = rownames(altExp(sce, alt_name))) %>%
  tibble::column_to_rownames("alt_id") %>%
  DataFrame()

test_that("Converting SCE to Seurat objects works as expected", {

  seurat_object <- sce_to_seurat(sce)

  # check that column names of Seurat object are derived from SCE object
  # they won't necessarily be equal if some cells with 0 counts were removed
  expect_true(colnames(seurat_object) %in% colnames(sce))

  # check that attached metadata/coldata/rowdata in SCE are present in seurat object
  coldata_sce <- as.data.frame(colData(sce))
  rowdata_sce <- as.data.frame(rowData(sce))
  expect_equal(coldata_sce, seurat_object@meta.data)
  expect_equal(rowdata_sce, seurat_object@var.features)
  expect_equal(metadata(sce), seurat_object@misc)

  # check that altExp data was converted and rowData
  expect_false(is.null(seurat_object[[alt_name]]))
  alt_rowdata_sce <- as.data.frame(rowData(altExp(sce, alt_name)))
  expect_equal(alt_rowdata_sce, seurat_object[[alt_name]]@var.features)

})
