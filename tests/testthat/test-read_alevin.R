dir <- system.file("extdata", package="scpcaData")
alevin_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevin_txome")
usa_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_usa_filtered")
intron_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_intron_filtered")

# expected alevin matrix size
sce_alevin_size <- c(60275, 10378)
# expected alevin-fry matrix size
sce_af_size <- c(60321, 614)
sample_ids = c("A1", "A2")

test_that("reading salmon alevin data works", {
  sce <- read_alevin(alevin_dir,
                     sample_id = sample_ids,
                     include_unspliced = FALSE)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_alevin_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  expect_equal(sce@metadata$mapping_tool, "alevin")
  expect_equal(sce@metadata$transcript_type, c("spliced"))
  expect_false(sce@metadata$include_unspliced)
  expect_false(is.null(sce@metadata$salmon_version))
  expect_false(is.null(sce@metadata$reference_index))
  expect_null(sce@metadata$alevinfry_version)
  expect_equal(sce@metadata$sample_id, sample_ids)
})

test_that("reading alevin-fry USA mode works", {
  sce <- read_alevin(usa_dir,
                     usa_mode = TRUE,
                     include_unspliced = TRUE,
                     sample_id = sample_ids)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_af_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  # check metadata
  expect_equal(sce@metadata$mapping_tool, "alevin-fry")
  expect_equal(sce@metadata$transcript_type, c("total", "spliced"))
  expect_true(sce@metadata$include_unspliced)
  expect_false(is.null(sce@metadata$salmon_version))
  expect_false(is.null(sce@metadata$reference_index))
  expect_false(is.null(sce@metadata$alevinfry_version))
  expect_equal(sce@metadata$sample_id, sample_ids)

  # no remaining unspliced
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)
})


test_that("reading alevin-fry intron mode works", {
  sce <- read_alevin(intron_dir,
                     include_unspliced = TRUE,
                     usa_mode = FALSE)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_af_size)
  # check that column names are barcodes
  col_barcode <- stringr::str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  # check metadata
  expect_equal(sce@metadata$mapping_tool, "alevin-fry")
  expect_equal(sce@metadata$transcript_type, c("total", "spliced"))
  expect_true(sce@metadata$include_unspliced)
  expect_false(is.null(sce@metadata$salmon_version))
  expect_false(is.null(sce@metadata$reference_index))
  expect_false(is.null(sce@metadata$alevinfry_version))

  # no remaining unspliced
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)
})
