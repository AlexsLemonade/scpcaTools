# Establish testing data
set.seed(27)
sce <- sim_sce(n_genes = 20, n_cells = 100, n_empty = 0)
reducedDim(sce, "PCA") <- matrix(runif(20 * 100, min = 0, max = 100), nrow = 100)
cluster_colname <- "clustered_results"
blus_param <- bluster::KmeansParam(centers = 10)


test_that("cluster_sce function works when correctly specified", {

  sce_clustered <- cluster_sce(sce,
                               BLUSPARAM = blus_param,
                               cluster_column_name = cluster_colname)

  # Does the output look as expected?
  expect_true(cluster_colname %in% names(colData(sce_clustered)))
  expect_equal(levels(sce_clustered[[cluster_colname]]),
               as.character(1:10))
})

test_that("cluster_sce function should fail when inputs incorrectly specified", {

  expect_error(
    cluster_sce("not an sce object",
                BLUSPARAM = blus_param,
                cluster_column_name = cluster_colname)
  )

  expect_error(
    cluster_sce(sce,
                pc_name = "definitely not the PCA name",
                BLUSPARAM = blus_param,
                cluster_column_name = cluster_colname)
  )

  expect_error(
    cluster_sce(sce,
                BLUSPARAM = "not a bluster param",
                cluster_column_name = cluster_colname)
  )

  expect_error(
    cluster_sce(sce,
                BLUSPARAM = blus_param)
  )

})
