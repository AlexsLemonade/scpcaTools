intron_mat <- Matrix::Matrix(
  data = c(0,1,2,100,10,
           1,1,100,0,10
           ),
  ncol = 2
)
rownames(intron_mat) <- c("a", "a-I", "b", "b-I", "c")

spliced_mat <- Matrix::Matrix(
  data = c(0,2,10,
           1,100,10
  ),
  ncol = 2
)
rownames(spliced_mat) <- c("a", "b", "c")

unspliced_mat <- Matrix::Matrix(
  data = c(1,102,10,
           2,100,10
  ),
  ncol = 2
)
rownames(unspliced_mat) <- c("a", "b", "c")

# We are checking only contents & names, because
# aggregate.Matrix used internally may add additional attributes

test_that("Check spliced collapse", {
  collapsed_mat <- collapse_intron_counts(intron_mat)
  # check the contents
  expect_equal(collapsed_mat@x,
               spliced_mat@x)
  # check the dimnames
  expect_equal(collapsed_mat@Dimnames,
               spliced_mat@Dimnames)
})

test_that("Check unspliced collapse", {
  collapsed_mat <- collapse_intron_counts(intron_mat, "unspliced")
  # check the contents
  expect_equal(collapsed_mat@x,
               unspliced_mat@x)
  # check the dimnames
  expect_equal(collapsed_mat@Dimnames,
               spliced_mat@Dimnames)
})
