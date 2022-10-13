intron_mat <- Matrix::Matrix(
  data = c(1.0,1,  3,100,10,
           2,1,100,  0,10
           ),
  ncol = 2,
  sparse = TRUE
)
rownames(intron_mat) <- c("a", "a-I", "b", "b-I", "c")

usa_mat <- Matrix::Matrix(
  data = c(0,1,1,  2,100, 1,10,
           1,1,1, 80,  0,20,10
  ),
  ncol = 2,
  sparse = TRUE
)
rownames(usa_mat) <- c("a", "a-U", "a-A", "b", "b-U","b-A", "c")

spliced_mat <- Matrix::Matrix(
  data = c(1,  3,10,
           2,100,10
  ),
  ncol = 2,
  sparse = TRUE
)
rownames(spliced_mat) <- c("a", "b", "c")

unspliced_mat <- Matrix::Matrix(
  data = c(2,103,10,
           3,100,10
  ),
  ncol = 2,
  sparse = TRUE
)
rownames(unspliced_mat) <- c("a", "b", "c")

# We are checking only contents & names, because
# aggregate.Matrix used internally may add additional attributes

test_that("Check spliced collapse", {
  collapsed_mat <- collapse_intron_counts(intron_mat)
  collapsed_usa <- collapse_intron_counts(usa_mat)
  # check the contents
  expect_equal(collapsed_mat@x,
               spliced_mat@x)
  # check the dimnames
  expect_equal(collapsed_mat@Dimnames,
               spliced_mat@Dimnames)
  # same checks for usa mode
  expect_equal(collapsed_usa@x,
               spliced_mat@x)
  expect_equal(collapsed_usa@Dimnames,
               spliced_mat@Dimnames)
})

test_that("Check unspliced collapse", {
  collapsed_mat <- collapse_intron_counts(intron_mat, "unspliced")
  collapsed_usa <- collapse_intron_counts(usa_mat, "unspliced")
  # check the contents
  expect_equal(collapsed_mat@x,
               unspliced_mat@x)
  # check the dimnames
  expect_equal(collapsed_mat@Dimnames,
               unspliced_mat@Dimnames)
  # same checks for usa mode
  expect_equal(collapsed_usa@x,
               unspliced_mat@x)
  expect_equal(collapsed_usa@Dimnames,
               unspliced_mat@Dimnames)

})
