library(maed)

test_that("FUNCTION add", {
  expect_equal(add(1, 2), 3)
})

test_that("FUNCTION hello", {
  expect_equal(hello(), "Hello, world!")
})

x <- runif(10)
meanC(x)
