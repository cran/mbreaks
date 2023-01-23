context("name of variables specification")

test_that("Misspecify y variable names", {
 y = rnorm(20)
 data = data.frame(y)
 expect_error(dosequa('z',data=data),'No matching dependent variable y. Please try again')
})


test_that("Misspecify z variable names", {
  y = rnorm(20)
  data = data.frame(y)
  expect_error(dosequa('y','z',data=data),'No z regressors found. Please try again')
})

test_that("Misspecify x variable names", {
  y = rnorm(20)
  data = data.frame(y)
  expect_error(dosequa('y',x_name = 'x',data=data),'No x regressors found. Please try again')
})


