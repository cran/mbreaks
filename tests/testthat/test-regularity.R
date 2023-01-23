context("regularity of matrix of regressors")

test_that("Duplicate regressors in z regressors", {
 y = rnorm(20)
 data = data.frame(y)
 data$c = rep(1,20)
 expect_error(dosequa('y','c',data=data),paste('There are duplicate regressors in both x and z regressors.'))
})


test_that("Same regressors in z and x regressors", {
  y = rnorm(20)
  data = data.frame(y)
  data$c = rep(1,20)
  expect_error(dosequa('y',x_name='c',data=data),paste('There are duplicate regressors in both x and z regressors.'))
})


test_that("Misspecified coefficient estimates of full sample regressors", {
  y = rnorm(100)
  data = data.frame(y)
  data$x = rnorm(100)
  fixb=1
  betaini=c(2,1)
  expect_error(dofix('y',x_name='x',fixb=1,betaini=c(2,1),data=data),paste('The initial estimates for full sample regressors are not in matrix form. Please specify beta0 in a p by 1 matrix'))
})
