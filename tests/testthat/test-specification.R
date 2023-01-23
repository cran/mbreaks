context("trimming options for models and tests")


#### try/catch errors about m (maximum breaks), eps1 (trimming level),
# and h (minimum segment length) for both model estimation and testing ####
## Model estimation
test_that("Maximum breaks cannot be larger than sample size allowed", {
  T = 50
  y = rnorm(T)
  y[c(25:50)] = y[c(25:50)] + 1
  data = data.frame(y)
  m = 10
  upd_h = floor(0.15*T)

  expect_warning(dosequa('y',data=data,m = m),
                 paste('Not enough observations for ',m+1,' segments with minimum length per segment =',
                       upd_h,'.The total required observations for such',m,
                       'breaks would be ',(m+1)*upd_h,'>T=',T,'\n'))
})

test_that("Maximum breaks cannot be negative", {
  T = 50
  y = rnorm(T)
  y[c(25:50)] = y[c(25:50)] + 1
  data = data.frame(y)
  m = -1
  expect_error(dosequa('y',data=data,m = m),'Maximum number of breaks cannot be less than 1')
})

test_that("Sample size too small", {
  T = 5
  y = rnorm(T)
  data = data.frame(y)
  m = 2
  max_eps1 = 0.25
  min_h = 5
  expect_error(dofix('y',data=data,fixn = m),
               paste('Not enough observations for any available trimming level'))
})

##Estimation without trimming value
test_that("0 trimming level with no specified h", {
  T = 50
  y = rnorm(T)
  data = data.frame(y)
  m = 2
  max_eps1 = 0.25
  min_h = 5
  expect_error(dofix('y',data=data,fixn = m,eps1 = 0),
               paste('Need to specify h when setting this option to estimation only'))
})

test_that("h is insufficient for number of regressors", {
  T = 50
  y = rnorm(T)
  data = data.frame(y)
  m = 2
  max_eps1 = 0.25
  min_h = 5
  expect_message(dofix('y',data=data,fixn = m,eps1 = 0, h=0),
                 paste('Not enough observations in minimum segment for number of regressors. Set to 15% sample size'))
})



## Testing
test_that("SupF test maximum breaks cannot be 0", {
  T = 10
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  m = 0
  expect_error(dotest('y',data=data,m = m),'Maximum number of breaks cannot be less than 1')
})



#### try/catch errors about eps1 (trimming level)
test_that("Invalid trimming level eps1", {
  T = 50
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  eps1 = 1 #must be either: 0, 0.05, 0.10, 0.15, 0.20 or 0.25
  expect_warning(doseqtests('y',data=data,eps1=eps1),'Invalid trimming level, set trimming level to 15%')
})


###Cannot specify 0 trimming level in testing
test_that("Test must have trimming level > 0", {
  T = 100
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  m = 5
  eps1 = 0
  expect_error(dotest('y',data=data,m = m,eps1 = 0),'The tests are undefined for 0% trimming level')
})


