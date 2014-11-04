context("PvarBreakTest")

test_that("Basic", {
  
  N = 1000
  x = rnorm(N) + rep(0:1, each=N/2)
  Test = PvarBreakTest(x)
  
  expect_that(Test <- PvarBreakTest(x), not(throws_error()))
  expect_that(summary(Test), not(throws_error()))
  expect_that(plot(Test), not(throws_error()))

})
  
  
  
  
  
