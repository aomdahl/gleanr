test_that("correlationAssessment gets averge frobnorm across matrices", {
  set.seed(22)
  m1 <- cbind(seq(0,10, length.out = 100),rnorm(100), seq(1000,0,length.out = 100))
  m2 <- cbind(seq(0,10, length.out = 100),rnorm(100), seq(1000,0,length.out = 100))

  mats <- list(m1,m2)
  expect_equal(correlationAssessment(mats), mean(c(norm(cor(m1), "F"), norm(cor(m2), "F"))))
})
