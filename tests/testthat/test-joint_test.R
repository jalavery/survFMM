# Unit tests for the joint_test function
test_that("joint_test function works correctly with a model", {
  # Create a sample model
  mod1 <- lm(mpg ~ wt + hp, data = mtcars)

  L <- matrix(c(0, 1, -1), nrow = 1)
  C <- matrix(0, nrow = 1)

  result <- joint_test(model = mod1, L = L, C = C)

  expect_equal(result$estimate, unname(mod1$coefficients[2])-unname(mod1$coefficients[3]))
  expect_type(result, "list")
  expect_named(result, c("estimate", "test_stat", "pval", "pval_formatted"))
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$test_stat))
  expect_true(!is.na(result$pval))
})

test_that("joint_test function works correctly with direct inputs", {
  coef_df_input <- matrix(c(0.5165, 1.2489), ncol = 1)
  cov_mtx_input <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

  L <- matrix(c(1, 0), nrow = 1)
  C <- matrix(0, nrow = 1)

  result <- joint_test(coef_df_input = coef_df_input,
                       cov_mtx_input = cov_mtx_input, L = L, C = C)

  expect_equal(result$estimate, coef_df_input[1])
  expect_type(result, "list")
  expect_named(result, c("estimate", "test_stat", "pval", "pval_formatted"))
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$test_stat))
  expect_true(!is.na(result$pval))
})

test_that("joint_test function handles NULL inputs correctly", {
  expect_error(joint_test(L = matrix(1), C = matrix(0)),
               "Either `model` or `coef_df_input` and `cov_mtx_input` are required.")
})

test_that("check against glht", {
  # check joint test script against glht
  m1 <- glm(response ~ age + marker + stage,
            family = binomial,
            data = gtsummary::trial)

  # a single test
  l <- matrix(
    c(0, 1, 0, 0, 0, 0),
    nrow = 1,
    byrow = TRUE
  )
  rownames(l) <- "Test age term"

  glht_results <- multcomp::glht(m1,
                                 linfct = l)

  s1 <- summary(glht_results, test = multcomp::adjusted("none"))

  # now try my function: match
  joint_test_single_est <- joint_test(coef_df_input = coef(m1),
             cov_mtx_input = unname(vcov(m1)),
             L = l,
             C = matrix(0))

  # test that estimates and p-value are equal
  expect_equal(unname(s1$test$coefficients), joint_test_single_est$estimate)
  expect_equal(unname(s1$test$pvalues)[1], joint_test_single_est$pval)

  # now test the difference in two covariates
  # a single test of equivalence of covariates: i match
  l2 <- matrix(
    c(0, 1, -1, 0, 0, 0),
    nrow = 1,
    byrow = TRUE
  )
  rownames(l2) <- "Test age vs marker term"

  glht_results2 <- multcomp::glht(m1, linfct = l2)

  s2 <- summary(glht_results2, test = multcomp::adjusted("none"))

  # now try my function: match
  joint_test_two_covars <- joint_test(coef_df_input = coef(m1),
             cov_mtx_input = unname(vcov(m1)),
             L = l2,
             C = matrix(0))

  # test that estimates and p-value are equal
  expect_equal(unname(s2$test$coefficients), joint_test_two_covars$estimate)
  expect_equal(unname(s2$test$pvalues)[1], joint_test_two_covars$pval)

  # joint test of equivalence of covariates: i match
  l3 <- matrix(
    c(0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0),
    nrow = 2,
    byrow = TRUE
  )
  rownames(l3) <- c("Test age term",
                    "Test marker term")

  # now try my function: matches SAS
  joint_test_joint <- joint_test(coef_df_input = coef(m1),
             cov_mtx_input = unname(vcov(m1)),
             L = l3,
             C = matrix(c(0, 0), nrow = 2))

  # test that p-value is equal to SAS output
  expect_equal(0.0709, round(joint_test_joint$pval, 4))

})
