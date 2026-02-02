#' Joint Test via Contrast Statement
#'
#' @param model
#' @param coef_df_input
#' @param cov_mtx_input
#' @param L
#' @param C
#'
#' @returns
#' @export
#'
#' @examples
joint_test <- function(model = NULL, coef_df_input = NULL,
                       cov_mtx_input = NULL,
                       L, C) {
  # if a model is supplied, pull out estimates and cov mtx
  if (!is.null(model)) {
    coef_df <- matrix(broom::tidy(model)$estimate)
    cov_mtx <- unname(vcov(model))
    # otherwise, use supplied ests and cov mtx
  } else {
    coef_df <- coef_df_input
    cov_mtx <- cov_mtx_input
  }

  # check input - if NULL then skip
  if (is.null(coef_df_input)) {
    tribble(
      ~estimate, ~std_error, ~lower,
      ~upper, ~test_stat, ~pval, ~pval_formatted,
      NA, NA, NA,
      NA, NA, NA, NA
    )
  } else {
    # browser()

    # compute test statistic and p-value for Wald test ~ chisq with df = nrow(c)
    # the [,1] just cleans up the variable name from test_stat[,1] to test_stat
    tibble(
      estimate = sum((L %*% coef_df)),
      # std error and conf int are not correct for joint test (estimate and test stat + pval are correct)
      # std_error = sqrt(L %*% cov_mtx %*% t(L))[[1]],
      # lower = estimate - qnorm(0.975) * std_error,
      # upper = estimate + qnorm(0.975) * std_error,
      test_stat = (t((L %*% coef_df) - C) %*% solve(L %*% cov_mtx %*% t(L)) %*% ((L %*% coef_df) - C))[, 1],
      pval = 1 - pchisq(q = test_stat, df = nrow(C)),
      pval_formatted = format.pval(pval, digits = 4)
    )
  }
}
