#' Joint Test via Contrast Statement
#'
#' @param model Model object that contrast statement applies to
#' @param coef_df_input Coefficients for contrast
#' @param cov_mtx_input Covariance matrix for contrast
#' @param L Vector or matrix corresponding to contrast coefficients
#' @param C Vector of contrasts
#'
#' @returns Tibble of contrast statement results containing the point estimate,
#'   test statistic, and p-value
#' @export
#'
#' @import
#' stats
#' dplyr
#' broom
#' tibble
joint_test <- function(model = NULL, coef_df_input = NULL,
                       cov_mtx_input = NULL,
                       L, C) {
  # if a model is supplied, pull out estimates and cov mtx
  if (!is.null(model)) {
    coef_df <- matrix(broom::tidy(model)$estimate)
    cov_mtx <- unname(stats::vcov(model))
    # otherwise, use supplied ests and cov mtx
  } else {
    coef_df <- coef_df_input
    cov_mtx <- cov_mtx_input
  }

  if (is.null(model) & is.null(coef_df_input) & is.null(cov_mtx_input)){
    stop("Either `model` or `coef_df_input` and `cov_mtx_input` are required.")
  }

  # check input - if model and coefficients are not supplied then skip
  if (is.null(model) & is.null(coef_df_input)) {
    tibble::tribble(
      ~estimate, ~std_error, ~lower,
      ~upper, ~test_stat, ~pval, ~pval_formatted,
      NA, NA, NA,
      NA, NA, NA, NA
    )
  } else {
    # browser()

    # compute test statistic and p-value for Wald test ~ chisq with df = nrow(c)
    # the [,1] just cleans up the variable name from test_stat[,1] to test_stat
    tibble::tibble(
      estimate = sum((L %*% coef_df)),
      # std error and conf int are not correct for joint test (estimate and test stat + pval are correct)
      # std_error = sqrt(L %*% cov_mtx %*% t(L))[[1]],
      # lower = estimate - qnorm(0.975) * std_error,
      # upper = estimate + qnorm(0.975) * std_error,
      test_stat = (t((L %*% coef_df) - C) %*% solve(L %*% cov_mtx %*% t(L)) %*% ((L %*% coef_df) - C))[, 1],
      pval = 1 - stats::pchisq(q = .data$test_stat, df = nrow(C)),
      pval_formatted = format.pval(.data$pval, digits = 4)
    )
  }
}
