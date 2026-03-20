# detcov()
# Determinant of the covariance matrix (robust or classical)
#
# Computes the determinant of the covariance matrix of a numeric matrix or data frame.
# Optionally uses a robust estimator (MASS::cov.trob) when available.
#
# Behaviour:
#   - removes rows with missing values
#   - returns NA if:
#       * not enough observations (n < p + 1)
#       * any variable has zero or non-finite variance
#       * covariance matrix cannot be computed or is not positive definite
#   - for univariate input, returns the variance
#
# Arguments:
#   x       numeric matrix or data frame
#   robust  logical; if TRUE (default), attempts robust covariance estimation
#
# Returns:
#   A numeric scalar (determinant of covariance matrix) or NA_real_
#
# Examples:
#
#   # Basic usage
#   detcov(matrix(rnorm(100), ncol = 2))
#
#   # With a data frame
#   df %>% summarise(d = detcov(pick(x, y)))
#
#   # Grouped computation
#   df %>%
#     group_by(group, condition) %>%
#     summarise(d = detcov(pick(x, y)), .groups = "drop")
#
#   # Compare raw / log / normalized distributions
#   df_det <- df %>%
#     summarise(d = detcov(pick(x, y)))
#
#   ggplot(df_det, aes(d)) + geom_histogram()
#   ggplot(df_det, aes(log(d))) + geom_histogram()
#
#   # Example with normalization
#   df_det %>%
#     mutate(d_norm = var_normalize(log(d))) %>%
#     ggplot(aes(d_norm)) + geom_histogram()
#
detcov <- function(x, robust = TRUE) {
  x <- as.matrix(x)
  x <- x[stats::complete.cases(x), , drop = FALSE]
  
  p <- ncol(x)
  if (nrow(x) < p + 1) return(NA_real_)
  
  sds <- apply(x, 2, stats::sd, na.rm = TRUE)
  if (any(!is.finite(sds)) || any(sds == 0)) return(NA_real_)
  
  S <- NULL
  
  if (isTRUE(robust) && requireNamespace("MASS", quietly = TRUE)) {
    fit <- tryCatch(
      suppressWarnings(MASS::cov.trob(x)),
      error = function(e) NULL
    )
    if (!is.null(fit) && !is.null(fit$cov)) {
      S <- fit$cov
    }
  }
  
  if (is.null(S)) {
    S <- tryCatch(stats::cov(x), error = function(e) NULL)
  }
  
  if (is.null(S) || anyNA(S) || any(!is.finite(S))) return(NA_real_)
  
  if (ncol(S) == 1) return(as.numeric(S[1, 1]))
  
  d <- tryCatch(as.numeric(det(S)), error = function(e) NA_real_)
  if (!is.finite(d) || d <= 0) return(NA_real_)
  
  d
}
