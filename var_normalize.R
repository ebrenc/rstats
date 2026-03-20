
# var_transform_modular.R
# Reversible variable transformations for tidyverse workflows.
#
# Core logic:
#   fit -> transform -> inverse
#
# Main functions:
#   - normal:  var_normal_fit(), var_normal_transform(), var_normal_inverse()
#   - uniform: var_uniform_fit(), var_uniform_transform(), var_uniform_inverse()
#
# Convenience wrappers:
#   - var_normalize()
#   - var_uniformize()
#   - var_denormalize()
#   - var_deuniformize()
#
# Metadata:
#   Transformation metadata are stored in:
#     attr(df, "var_transform_meta")
#   and are used for inverse transformations and reporting.
#
# See minimal examples at the end of this file.

.ensure_bestnormalize <- function() {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    stop(
      "Package 'bestNormalize' is required. Install it with install.packages('bestNormalize').",
      call. = FALSE
    )
  }
}

.ensure_scales <- function() {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop(
      "Package 'scales' is required for ggplot2 scale transforms. Install it with install.packages('scales').",
      call. = FALSE
    )
  }
}

.safe_shapiro_p <- function(x, max_n = 5000, seed = 123) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3 || length(unique(x)) < 3) {
    return(NA_real_)
  }
  
  if (length(x) > max_n) {
    set.seed(seed)
    x <- sample(x, max_n)
  }
  
  suppressWarnings(stats::shapiro.test(x)$p.value)
}

.apply_pretransform <- function(x, method) {
  if (method == "raw") return(x)
  if (method == "log") {
    if (any(x <= 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(log(x))
  }
  if (method == "sqrt") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(sqrt(x))
  }
  stop("Unknown pretransform: ", method, call. = FALSE)
}

.invert_pretransform <- function(x, method) {
  if (method == "raw") return(x)
  if (method == "log") return(exp(x))
  if (method == "sqrt") return(x^2)
  stop("Unknown pretransform: ", method, call. = FALSE)
}

.format_poly_equation <- function(coef_vec, pretransform = "raw", x_label = "x") {
  input_name <- if (pretransform == "raw") x_label else paste0(pretransform, "(", x_label, ")")
  terms <- c("1", input_name)
  if (length(coef_vec) >= 3) {
    terms <- c(terms, paste0(input_name, "^", 2:(length(coef_vec) - 1)))
  }
  terms <- terms[seq_along(coef_vec)]
  paste0("y = ", paste0(signif(coef_vec, 6), " * ", terms, collapse = " + "))
}

.eval_poly <- function(z, coef_vec) {
  # Horner would be faster, but this is readable and enough for moderate degree.
  sum(coef_vec * z^(0:(length(coef_vec) - 1)))
}

.eval_poly_vec <- function(z, coef_vec) {
  vapply(z, .eval_poly, numeric(1), coef_vec = coef_vec)
}

.get_meta <- function(x, fit = NULL, expected_method = NULL) {
  if (is.null(fit)) {
    fit <- attr(x, "transform_meta")
  }
  if (is.null(fit) || !is.list(fit) || is.null(fit$method)) {
    stop("No valid transformation metadata found. Supply 'fit' explicitly or pass a vector with attribute 'transform_meta'.", call. = FALSE)
  }
  if (!is.null(expected_method) && !identical(fit$method, expected_method)) {
    stop("Transformation metadata has method '", fit$method, "', expected '", expected_method, "'.", call. = FALSE)
  }
  fit
}

.attach_meta_if_requested <- function(x, fit, attach_meta = TRUE) {
  if (isTRUE(attach_meta)) attr(x, "transform_meta") <- fit
  x
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------------
# Metadata store helpers
# --------------------------------

append_transform_meta <- function(meta_store = NULL, name, fit) {
  if (is.null(meta_store)) meta_store <- list()
  meta_store[[name]] <- fit
  meta_store
}

collect_transform_meta <- function(df) {
  metas <- lapply(df, attr, which = "transform_meta", exact = TRUE)
  metas[!vapply(metas, is.null, logical(1))]
}

# --------------------------------
# Normal transformation: fit / transform / inverse / ggplot scale
# --------------------------------

var_normal_fit <- function(x, poly_n = 3, pretransform = c("auto", "raw", "log", "sqrt")) {
  .ensure_bestnormalize()
  pretransform <- match.arg(pretransform)
  
  x <- as.numeric(x)
  ok <- is.finite(x)
  x_ok <- x[ok]
  
  if (!any(ok)) {
    stop("x has no finite values.", call. = FALSE)
  }
  
  if (length(unique(x_ok)) <= 1) {
    coef_vec <- c(`(Intercept)` = 0)
    fit <- list(
      method = "normal",
      preprocess = "raw",
      poly_n = 0L,
      coef = coef_vec,
      shapiro_p = NA_real_,
      equation = "y = 0 * 1",
      x_pre_range = range(x_ok, na.rm = TRUE),
      x_orig_range = range(x_ok, na.rm = TRUE)
    )
    class(fit) <- c("var_normal_fit", "var_transform_fit", "list")
    return(fit)
  }
  
  candidates <- if (pretransform == "auto") c("raw", "log", "sqrt") else pretransform
  results <- list()
  
  for (cand in candidates) {
    x_pre <- .apply_pretransform(x_ok, cand)
    if (!all(is.finite(x_pre))) next
    if (length(unique(x_pre)) <= 1) next
    
    ord_obj <- suppressWarnings(bestNormalize::orderNorm(x_pre, warn = FALSE))
    y_target <- ord_obj$x.t
    
    fit_lm <- stats::lm(y_target ~ stats::poly(x_pre, poly_n, raw = TRUE))
    pred <- as.numeric(stats::predict(fit_lm, newdata = data.frame(x_pre = x_pre)))
    p_val <- .safe_shapiro_p(pred)
    
    results[[length(results) + 1]] <- list(
      preprocess = cand,
      poly_n = poly_n,
      coef = stats::coef(fit_lm),
      shapiro_p = p_val,
      score = p_val,
      equation = .format_poly_equation(stats::coef(fit_lm), pretransform = cand),
      x_pre_range = range(x_pre, na.rm = TRUE),
      x_orig_range = range(x_ok, na.rm = TRUE)
    )
  }
  
  if (length(results) == 0) {
    stop("No valid normal transformation candidate could be fitted.", call. = FALSE)
  }
  
  scores <- purrr::map_dbl(results, "score")
  
  if (all(is.na(scores)) || all(is.infinite(scores))) {
    best <- results[[1]]
  } else {
    best <- results[[which.max(scores)]]
  }
  
  best$score <- NULL
  
  fit <- c(list(method = "normal"), best)
  class(fit) <- c("var_normal_fit", "var_transform_fit", "list")
  fit
}

var_normal_transform <- function(x, fit, attach_meta = TRUE) {
  fit <- .get_meta(x = NULL, fit = fit, expected_method = "normal")
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(.attach_meta_if_requested(out, fit, attach_meta))

  x_pre <- .apply_pretransform(x[ok], fit$preprocess)
  keep <- is.finite(x_pre)

  if (length(fit$coef) == 1L) {
    out[ok][keep] <- fit$coef[1]
  } else {
    out_ok <- rep(NA_real_, sum(ok))
    out_ok[keep] <- .eval_poly_vec(x_pre[keep], as.numeric(fit$coef))
    out[ok] <- out_ok
  }

  .attach_meta_if_requested(out, fit, attach_meta)
}

var_normal_inverse <- function(x, fit = NULL, grid_n = 4000, fallback_uniroot = TRUE,
                               lower = -1e6, upper = 1e6, attach_meta = FALSE) {
  fit <- .get_meta(x, fit, expected_method = "normal")
  x_in <- as.numeric(x)
  out <- rep(NA_real_, length(x_in))
  ok <- is.finite(x_in)
  if (!any(ok)) return(.attach_meta_if_requested(out, fit, attach_meta))

  coef_vec <- as.numeric(fit$coef)
  if (length(coef_vec) == 1L) {
    out[ok] <- .invert_pretransform(rep(0, sum(ok)), fit$preprocess)
    return(.attach_meta_if_requested(out, fit, attach_meta))
  }

  # Fast approximate inverse by monotonic interpolation over observed x_pre range.
  z_grid <- seq(fit$x_pre_range[1], fit$x_pre_range[2], length.out = grid_n)
  y_grid <- .eval_poly_vec(z_grid, coef_vec)
  ord <- order(y_grid)
  y_grid <- y_grid[ord]
  z_grid <- z_grid[ord]

  # Remove duplicated y values so approx() is well behaved.
  keep <- !duplicated(y_grid)
  y_grid <- y_grid[keep]
  z_grid <- z_grid[keep]

  y_ok <- x_in[ok]
  z_hat <- approx(x = y_grid, y = z_grid, xout = y_ok, rule = 1, ties = "ordered")$y

  # Optional fallback with uniroot for out-of-range values.
  need_root <- is.na(z_hat) & isTRUE(fallback_uniroot)
  if (any(need_root)) {
    root_fun <- function(target) {
      tryCatch(
        stats::uniroot(function(z) .eval_poly(z, coef_vec) - target, lower = lower, upper = upper)$root,
        error = function(e) NA_real_
      )
    }
    z_hat[need_root] <- vapply(y_ok[need_root], root_fun, numeric(1))
  }

  out_ok <- rep(NA_real_, length(y_ok))
  ok2 <- is.finite(z_hat)
  out_ok[ok2] <- .invert_pretransform(z_hat[ok2], fit$preprocess)
  out[ok] <- out_ok

  .attach_meta_if_requested(out, fit, attach_meta)
}

var_normal_trans <- function(fit, n = 6, digits = 2) {
  scales::trans_new(
    name = paste0("var_normal_", fit$preprocess),
    transform = function(x) var_normal_transform(x, fit),
    inverse = function(x) var_normal_inverse(x, fit),
    breaks = function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(numeric(0))
      br <- pretty(range(x, na.rm = TRUE), n = n)
      if (fit$preprocess %in% c("log", "sqrt")) br <- br[br > 0]
      br[is.finite(br)]
    },
    format = function(x) format(round(x, digits), trim = TRUE)
  )
}

# Convenience wrappers
var_normalize <- function(x, poly_n = 3, pretransform = c("auto", "raw", "log", "sqrt"),
                          return_fit = FALSE, attach_meta = TRUE) {
  fit <- var_normal_fit(x = x, poly_n = poly_n, pretransform = pretransform)
  x_t <- var_normal_transform(x = x, fit = fit, attach_meta = attach_meta)
  if (isTRUE(return_fit)) return(list(x = x_t, fit = fit))
  x_t
}

var_denormalize <- function(x, fit = NULL, ...) {
  var_normal_inverse(x = x, fit = fit, ...)
}

# --------------------------------
# Uniform transformation: fit / transform / inverse / ggplot scale
# --------------------------------

var_uniform_fit <- function(x, ties_method = "average") {
  x <- as.numeric(x)
  ok <- is.finite(x)
  x_ok <- x[ok]

  if (!any(ok)) {
    stop("x has no finite values.", call. = FALSE)
  }

  fit <- list(
    method = "uniform",
    ties_method = ties_method,
    reference = sort(x_ok),
    valid_n = length(x_ok),
    equation = paste0("u = (rank(x, ties.method = '", ties_method, "') - 0.5) / n")
  )
  class(fit) <- c("var_uniform_fit", "var_transform_fit", "list")
  fit
}

var_uniform_transform <- function(x, fit, attach_meta = TRUE) {
  fit <- .get_meta(x = NULL, fit = fit, expected_method = "uniform")
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(.attach_meta_if_requested(out, fit, attach_meta))

  ranks_ok <- rank(x[ok], ties.method = fit$ties_method, na.last = "keep")
  out[ok] <- (ranks_ok - 0.5) / sum(ok)

  .attach_meta_if_requested(out, fit, attach_meta)
}

var_uniform_inverse <- function(x, fit = NULL, quantile_type = 8, attach_meta = FALSE) {
  fit <- .get_meta(x, fit, expected_method = "uniform")
  u <- as.numeric(x)
  out <- rep(NA_real_, length(u))
  ok <- is.finite(u)
  if (!any(ok)) return(.attach_meta_if_requested(out, fit, attach_meta))

  probs <- pmin(pmax(u[ok], 0), 1)
  out[ok] <- as.numeric(stats::quantile(fit$reference, probs = probs, type = quantile_type, na.rm = TRUE))

  .attach_meta_if_requested(out, fit, attach_meta)
}

var_uniform_trans <- function(fit, quantile_type = 8) {
  .ensure_scales()
  fit <- .get_meta(x = NULL, fit = fit, expected_method = "uniform")
  scales::trans_new(
    name = paste0("var_uniform_", fit$ties_method),
    transform = function(x) var_uniform_transform(x, fit = fit, attach_meta = FALSE),
    inverse   = function(x) var_uniform_inverse(x, fit = fit, quantile_type = quantile_type, attach_meta = FALSE)
  )
}

# Convenience wrappers
var_uniformize <- function(x, ties_method = "average", return_fit = FALSE, attach_meta = TRUE) {
  fit <- var_uniform_fit(x = x, ties_method = ties_method)
  x_t <- var_uniform_transform(x = x, fit = fit, attach_meta = attach_meta)
  if (isTRUE(return_fit)) return(list(x = x_t, fit = fit))
  x_t
}

var_deuniformize <- function(x, fit = NULL, ...) {
  var_uniform_inverse(x = x, fit = fit, ...)
}

# --------------------------------
# Data-frame helper
# --------------------------------

var_transform_df <- function(df, vars, method = c("normal", "uniform", "both"), meta_store = TRUE) {
  
  method <- match.arg(method)
  vars_quo <- rlang::enquo(vars)
  vars_expr <- rlang::get_expr(vars_quo)
  
  if (is.character(vars_expr)) {
    vars <- intersect(vars_expr, names(df))
  } else {
    vars <- names(tidyselect::eval_select(vars_quo, df))
  }
  
  meta_list <- list()
  
  for (v in vars) {
    
    if (method %in% c("normal", "both")) {
      fit_n <- var_normal_fit(df[[v]])
      new_name <- paste0(v, "_normal")
      df[[new_name]] <- var_normal_transform(df[[v]], fit_n)
      meta_list[[new_name]] <- fit_n
    }
    
    if (method %in% c("uniform", "both")) {
      fit_u <- var_uniform_fit(df[[v]])
      new_name <- paste0(v, "_uniform")
      df[[new_name]] <- var_uniform_transform(df[[v]], fit_u)
      meta_list[[new_name]] <- fit_u
    }
  }
  
  if (length(meta_list) > 0) {
    attr(df, "var_transform_meta") <- meta_list
  }
  
  if (meta_store) {
    list(df = df, meta_store = meta_list)
  } else {
    df
  }
}

# --------------------------------
# Printing helpers
# --------------------------------

print.var_normal_fit <- function(x, ...) {
  cat("<var_normal_fit>\n")
  cat("  preprocess:", x$preprocess, "\n")
  cat("  poly_n:", x$poly_n, "\n")
  cat("  shapiro_p:", x$shapiro_p, "\n")
  cat("  equation:", x$equation, "\n")
  invisible(x)
}

print.var_uniform_fit <- function(x, ...) {
  cat("<var_uniform_fit>\n")
  cat("  ties_method:", x$ties_method, "\n")
  cat("  valid_n:", x$valid_n, "\n")
  cat("  equation:", x$equation, "\n")
  invisible(x)
}

# --------------------------------
# Minimal examples
# --------------------------------

# TRANSFORMAR

# df <- df %>% 
#   mutate(
#     across(c(euc, mah, pil, bha), var_normalize, .names = "{.col}_normal"),
#     across(c(euc, mah, pil, bha), var_uniformize, .names = "{.col}_uniform")
#   )

# DESTRANSFORMAR

# meta_store <- attr(df, "var_transform_meta")
# df <- df %>% 
#   mutate(
#     across(
#       ends_with("_normal"),
#       ~ var_normal_inverse(.x, meta_store[[cur_column()]]),
#       .names = "{.col}_back"
#     ),
#     across(
#       ends_with("_uniform"),
#       ~ var_uniform_inverse(.x, meta_store[[cur_column()]]),
#       .names = "{.col}_back"
#     )
#   )
