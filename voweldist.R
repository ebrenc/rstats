# voweldist(): compute within- and between-group vowel distance metrics
#
# Computes vowel-distance and category-separation measures for a binary
# contrast (`contrast_var`) using one or more acoustic dimensions
# (`dependent_vars`).
#
# The function can return:
# - token-level distances:
#     * Euclidean
#     * Mahalanobis
# - category-level separation measures:
#     * Pillai trace
#     * Bhattacharyya distance
#
# Stratification:
# `condition_vars` define the initial analysis strata. When multiple
# variables are supplied, they are internally combined into `.condition_id`.
# If `condition_vars = NULL`, a single global stratum is used.
#
# Speaker handling:
# `speaker_var` is optional.
# - If supplied, the most specific clouds may initially be speaker-specific.
# - If `speaker_var = NULL`, computations are performed on pooled data.
#
# Within-group mode (`inter_group = FALSE`):
# Metrics are computed within each group/stratum.
# Depending on the metric, the result may be:
# - token-level (e.g. Euclidean, Mahalanobis), or
# - stratum-level and then joined back to all tokens in that stratum
#   (e.g. Pillai, Bhattacharyya).
#
# Between-group mode (`inter_group = TRUE`):
# Metrics are computed between each focal group and `reference_group`,
# within the requested stratification, and then joined back to token-level
# rows.
#
# Cloud resolution / fallback:
# If the requested stratification is not feasible, the function
# automatically falls back to simpler clouds by progressively removing
# `speaker_var` and then `condition_vars` from right to left until a usable
# cloud is found, or until global pooling is reached.
#
# The effective cloud actually used is reported via messages and stored in:
# - `cloud_wit` for within-group computations
# - `cloud_bet` for between-group computations
#
# Metric-specific notes:
# - `contrast_var` must have exactly two levels.
# - Euclidean uses the centroid of the resolved reference cloud.
#   In between-group mode, Euclidean may still be computed when the
#   reference cloud contains only one complete observation, because a
#   centroid can still be defined.
# - Mahalanobis requires a usable covariance matrix in the resolved
#   reference cloud. If covariance cannot be estimated or inverted,
#   the result is `NA`.
# - Bhattacharyya and Pillai require sufficiently populated clouds for both
#   categories; otherwise fallback or `NA` may occur.
#
# Output:
# Depending on the requested metrics and mode, the function may return:
# - `dist_wit_euc`, `dist_wit_mah`, `dist_wit_pil`, `dist_wit_bha`
# - `dist_bet_euc`, `dist_bet_mah`, `dist_bet_pil`, `dist_bet_bha`
# - `cloud_wit` and/or `cloud_bet`
#
# Diagnostics:
# If `diagnostics = TRUE`, additional metric-specific cloud and failure
# columns may be returned. If `diagnostics = FALSE`, these auxiliary
# diagnostic columns are removed from the final output.
#
# Example 1: within-group distances with stratification
# df %>%
#   voweldist(
#     dependent_vars = c(B1B0, B2B1),
#     contrast_var = Vowel,
#     condition_vars = c(Time, WordAge),
#     speaker_var = Speaker,
#     compute_euclidean = TRUE,
#     compute_mahalanobis = TRUE
#   )
#
# Example 2: between-group comparison with reference group
# df %>%
#   voweldist(
#     dependent_vars = c(B1B0, B2B1),
#     contrast_var = Vowel,
#     condition_vars = c(Time, WordAge),
#     speaker_var = Speaker,
#     group_var = Condition,
#     reference_group = "Native",
#     inter_group = TRUE
#   )
#
# Example 3: global pooled analysis
# df %>%
#   voweldist(
#     dependent_vars = c(B1B0, B2B1),
#     contrast_var = Vowel,
#     condition_vars = NULL,
#     speaker_var = NULL
#   )

voweldist = function(
    x,                                  # data frame of application
    dependent_vars,                     # e.g. c(B1B0, B2B1)
    inter_group = FALSE,                # e.g. FALSE (within ES), TRUE (ES vs NS)
    contrast_var,                       # e.g. Vowel
    condition_vars,                     # e.g. TestingTime OR c(Task, TestingTime)
    speaker_var,                        # e.g. Participant
    group_var = NULL,                   # e.g. NativeParticipant (or Group)
    reference_group = NULL,             # e.g. "NS"
    focal_group = NULL,                 # e.g. "ES" (or "Control")
    within_ref = "self_speaker",        # "self_speaker", "inter_speaker_mean", "population_pooled"
    compute_euclidean = FALSE,
    compute_mahalanobis = TRUE,
    compute_pillai = FALSE,
    compute_bhatt = TRUE,
    diagnostics = FALSE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(tibble)
    library(rlang)
    library(tidyselect)
  })
  
  xdata = x %>% tibble::as_tibble()
  
  measure_nm  = names(tidyselect::eval_select(enquo(dependent_vars), xdata))
  n_measures  = length(measure_nm)
  min_n_cloud = n_measures + 1
  
  if (is.null(compute_bhatt))       compute_bhatt       <- !inter_group
  if (is.null(compute_mahalanobis)) compute_mahalanobis <-  inter_group
  
  options(scipen = 999)
  
  xdata = xdata %>% mutate(.token_rowid = as.character(dplyr::row_number()))
  
  # ---- capture columns ----
  contrast_q  = enquo(contrast_var)
  speaker_q   = enquo(speaker_var)
  group_q     = enquo(group_var)
  
  contrast_nm = as_name(contrast_q)
  speaker_user_supplied = !quo_is_null(speaker_q)
  if (speaker_user_supplied) {
    speaker_nm = as_name(speaker_q)
  } else {
    speaker_nm = ".speaker_dummy_internal"
    xdata = xdata %>% mutate(!!speaker_nm := "ALL")
    speaker_q = rlang::sym(speaker_nm)
  }
  speaker_spec_nm = if (speaker_user_supplied) speaker_nm else NULL
  group_nm    = if (quo_is_null(group_q)) NA_character_ else as_name(group_q)
  token_id_nm = ".token_rowid"
  
  # ---- condition vars -> .condition_id ----
  cond_sel = tidyselect::eval_select(enquo(condition_vars), xdata)
  
  if (length(cond_sel) == 0) {
    xdata = xdata %>% mutate(.condition_id = "ALL")
    cond_nm = character(0)
  } else {
    cond_nm = names(cond_sel)
    xdata = xdata %>% tidyr::unite(".condition_id", all_of(cond_nm), sep = "×", remove = FALSE)
  }
  
  condition_nm = ".condition_id"
  condition_q  = rlang::sym(condition_nm)
  condition_vars_nm = cond_nm
  
  vars_to_char = c(
    contrast_nm, condition_nm, speaker_nm, token_id_nm,
    if (!is.na(group_nm)) group_nm else NULL
  ) %>% unique() %>% na.omit()
  
  xdata = xdata %>% mutate(across(all_of(vars_to_char), as.character))
  
  rm(vars_to_char, cond_sel)
  
  if (!inter_group) {
    within_ref <- match.arg(within_ref, c("self_speaker", "inter_speaker_mean", "population_pooled"))
  }
  
  contrast_levels = xdata %>%
    distinct(!!contrast_q) %>%
    pull(!!contrast_q) %>%
    as.character()
  
  if (length(contrast_levels) != 2) {
    stop("`contrast_var` must have exactly 2 levels in the data passed to voweldist().")
  }
  
  lvl1 = contrast_levels[1]
  lvl2 = contrast_levels[2]
  
  # ---- sanity checks for BETWEEN ----
  auto_two_group = FALSE
  two_groups = NULL
  
  if (inter_group) {
    if (is.na(group_nm)) {
      stop("inter_group=TRUE requires `group_var`.")
    }
    
    present_groups = xdata %>%
      distinct(!!group_q) %>%
      pull(!!group_q) %>%
      as.character()
    
    if (is.null(reference_group)) {
      if (!is.null(focal_group)) {
        ref_guess = setdiff(present_groups, focal_group)
        if (length(ref_guess) != 1) {
          stop("Between-group: reference_group is NULL and focal_group does not uniquely determine the reference group.")
        }
        reference_group = ref_guess
        focal_groups = focal_group
      } else {
        if (length(present_groups) == 2) {
          auto_two_group = TRUE
          two_groups = present_groups
        } else {
          stop("Between-group: reference_group is NULL and group_var has != 2 levels. Provide reference_group.")
        }
      }
    } else {
      if (!reference_group %in% present_groups) {
        stop(paste0("Between-group: reference_group='", reference_group, "' not found in `", group_nm, "`."))
      }
      if (is.null(focal_group)) {
        focal_groups = setdiff(present_groups, reference_group)
        if (length(focal_groups) == 0) {
          stop("Between-group: no focal groups found.")
        }
      } else {
        if (!focal_group %in% present_groups) {
          stop(paste0("Between-group: focal_group='", focal_group, "' not found in `", group_nm, "`."))
        }
        focal_groups = focal_group
      }
    }
  }
  
  # =========================
  # Helpers
  # =========================
  

  .cov_matrix_is_usable <- function(S, require_invertible = FALSE) {
    if (is.null(S)) return(FALSE)
    if (is.null(dim(S))) S <- matrix(S, 1, 1)
    if (anyNA(S) || any(!is.finite(S))) return(FALSE)
    
    if (isTRUE(require_invertible)) {
      return(tryCatch({
        chol(S)
        TRUE
      }, error = function(e) FALSE))
    }
    
    TRUE
  }
  
  .cloud_cache <- new.env(parent = emptyenv())
  .resolved_cloud_cache <- new.env(parent = emptyenv())
  
  .short_hash_chr <- function(x) {
    if (length(x) == 0) return("0")
    raw_x <- serialize(x, NULL, ascii = FALSE)
    ints <- as.integer(raw_x)
    h1 <- 2166136261
    h2 <- 16777619
    mod <- 2147483647
    for (b in ints) {
      h1 <- (h1 + b + 1) %% mod
      h1 <- (h1 * h2) %% mod
    }
    paste0(format(length(ints), scientific = FALSE), "_", sprintf("%08x", as.integer(h1)))
  }

  .cloud_cache_key <- function(df_cloud, require_invertible = FALSE, robust = TRUE) {
    if (is.null(df_cloud) || !is.data.frame(df_cloud)) return(NULL)
    id_cols <- intersect(c(token_id_nm, ".token_rowid"), names(df_cloud))
    
    if (length(id_cols) > 0) {
      ids <- sort(unique(as.character(df_cloud[[id_cols[1]]])))
      id_part <- .short_hash_chr(ids)
    } else {
      mat <- as.matrix(dplyr::select(df_cloud, dplyr::all_of(measure_nm)))
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
      if (nrow(mat) == 0) return(NULL)
      rows_sig <- apply(signif(mat, 12), 1, paste, collapse = ",")
      id_part <- .short_hash_chr(rows_sig)
    }
    
    paste0(
      "inv=", require_invertible,
      ";rob=", robust,
      ";meas=", paste(measure_nm, collapse = ","),
      ";ids_hash=", id_part
    )
  }
  
  .compute_cloud_stats <- function(df_cloud, measure_nm, min_n = min_n_cloud, robust = TRUE, require_invertible = FALSE) {
    res <- list(
      ok = FALSE,
      reason = NA_character_,
      center = NULL,
      cov = NULL,
      chol = NULL,
      cov_inv = NULL,
      logdet = NA_real_,
      det = NA_real_,
      n = 0L,
      data = NULL
    )
    
    if (is.null(df_cloud) || !is.data.frame(df_cloud)) {
      res$reason <- "null_cloud"
      return(res)
    }
    if (!all(measure_nm %in% names(df_cloud))) {
      res$reason <- "missing_measure_columns"
      return(res)
    }
    
    mat <- as.matrix(dplyr::select(df_cloud, dplyr::all_of(measure_nm)))
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    res$n <- nrow(mat)
    res$data <- as.data.frame(mat)
    
    if (nrow(mat) == 0) {
      res$reason <- "no_complete_cases"
      return(res)
    }
    if (nrow(mat) < min_n) {
      res$reason <- paste0("too_few_observations(n=", nrow(mat), ", min=", min_n, ")")
      return(res)
    }
    if (nrow(mat) <= 1) {
      res$reason <- "n_le_1"
      return(res)
    }
    
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
    if (any(!is.finite(sds))) {
      res$reason <- "sd_not_finite"
      return(res)
    }
    if (any(sds == 0)) {
      res$reason <- "zero_variance"
      return(res)
    }
    
    fit <- NULL
    if (isTRUE(robust) && requireNamespace("MASS", quietly = TRUE)) {
      fit <- tryCatch(
        suppressWarnings(MASS::cov.trob(mat)),
        error = function(e) NULL
      )
    }
    
    if (!is.null(fit) && !is.null(fit$cov) && !is.null(fit$center) &&
        all(is.finite(fit$center)) &&
        .cov_matrix_is_usable(fit$cov, require_invertible = FALSE)) {
      center <- fit$center
      S <- fit$cov
    } else {
      center <- colMeans(mat)
      S <- tryCatch(stats::cov(mat), error = function(e) NULL)
    }
    
    if (is.null(S) || anyNA(S) || any(!is.finite(S))) {
      res$reason <- "covariance_not_finite"
      return(res)
    }
    if (is.null(dim(S))) S <- matrix(S, 1, 1)
    
    cholS <- tryCatch(chol(S), error = function(e) NULL)
    if (isTRUE(require_invertible) && is.null(cholS)) {
      res$reason <- "covariance_not_invertible"
      return(res)
    }
    
    if (!is.null(cholS)) {
      logdet <- 2 * sum(log(diag(cholS)))
      cov_inv <- tryCatch(chol2inv(cholS), error = function(e) NULL)
      detS <- exp(logdet)
    } else {
      logdet <- tryCatch(as.numeric(determinant(S, logarithm = TRUE)$modulus), error = function(e) NA_real_)
      cov_inv <- tryCatch(solve(S), error = function(e) NULL)
      detS <- tryCatch(as.numeric(det(S)), error = function(e) NA_real_)
    }
    
    res$ok <- TRUE
    res$reason <- NA_character_
    res$center <- center
    res$cov <- S
    res$chol <- cholS
    res$cov_inv <- cov_inv
    res$logdet <- logdet
    res$det <- detS
    res
  }
  
  .compute_center_only_stats <- function(df_cloud, measure_nm) {
    res <- list(
      ok = FALSE,
      reason = NA_character_,
      center = NULL,
      n = 0L,
      data = NULL
    )
    
    if (is.null(df_cloud) || !is.data.frame(df_cloud)) {
      res$reason <- "null_cloud"
      return(res)
    }
    if (!all(measure_nm %in% names(df_cloud))) {
      res$reason <- "missing_measure_columns"
      return(res)
    }
    
    mat <- as.matrix(dplyr::select(df_cloud, dplyr::all_of(measure_nm)))
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    res$n <- nrow(mat)
    res$data <- as.data.frame(mat)
    
    if (nrow(mat) == 0) {
      res$reason <- "no_complete_cases"
      return(res)
    }
    
    center <- colMeans(mat)
    if (any(!is.finite(center))) {
      res$reason <- "center_not_finite"
      return(res)
    }
    
    res$ok <- TRUE
    res$reason <- NA_character_
    res$center <- center
    res
  }
  
  .resolve_center_only_cloud <- function(pool, row_meta, spec_list, measure_nm,
                                         contrast_value = NULL,
                                         exclude_same_speaker = FALSE) {
    pool0 <- pool
    
    if (!is.null(contrast_value) && contrast_nm %in% names(pool0)) {
      if (is.na(contrast_value)) {
        pool0 <- pool0 %>% dplyr::filter(is.na(.data[[contrast_nm]]))
      } else {
        pool0 <- pool0 %>% dplyr::filter(.data[[contrast_nm]] == contrast_value)
      }
    }
    
    failure_log <- character(0)
    
    for (spec_name in names(spec_list)) {
      spec_cols <- spec_list[[spec_name]]
      cand <- .subset_pool_by_spec(pool0, row_meta, spec_cols)
      
      if (isTRUE(exclude_same_speaker) &&
          !is.null(speaker_nm) &&
          speaker_nm %in% names(cand) &&
          speaker_nm %in% names(row_meta)) {
        cand <- cand %>% dplyr::filter(.data[[speaker_nm]] != row_meta[[speaker_nm]][1])
      }
      
      stats_obj <- .compute_center_only_stats(cand, measure_nm)
      
      if (isTRUE(stats_obj$ok)) {
        return(list(
          data = cand,
          spec_name = spec_name,
          spec_cols = spec_cols,
          effective_spec_name = .effective_spec_name(spec_cols, row_meta),
          failure_reason = NA_character_,
          stats = stats_obj
        ))
      }
      
      failure_log <- c(failure_log, paste0(spec_name, ": ", stats_obj$reason))
    }
    
    list(
      data = NULL,
      spec_name = NA_character_,
      spec_cols = NULL,
      effective_spec_name = NA_character_,
      failure_reason = paste(failure_log, collapse = " || "),
      stats = NULL
    )
  }
  
  .get_cloud_stats <- function(df_cloud, measure_nm, min_n = min_n_cloud, robust = TRUE, require_invertible = FALSE) {
    key <- .cloud_cache_key(df_cloud, require_invertible = require_invertible, robust = robust)
    if (!is.null(key) && exists(key, envir = .cloud_cache, inherits = FALSE)) {
      return(get(key, envir = .cloud_cache, inherits = FALSE))
    }
    
    stats <- .compute_cloud_stats(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n,
      robust = robust,
      require_invertible = require_invertible
    )
    
    if (!is.null(key)) assign(key, stats, envir = .cloud_cache)
    stats
  }

  .cache_value_string <- function(x) {
    if (length(x) == 0) return("")
    if (all(is.na(x))) return("<NA>")
    paste(ifelse(is.na(x), "<NA>", as.character(x)), collapse = "|")
  }
  
  .row_cache_key <- function(row_meta, cols) {
    cols <- unique(cols)
    cols <- cols[!is.na(cols)]
    cols <- cols[cols %in% names(row_meta)]
    if (length(cols) == 0) return("row=<global>")
    pieces <- vapply(cols, function(col) {
      paste0(col, "=", .cache_value_string(row_meta[[col]]))
    }, FUN.VALUE = character(1))
    paste(c("row", pieces), collapse = ";")
  }
  
  .resolve_cloud_detail <- function(pool, row_meta, spec_list, measure_nm, contrast_value = NULL, exclude_same_speaker = FALSE, require_invertible = FALSE) {
    pool0 <- pool
    
    if (!is.null(contrast_value) && contrast_nm %in% names(pool0)) {
      if (is.na(contrast_value)) {
        pool0 <- pool0 %>% dplyr::filter(is.na(.data[[contrast_nm]]))
      } else {
        pool0 <- pool0 %>% dplyr::filter(.data[[contrast_nm]] == contrast_value)
      }
    }
    
    failure_log <- character(0)
    
    for (spec_name in names(spec_list)) {
      spec_cols <- spec_list[[spec_name]]
      cand <- .subset_pool_by_spec(pool0, row_meta, spec_cols)
      
      if (isTRUE(exclude_same_speaker) &&
          !is.null(speaker_nm) &&
          speaker_nm %in% names(cand) &&
          speaker_nm %in% names(row_meta)) {
        cand <- cand %>% dplyr::filter(.data[[speaker_nm]] != row_meta[[speaker_nm]][1])
      }
      
      stats_obj <- .get_cloud_stats(
        df_cloud = cand,
        measure_nm = measure_nm,
        min_n = min_n_cloud,
        robust = TRUE,
        require_invertible = require_invertible
      )
      
      if (isTRUE(stats_obj$ok)) {
        return(list(
          data = cand,
          spec_name = spec_name,
          spec_cols = spec_cols,
          effective_spec_name = .effective_spec_name(spec_cols, row_meta),
          failure_reason = NA_character_,
          stats = stats_obj
        ))
      }
      
      failure_log <- c(failure_log, paste0(spec_name, ": ", stats_obj$reason))
    }
    
    list(
      data = NULL,
      spec_name = NA_character_,
      spec_cols = NULL,
      effective_spec_name = NA_character_,
      failure_reason = paste(failure_log, collapse = " || "),
      stats = NULL
    )
  }
  
  .resolve_within_cloud_cached <- function(base, row_meta, contrast_value, within_ref, require_invertible = FALSE, exclude_same_speaker = FALSE) {
    row_group <- if (!is.na(group_nm) && group_nm %in% names(row_meta)) {
      as.character(row_meta[[group_nm]][1])
    } else {
      NA_character_
    }
    
    specs_here <- .make_within_reference_specs(row_group, within_ref)
    cache_cols <- unique(c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm, condition_nm, condition_vars_nm))
    key <- paste(
      "within",
      paste0("contrast=", ifelse(is.null(contrast_value), "<NULL>", as.character(contrast_value))),
      paste0("within_ref=", within_ref),
      paste0("require_invertible=", require_invertible),
      paste0("exclude_same_speaker=", exclude_same_speaker),
      paste0("specs=", paste(names(specs_here), collapse = ">")),
      .row_cache_key(row_meta, cache_cols),
      sep = ";"
    )
    
    if (exists(key, envir = .resolved_cloud_cache, inherits = FALSE)) {
      return(get(key, envir = .resolved_cloud_cache, inherits = FALSE))
    }
    
    pool_base <- base
    if (!is.na(group_nm) && group_nm %in% names(base) && !is.na(row_group)) {
      pool_base <- pool_base %>% dplyr::filter(.data[[group_nm]] == row_group)
    }
    
    out <- .resolve_cloud_detail(
      pool = pool_base,
      row_meta = row_meta,
      spec_list = specs_here,
      measure_nm = measure_nm,
      contrast_value = contrast_value,
      exclude_same_speaker = exclude_same_speaker,
      require_invertible = require_invertible
    )
    
    assign(key, out, envir = .resolved_cloud_cache)
    out
  }
  
  .cov_estimation_reason <- function(df_cloud, measure_nm, min_n = min_n_cloud, robust = TRUE, require_invertible = FALSE) {
    .get_cloud_stats(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n,
      robust = robust,
      require_invertible = require_invertible
    )$reason
  }
  
  .cloud_center_cov <- function(df_cloud, measure_nm, min_n = min_n_cloud, robust = TRUE, require_invertible = FALSE, center_name = c("center", "mu"), cov_name = c("cov", "S")) {
    
    center_name <- match.arg(center_name)
    cov_name <- match.arg(cov_name)
    
    stats_obj <- .get_cloud_stats(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n,
      robust = robust,
      require_invertible = require_invertible
    )
    
    if (!isTRUE(stats_obj$ok)) return(NULL)
    
    out <- list()
    out[[center_name]] <- stats_obj$center
    out[[cov_name]] <- stats_obj$cov
    out$n <- stats_obj$n
    out$data <- stats_obj$data
    out$chol <- stats_obj$chol
    out$cov_inv <- stats_obj$cov_inv
    out$logdet <- stats_obj$logdet
    out$det <- stats_obj$det
    out
  }
  
  .mah_cov <- function(df_cloud) {
    .cloud_center_cov(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n_cloud,
      robust = TRUE,
      require_invertible = TRUE,
      center_name = "center",
      cov_name = "cov"
    )
  }
  
  .mu_cov <- function(df_cloud) {
    .cloud_center_cov(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n_cloud,
      robust = TRUE,
      require_invertible = FALSE,
      center_name = "mu",
      cov_name = "S"
    )
  }
  
  .resolve_cloud_generic <- function(pool, row_meta, spec_list, measure_nm, contrast_value = NULL, exclude_same_speaker = FALSE, validator, postprocess = identity) {
    
    pool0 <- pool
    
    if (!is.null(contrast_value) && contrast_nm %in% names(pool0)) {
      if (is.na(contrast_value)) {
        pool0 <- pool0 %>% dplyr::filter(is.na(.data[[contrast_nm]]))
      } else {
        pool0 <- pool0 %>% dplyr::filter(.data[[contrast_nm]] == contrast_value)
      }
    }
    
    for (spec_name in names(spec_list)) {
      spec_cols <- spec_list[[spec_name]]
      
      cand <- .subset_pool_by_spec(pool0, row_meta, spec_cols)
      
      if (isTRUE(exclude_same_speaker) &&
          !is.null(speaker_nm) &&
          speaker_nm %in% names(cand) &&
          speaker_nm %in% names(row_meta)) {
        cand <- cand %>% dplyr::filter(.data[[speaker_nm]] != row_meta[[speaker_nm]][1])
      }
      
      ok <- validator(cand)
      if (isTRUE(ok)) {
        return(postprocess(list(
          data = cand,
          spec_name = spec_name,
          spec_cols = spec_cols,
          effective_spec_name = .effective_spec_name(spec_cols, row_meta)
        )))
      }
    }
    
    NULL
  }
  
  .resolve_cloud <- function(pool, row_meta, spec_list, measure_nm, contrast_value = NULL, exclude_same_speaker = FALSE, require_invertible = FALSE) {
    
    .resolve_cloud_generic(
      pool = pool,
      row_meta = row_meta,
      spec_list = spec_list,
      measure_nm = measure_nm,
      contrast_value = contrast_value,
      exclude_same_speaker = exclude_same_speaker,
      validator = function(cand) {
        .cloud_viable(
          cand,
          measure_nm = measure_nm,
          min_n = min_n_cloud,
          require_invertible = require_invertible
        )
      }
    )
  }
  
  .resolve_mah_cloud <- function(pool, row_meta, spec_list, measure_nm, contrast_value = NULL, exclude_same_speaker = FALSE) {
    
    .resolve_cloud_generic(
      pool = pool,
      row_meta = row_meta,
      spec_list = spec_list,
      measure_nm = measure_nm,
      contrast_value = contrast_value,
      exclude_same_speaker = exclude_same_speaker,
      validator = function(cand) {
        !is.null(.mah_cov(cand))
      }
    )
  }
  
  .spec_label <- function(cols) {
    if (length(cols) == 0) return("global")
    paste(cols, collapse = " * ")
  }
  
  .make_specs <- function(condition_vars_nm, speaker_nm = NULL, include_speaker = FALSE) {
    specs <- list()
    
    if (length(condition_vars_nm) > 0) {
      if (isTRUE(include_speaker) && !is.null(speaker_nm)) {
        specs[[length(specs) + 1]] <- c(condition_vars_nm, speaker_spec_nm)
      }
      specs[[length(specs) + 1]] <- condition_vars_nm
      
      if (length(condition_vars_nm) > 1) {
        for (i in seq(length(condition_vars_nm) - 1, 1, by = -1)) {
          specs[[length(specs) + 1]] <- condition_vars_nm[1:i]
        }
      }
    } else {
      if (isTRUE(include_speaker) && !is.null(speaker_nm)) {
        specs[[length(specs) + 1]] <- speaker_spec_nm
      }
    }
    
    specs[[length(specs) + 1]] <- character(0)
    keep <- !duplicated(vapply(specs, paste, collapse = "|", FUN.VALUE = character(1)))
    specs <- specs[keep]
    names(specs) <- vapply(specs, .spec_label, FUN.VALUE = character(1))
    specs
  }
  
  .is_reference_group <- function(row_group) {
    !is.na(group_nm) &&
      !is.null(reference_group) &&
      !is.na(row_group) &&
      identical(as.character(row_group), as.character(reference_group))
  }
  
  .make_within_reference_specs <- function(row_group, within_ref) {
    if (within_ref == "self_speaker") {
      if (.is_reference_group(row_group)) {
        specs <- list()
        
        if (length(condition_vars_nm) > 0) {
          specs[[length(specs) + 1]] <- c(condition_vars_nm, speaker_spec_nm)
          specs[[length(specs) + 1]] <- condition_vars_nm
          
          if (length(condition_vars_nm) > 1) {
            for (i in seq(length(condition_vars_nm) - 1, 1, by = -1)) {
              specs[[length(specs) + 1]] <- condition_vars_nm[1:i]
            }
          }
        } else {
          specs[[length(specs) + 1]] <- speaker_spec_nm
        }
        
        specs[[length(specs) + 1]] <- character(0)
        
        keep <- !duplicated(vapply(specs, paste, collapse = "|", FUN.VALUE = character(1)))
        specs <- specs[keep]
        names(specs) <- vapply(specs, .spec_label, FUN.VALUE = character(1))
        return(specs)
        
      } else {
        return(.make_specs(condition_vars_nm, speaker_nm = speaker_spec_nm, include_speaker = TRUE))
      }
    } else {
      return(.make_specs(condition_vars_nm, speaker_nm = speaker_spec_nm, include_speaker = FALSE))
    }
  }
  
  .make_between_specs <- function(condition_vars_nm) {
    .make_specs(condition_vars_nm, speaker_nm = NULL, include_speaker = FALSE)
  }
  
  .make_within_category_specs <- function(row_group) {
    if (.is_reference_group(row_group)) {
      .make_specs(condition_vars_nm, speaker_nm = speaker_spec_nm, include_speaker = FALSE)
    } else {
      .make_specs(condition_vars_nm, speaker_nm = speaker_spec_nm, include_speaker = TRUE)
    }
  }
  
  .effective_spec_name <- function(spec_cols, row_meta) {
    if (length(spec_cols) == 0) return("global")
    
    spec_cols <- spec_cols[spec_cols %in% names(row_meta)]
    if (length(spec_cols) == 0) return("global")
    
    keep <- spec_cols[!purrr::map_lgl(spec_cols, ~ is.na(row_meta[[.x]][1]))]
    if (length(keep) == 0) return("global")
    
    paste(keep, collapse = " * ")
  }
  
  .subset_pool_by_spec <- function(pool, row_meta, spec_cols) {
    out <- pool
    if (length(spec_cols) == 0) return(out)
    
    spec_cols <- spec_cols[spec_cols %in% names(pool) & spec_cols %in% names(row_meta)]
    if (length(spec_cols) == 0) return(out)
    
    for (v in spec_cols) {
      val <- row_meta[[v]][1]
      if (is.na(val)) {
        out <- out %>% dplyr::filter(is.na(.data[[v]]))
      } else {
        out <- out %>% dplyr::filter(.data[[v]] == val)
      }
    }
    out
  }
  
  .emit_effective_clouds <- function(prefix, values, sep = " || ") {
    
    values <- values[!is.na(values)]
    values <- unique(values)
    values <- stringr::str_replace_all(values, fixed("*"), "×")
    
    if (length(values) > 0) {
      message(prefix, paste(values, collapse = sep))
    }
    
    invisible(values)
  }
  
  .emit_clouds_used <- function(prefix, values, sep = " || ") {
    .emit_effective_clouds(prefix, values, sep = sep)
  }


  .expected_cloud_name <- function(include_speaker = TRUE) {
    cols <- condition_vars_nm
    if (isTRUE(include_speaker) && !is.null(speaker_nm)) {
      cols <- c(cols, speaker_nm)
    }
    .spec_label(cols)
  }

  .clean_cloud_label <- function(x) {
    x <- as.character(x)
    x <- gsub("^[^:]+: ", "", x)
    x
  }

  .emit_cloud_fallbacks <- function(prefix, used_clouds, expected_cloud) {
    used_clouds <- used_clouds[!is.na(used_clouds)]
    if (length(used_clouds) == 0) return(invisible(NULL))
    used_clouds <- unique(.clean_cloud_label(used_clouds))
    fallback <- used_clouds[used_clouds != expected_cloud]
    if (length(fallback) > 0) {
      message(
        prefix,
        "Expected cloud: ", expected_cloud,
        ". Effective cloud(s): ",
        paste(fallback, collapse = " || ")
      )
    }
    invisible(fallback)
  }

  .cov_is_invertible <- function(df_cloud, measure_nm, min_n = min_n_cloud) {
    stats_obj <- .get_cloud_stats(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n,
      robust = TRUE,
      require_invertible = TRUE
    )
    isTRUE(stats_obj$ok)
  }
  
  .cloud_viable <- function(df_cloud, measure_nm, min_n = min_n_cloud, require_invertible = FALSE) {
    stats_obj <- .get_cloud_stats(
      df_cloud = df_cloud,
      measure_nm = measure_nm,
      min_n = min_n,
      robust = TRUE,
      require_invertible = require_invertible
    )
    isTRUE(stats_obj$ok)
  }
  
  .emit_failure_reasons <- function(prefix, reasons) {
    reasons <- reasons[!is.na(reasons)]
    if (length(reasons) == 0) return(invisible(NULL))
    
    tb <- sort(table(reasons), decreasing = TRUE)
    msg <- paste(names(tb), as.integer(tb), sep = "=", collapse = "; ")
    message(prefix, msg)
    invisible(tb)
  }
  
  .emit_subject_failure_summary <- function(df_fail, speaker_col = speaker_nm, group_col = group_nm, max_show = 8) {
    if (is.null(df_fail) || !is.data.frame(df_fail) || nrow(df_fail) == 0) return(invisible(NULL))
    
    df_show <- df_fail %>% dplyr::slice_head(n = max_show)
    detail_lines <- apply(df_show, 1, function(row) {
      grp_txt <- if (!is.na(group_col) && group_col %in% names(df_show)) paste0(row[[group_col]], " / ") else ""
      paste0("- ", grp_txt, row[[speaker_col]], ": ", row[["failure_reason"]])
    })
    
    message(
      "Mahalanobis: impossible for ", nrow(df_fail), " subject(s). First ", nrow(df_show),
      ":
", paste(detail_lines, collapse = "
"),
      "
Suggestion: compute these subject(s) separately, or set their group as the reference via `group_var` and `reference_group`."
    )
    
    invisible(df_fail)
  }
  
  .resolve_mah_cloud_detail <- function(pool, row_meta, spec_list, measure_nm, contrast_value = NULL, exclude_same_speaker = FALSE) {
    .resolve_cloud_detail(
      pool = pool,
      row_meta = row_meta,
      spec_list = spec_list,
      measure_nm = measure_nm,
      contrast_value = contrast_value,
      exclude_same_speaker = exclude_same_speaker,
      require_invertible = TRUE
    )
  }
  
  .mah_one <- function(data.x, data.y = NULL, stats_y = NULL) {
    fail <- function(reason) {
      res <- NA_real_
      attr(res, "reason") <- reason
      res
    }
    
    if (is.null(data.x)) return(fail("null_input"))
    if (!is.data.frame(data.x)) return(fail("non_dataframe_input"))
    
    dx <- data.x %>%
      dplyr::select(dplyr::all_of(measure_nm)) %>%
      stats::na.omit()
    
    if (nrow(dx) != 1) return(fail("token_not_single_row"))
    
    if (is.null(stats_y)) {
      if (is.null(data.y) || !is.data.frame(data.y)) return(fail("non_dataframe_input"))
      stats_y <- .get_cloud_stats(
        df_cloud = data.y,
        measure_nm = measure_nm,
        min_n = min_n_cloud,
        robust = TRUE,
        require_invertible = TRUE
      )
    }
    
    if (!isTRUE(stats_y$ok)) return(fail(stats_y$reason))
    
    z <- as.numeric(as.matrix(dx) - matrix(stats_y$center, nrow = 1))
    
    out <- tryCatch({
      val <- z %*% stats_y$cov_inv %*% z
      as.numeric(val[1, 1])
    }, error = function(e) fail(paste0("mahalanobis_error: ", conditionMessage(e))))
    
    if (length(out) == 1 && is.numeric(out) && !is.na(out) && !is.nan(out)) {
      out <- as.numeric(out)
      if (!is.finite(out) || out < 0) return(fail("mahalanobis_non_finite"))
      return(sqrt(out))
    }
    
    out
  }
  
  .cloud_clean <- function(df) {
    df %>%
      dplyr::select(dplyr::all_of(measure_nm)) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x)))
  }
  
  .bhatt <- function(cloudA, cloudB, statsA = NULL, statsB = NULL) {
    A <- .cloud_clean(cloudA)
    B <- .cloud_clean(cloudB)
    if (nrow(A) == 0 || nrow(B) == 0) return(NA_real_)
    
    .as_bhatt_stats <- function(stats_obj, df_cloud) {
      if (is.null(stats_obj)) {
        return(.mu_cov(df_cloud))
      }
      if (!is.null(stats_obj$mu) && !is.null(stats_obj$S)) {
        return(stats_obj)
      }
      if (!is.null(stats_obj$center) && !is.null(stats_obj$cov)) {
        return(list(
          mu = stats_obj$center,
          S = stats_obj$cov,
          logdet = stats_obj$logdet,
          det = stats_obj$det,
          n = stats_obj$n,
          data = stats_obj$data,
          chol = stats_obj$chol,
          cov_inv = stats_obj$cov_inv
        ))
      }
      NULL
    }
    
    statsA <- .as_bhatt_stats(statsA, A)
    statsB <- .as_bhatt_stats(statsB, B)
    if (is.null(statsA) || is.null(statsB)) return(NA_real_)
    
    mu1 <- statsA$mu; S1 <- statsA$S
    mu2 <- statsB$mu; S2 <- statsB$S
    S    <- (S1 + S2) / 2
    
    cholS <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(cholS)) return(NA_real_)
    
    invS <- tryCatch(chol2inv(cholS), error = function(e) NULL)
    if (is.null(invS)) return(NA_real_)
    
    dmu <- matrix(mu2 - mu1, ncol = 1)
    term1 <- as.numeric((t(dmu) %*% invS %*% dmu) / 8)
    
    logdetS  <- 2 * sum(log(diag(cholS)))
    logdetS1 <- statsA$logdet
    logdetS2 <- statsB$logdet
    
    if (any(!is.finite(c(logdetS, logdetS1, logdetS2)))) return(NA_real_)
    
    term2 <- 0.5 * (logdetS - 0.5 * (logdetS1 + logdetS2))
    out <- term1 + term2
    if (!is.finite(out)) NA_real_ else out
  }
  
  .pillai_reason <- function(df_ab, class_col) {
    df_ab <- df_ab %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
      dplyr::filter(!is.na(.data[[class_col]])) %>%
      dplyr::mutate(.class_tmp = as.character(.data[[class_col]]))
    
    levs <- unique(df_ab$.class_tmp)
    if (length(levs) != 2) return("n_classes_ne_2")
    
    df1 <- df_ab %>%
      dplyr::filter(.class_tmp == levs[1]) %>%
      dplyr::select(dplyr::all_of(measure_nm))
    
    df2 <- df_ab %>%
      dplyr::filter(.class_tmp == levs[2]) %>%
      dplyr::select(dplyr::all_of(measure_nm))
    
    if (nrow(df1) == 0 || nrow(df2) == 0) return("empty_class")
    
    X1 <- as.matrix(df1)
    X2 <- as.matrix(df2)
    
    if (is.null(dim(X1))) X1 <- matrix(X1, ncol = 1)
    if (is.null(dim(X2))) X2 <- matrix(X2, ncol = 1)
    
    p  <- ncol(X1)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    
    if (n1 < 2 || n2 < 2) return("class_with_lt_2_observations")
    if ((n1 + n2) < (p + 2)) return("insufficient_observations")
    if (any(!is.finite(X1)) || any(!is.finite(X2))) return("non_finite_values")
    
    sds1 <- apply(X1, 2, stats::sd, na.rm = TRUE)
    sds2 <- apply(X2, 2, stats::sd, na.rm = TRUE)
    
    if (any(!is.finite(sds1)) || any(!is.finite(sds2))) return("sd_not_finite")
    if (all(sds1 == 0) && all(sds2 == 0)) return("both_groups_constant")
    
    S1 <- tryCatch(stats::cov(X1), error = function(e) NULL)
    S2 <- tryCatch(stats::cov(X2), error = function(e) NULL)
    
    if (is.null(S1) || is.null(S2)) return("covariance_error")
    
    if (is.null(dim(S1))) {
      S1 <- matrix(S1, 1, 1)
      S2 <- matrix(S2, 1, 1)
    }
    
    if (anyNA(S1) || anyNA(S2) || any(!is.finite(S1)) || any(!is.finite(S2))) {
      return("covariance_not_finite")
    }
    
    NA_character_
  }
  
  .pillai_two_clouds <- function(df_ab, class_col, debug = FALSE, debug_id = NULL) {
    fail <- function(reason) {
      res <- NA_real_
      attr(res, "reason") <- reason
      res
    }
    
    reason <- .pillai_reason(df_ab, class_col)
    if (!is.na(reason)) {
      return(fail(reason))
    }
    
    df_ab <- df_ab %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
      dplyr::filter(!is.na(.data[[class_col]])) %>%
      dplyr::mutate(.class_tmp = as.character(.data[[class_col]]))
    
    levs <- unique(df_ab$.class_tmp)
    if (length(levs) != 2) {
      return(fail("n_classes_ne_2"))
    }
    
    X1 <- as.matrix(df_ab[df_ab$.class_tmp == levs[1], measure_nm, drop = FALSE])
    X2 <- as.matrix(df_ab[df_ab$.class_tmp == levs[2], measure_nm, drop = FALSE])
    
    if (is.null(dim(X1))) X1 <- matrix(X1, ncol = 1)
    if (is.null(dim(X2))) X2 <- matrix(X2, ncol = 1)
    
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    
    S1 <- tryCatch(stats::cov(X1), error = function(e) NULL)
    S2 <- tryCatch(stats::cov(X2), error = function(e) NULL)
    
    if (is.null(S1) || is.null(S2)) {
      return(fail("covariance_error"))
    }
    
    if (is.null(dim(S1))) {
      S1 <- matrix(S1, 1, 1)
      S2 <- matrix(S2, 1, 1)
    }
    
    if (anyNA(S1) || anyNA(S2) || any(!is.finite(S1)) || any(!is.finite(S2))) {
      return(fail("covariance_not_finite"))
    }
    
    E <- (n1 - 1) * S1 + (n2 - 1) * S2
    d <- matrix(colMeans(X1) - colMeans(X2), ncol = 1)
    H <- (n1 * n2) / (n1 + n2) * (d %*% t(d))
    
    if (anyNA(E) || any(!is.finite(E))) {
      return(fail("E_not_finite"))
    }
    if (anyNA(H) || any(!is.finite(H))) {
      return(fail("H_not_finite"))
    }
    
    Tmat <- E + H
    
    if (anyNA(Tmat) || any(!is.finite(Tmat))) {
      return(fail("Tmat_not_finite"))
    }
    
    eps <- 1e-8
    Tmat <- Tmat + diag(eps, nrow(Tmat))
    
    if (anyNA(Tmat) || any(!is.finite(Tmat))) {
      return(fail("Tmat_not_finite_after_regularization"))
    }
    
    s <- tryCatch(svd(Tmat), error = function(e) NULL)
    if (is.null(s)) {
      return(fail("svd_error"))
    }
    
    if (anyNA(s$d) || any(!is.finite(s$d))) {
      return(fail("svd_values_not_finite"))
    }
    
    Tinv <- tryCatch(
      s$v %*% (diag(ifelse(s$d > 1e-9, 1 / s$d, 0), nrow = length(s$d))) %*% t(s$u),
      error = function(e) NULL
    )
    
    if (is.null(Tinv) || anyNA(Tinv) || any(!is.finite(Tinv))) {
      return(fail("Tinv_not_finite"))
    }
    
    V <- tryCatch(sum(diag(H %*% Tinv)), error = function(e) NA_real_)
    
    if (is.na(V) || is.nan(V) || !is.finite(V)) {
      return(fail("pillai_not_finite"))
    }
    
    as.numeric(max(0, min(V, 1)))
  }


  .pillai_cloud_reason <- function(df_cloud, contrast_col) {
    if (is.null(df_cloud) || !is.data.frame(df_cloud)) return("null_cloud")
    if (!all(c(measure_nm, contrast_col) %in% names(df_cloud))) return("missing_columns")

    df_use <- df_cloud %>%
      dplyr::select(dplyr::all_of(c(measure_nm, contrast_col))) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
      dplyr::filter(!is.na(.data[[contrast_col]]))

    .pillai_reason(df_use, contrast_col)
  }

  .resolve_pillai_cloud_detail <- function(pool, row_meta, spec_list, contrast_col,
                                           exclude_same_speaker = FALSE) {
    failure_log <- character(0)

    for (spec_name in names(spec_list)) {
      spec_cols <- spec_list[[spec_name]]
      cand <- .subset_pool_by_spec(pool, row_meta, spec_cols)

      if (isTRUE(exclude_same_speaker) &&
          !is.null(speaker_nm) &&
          !is.na(speaker_nm) &&
          speaker_nm %in% names(cand) &&
          speaker_nm %in% names(row_meta)) {
        cand <- cand %>% dplyr::filter(.data[[speaker_nm]] != row_meta[[speaker_nm]][1])
      }

      reason <- .pillai_cloud_reason(cand, contrast_col)
      if (is.na(reason)) {
        return(list(
          data = cand,
          spec_name = spec_name,
          spec_cols = spec_cols,
          effective_spec_name = .effective_spec_name(spec_cols, row_meta),
          failure_reason = NA_character_
        ))
      }
      failure_log <- c(failure_log, paste0(spec_name, ": ", reason))
    }

    list(
      data = NULL,
      spec_name = NA_character_,
      spec_cols = NULL,
      effective_spec_name = NA_character_,
      failure_reason = paste(failure_log, collapse = " || ")
    )
  }

  .resolve_pillai_cloud_cached <- function(base, row_meta, within_ref,
                                           exclude_same_speaker = FALSE) {
    row_group <- if (!is.na(group_nm) && group_nm %in% names(row_meta)) {
      as.character(row_meta[[group_nm]][1])
    } else {
      NA_character_
    }

    specs_here <- .make_within_reference_specs(row_group, within_ref)
    cache_cols <- unique(c(
      if (!is.na(group_nm)) group_nm else NULL,
      if (!is.na(speaker_nm)) speaker_nm else NULL,
      if (!is.na(condition_nm)) condition_nm else NULL,
      condition_vars_nm
    ))

    key <- paste(
      "pillai",
      paste0("within_ref=", within_ref),
      paste0("exclude_same_speaker=", exclude_same_speaker),
      paste0("specs=", paste(names(specs_here), collapse = ">")),
      .row_cache_key(row_meta, cache_cols),
      sep = ";"
    )

    if (exists(key, envir = .resolved_cloud_cache, inherits = FALSE)) {
      return(get(key, envir = .resolved_cloud_cache, inherits = FALSE))
    }

    pool_base <- base
    if (!is.na(group_nm) && group_nm %in% names(base) && !is.na(row_group)) {
      pool_base <- pool_base %>% dplyr::filter(.data[[group_nm]] == row_group)
    }

    out <- .resolve_pillai_cloud_detail(
      pool = pool_base,
      row_meta = row_meta,
      spec_list = specs_here,
      contrast_col = "contrast_tmp",
      exclude_same_speaker = exclude_same_speaker
    )

    assign(key, out, envir = .resolved_cloud_cache)
    out
  }

  .resolve_pillai_between_cloud_detail <- function(pool, row_meta, spec_list, contrast_value) {
    failure_log <- character(0)

    for (spec_name in names(spec_list)) {
      spec_cols <- spec_list[[spec_name]]
      cand <- .subset_pool_by_spec(pool, row_meta, spec_cols)
      cand <- cand %>% dplyr::filter(.data[[contrast_nm]] == contrast_value)

      if (all(measure_nm %in% names(cand))) {
        cand_meas <- cand %>%
          dplyr::select(dplyr::all_of(measure_nm)) %>%
          dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x)))
      } else {
        cand_meas <- tibble::tibble()
      }

      if (nrow(cand_meas) >= 2) {
        return(list(
          data = cand_meas,
          spec_name = spec_name,
          spec_cols = spec_cols,
          effective_spec_name = .effective_spec_name(spec_cols, row_meta),
          failure_reason = NA_character_
        ))
      }

      failure_log <- c(failure_log, paste0(spec_name, ": insufficient_reference_points"))
    }

    list(
      data = NULL,
      spec_name = NA_character_,
      spec_cols = NULL,
      effective_spec_name = NA_character_,
      failure_reason = paste(failure_log, collapse = " || ")
    )
  }

  .upsert_cloud_col <- function(df, cloud_src, target_col = "cloud") {
    target_sym <- rlang::sym(target_col)
    
    if (target_col %in% names(df)) {
      df %>%
        dplyr::mutate(
          !!target_sym := dplyr::coalesce(
            !!target_sym,
            .clean_cloud_label({{ cloud_src }})
          )
        )
    } else {
      df %>%
        dplyr::mutate(
          !!target_sym := .clean_cloud_label({{ cloud_src }})
        )
    }
  }

  .resolve_between_pair_table <- function(target_index, pool_all, ref_g) {
    
    target_index %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        tmp = list({
          row_meta <- dplyr::pick(dplyr::everything())
          fg <- row_meta[[group_nm]][1]
          
          focal_pool <- pool_all %>% dplyr::filter(.data[[group_nm]] == fg)
          ref_pool   <- pool_all %>% dplyr::filter(.data[[group_nm]] == ref_g)
          
          # focal: mateix comportament que within-category
          focal_specs <- .make_within_category_specs(fg)
          
          # reference: sempre pooled per condition_vars, sense speaker
          ref_specs <- .make_between_specs(condition_vars_nm)
          
          focal_cloud <- .resolve_cloud(
            pool = focal_pool,
            row_meta = row_meta,
            spec_list = focal_specs,
            measure_nm = measure_nm,
            contrast_value = row_meta[[contrast_nm]][1],
            require_invertible = FALSE
          )
          
          ref_cloud <- .resolve_cloud(
            pool = ref_pool,
            row_meta = row_meta,
            spec_list = ref_specs,
            measure_nm = measure_nm,
            contrast_value = row_meta[[contrast_nm]][1],
            require_invertible = FALSE
          )
          
          tibble::tibble(
            focal_cloud_data = list(if (is.null(focal_cloud)) NULL else focal_cloud$data),
            focal_cloud_name = if (is.null(focal_cloud)) NA_character_ else focal_cloud$effective_spec_name,
            reference_cloud_data = list(if (is.null(ref_cloud)) NULL else ref_cloud$data),
            reference_cloud_name = if (is.null(ref_cloud)) NA_character_ else paste0(ref_g, ": ", ref_cloud$effective_spec_name)
          )
        })
      ) %>%
      tidyr::unnest_wider(tmp) %>%
      dplyr::ungroup()
  }
  
  .keep_existing <- function(cols, df) {
    cols[cols %in% names(df)]
  }
  
  .mah_debug <- FALSE
  
  # =========================
  # EUCLIDEAN
  # =========================

  if (isTRUE(compute_euclidean)) {

    .euc_vec = function(x_mat, mu_mat) {
      d = sqrt(rowSums((x_mat - mu_mat)^2))
      d[!is.finite(d)] = NA_real_
      d
    }

    if (!inter_group) {
      base_cols <- unique(c(
        if (!is.na(group_nm)) group_nm else NULL,
        if (!is.na(speaker_nm)) speaker_nm else NULL,
        if (!is.na(condition_nm)) condition_nm else NULL,
        condition_vars_nm,
        contrast_nm,
        token_id_nm,
        measure_nm
      ))
      base_cols_present <- base_cols[base_cols %in% names(xdata)]

      base_euc <- xdata %>%
        dplyr::select(dplyr::all_of(base_cols_present)) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct()

      target_keys <- unique(c(
        if (!is.na(group_nm)) group_nm else NULL,
        if (!is.na(speaker_nm)) speaker_nm else NULL,
        if (!is.na(condition_nm)) condition_nm else NULL,
        condition_vars_nm
      ))
      target_keys <- target_keys[target_keys %in% names(xdata)]

      target_index <- xdata %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(target_keys)))

      exclude_same_euc <- within_ref %in% c("population_pooled", "inter_speaker_mean")

      within_euc_refs <- target_index %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmp = list({
            row_meta <- dplyr::pick(dplyr::everything())

            ref_lvl1 <- .resolve_within_cloud_cached(
              base = base_euc,
              row_meta = row_meta,
              contrast_value = lvl2,
              within_ref = within_ref,
              require_invertible = FALSE,
              exclude_same_speaker = exclude_same_euc
            )

            ref_lvl2 <- .resolve_within_cloud_cached(
              base = base_euc,
              row_meta = row_meta,
              contrast_value = lvl1,
              within_ref = within_ref,
              require_invertible = FALSE,
              exclude_same_speaker = exclude_same_euc
            )

            mu_lvl1 <- if (!is.null(ref_lvl1$data) && nrow(ref_lvl1$data) > 0) {
              as.numeric(colMeans(ref_lvl1$data[, measure_nm, drop = FALSE], na.rm = TRUE))
            } else {
              rep(NA_real_, length(measure_nm))
            }
            mu_lvl2 <- if (!is.null(ref_lvl2$data) && nrow(ref_lvl2$data) > 0) {
              as.numeric(colMeans(ref_lvl2$data[, measure_nm, drop = FALSE], na.rm = TRUE))
            } else {
              rep(NA_real_, length(measure_nm))
            }

            names(mu_lvl1) <- paste0("mu_lvl1_", measure_nm)
            names(mu_lvl2) <- paste0("mu_lvl2_", measure_nm)
            
            out_row <- c(
              as.list(stats::setNames(as.numeric(mu_lvl1), paste0("mu_lvl1_", measure_nm))),
              as.list(stats::setNames(as.numeric(mu_lvl2), paste0("mu_lvl2_", measure_nm))),
              list(
                euc_ref_cloud_lvl1 = if (is.null(ref_lvl1$data)) NA_character_ else ref_lvl1$effective_spec_name,
                euc_ref_cloud_lvl2 = if (is.null(ref_lvl2$data)) NA_character_ else ref_lvl2$effective_spec_name,
                euc_fail_reason_lvl1 = if (is.null(ref_lvl1$data)) ref_lvl1$failure_reason else NA_character_,
                euc_fail_reason_lvl2 = if (is.null(ref_lvl2$data)) ref_lvl2$failure_reason else NA_character_
              )
            )
            
            tibble::as_tibble_row(out_row)
          })
        ) %>%
        tidyr::unnest_wider(tmp) %>%
        dplyr::ungroup()

      xdata <- xdata %>%
        dplyr::left_join(
          within_euc_refs %>% dplyr::select(-dplyr::any_of(c("euc_fail_reason_lvl1", "euc_fail_reason_lvl2"))),
          by = target_keys
        )

      mu_lvl1_cols <- paste0("mu_lvl1_", measure_nm)
      mu_lvl2_cols <- paste0("mu_lvl2_", measure_nm)

      euc_distances <- xdata %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          dist_wit_euc = {
            x_vec <- c_across(dplyr::all_of(measure_nm))
            mu_vec <- if (.data[[contrast_nm]] == lvl1) {
              c_across(dplyr::all_of(mu_lvl1_cols))
            } else if (.data[[contrast_nm]] == lvl2) {
              c_across(dplyr::all_of(mu_lvl2_cols))
            } else {
              rep(NA_real_, length(measure_nm))
            }
            .euc_vec(matrix(as.numeric(x_vec), nrow = 1), matrix(as.numeric(mu_vec), nrow = 1))[1]
          },
          euc_ref_cloud_wit = dplyr::case_when(
            .data[[contrast_nm]] == lvl1 ~ .data$euc_ref_cloud_lvl1,
            .data[[contrast_nm]] == lvl2 ~ .data$euc_ref_cloud_lvl2,
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-dplyr::all_of(c(mu_lvl1_cols, mu_lvl2_cols, "euc_ref_cloud_lvl1", "euc_ref_cloud_lvl2")))

      xdata <- .upsert_cloud_col(euc_distances, .data$euc_ref_cloud_wit, "cloud_wit")

      rm(base_euc, target_index, within_euc_refs, euc_distances)
    }

    if (inter_group) {
      .resolve_euc_between_center <- function(pool, row_meta, ref_g, contrast_value) {
        ref_specs <- .make_specs(
          condition_vars_nm,
          speaker_nm = NULL,
          include_speaker = FALSE
        )
        
        ref_pool <- pool %>% dplyr::filter(.data[[group_nm]] == ref_g)
        
        ref_cloud <- .resolve_center_only_cloud(
          pool = ref_pool,
          row_meta = row_meta,
          spec_list = ref_specs,
          measure_nm = measure_nm,
          contrast_value = contrast_value,
          exclude_same_speaker = FALSE
        )
        
        if (is.null(ref_cloud$data)) {
          return(list(
            mu = rep(NA_real_, length(measure_nm)),
            effective_spec_name = NA_character_,
            failure_reason = ref_cloud$failure_reason
          ))
        }
        
        list(
          mu = as.numeric(ref_cloud$stats$center),
          effective_spec_name = paste0(ref_g, ": ", ref_cloud$effective_spec_name),
          failure_reason = NA_character_
        )
      }

      .euc_between_one_direction <- function(ref_g, focal_gs) {
        purrr::map_dfr(focal_gs, function(fg) {
          focal_dat <- xdata %>%
            dplyr::filter(.data[[group_nm]] == fg) %>%
            dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
            dplyr::select(dplyr::all_of(unique(c(
              group_nm, speaker_nm, contrast_nm, token_id_nm, condition_nm, condition_vars_nm, measure_nm
            ))))

          if (nrow(focal_dat) == 0) return(tibble::tibble())

          focal_dat %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              tmp = list({
                row_meta <- dplyr::pick(dplyr::everything())

                ref_center <- .resolve_euc_between_center(
                  pool = xdata,
                  row_meta = row_meta,
                  ref_g = ref_g,
                  contrast_value = row_meta[[contrast_nm]][1]
                )

                x_vec <- as.numeric(row_meta[measure_nm])
                mu_vec <- ref_center$mu

                val <- .euc_vec(matrix(x_vec, nrow = 1), matrix(mu_vec, nrow = 1))[1]
                rr <- attr(val, "reason")
                if (is.null(rr)) rr <- NA_character_
                rr <- dplyr::coalesce(rr, ref_center$failure_reason)

                tibble::tibble(
                  dist_bet_euc = as.numeric(val),
                  euc_ref_cloud_bet = ref_center$effective_spec_name,
                  euc_fail_reason_bet = rr
                )
              })
            ) %>%
            tidyr::unnest_wider(tmp) %>%
            dplyr::ungroup() %>%
            dplyr::select(
              dplyr::all_of(unique(c(group_nm, speaker_nm, contrast_nm, token_id_nm, condition_nm))),
              dist_bet_euc,
              euc_ref_cloud_bet,
              euc_fail_reason_bet
            )
        })
      }

      if (isTRUE(auto_two_group)) {
        euc_segment_all <- dplyr::bind_rows(
          .euc_between_one_direction(ref_g = two_groups[2], focal_gs = two_groups[1]),
          .euc_between_one_direction(ref_g = two_groups[1], focal_gs = two_groups[2])
        )
      } else {
        euc_segment_all <- .euc_between_one_direction(ref_g = reference_group, focal_gs = focal_groups)
      }

      xdata <- xdata %>%
        dplyr::left_join(
          euc_segment_all %>% dplyr::select(-dplyr::any_of("euc_fail_reason_bet")),
          by = c(group_nm, speaker_nm, contrast_nm, token_id_nm, condition_nm) %>% unique()
        )

      xdata <- .upsert_cloud_col(xdata, .data$euc_ref_cloud_bet, "cloud_bet")

      rm(euc_segment_all, .euc_between_one_direction, .resolve_euc_between_center)
    }

    rm(.euc_vec)
  }

  # =========================
  # MAHALANOBIS
  # =========================
  
  if (isTRUE(compute_mahalanobis)) {
    
    if (!inter_group) {
      base_cols <- unique(c(
        if (!is.na(group_nm)) group_nm else NULL,
        speaker_nm, condition_nm, condition_vars_nm, contrast_nm, token_id_nm, measure_nm
      ))
      
      base_cols_present <- base_cols[base_cols %in% names(xdata)]
      
      base <- xdata %>%
        dplyr::select(dplyr::all_of(base_cols_present)) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct()
      
      target_index <- xdata %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(c(
          if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm, condition_nm, condition_vars_nm
        ))))
      
      exclude_same_wit <- within_ref %in% c("population_pooled", "inter_speaker_mean")
      
      within_cloud_refs <- target_index %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmp = list({
            row_meta <- dplyr::pick(dplyr::everything())
            ref_lvl1_inv <- .resolve_within_cloud_cached(
              base = base,
              row_meta = row_meta,
              contrast_value = lvl1,
              within_ref = within_ref,
              require_invertible = TRUE,
              exclude_same_speaker = exclude_same_wit
            )
            ref_lvl2_inv <- .resolve_within_cloud_cached(
              base = base,
              row_meta = row_meta,
              contrast_value = lvl2,
              within_ref = within_ref,
              require_invertible = TRUE,
              exclude_same_speaker = exclude_same_wit
            )
            ref_lvl1_bha <- if (isTRUE(compute_bhatt)) {
              .resolve_within_cloud_cached(
                base = base,
                row_meta = row_meta,
                contrast_value = lvl1,
                within_ref = within_ref,
                require_invertible = FALSE,
                exclude_same_speaker = FALSE
              )
            } else {
              NULL
            }
            ref_lvl2_bha <- if (isTRUE(compute_bhatt)) {
              .resolve_within_cloud_cached(
                base = base,
                row_meta = row_meta,
                contrast_value = lvl2,
                within_ref = within_ref,
                require_invertible = FALSE,
                exclude_same_speaker = FALSE
              )
            } else {
              NULL
            }
            tibble::tibble(
              mah_ref_lvl1 = list(ref_lvl1_inv),
              mah_ref_lvl2 = list(ref_lvl2_inv),
              bha_ref_lvl1 = list(ref_lvl1_bha),
              bha_ref_lvl2 = list(ref_lvl2_bha),
              mah_ref_cloud_lvl1 = if (is.null(ref_lvl1_inv$data)) NA_character_ else ref_lvl1_inv$effective_spec_name,
              mah_ref_cloud_lvl2 = if (is.null(ref_lvl2_inv$data)) NA_character_ else ref_lvl2_inv$effective_spec_name,
              mah_fail_reason_lvl1 = if (is.null(ref_lvl1_inv$data)) ref_lvl1_inv$failure_reason else NA_character_,
              mah_fail_reason_lvl2 = if (is.null(ref_lvl2_inv$data)) ref_lvl2_inv$failure_reason else NA_character_,
              bha_cloud_1 = if (is.null(ref_lvl1_bha) || is.null(ref_lvl1_bha$data)) NA_character_ else ref_lvl1_bha$effective_spec_name,
              bha_cloud_2 = if (is.null(ref_lvl2_bha) || is.null(ref_lvl2_bha$data)) NA_character_ else ref_lvl2_bha$effective_spec_name
            )
          })
        ) %>%
        tidyr::unnest_wider(tmp) %>%
        dplyr::ungroup()
      
      tokens_grp <- .keep_existing(
        c(if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm,
          condition_nm,
          condition_vars_nm,
          token_id_nm),
        base
      )
      
      token_cloud_cols <- .keep_existing(
        c(if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm,
          condition_nm,
          condition_vars_nm),
        within_cloud_refs
      )
      
      tokens_nested <- base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(tokens_grp, contrast_nm)))) %>%
        tidyr::nest() %>%
        dplyr::rename(token_data = data) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
          within_cloud_refs,
          by = token_cloud_cols
        )
      
      mah_distances <- tokens_nested %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmp = list({
            ref_cloud <- if (.data[[contrast_nm]] == lvl1) mah_ref_lvl2[[1]] else mah_ref_lvl1[[1]]
            ref_name  <- if (.data[[contrast_nm]] == lvl1) mah_ref_cloud_lvl2 else mah_ref_cloud_lvl1
            fail_reason <- if (.data[[contrast_nm]] == lvl1) mah_fail_reason_lvl2 else mah_fail_reason_lvl1
            
            mah_value <- if (is.null(ref_cloud$data)) {
              tmp_val <- NA_real_
              attr(tmp_val, "reason") <- fail_reason
              tmp_val
            } else {
              .mah_one(token_data, stats_y = ref_cloud$stats)
            }
            
            reason <- attr(mah_value, "reason")
            if (is.null(reason)) reason <- NA_character_
            
            tibble::tibble(
              dist_wit_mah = as.numeric(mah_value),
              mah_ref_cloud_wit = ref_name,
              mah_fail_reason_wit = reason
            )
          })
        ) %>%
        tidyr::unnest_wider(tmp) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(
          !!!rlang::syms(c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm, condition_nm, token_id_nm)),
          !!contrast_nm := .data[[contrast_nm]],
          dist_wit_mah,
          mah_ref_cloud_wit,
          mah_fail_reason_wit
        )
      
      # used_refs_wit <- mah_distances %>%
      #   dplyr::filter(!is.na(dist_wit_mah)) %>%
      #   dplyr::distinct(mah_ref_cloud_wit) %>%
      #   dplyr::pull(mah_ref_cloud_wit)
      # 
      # .emit_effective_clouds(
      #   "Mahalanobis (within):\n    ",
      #   used_refs_wit,
      #   sep = " || "
      # )
      
      fail_subjects_wit <- mah_distances %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm)))) %>%
        dplyr::summarise(
          n_tokens = dplyr::n(),
          n_success = sum(!is.na(dist_wit_mah)),
          failure_reason = if (all(is.na(mah_fail_reason_wit))) NA_character_ else dplyr::first(mah_fail_reason_wit[!is.na(mah_fail_reason_wit)]),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_success == 0)
      
      .emit_subject_failure_summary(fail_subjects_wit)
      xdata <- xdata %>%
        dplyr::left_join(
          mah_distances %>% dplyr::select(-dplyr::any_of("mah_fail_reason_wit")),
          by = c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm, condition_nm, token_id_nm, contrast_nm)
        )
      
      xdata <- .upsert_cloud_col(xdata, .data$mah_ref_cloud_wit, "cloud_wit")
      
      rm(base, target_index, within_cloud_refs, tokens_nested, mah_distances, fail_subjects_wit)
    }
    
    if (inter_group) {
      pool_all <- xdata %>%
        dplyr::select(!!!rlang::syms(unique(c(group_nm, speaker_nm, condition_nm, condition_vars_nm, contrast_nm, token_id_nm, measure_nm)))) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct()
      
      .mah_between_one_direction <- function(ref_g, focal_gs) {
        token_index_cols <- .keep_existing(
          c(group_nm, speaker_nm, condition_nm, condition_vars_nm, contrast_nm, token_id_nm),
          pool_all
        )
        
        token_index <- pool_all %>%
          dplyr::filter(.data[[group_nm]] %in% focal_gs) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(token_index_cols))) %>%
          tidyr::nest() %>%
          rename(token_data = data) %>%
          dplyr::ungroup()
        
        ref_specs <- .make_between_specs(condition_vars_nm)
        ref_pool_all <- pool_all %>% dplyr::filter(.data[[group_nm]] == ref_g)
        
        token_index %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            tmp = list({
              row_meta <- dplyr::pick(dplyr::everything())
              
              ref_cloud <- .resolve_mah_cloud_detail(
                pool = ref_pool_all,
                row_meta = row_meta,
                spec_list = ref_specs,
                measure_nm = measure_nm,
                contrast_value = row_meta[[contrast_nm]][1],
                exclude_same_speaker = FALSE
              )
              
              mah_value <- if (is.null(ref_cloud$data)) {
                tmp_val <- NA_real_
                attr(tmp_val, "reason") <- ref_cloud$failure_reason
                tmp_val
              } else {
                .mah_one(token_data, stats_y = ref_cloud$stats)
              }
              
              reason <- attr(mah_value, "reason")
              if (is.null(reason)) reason <- NA_character_
              
              tibble::tibble(
                dist_bet_mah = as.numeric(mah_value),
                mah_ref_cloud_bet = if (is.null(ref_cloud$data)) NA_character_ else paste0(ref_g, ": ", ref_cloud$effective_spec_name),
                mah_fail_reason_bet = reason
              )
            })
          ) %>%
          tidyr::unnest_wider(tmp) %>%
          dplyr::ungroup() %>%
          dplyr::transmute(
            !!group_nm := .data[[group_nm]],
            !!speaker_nm := .data[[speaker_nm]],
            !!condition_nm := .data[[condition_nm]],
            !!contrast_nm := .data[[contrast_nm]],
            !!token_id_nm := .data[[token_id_nm]],
            dist_bet_mah,
            mah_ref_cloud_bet,
            mah_fail_reason_bet
          )
      }
      
      if (isTRUE(auto_two_group)) {
        mah_between_all <- dplyr::bind_rows(
          .mah_between_one_direction(ref_g = two_groups[2], focal_gs = two_groups[1]),
          .mah_between_one_direction(ref_g = two_groups[1], focal_gs = two_groups[2])
        )
      } else {
        mah_between_all <- .mah_between_one_direction(ref_g = reference_group, focal_gs = focal_groups)
      }
      
      # used_refs_bet <- mah_between_all %>%
      #   dplyr::filter(!is.na(dist_bet_mah)) %>%
      #   dplyr::distinct(mah_ref_cloud_bet) %>%
      #   dplyr::pull(mah_ref_cloud_bet)
      # 
      # .emit_effective_clouds(
      #   "Mahalanobis (between):\n    ",
      #   used_refs_bet,
      #   sep = " <-> "
      # )
      
      xdata <- xdata %>%
        dplyr::left_join(
          mah_between_all %>% dplyr::select(-dplyr::any_of("mah_fail_reason_bet")),
          by = c(group_nm, speaker_nm, condition_nm, contrast_nm, token_id_nm)
        )
      
      xdata <- .upsert_cloud_col(xdata, .data$mah_ref_cloud_bet, "cloud_bet")
      
      rm(pool_all, mah_between_all, .mah_between_one_direction)
    }
  }
  
  # =========================
  # PILLAI
  # =========================

  if (isTRUE(compute_pillai)) {

    if (!inter_group) {
      base_cols <- unique(c(
        if (!is.na(group_nm)) group_nm else NULL,
        if (!is.na(speaker_nm)) speaker_nm else NULL,
        if (!is.na(condition_nm)) condition_nm else NULL,
        condition_vars_nm,
        contrast_nm,
        token_id_nm,
        measure_nm
      ))
      base_cols_present <- base_cols[base_cols %in% names(xdata)]

      base_pillai <- xdata %>%
        dplyr::select(dplyr::all_of(base_cols_present)) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct() %>%
        dplyr::rename(contrast_tmp = !!contrast_q)

      target_keys <- unique(c(
        if (!is.na(group_nm)) group_nm else NULL,
        if (!is.na(speaker_nm)) speaker_nm else NULL,
        if (!is.na(condition_nm)) condition_nm else NULL,
        condition_vars_nm
      ))
      target_keys <- target_keys[target_keys %in% names(xdata)]

      target_index <- xdata %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(target_keys)))

      exclude_same_pillai <- within_ref %in% c("population_pooled", "inter_speaker_mean")

      pillai_cloud_refs <- target_index %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmp = list({
            row_meta <- dplyr::pick(dplyr::everything())

            ref_cloud <- .resolve_pillai_cloud_cached(
              base = base_pillai,
              row_meta = row_meta,
              within_ref = within_ref,
              exclude_same_speaker = exclude_same_pillai
            )

            tibble::tibble(
              pil_cloud_data = list(ref_cloud$data),
              pil_ref_cloud_wit = if (is.null(ref_cloud$data)) NA_character_ else ref_cloud$effective_spec_name,
              pil_fail_reason_wit = if (is.null(ref_cloud$data)) ref_cloud$failure_reason else NA_character_
            )
          })
        ) %>%
        tidyr::unnest_wider(tmp) %>%
        dplyr::ungroup()

      pillai <- pillai_cloud_refs %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          dist_wit_pil = {
            if (is.null(pil_cloud_data[[1]])) {
              tmp_val <- NA_real_
              attr(tmp_val, "reason") <- pil_fail_reason_wit
              tmp_val
            } else {
              .pillai_two_clouds(
                pil_cloud_data[[1]],
                "contrast_tmp",
                debug = diagnostics,
                debug_id = paste(
                  if (!is.na(speaker_nm) && speaker_nm %in% names(dplyr::pick(dplyr::everything()))) .data[[speaker_nm]] else "ALL",
                  if (!is.na(condition_nm) && condition_nm %in% names(dplyr::pick(dplyr::everything()))) .data[[condition_nm]] else "GLOBAL",
                  sep = " || "
                )
              )
            }
          },
          pil_num_reason = {
            rr <- attr(dist_wit_pil, "reason")
            if (is.null(rr)) NA_character_ else rr
          }
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          pil_fail_reason_wit = dplyr::coalesce(.data$pil_fail_reason_wit, .data$pil_num_reason),
          dist_wit_pil = as.numeric(.data$dist_wit_pil)
        ) %>%
        dplyr::select(
          dplyr::all_of(target_keys),
          dist_wit_pil,
          pil_ref_cloud_wit,
          pil_fail_reason_wit
        )

      xdata <- xdata %>%
        dplyr::left_join(
          pillai %>% dplyr::select(-dplyr::any_of("pil_fail_reason_wit")),
          by = target_keys
        )
      xdata <- .upsert_cloud_col(xdata, .data$pil_ref_cloud_wit, "cloud_wit")

      if (isTRUE(diagnostics)) {
        total_n  <- nrow(target_index)
        success_n <- sum(!is.na(pillai$dist_wit_pil))
        message("\n--- PILLAI DIAGNOSTIC SUMMARY (WITHIN) ---")
        message("Total cases attempted: ", total_n)
        message("Success (Valid numbers): ", success_n, " (", round(100 * success_n / max(total_n, 1), 1), "%)")
        fail_reasons <- pillai$pil_fail_reason_wit
        fail_reasons <- fail_reasons[!is.na(fail_reasons)]
        if (length(fail_reasons) > 0) {
          tb <- sort(table(fail_reasons), decreasing = TRUE)
          message("\nBreakdown of failures:\n")
          print(tb)
        } else if (success_n == 0) {
          message("![WARNING] All calculations resulted in NaN or were skipped.")
        }
        message("----------------------------------------------")
      }

      rm(base_pillai, target_index, pillai_cloud_refs, pillai)
    }

    if (inter_group) {
      .pillai_group_tmp <- function(dat) {
        .pillai_two_clouds(
          df_ab     = dat,
          class_col = "group_tmp",
          debug     = diagnostics,
          debug_id  = "Inter-Group"
        )
      }

      .pillai_between_one_direction <- function(ref_g, focal_gs) {
        ref_pool_all <- xdata %>%
          dplyr::filter(.data[[group_nm]] == ref_g) %>%
          dplyr::select(
            dplyr::all_of(unique(c(
              group_nm, speaker_nm, condition_nm, condition_vars_nm, contrast_nm, measure_nm
            )))
          ) %>%
          dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x)))

        purrr::map_dfr(focal_gs, function(fg) {
          focal_dat <- xdata %>%
            dplyr::filter(.data[[group_nm]] == fg) %>%
            dplyr::select(
              dplyr::all_of(unique(c(
                group_nm, speaker_nm, condition_nm, condition_vars_nm, contrast_nm, token_id_nm, measure_nm
              )))
            ) %>%
            dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x)))

          if (nrow(focal_dat) == 0) return(tibble::tibble())

          focal_dat %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              tmp = list({
                row_meta <- dplyr::pick(dplyr::everything())

                ref_specs <- .make_specs(
                  include_speaker = FALSE,
                  include_condition = TRUE,
                  condition_vars = condition_vars_nm,
                  speaker_nm = speaker_nm,
                  condition_nm = condition_nm
                )

                ref_cloud <- .resolve_pillai_between_cloud_detail(
                  pool = ref_pool_all,
                  row_meta = row_meta,
                  spec_list = ref_specs,
                  contrast_value = row_meta[[contrast_nm]][1]
                )

                token_data <- tibble::as_tibble(row_meta) %>%
                  dplyr::select(dplyr::all_of(measure_nm))

                pil_value <- if (is.null(ref_cloud$data)) {
                  tmp_val <- NA_real_
                  attr(tmp_val, "reason") <- ref_cloud$failure_reason
                  tmp_val
                } else {
                  dat_ab <- dplyr::bind_rows(
                    dplyr::mutate(token_data, group_tmp = "focal"),
                    dplyr::mutate(ref_cloud$data, group_tmp = "reference")
                  )
                  .pillai_group_tmp(dat_ab)
                }

                rr <- attr(pil_value, "reason")
                if (is.null(rr)) rr <- NA_character_
                rr <- dplyr::coalesce(rr, ref_cloud$failure_reason)

                tibble::tibble(
                  dist_bet_pil = as.numeric(pil_value),
                  pil_ref_cloud_bet = if (is.null(ref_cloud$data)) NA_character_ else paste0(ref_g, ": ", ref_cloud$effective_spec_name),
                  pil_fail_reason_bet = rr
                )
              })
            ) %>%
            tidyr::unnest_wider(tmp) %>%
            dplyr::ungroup() %>%
            dplyr::transmute(
              !!group_nm := .data[[group_nm]],
              !!speaker_nm := .data[[speaker_nm]],
              !!condition_nm := .data[[condition_nm]],
              !!contrast_nm := .data[[contrast_nm]],
              !!token_id_nm := .data[[token_id_nm]],
              dist_bet_pil,
              pil_ref_cloud_bet,
              pil_fail_reason_bet
            )
        })
      }

      if (isTRUE(auto_two_group)) {
        pillai_all <- dplyr::bind_rows(
          .pillai_between_one_direction(ref_g = two_groups[2], focal_gs = two_groups[1]),
          .pillai_between_one_direction(ref_g = two_groups[1], focal_gs = two_groups[2])
        )
      } else {
        pillai_all <- .pillai_between_one_direction(ref_g = reference_group, focal_gs = focal_groups)
      }

      xdata <- xdata %>%
        dplyr::left_join(
          pillai_all %>% dplyr::select(-dplyr::any_of("pil_fail_reason_bet")),
          by = c(group_nm, speaker_nm, condition_nm, contrast_nm, token_id_nm)
        )
      xdata <- .upsert_cloud_col(xdata, .data$pil_ref_cloud_bet, "cloud_bet")
      xdata <- xdata %>%
        {
          if (isTRUE(auto_two_group)) . else dplyr::mutate(., dist_bet_pil = dplyr::if_else(!!group_q == reference_group, NA_real_, dist_bet_pil))
        }

      if (isTRUE(diagnostics)) {
        total_n <- nrow(pillai_all)
        success_n <- sum(!is.na(pillai_all$dist_bet_pil))
        message("\n--- PILLAI DIAGNOSTIC SUMMARY (BETWEEN) ---")
        message("Total cases attempted: ", total_n)
        message("Success (Valid numbers): ", success_n, " (", round(100 * success_n / max(total_n, 1), 1), "%)")
        fail_reasons <- pillai_all$pil_fail_reason_bet
        fail_reasons <- fail_reasons[!is.na(fail_reasons)]
        if (length(fail_reasons) > 0) {
          tb <- sort(table(fail_reasons), decreasing = TRUE)
          message("\nBreakdown of failures:\n")
          print(tb)
        } else if (success_n == 0) {
          message("![WARNING] All calculations resulted in NaN or were skipped.")
        }
        message("----------------------------------------------")
      }

      rm(pillai_all, .pillai_between_one_direction, .pillai_group_tmp)
    }
  }

  # --- Lógica de Limpieza y Diagnóstico Post-Cálculo ---

  if (isTRUE(compute_pillai)) {
    target_col <- if (inter_group) "dist_bet_pil" else "dist_wit_pil"

    if (target_col %in% names(xdata)) {
      all_failed <- all(is.na(xdata[[target_col]]) | is.nan(xdata[[target_col]]))
      if (isTRUE(all_failed) && nrow(xdata) > 0) {
        xdata <- xdata %>% dplyr::select(-dplyr::all_of(target_col))
      }
    }
  }

  # =========================
  # BHATTACHARYYA
  # =========================
  
  if (isTRUE(compute_bhatt)) {
    
    if (!inter_group) {
      base <- xdata %>%
        dplyr::select(!!!rlang::syms(unique(c(
          if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm, condition_nm, condition_vars_nm, contrast_nm, measure_nm
        )))) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct()
      
      target_index <- xdata %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(c(
          if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm, condition_nm, condition_vars_nm
        ))))
      
      bha_ref_cols <- .keep_existing(
        c(if (!is.na(group_nm)) group_nm else NULL,
          speaker_nm, condition_nm, condition_vars_nm,
          "bha_ref_lvl1", "bha_ref_lvl2", "bha_cloud_1", "bha_cloud_2"),
        xdata
      )
      
      if (all(c("bha_ref_lvl1", "bha_ref_lvl2", "bha_cloud_1", "bha_cloud_2") %in% names(xdata))) {
        bh_within <- xdata %>%
          dplyr::distinct(dplyr::across(dplyr::all_of(bha_ref_cols)))
      } else {
        bh_within <- target_index %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            tmp = list({
              row_meta <- dplyr::pick(dplyr::everything())
              cloud1 <- .resolve_within_cloud_cached(
                base = base,
                row_meta = row_meta,
                contrast_value = lvl1,
                within_ref = within_ref,
                require_invertible = FALSE,
                exclude_same_speaker = FALSE
              )
              cloud2 <- .resolve_within_cloud_cached(
                base = base,
                row_meta = row_meta,
                contrast_value = lvl2,
                within_ref = within_ref,
                require_invertible = FALSE,
                exclude_same_speaker = FALSE
              )
              tibble::tibble(
                bha_ref_lvl1 = list(cloud1),
                bha_ref_lvl2 = list(cloud2),
                bha_cloud_1 = if (is.null(cloud1$data)) NA_character_ else cloud1$effective_spec_name,
                bha_cloud_2 = if (is.null(cloud2$data)) NA_character_ else cloud2$effective_spec_name
              )
            })
          ) %>%
          tidyr::unnest_wider(tmp) %>%
          dplyr::ungroup()
      }
      
      bh_within <- bh_within %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          tmp = list({
            cloud1 <- bha_ref_lvl1[[1]]
            cloud2 <- bha_ref_lvl2[[1]]
            if (is.null(cloud1$data) || is.null(cloud2$data)) {
              tibble::tibble(
                dist_wit_bha = NA_real_,
                cloud_1 = NA_character_,
                cloud_2 = NA_character_,
                cloud = NA_character_
              )
            } else {
              tibble::tibble(
                dist_wit_bha = .bhatt(cloud1$data, cloud2$data, statsA = cloud1$stats, statsB = cloud2$stats),
                cloud_1 = bha_cloud_1,
                cloud_2 = bha_cloud_2,
                cloud = ifelse(
                  identical(bha_cloud_1, bha_cloud_2),
                  bha_cloud_1,
                  paste0(
                    "[", bha_cloud_1, "] (", lvl1, ") vs ",
                    "[", bha_cloud_2, "] (", lvl2, ")"
                  )
                )
              )
            }
          })
        ) %>%
        tidyr::unnest_wider(tmp) %>%
        dplyr::ungroup()
      
      used_pairs <- bh_within %>%
        dplyr::filter(!is.na(dist_wit_bha)) %>%
        dplyr::distinct(cloud) %>%
        dplyr::pull(cloud)
      
      .emit_clouds_used(
        "clouds used: ",
        used_pairs
      )
      
      bh_wit_join_cols <- c(
        if (!is.na(group_nm)) group_nm else NULL,
        speaker_nm, condition_nm
      )
      
      xdata <- xdata %>%
        dplyr::left_join(
          bh_within %>% dplyr::select(dplyr::all_of(bh_wit_join_cols), dist_wit_bha, bha_ref_cloud_wit = cloud),
          by = c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm, condition_nm)
        ) %>%
        dplyr::select(-dplyr::any_of("cloud_bha"))
      
      xdata <- .upsert_cloud_col(xdata, .data$bha_ref_cloud_wit, "cloud_wit")
      
      rm(base, target_index, bh_within, used_pairs, bha_ref_cols)
    }
    
    if (inter_group) {
      pool_all <- xdata %>%
        dplyr::select(!!!rlang::syms(unique(c(
          group_nm, speaker_nm, condition_nm, condition_vars_nm, contrast_nm, measure_nm
        )))) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(measure_nm), ~ !is.na(.x))) %>%
        dplyr::distinct()
      
      if (isTRUE(auto_two_group)) {
        bh_between_all <- dplyr::bind_rows(
          .resolve_between_pair_table(
            target_index = xdata %>%
              dplyr::filter(.data[[group_nm]] == two_groups[1]) %>%
              dplyr::distinct(
                !!group_q, !!speaker_q, !!condition_q,
                !!!rlang::syms(condition_vars_nm), !!contrast_q
              ),
            pool_all = pool_all,
            ref_g = two_groups[2]
          ) %>%
            dplyr::mutate(.ref_group = two_groups[2]),
          
          .resolve_between_pair_table(
            target_index = xdata %>%
              dplyr::filter(.data[[group_nm]] == two_groups[2]) %>%
              dplyr::distinct(
                !!group_q, !!speaker_q, !!condition_q,
                !!!rlang::syms(condition_vars_nm), !!contrast_q
              ),
            pool_all = pool_all,
            ref_g = two_groups[1]
          ) %>%
            dplyr::mutate(.ref_group = two_groups[1])
        ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            dist_bet_bha = if (is.null(focal_cloud_data[[1]]) || is.null(reference_cloud_data[[1]])) {
              NA_real_
            } else {
              .bhatt(focal_cloud_data[[1]], reference_cloud_data[[1]])
            },
            cloud_focal = focal_cloud_name,
            cloud_ref = reference_cloud_name
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(
            -.ref_group,
            dplyr::everything()
          ) %>%
          dplyr::mutate(ref_group = .data[[group_nm]]) %>%
          dplyr::mutate(ref_group = dplyr::if_else(
            !is.na(cloud_ref),
            sub(":.*$", "", cloud_ref),
            ref_group
          )) %>%
          dplyr::select(
            -focal_cloud_data, -reference_cloud_data,
            -focal_cloud_name, -reference_cloud_name
          )
      } else {
        target_index <- xdata %>%
          dplyr::filter(.data[[group_nm]] %in% focal_groups) %>%
          dplyr::distinct(
            !!group_q, !!speaker_q, !!condition_q,
            !!!rlang::syms(condition_vars_nm), !!contrast_q
          )
        
        bh_between_all <- .resolve_between_pair_table(
          target_index = target_index,
          pool_all = pool_all,
          ref_g = reference_group
        ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            dist_bet_bha = if (is.null(focal_cloud_data[[1]]) || is.null(reference_cloud_data[[1]])) {
              NA_real_
            } else {
              .bhatt(focal_cloud_data[[1]], reference_cloud_data[[1]])
            },
            cloud_focal = focal_cloud_name,
            cloud_ref = reference_cloud_name,
            ref_group = reference_group
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(
            -focal_cloud_data, -reference_cloud_data,
            -focal_cloud_name, -reference_cloud_name
          )
      }
      
      used_pairs_bet <- bh_between_all %>%
        dplyr::filter(!is.na(dist_bet_bha)) %>%
        dplyr::distinct(cloud_focal, cloud_ref) %>%
        dplyr::mutate(
          cloud_ref_clean = sub("^[^:]+: ", "", cloud_ref),
          lbl = ifelse(
            cloud_focal == cloud_ref_clean,
            cloud_focal,
            paste0(cloud_focal, " <-> ", cloud_ref_clean)
          )
        ) %>%
        dplyr::pull(lbl)
      
      .emit_clouds_used(
        "clouds used: ",
        used_pairs_bet
      )
      
      ref_cloud_lookup <- bh_between_all %>%
        dplyr::filter(!is.na(cloud_ref)) %>%
        dplyr::distinct(ref_group, !!rlang::sym(contrast_nm), cloud_ref) %>%
        dplyr::rename(
          .ref_group_lookup = ref_group,
          .contrast_lookup = !!rlang::sym(contrast_nm),
          .cloud_ref_lookup = cloud_ref
        )
      
      xdata <- xdata %>%
        dplyr::left_join(
          bh_between_all %>% 
            dplyr::select(
              !!group_q, !!speaker_q, !!condition_q, !!contrast_q,
              dist_bet_bha, cloud_focal
            ),
          by = c(group_nm, speaker_nm, condition_nm, contrast_nm)
        ) %>%
        dplyr::left_join(
          ref_cloud_lookup,
          by = stats::setNames(
            c(".ref_group_lookup", ".contrast_lookup"),
            c(group_nm, contrast_nm)
          )
        ) %>%
        dplyr::mutate(
          bha_ref_cloud_bet = dplyr::if_else(
            .data[[group_nm]] == reference_group,
            sub("^[^:]+: ", "", .cloud_ref_lookup),
            cloud_focal
          )
        ) %>%
        dplyr::mutate(
          dist_bet_bha = dplyr::if_else(
            .data[[group_nm]] == reference_group,
            NA_real_,
            dist_bet_bha
          )
        ) %>%
        dplyr::select(-cloud_focal, -.cloud_ref_lookup)
      
      xdata <- .upsert_cloud_col(xdata, .data$bha_ref_cloud_bet, "cloud_bet")
      
      rm(pool_all, bh_between_all, used_pairs_bet, ref_cloud_lookup)
      if (exists("target_index")) rm(target_index)
    }
  }
  
  # =========================
  # Final
  # =========================
  
  if (isTRUE(compute_mahalanobis)) {
    if (inter_group) {
      xdata <- xdata %>% mutate(dist_bet_mah = ifelse(is.nan(dist_bet_mah), NA_real_, dist_bet_mah))
    } else {
      xdata <- xdata %>% mutate(dist_wit_mah = ifelse(is.nan(dist_wit_mah), NA_real_, dist_wit_mah))
    }
  }
  
  if ("dist_wit_bha" %in% names(xdata)) {
    xdata <- xdata %>% mutate(dist_wit_bha = ifelse(is.nan(dist_wit_bha), NA_real_, dist_wit_bha))
  }
  if ("dist_bet_bha" %in% names(xdata)) {
    xdata <- xdata %>% mutate(dist_bet_bha = ifelse(is.nan(dist_bet_bha), NA_real_, dist_bet_bha))
  }
  
  xdata <- xdata %>%
    arrange(as.integer(.token_rowid)) %>%
    select(-matches("^n_groups(\\.x|\\.y)?$")) %>%
    select(-any_of(c(".token_rowid", ".condition_id", if (!speaker_user_supplied) ".speaker_dummy_internal" else NULL)))
  
  if (!isTRUE(diagnostics)) {
    xdata <- xdata %>%
      dplyr::select(
        -dplyr::any_of(c(
          "euc_ref_cloud_wit", "mah_ref_cloud_wit", "pil_ref_cloud_wit", "bha_ref_cloud_wit",
          "euc_ref_cloud_bet", "mah_ref_cloud_bet", "pil_ref_cloud_bet", "bha_ref_cloud_bet",
          "euc_fail_reason_lvl1", "euc_fail_reason_lvl2",
          "mah_fail_reason_wit", "mah_fail_reason_bet",
          "pil_fail_reason_wit", "pil_fail_reason_bet"
        ))
      )
  }
  
  xdata
}
