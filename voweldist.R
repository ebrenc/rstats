# voweldist(): compute within- and between-group vowel distance metrics
#
# Computes token-level distances (Euclidean, Mahalanobis) and
# category-level separation measures (Pillai trace, Bhattacharyya)
# between levels of a binary contrast (contrast_var).
#
# Stratification:
# `condition_vars` (if provided) define analysis strata. When multiple
# variables are supplied, they are internally combined into `.condition_id`.
# If NULL, a single global stratum is used.
#
# Speaker handling:
# `speaker_var` is optional. If NULL, computations are performed on pooled
# data (no speaker-level stratification).
#
# Within-group (inter_group = FALSE):
# Distances are computed for each token relative to a reference distribution
# defined by `within_ref` (e.g., inter-speaker mean, pooled, etc.), within
# each stratum.
#
# Between-group (inter_group = TRUE):
# Separation metrics are computed between each focal group and a
# `reference_group`, within each stratum, and joined back to token-level data.
#
# Cloud resolution:
# If the requested stratification is not feasible (e.g., insufficient data or
# non-invertible covariance), the function automatically falls back to simpler
# clouds (e.g., removing speaker or condition variables, or global).
# The effective clouds used are reported via messages and stored in `cloud`.
#
# Notes:
# - `contrast_var` must have exactly two levels.
# - Mahalanobis and Bhattacharyya require sufficient data to estimate
#   covariance matrices; otherwise fallback or NA may occur.
# 
# # Example 1: within-group distances with stratification
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
# # Example 2: between-group comparison with reference group
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
# # Example 3: global analysis (no conditions, no speaker)
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
    
    if (!inter_group) {
      .euc_vec = function(x_mat, mu_mat) {
        d = sqrt(rowSums((x_mat - mu_mat)^2))
        d[d == 0] = NA_real_
        d
      }
      
      d1 = xdata %>%
        select(!!speaker_q, !!contrast_q, all_of(token_id_nm), !!condition_q, all_of(measure_nm)) %>%
        filter(!!contrast_q == lvl1) %>%
        select(!all_of(contrast_nm))
      
      d2 = xdata %>%
        select(!!speaker_q, !!contrast_q, !!condition_q, all_of(measure_nm)) %>%
        filter(!!contrast_q == lvl2) %>%
        select(!all_of(contrast_nm))
      
      if (nrow(d2) == 0) stop("Within-group Euclidean: opposite contrast level has 0 rows.")
      
      if (within_ref == "self_speaker") {
        ref_lvl1 = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        joined = left_join(d1, ref_lvl1, by = c(speaker_nm, condition_nm), suffix = c(".x", ".y"))
        x_mat  = as.matrix(joined[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(joined[, paste0(measure_nm, ".y"), drop = FALSE])
        
        euc_level_1 = tibble(distance = .euc_vec(x_mat, mu_mat)) %>%
          bind_cols(d1, .) %>%
          select(!all_of(measure_nm)) %>%
          mutate(!!contrast_nm := lvl1)
        
      } else if (within_ref == "inter_speaker_mean") {
        ref_by_speaker = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
          rename(ref_speaker = !!speaker_q)
        
        tokens = d1 %>% rename(target_speaker = !!speaker_q)
        
        pairs = left_join(tokens, ref_by_speaker, by = condition_nm,
                          relationship = "many-to-many", suffix = c(".x", ".y")) %>%
          filter(target_speaker != ref_speaker)
        
        x_mat  = as.matrix(pairs[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(pairs[, paste0(measure_nm, ".y"), drop = FALSE])
        pairs$distance = .euc_vec(x_mat, mu_mat)
        
        euc_level_1 = pairs %>%
          group_by(target_speaker, across(all_of(token_id_nm)), !!condition_q) %>%
          summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
          rename(!!speaker_nm := target_speaker) %>%
          mutate(!!contrast_nm := lvl1)
        
      } else {
        ref_means = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
          rename(ref_speaker = !!speaker_q)
        
        speakers_by_cond = xdata %>%
          distinct(!!speaker_q, !!condition_q) %>%
          rename(target_speaker = !!speaker_q)
        
        pooled_ref = left_join(speakers_by_cond, ref_means, by = condition_nm, relationship = "many-to-many") %>%
          filter(target_speaker != ref_speaker) %>%
          group_by(target_speaker, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        joined = left_join(d1 %>% rename(target_speaker = !!speaker_q),
                           pooled_ref,
                           by = c("target_speaker", condition_nm),
                           suffix = c(".x", ".y"))
        
        x_mat  = as.matrix(joined[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(joined[, paste0(measure_nm, ".y"), drop = FALSE])
        
        euc_level_1 = tibble(distance = .euc_vec(x_mat, mu_mat)) %>%
          bind_cols(d1, .) %>%
          select(!all_of(measure_nm)) %>%
          mutate(!!contrast_nm := lvl1)
      }
      
      d1 = xdata %>%
        select(!!speaker_q, !!contrast_q, all_of(token_id_nm), !!condition_q, all_of(measure_nm)) %>%
        filter(!!contrast_q == lvl2) %>%
        select(!all_of(contrast_nm))
      
      d2 = xdata %>%
        select(!!speaker_q, !!contrast_q, !!condition_q, all_of(measure_nm)) %>%
        filter(!!contrast_q == lvl1) %>%
        select(!all_of(contrast_nm))
      
      if (nrow(d2) == 0) stop("Within-group Euclidean: opposite contrast level has 0 rows.")
      
      if (within_ref == "self_speaker") {
        ref_lvl2 = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        joined = left_join(d1, ref_lvl2, by = c(speaker_nm, condition_nm), suffix = c(".x", ".y"))
        x_mat  = as.matrix(joined[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(joined[, paste0(measure_nm, ".y"), drop = FALSE])
        
        euc_level_2 = tibble(distance = .euc_vec(x_mat, mu_mat)) %>%
          bind_cols(d1, .) %>%
          select(!all_of(measure_nm)) %>%
          mutate(!!contrast_nm := lvl2)
        
      } else if (within_ref == "inter_speaker_mean") {
        ref_by_speaker = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
          rename(ref_speaker = !!speaker_q)
        
        tokens = d1 %>% rename(target_speaker = !!speaker_q)
        
        pairs = left_join(tokens, ref_by_speaker, by = condition_nm,
                          relationship = "many-to-many", suffix = c(".x", ".y")) %>%
          filter(target_speaker != ref_speaker)
        
        x_mat  = as.matrix(pairs[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(pairs[, paste0(measure_nm, ".y"), drop = FALSE])
        pairs$distance = .euc_vec(x_mat, mu_mat)
        
        euc_level_2 = pairs %>%
          group_by(target_speaker, across(all_of(token_id_nm)), !!condition_q) %>%
          summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
          rename(!!speaker_nm := target_speaker) %>%
          mutate(!!contrast_nm := lvl2)
        
      } else {
        ref_means = d2 %>%
          group_by(!!speaker_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
          rename(ref_speaker = !!speaker_q)
        
        speakers_by_cond = xdata %>%
          distinct(!!speaker_q, !!condition_q) %>%
          rename(target_speaker = !!speaker_q)
        
        pooled_ref = left_join(speakers_by_cond, ref_means, by = condition_nm,
                               relationship = "many-to-many") %>%
          filter(target_speaker != ref_speaker) %>%
          group_by(target_speaker, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        joined = left_join(d1 %>% rename(target_speaker = !!speaker_q),
                           pooled_ref,
                           by = c("target_speaker", condition_nm),
                           suffix = c(".x", ".y"))
        
        x_mat  = as.matrix(joined[, paste0(measure_nm, ".x"), drop = FALSE])
        mu_mat = as.matrix(joined[, paste0(measure_nm, ".y"), drop = FALSE])
        
        euc_level_2 = tibble(distance = .euc_vec(x_mat, mu_mat)) %>%
          bind_cols(d1, .) %>%
          select(!all_of(measure_nm)) %>%
          mutate(!!contrast_nm := lvl2)
      }
      
      euc_distances = bind_rows(euc_level_1, euc_level_2) %>%
        rename(dist_wit_euc = distance)
      
      xdata = xdata %>%
        left_join(euc_distances, by = c(speaker_nm, contrast_nm, token_id_nm, condition_nm) %>% unique())
      
      rm(euc_level_1, euc_level_2, euc_distances, d1, d2, .euc_vec)
    }
    
    if (inter_group) {
      .euc_between_one_direction = function(ref_g, focal_gs) {
        d2_stratum = xdata %>%
          select(!!group_q, !!contrast_q, !!condition_q, all_of(measure_nm)) %>%
          filter(!!group_q == ref_g) %>%
          group_by(!!contrast_q, !!condition_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        d2_pooled = xdata %>%
          select(!!group_q, !!contrast_q, all_of(measure_nm)) %>%
          filter(!!group_q == ref_g) %>%
          group_by(!!contrast_q) %>%
          summarise(across(all_of(measure_nm), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
        
        purrr::map_dfr(focal_gs, function(fg) {
          d1 = xdata %>%
            select(!!group_q, !!speaker_q, !!contrast_q, all_of(token_id_nm), !!condition_q, all_of(measure_nm)) %>%
            filter(!!group_q == fg)
          
          if (nrow(d1) == 0) return(tibble())
          
          euc_segment = d1 %>%
            left_join(d2_stratum, by = c(contrast_nm, condition_nm), suffix = c(".x", ".y"))
          
          euc_segment_fb = d1 %>%
            left_join(d2_pooled, by = contrast_nm, suffix = c(".x", ".fb"))
          
          for (m in measure_nm) {
            y_col  = paste0(m, ".y")
            fb_col = paste0(m, ".fb")
            if (y_col %in% names(euc_segment) && fb_col %in% names(euc_segment_fb)) {
              euc_segment[[y_col]] = dplyr::coalesce(euc_segment[[y_col]], euc_segment_fb[[fb_col]])
            }
          }
          
          x_mat  = as.matrix(euc_segment[, paste0(measure_nm, ".x"), drop = FALSE])
          mu_mat = as.matrix(euc_segment[, paste0(measure_nm, ".y"), drop = FALSE])
          dist = sqrt(rowSums((x_mat - mu_mat)^2))
          dist[dist == 0] = NA_real_
          
          tibble(dist_bet_euc = dist) %>%
            bind_cols(d1 %>% select(!!group_q, !!speaker_q, !!contrast_q, all_of(token_id_nm), !!condition_q), .)
        })
      }
      
      if (isTRUE(auto_two_group)) {
        euc_segment_all = bind_rows(
          .euc_between_one_direction(ref_g = two_groups[2], focal_gs = two_groups[1]),
          .euc_between_one_direction(ref_g = two_groups[1], focal_gs = two_groups[2])
        )
      } else {
        euc_segment_all = .euc_between_one_direction(ref_g = reference_group, focal_gs = focal_groups)
      }
      
      xdata = xdata %>%
        left_join(euc_segment_all, by = c(group_nm, speaker_nm, contrast_nm, token_id_nm, condition_nm) %>% unique())
      
      rm(euc_segment_all, .euc_between_one_direction)
    }
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
        ) %>%
        dplyr::mutate(cloud = .clean_cloud_label(.data$mah_ref_cloud_wit))
      
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
        ) %>%
        dplyr::mutate(cloud = .clean_cloud_label(.data$mah_ref_cloud_bet))
      
      rm(pool_all, mah_between_all, .mah_between_one_direction)
    }
  }
  
  # =========================
  # PILLAI
  # =========================
  
  if (isTRUE(compute_pillai)) {
    
    if (!inter_group) {
      base_pillai = xdata %>%
        rename(contrast_raw = !!contrast_q) %>%
        select(!!speaker_q, !!condition_q, all_of(measure_nm), contrast_raw)
      
      tok_lvl1 = base_pillai %>% filter(contrast_raw == lvl1) %>% select(-contrast_raw)
      tok_lvl2 = base_pillai %>% filter(contrast_raw == lvl2) %>% select(-contrast_raw)
      
      pillai_errors <- list()
      
      if (within_ref == "self_speaker") {
        pillai = base_pillai %>%
          mutate(contrast_tmp = contrast_raw) %>%
          select(-contrast_raw) %>%
          group_by(!!condition_q, !!speaker_q) %>%
          nest() %>%
          rename(token_data = data) %>%
          mutate(
            dist_wit_pil = purrr::map2_dbl(
              token_data,
              paste(.data[[speaker_nm]], .data[[condition_nm]], sep = " || "),
              ~ .pillai_two_clouds(.x, "contrast_tmp", debug = diagnostics, debug_id = .y)
            )
          ) %>%
          ungroup() %>%
          select(!!speaker_q, !!condition_q, dist_wit_pil)
        
      } else if (within_ref == "inter_speaker_mean") {
        targets = base_pillai %>%
          distinct(!!speaker_q, !!condition_q) %>%
          rename(target_speaker = !!speaker_q)
        
        ref_index = base_pillai %>%
          distinct(!!speaker_q, !!condition_q) %>%
          rename(ref_speaker = !!speaker_q)
        
        pairs = targets %>%
          left_join(ref_index, by = condition_nm, relationship = "many-to-many") %>%
          filter(.data[["target_speaker"]] != .data[["ref_speaker"]]) %>%
          left_join(tok_lvl1 %>% rename(target_speaker = !!speaker_q),
                    by = c("target_speaker", condition_nm),
                    relationship = "many-to-many") %>%
          left_join(tok_lvl2 %>% rename(ref_speaker = !!speaker_q),
                    by = c("ref_speaker", condition_nm),
                    relationship = "many-to-many",
                    suffix = c(".A", ".B"))
        
        pillai_pairs = pairs %>%
          group_by(target_speaker, ref_speaker, !!condition_q) %>%
          nest() %>%
          rename(token_data = data) %>%
          mutate(pil = purrr::map_dbl(token_data, ~ {
            dfA = dplyr::select(.x, dplyr::ends_with(".A")) %>% rlang::set_names(measure_nm)
            dfB = dplyr::select(.x, dplyr::ends_with(".B")) %>% rlang::set_names(measure_nm)
            if (nrow(dfA) == 0 || nrow(dfB) == 0) return(NA_real_)
            df_ab = bind_rows(
              mutate(dfA, contrast_tmp = "A"),
              mutate(dfB, contrast_tmp = "B")
            )
            .pillai_two_clouds(df_ab, "contrast_tmp")
          })) %>%
          ungroup() %>%
          select(target_speaker, !!condition_q, pil)
        
        pillai = pillai_pairs %>%
          group_by(target_speaker, !!condition_q) %>%
          summarise(dist_wit_pil = mean(pil, na.rm = TRUE), .groups = "drop") %>%
          mutate(dist_wit_pil = ifelse(is.nan(dist_wit_pil), NA_real_, dist_wit_pil)) %>%
          rename(!!speaker_nm := target_speaker)
        
        rm(targets, ref_index, pairs, pillai_pairs)
        
      } else {
        targets = base_pillai %>%
          distinct(!!speaker_q, !!condition_q) %>%
          rename(target_speaker = !!speaker_q)
        
        pooled_lvl2 = targets %>%
          left_join(tok_lvl2 %>% rename(ref_speaker = !!speaker_q),
                    by = condition_nm,
                    relationship = "many-to-many") %>%
          filter(.data[["target_speaker"]] != .data[["ref_speaker"]]) %>%
          group_by(target_speaker, !!condition_q) %>%
          summarise(dataB = list(select(pick(everything()), all_of(measure_nm))), .groups = "drop")
        
        pillai = targets %>%
          left_join(tok_lvl1 %>% rename(target_speaker = !!speaker_q),
                    by = c("target_speaker", condition_nm),
                    relationship = "many-to-many") %>%
          group_by(target_speaker, !!condition_q) %>%
          nest() %>%
          rename(token_data = data) %>%
          left_join(pooled_lvl2, by = c("target_speaker", condition_nm)) %>%
          mutate(dist_wit_pil = purrr::map2_dbl(token_data, dataB, ~ {
            dfA = select(.x, all_of(measure_nm))
            dfB = bind_rows(.y)
            if (nrow(dfA) == 0 || nrow(dfB) == 0) return(NA_real_)
            df_ab = bind_rows(
              mutate(dfA, contrast_tmp = "A"),
              mutate(dfB, contrast_tmp = "B")
            )
            .pillai_two_clouds(df_ab, "contrast_tmp")
          })) %>%
          ungroup() %>%
          select(target_speaker, !!condition_q, dist_wit_pil) %>%
          rename(!!speaker_nm := target_speaker)
        
        rm(targets, pooled_lvl2)
      }
      
      xdata = xdata %>% left_join(pillai, by = c(speaker_nm, condition_nm))
      rm(base_pillai, tok_lvl1, tok_lvl2, pillai)
    }
    
    if (inter_group) {
      .pillai_group_tmp = function(dat) {
        .pillai_two_clouds(dat, "group_tmp", debug = diagnostics, debug_id = "Inter-Group")
        dat = dat %>% filter(if_all(all_of(measure_nm), ~ !is.na(.x)))
        if (nrow(dat) <= length(measure_nm)) return(NA_real_)
        if (n_distinct(dat$group_tmp) < 2) return(NA_real_)
        tryCatch({
          summary(
            manova(as.matrix(select(dat, all_of(measure_nm))) ~ group_tmp, data = dat)
          )$stats["group_tmp", "Pillai"] %>% as.numeric()
        }, error = function(e) NA_real_)
      }
      
      .pillai_group_tmp = function(dat) {
        # Simplemente llamamos a la función robusta y que ella se encargue de todo
        .pillai_two_clouds(
          df_ab    = dat, 
          class_col = "group_tmp", 
          debug     = diagnostics, 
          debug_id  = "Inter-Group"
        )
      }
      
      .pillai_between_one_direction = function(ref_g, focal_gs) {
        ref_has_var = purrr::map_lgl(condition_vars_nm, function(vnm) {
          any(!is.na(xdata[[vnm]][xdata[[group_nm]] == ref_g]))
        })
        names(ref_has_var) = condition_vars_nm
        
        ref_pool = xdata %>%
          filter(.data[[group_nm]] == ref_g) %>%
          select(all_of(condition_vars_nm), !!contrast_q, all_of(measure_nm)) %>%
          filter(if_all(all_of(measure_nm), ~ !is.na(.x)))
        
        purrr::map_dfr(focal_gs, function(fg) {
          focal_dat = xdata %>%
            filter(.data[[group_nm]] == fg) %>%
            select(!!speaker_q, !!condition_q, all_of(condition_vars_nm), !!contrast_q, all_of(measure_nm)) %>%
            filter(if_all(all_of(measure_nm), ~ !is.na(.x)))
          
          if (nrow(focal_dat) == 0) return(tibble())
          
          pillai_stratum_spk = focal_dat %>%
            group_by(!!speaker_q, !!condition_q, across(all_of(condition_vars_nm)), !!contrast_q) %>%
            group_modify(~{
              focal_meas = select(.x, all_of(measure_nm))
              ref_sub = ref_pool %>% filter(.data[[contrast_nm]] == .y[[contrast_nm]])
              
              for (vnm in condition_vars_nm) {
                if (isTRUE(ref_has_var[[vnm]])) {
                  val = .y[[vnm]]
                  if (is.na(val)) ref_sub = ref_sub %>% filter(is.na(.data[[vnm]]))
                  else ref_sub = ref_sub %>% filter(.data[[vnm]] == val)
                }
              }
              
              dat_ab = bind_rows(
                mutate(focal_meas, group_tmp = "focal"),
                mutate(select(ref_sub, all_of(measure_nm)), group_tmp = "reference")
              )
              
              tibble(pillai = .pillai_group_tmp(dat_ab))
            }) %>%
            ungroup() %>%
            mutate(!!rlang::sym(group_nm) := fg) %>%
            select(!!rlang::sym(group_nm), !!speaker_q, !!condition_q, !!contrast_q, pillai)
          
          pillai_fb_spk = focal_dat %>%
            group_by(!!speaker_q, !!contrast_q) %>%
            nest() %>%
            rename(token_data = data) %>%
            ungroup() %>%
            mutate(
              pillai_fb = purrr::map2_dbl(
                token_data, .data[[contrast_nm]],
                function(dat_focal, contr_val) {
                  focal_meas = select(dat_focal, all_of(measure_nm))
                  ref_sub = ref_pool %>% filter(.data[[contrast_nm]] == contr_val)
                  dat_ab = bind_rows(
                    mutate(focal_meas, group_tmp = "focal"),
                    mutate(select(ref_sub, all_of(measure_nm)), group_tmp = "reference")
                  )
                  .pillai_group_tmp(dat_ab)
                }
              )
            ) %>%
            mutate(!!rlang::sym(group_nm) := fg) %>%
            select(!!rlang::sym(group_nm), !!speaker_q, !!contrast_q, pillai_fb)
          
          pillai_stratum_spk %>%
            left_join(pillai_fb_spk, by = c(group_nm, speaker_nm, contrast_nm)) %>%
            mutate(pillai = coalesce(pillai, pillai_fb)) %>%
            select(-pillai_fb)
        })
      }
      
      if (isTRUE(auto_two_group)) {
        pillai_all = bind_rows(
          .pillai_between_one_direction(ref_g = two_groups[2], focal_gs = two_groups[1]),
          .pillai_between_one_direction(ref_g = two_groups[1], focal_gs = two_groups[2])
        )
      } else {
        pillai_all = .pillai_between_one_direction(ref_g = reference_group, focal_gs = focal_groups)
      }
      
      xdata = xdata %>%
        left_join(pillai_all, by = c(group_nm, speaker_nm, condition_nm, contrast_nm)) %>%
        rename(dist_bet_pil = pillai) %>%
        {
          if (isTRUE(auto_two_group)) . else mutate(., dist_bet_pil = if_else(!!group_q == reference_group, NA_real_, dist_bet_pil))
        } %>%
        select(-any_of("pillai"))
      
      rm(pillai_all, .pillai_between_one_direction, .pillai_group_tmp)
    }
  }
  
  # --- Lógica de Limpieza y Diagnóstico Post-Cálculo ---
  
  if (isTRUE(compute_pillai)) {
    # 1. Identificar qué columna de Pillai existe en xdata
    # Buscamos 'dist_bet_pil' si es inter_group, o 'dist_wit_pil' si no.
    target_col <- if (inter_group) "dist_bet_pil" else "dist_wit_pil"
    
    if (target_col %in% names(xdata)) {
      
      # 2. Extraer motivos y verificar éxitos
      reasons <- purrr::map_chr(xdata[[target_col]], ~ {
        r <- attr(.x, "reason")
        if (!is.null(r)) return(r)
        if (is.na(.x) || is.nan(.x)) return("numerical_failure/NA")
        return("Success")
      })
      
      n_total <- length(reasons)
      n_success <- sum(reasons == "Success")
      all_failed <- (n_success == 0)
      
      # 3. Reporte de Diagnóstico
      if (isTRUE(diagnostics)) {
        type_label <- if (inter_group) "BETWEEN" else "WITHIN"
        cat(sprintf("\n--- PILLAI DIAGNOSTIC SUMMARY (%s) ---\n", type_label))
        cat(sprintf("Total cases attempted: %d\n", n_total))
        
        # Evitar NaN% si n_total es 0
        perc_success <- if(n_total > 0) (n_success/n_total)*100 else 0
        cat(sprintf("Success (Valid numbers): %d (%.1f%%)\n", n_success, perc_success))
        
        if (all_failed && n_total > 0) {
          cat(crayon::red("![WARNING] All calculations resulted in NaN or were skipped.\n"))
        } else if (n_success < n_total) {
          cat("\nBreakdown of failures:\n")
          print(table(reasons[reasons != "Success"]))
        }
        cat("----------------------------------------------\n")
      }
      
      # 4. Eliminación de columna si el 100% ha fallado
      if (all_failed && n_total > 0) {
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
          bh_within %>% dplyr::select(dplyr::all_of(bh_wit_join_cols), dist_wit_bha, cloud_bha = cloud),
          by = c(if (!is.na(group_nm)) group_nm else NULL, speaker_nm, condition_nm)
        ) %>%
        dplyr::mutate(
          cloud = dplyr::if_else(!is.na(.data$cloud_bha), .data$cloud_bha, .data$cloud)
        ) %>%
        dplyr::select(-dplyr::any_of("cloud_bha"))
      
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
          cloud_bha = dplyr::if_else(
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
        dplyr::mutate(
          cloud = cloud_bha
        ) %>%
        dplyr::select(-cloud_focal, -.cloud_ref_lookup, -dplyr::any_of("cloud_bha"))
      
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
  
  if ("cloud" %in% names(xdata) && !isTRUE(compute_bhatt)) {
    if (inter_group) {
      used_clouds_msg <- xdata %>%
        dplyr::filter(!is.na(.data$cloud)) %>%
        dplyr::distinct(.data$cloud) %>%
        dplyr::pull(.data$cloud)
      .emit_clouds_used("clouds used: ", used_clouds_msg)
    } else {
      used_clouds_msg <- xdata %>%
        dplyr::filter(!is.na(.data$cloud)) %>%
        dplyr::distinct(.data$cloud) %>%
        dplyr::pull(.data$cloud)
      .emit_clouds_used("clouds used: ", used_clouds_msg)
    }
  }
  
  xdata <- xdata %>%
    arrange(as.integer(.token_rowid)) %>%
    select(-matches("^n_groups(\\.x|\\.y)?$")) %>%
    select(-any_of(c(".token_rowid", ".condition_id", if (!speaker_user_supplied) ".speaker_dummy_internal" else NULL)))
  
  if (!diagnostics) {
    diagnostic_cols <- grep("_cloud_", names(xdata), value = TRUE)
    diagnostic_cols <- setdiff(diagnostic_cols, c("cloud"))
    xdata <- xdata %>% select(-all_of(diagnostic_cols))
  }
  
  xdata
}
