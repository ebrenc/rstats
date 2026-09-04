# vowelpillai(): dedicated Pillai-trace analysis for vowel overlap/separation
#
# Purpose
# -------
# `vowelpillai()` is a stand-alone companion to `voweldist()`. It focuses on
# Pillai trace and returns the information needed to interpret it responsibly:
# sample sizes, MANOVA statistics and p-values, the Stanley & Sneller (2023)
# sample-size-dependent null threshold, and (optionally) an empirical
# permutation calibration.
#
# The function deliberately does NOT return a binary "merged/distinct"
# classification. Pillai, p-values, sample size, plots, and substantive
# phonetic evidence should be interpreted together.
#
# Two modes
# ---------
# WITHIN (`inter_group = FALSE`)
#   `contrast_var` is the two-class factor whose separation is measured:
#       /ae/ vs /ʌ/
#   If `group_var` is supplied, analyses are kept within each group.
#
# BETWEEN (`inter_group = TRUE`)
#   `group_var` is the two-class factor tested by the MANOVA, while
#   `contrast_var` STRATIFIES the analysis:
#       Learner /u/ vs Native /u/
#       Learner /i/ vs Native /i/
#   Thus contrast levels are never pooled across groups.
#
# Cloud resolution
# ----------------
# Within mode follows the same semantics adopted for `voweldist()`:
#
#   initial_speaker_stratification = TRUE
#       condition_vars × speaker -> condition_vars -> ... -> group pool
#
#   initial_speaker_stratification = FALSE
#       condition_vars -> ... -> group pool
#
#   exclude_focal_speaker_from_ref = TRUE
#       constructs a leave-one-speaker-out within reference; therefore a
#       speaker-specific reference cloud is skipped automatically.
#
# Between mode has its own hierarchy:
#   - focal cloud: condition_vars × speaker -> condition_vars -> ... -> focal pool
#   - reference cloud: condition_vars -> ... -> reference pool
# This allows, for example, a learner's /u/ to be compared with pooled Native
# /u/, while retaining matching conditions whenever the reference data permit.
#
# Stanley & Sneller (2023)
# ------------------------
# For a bivariate F1-F2 Pillai analysis with two classes and total sample size n,
# Stanley & Sneller recommend the 95th-percentile null threshold
#
#     p95 = exp(1) / m
#
# where m = n / 2 is the mean class size. Equivalently,
#
#     p95 = 2 * exp(1) / n.
#
# Their simulations used two dependent variables and class sizes from 5 to 100.
# `vowelpillai()` therefore:
#   - reports the threshold only when there are exactly 2 dependent variables;
#   - reports whether the observed class sizes are inside or outside the
#     simulated 5--100-per-class range;
#   - labels the result only as above/at-or-below the 95% null threshold,
#     rather than automatically declaring a merger or distinction.
#
# Their simulations also showed that TOTAL n matters much more than imbalance
# between the two classes. For this reason, `vowelpillai()` always reports
# n_1, n_2, n_total, mean_n_group and n_balance.
#
# Permutation calibration
# -----------------------
# `calibration = "permutation"` or `"both"` additionally permutes the two class
# labels within each resolved comparison, preserving the observed class sizes.
# It returns an empirical 95th-percentile null Pillai and an empirical p-value.
# This is useful as a data-specific complement to the Stanley-Sneller formula
# and can also be used when the analysis has more than two dependent variables.
#
# Important scope note
# --------------------
# The MANOVA fitted here has ONE predictor: the tested two-class factor.
# `condition_vars` are stratifiers, not covariates in the MANOVA. This mirrors
# the models used in the Stanley & Sneller sample-size simulations. Model-based
# Pillai with additional fixed effects is intentionally left for a later
# extension rather than silently applying the Stanley-Sneller threshold outside
# the setting in which it was derived.
#
# Output
# ------
# One row is returned per requested analysis unit. Important columns include:
#
#   mode
#   factor_tested
#   level_1, level_2
#   requested_cloud, cloud
#   focal_cloud, reference_cloud        # between mode
#   n_input, n_excluded
#   n_1, n_2, n_total, mean_n_group, n_balance
#   n_dimensions
#   pillai
#   approx_F, num_df, den_df, manova_p
#   ss_null95, ss_above_null95, ss_relation, ss_scope
#   perm_null95, perm_p, perm_n_valid
#   status
#
# If `diagnostics = TRUE`, `failure_reason` and `fallback_log` are retained.
# If `keep_resamples = TRUE`, permutation Pillai values are retained in the
# list-column `perm_values`.
#
# Example: within-speaker vowel contrast
# --------------------------------------
# df %>%
#   vowelpillai(
#     dependent_vars = c(F1, F2),
#     contrast_var = Vowel,
#     contrast_levels = c("æ", "ʌ"),
#     condition_vars = c(Time, WordAge),
#     speaker_var = Speaker,
#     group_var = Condition,
#     initial_speaker_stratification = TRUE,
#     calibration = "both",
#     n_perm = 999
#   )
#
# Example: same vowel across groups
# ---------------------------------
# df %>%
#   vowelpillai(
#     dependent_vars = c(F1, F2),
#     inter_group = TRUE,
#     contrast_var = Vowel,
#     condition_vars = Time,
#     speaker_var = Speaker,
#     group_var = Condition,
#     reference_group = "Native",
#     focal_group = "Learner",
#     calibration = "both"
#   )
#
vowelpillai <- function(
    x,
    dependent_vars,
    inter_group = FALSE,
    contrast_var,
    contrast_levels = NULL,
    condition_vars = NULL,
    speaker_var = NULL,
    group_var = NULL,
    reference_group = NULL,
    focal_group = NULL,
    initial_speaker_stratification = TRUE,
    exclude_focal_speaker_from_ref = FALSE,
    calibration = c("stanley_sneller", "permutation", "both", "none"),
    n_perm = 999,
    keep_resamples = FALSE,
    seed = NULL,
    diagnostics = FALSE
) {

  # ---------------------------------------------------------------------------
  # Dependencies and basic validation
  # ---------------------------------------------------------------------------

  needed_pkgs <- c("dplyr", "tibble", "tidyr", "purrr", "rlang", "tidyselect")
  missing_pkgs <- needed_pkgs[
    !vapply(needed_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (length(missing_pkgs) > 0L) {
    stop(
      "Missing required package(s): ",
      paste(missing_pkgs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  calibration <- match.arg(calibration)

  .check_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  .check_flag(inter_group, "inter_group")
  .check_flag(initial_speaker_stratification, "initial_speaker_stratification")
  .check_flag(exclude_focal_speaker_from_ref, "exclude_focal_speaker_from_ref")
  .check_flag(keep_resamples, "keep_resamples")
  .check_flag(diagnostics, "diagnostics")

  if (calibration %in% c("permutation", "both")) {
    if (
      length(n_perm) != 1L ||
      !is.numeric(n_perm) ||
      is.na(n_perm) ||
      !is.finite(n_perm) ||
      n_perm < 1 ||
      n_perm != as.integer(n_perm)
    ) {
      stop(
        "`n_perm` must be a positive integer when permutation calibration is used.",
        call. = FALSE
      )
    }
    n_perm <- as.integer(n_perm)
  }

  if (!is.null(seed)) {
    if (
      length(seed) != 1L ||
      !is.numeric(seed) ||
      is.na(seed) ||
      !is.finite(seed)
    ) {
      stop("`seed` must be NULL or one finite number.", call. = FALSE)
    }

    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_random_seed) {
      old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }

    on.exit({
      if (had_random_seed) {
        assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(seed)
  }

  xdata <- tibble::as_tibble(x)

  if (nrow(xdata) == 0L) {
    stop("`x` has no rows.", call. = FALSE)
  }

  dep_q <- rlang::enquo(dependent_vars)
  dep_sel <- tidyselect::eval_select(dep_q, xdata)
  measure_nm <- names(dep_sel)

  if (length(measure_nm) == 0L) {
    stop("`dependent_vars` selected no columns.", call. = FALSE)
  }

  non_numeric <- measure_nm[
    !vapply(xdata[measure_nm], is.numeric, FUN.VALUE = logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "All `dependent_vars` must be numeric. Non-numeric: ",
      paste(non_numeric, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  contrast_q <- rlang::enquo(contrast_var)
  if (rlang::quo_is_missing(contrast_q) || rlang::quo_is_null(contrast_q)) {
    stop("`contrast_var` is required.", call. = FALSE)
  }
  contrast_nm <- rlang::as_name(contrast_q)

  speaker_q <- rlang::enquo(speaker_var)
  speaker_user_supplied <- !rlang::quo_is_null(speaker_q)
  speaker_nm <- if (speaker_user_supplied) rlang::as_name(speaker_q) else NULL

  group_q <- rlang::enquo(group_var)
  group_user_supplied <- !rlang::quo_is_null(group_q)
  group_nm <- if (group_user_supplied) rlang::as_name(group_q) else NULL

  cond_q <- rlang::enquo(condition_vars)
  cond_sel <- tidyselect::eval_select(cond_q, xdata)
  condition_nm <- names(cond_sel)

  required_cols <- unique(c(
    measure_nm,
    contrast_nm,
    speaker_nm,
    group_nm,
    condition_nm
  ))
  absent_cols <- setdiff(required_cols, names(xdata))
  if (length(absent_cols) > 0L) {
    stop(
      "Column(s) not found in `x`: ",
      paste(absent_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Keep class-like variables stable and easy to compare.
  xdata[[contrast_nm]] <- as.character(xdata[[contrast_nm]])
  if (speaker_user_supplied) {
    xdata[[speaker_nm]] <- as.character(xdata[[speaker_nm]])
  }
  if (group_user_supplied) {
    xdata[[group_nm]] <- as.character(xdata[[group_nm]])
  }

  xdata$.vp_rowid <- seq_len(nrow(xdata))

  if (isTRUE(exclude_focal_speaker_from_ref) && !speaker_user_supplied) {
    stop(
      "`exclude_focal_speaker_from_ref = TRUE` requires `speaker_var`.",
      call. = FALSE
    )
  }

  if (isTRUE(inter_group) && isTRUE(exclude_focal_speaker_from_ref)) {
    warning(
      "`exclude_focal_speaker_from_ref` applies only to within-group analyses ",
      "and is ignored when `inter_group = TRUE`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Contrast levels
  # ---------------------------------------------------------------------------

  present_contrast_levels <- unique(xdata[[contrast_nm]])
  present_contrast_levels <- present_contrast_levels[!is.na(present_contrast_levels)]

  if (!isTRUE(inter_group)) {
    if (is.null(contrast_levels)) {
      if (length(present_contrast_levels) != 2L) {
        stop(
          "Within mode requires exactly two levels of `contrast_var`, or ",
          "explicit `contrast_levels = c(level1, level2)`. Present levels: ",
          paste(present_contrast_levels, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      contrast_levels_use <- present_contrast_levels
    } else {
      contrast_levels_use <- as.character(contrast_levels)
      if (length(contrast_levels_use) != 2L || anyNA(contrast_levels_use)) {
        stop(
          "Within mode requires `contrast_levels` to contain exactly two non-NA levels.",
          call. = FALSE
        )
      }
      if (anyDuplicated(contrast_levels_use)) {
        stop("The two `contrast_levels` must be different.", call. = FALSE)
      }
      absent_levels <- setdiff(contrast_levels_use, present_contrast_levels)
      if (length(absent_levels) > 0L) {
        stop(
          "`contrast_levels` not found in `contrast_var`: ",
          paste(absent_levels, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }

    xdata <- xdata[
      !is.na(xdata[[contrast_nm]]) &
        xdata[[contrast_nm]] %in% contrast_levels_use,
      ,
      drop = FALSE
    ]

  } else {

    if (is.null(contrast_levels)) {
      contrast_levels_use <- present_contrast_levels
    } else {
      contrast_levels_use <- as.character(contrast_levels)
      if (length(contrast_levels_use) < 1L || anyNA(contrast_levels_use)) {
        stop(
          "In between mode, `contrast_levels` must contain at least one non-NA level.",
          call. = FALSE
        )
      }
      contrast_levels_use <- unique(contrast_levels_use)
      absent_levels <- setdiff(contrast_levels_use, present_contrast_levels)
      if (length(absent_levels) > 0L) {
        stop(
          "`contrast_levels` not found in `contrast_var`: ",
          paste(absent_levels, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }

    xdata <- xdata[
      !is.na(xdata[[contrast_nm]]) &
        xdata[[contrast_nm]] %in% contrast_levels_use,
      ,
      drop = FALSE
    ]
  }

  # ---------------------------------------------------------------------------
  # Group handling for between mode
  # ---------------------------------------------------------------------------

  directions <- NULL

  if (isTRUE(inter_group)) {

    if (!group_user_supplied) {
      stop("`inter_group = TRUE` requires `group_var`.", call. = FALSE)
    }

    present_groups <- unique(xdata[[group_nm]])
    present_groups <- present_groups[!is.na(present_groups)]

    if (length(present_groups) < 2L) {
      stop("Between mode requires at least two non-NA groups.", call. = FALSE)
    }

    if (!is.null(reference_group)) {
      reference_group <- as.character(reference_group)
      if (length(reference_group) != 1L || is.na(reference_group)) {
        stop("`reference_group` must be NULL or one non-NA value.", call. = FALSE)
      }
      if (!reference_group %in% present_groups) {
        stop(
          "`reference_group = \"", reference_group,
          "\"` is not present in `", group_nm, "`.",
          call. = FALSE
        )
      }
    }

    if (!is.null(focal_group)) {
      focal_group <- as.character(focal_group)
      if (anyNA(focal_group) || length(focal_group) < 1L) {
        stop("`focal_group` must contain non-NA group value(s).", call. = FALSE)
      }
      absent_focal <- setdiff(focal_group, present_groups)
      if (length(absent_focal) > 0L) {
        stop(
          "`focal_group` value(s) not present: ",
          paste(absent_focal, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }

    if (is.null(reference_group) && is.null(focal_group)) {
      if (length(present_groups) != 2L) {
        stop(
          "If `reference_group` and `focal_group` are both NULL, ",
          "`group_var` must have exactly two levels.",
          call. = FALSE
        )
      }

      # Same auto-two-group philosophy as voweldist(): both directions.
      directions <- tibble::tibble(
        focal_group = c(present_groups[1], present_groups[2]),
        reference_group = c(present_groups[2], present_groups[1])
      )

    } else if (is.null(reference_group) && !is.null(focal_group)) {

      ref_guess <- setdiff(present_groups, focal_group)
      if (length(ref_guess) != 1L) {
        stop(
          "`focal_group` does not uniquely determine one reference group. ",
          "Provide `reference_group`.",
          call. = FALSE
        )
      }

      directions <- tibble::tibble(
        focal_group = focal_group,
        reference_group = ref_guess
      )

    } else {

      focal_groups_use <- if (is.null(focal_group)) {
        setdiff(present_groups, reference_group)
      } else {
        focal_group
      }

      focal_groups_use <- setdiff(focal_groups_use, reference_group)

      if (length(focal_groups_use) == 0L) {
        stop("No focal groups remain after excluding the reference group.", call. = FALSE)
      }

      directions <- tibble::tibble(
        focal_group = focal_groups_use,
        reference_group = rep(reference_group, length(focal_groups_use))
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Cloud helpers
  # ---------------------------------------------------------------------------

  .value_is_na <- function(x) {
    length(x) == 0L || all(is.na(x))
  }

  .spec_label <- function(cols) {
    if (length(cols) == 0L) return("global")
    paste(cols, collapse = " * ")
  }

  .make_specs <- function(
      condition_cols,
      speaker_col = NULL,
      include_speaker = FALSE
  ) {

    specs <- list()

    if (length(condition_cols) > 0L) {

      if (isTRUE(include_speaker) && !is.null(speaker_col)) {
        specs[[length(specs) + 1L]] <- c(condition_cols, speaker_col)
      }

      specs[[length(specs) + 1L]] <- condition_cols

      if (length(condition_cols) > 1L) {
        for (i in seq.int(length(condition_cols) - 1L, 1L)) {
          specs[[length(specs) + 1L]] <- condition_cols[seq_len(i)]
        }
      }

    } else if (isTRUE(include_speaker) && !is.null(speaker_col)) {

      specs[[length(specs) + 1L]] <- speaker_col
    }

    specs[[length(specs) + 1L]] <- character(0)

    signatures <- vapply(
      specs,
      function(z) paste(z, collapse = "\r"),
      FUN.VALUE = character(1)
    )
    specs <- specs[!duplicated(signatures)]
    names(specs) <- vapply(specs, .spec_label, FUN.VALUE = character(1))
    specs
  }

  .subset_by_spec <- function(pool, row_meta, spec_cols) {

    out <- pool

    if (length(spec_cols) == 0L) {
      return(out)
    }

    valid_cols <- spec_cols[
      spec_cols %in% names(pool) &
        spec_cols %in% names(row_meta)
    ]

    if (length(valid_cols) == 0L) {
      return(out)
    }

    for (v in valid_cols) {
      val <- row_meta[[v]][1]

      if (is.na(val)) {
        out <- out[is.na(out[[v]]), , drop = FALSE]
      } else {
        out <- out[!is.na(out[[v]]) & out[[v]] == val, , drop = FALSE]
      }
    }

    out
  }

  .effective_spec_label <- function(spec_cols, row_meta) {

    if (length(spec_cols) == 0L) return("global")

    present <- spec_cols[
      spec_cols %in% names(row_meta)
    ]

    if (length(present) == 0L) return("global")

    keep <- present[
      vapply(
        present,
        function(v) !.value_is_na(row_meta[[v]][1]),
        FUN.VALUE = logical(1)
      )
    ]

    if (length(keep) == 0L) "global" else paste(keep, collapse = " * ")
  }

  .requested_cloud_label <- function(
      row_meta,
      include_speaker = FALSE
  ) {
    cols <- condition_nm
    if (isTRUE(include_speaker) && speaker_user_supplied) {
      cols <- c(cols, speaker_nm)
    }
    .effective_spec_label(cols, row_meta)
  }

  .complete_cloud <- function(df) {
    if (nrow(df) == 0L) return(df)
    cc <- stats::complete.cases(df[, measure_nm, drop = FALSE])
    df[cc, , drop = FALSE]
  }

  # ---------------------------------------------------------------------------
  # Pillai and MANOVA helpers
  # ---------------------------------------------------------------------------

  .pillai_direct <- function(Y, cls) {

    Y <- as.matrix(Y)
    cls <- as.character(cls)

    levs <- unique(cls)
    levs <- levs[!is.na(levs)]

    if (length(levs) != 2L) return(NA_real_)

    X1 <- Y[cls == levs[1], , drop = FALSE]
    X2 <- Y[cls == levs[2], , drop = FALSE]

    n1 <- nrow(X1)
    n2 <- nrow(X2)

    if (n1 < 2L || n2 < 2L) return(NA_real_)

    S1 <- tryCatch(stats::cov(X1), error = function(e) NULL)
    S2 <- tryCatch(stats::cov(X2), error = function(e) NULL)

    if (is.null(S1) || is.null(S2)) return(NA_real_)

    if (is.null(dim(S1))) {
      S1 <- matrix(S1, 1L, 1L)
      S2 <- matrix(S2, 1L, 1L)
    }

    if (
      anyNA(S1) || anyNA(S2) ||
      any(!is.finite(S1)) || any(!is.finite(S2))
    ) {
      return(NA_real_)
    }

    E <- (n1 - 1) * S1 + (n2 - 1) * S2
    d <- matrix(colMeans(X1) - colMeans(X2), ncol = 1L)
    H <- (n1 * n2) / (n1 + n2) * (d %*% t(d))
    Tmat <- E + H

    if (
      anyNA(Tmat) ||
      any(!is.finite(Tmat))
    ) {
      return(NA_real_)
    }

    s <- tryCatch(svd(Tmat), error = function(e) NULL)
    if (is.null(s) || anyNA(s$d) || any(!is.finite(s$d))) {
      return(NA_real_)
    }

    tol <- max(dim(Tmat)) * max(s$d, 1) * .Machine$double.eps
    inv_d <- ifelse(s$d > tol, 1 / s$d, 0)

    Tinv <- tryCatch(
      s$v %*% (diag(inv_d, nrow = length(inv_d))) %*% t(s$u),
      error = function(e) NULL
    )

    if (
      is.null(Tinv) ||
      anyNA(Tinv) ||
      any(!is.finite(Tinv))
    ) {
      return(NA_real_)
    }

    out <- tryCatch(
      sum(diag(H %*% Tinv)),
      error = function(e) NA_real_
    )

    if (!is.finite(out)) return(NA_real_)
    as.numeric(max(0, min(out, 1)))
  }

  .stat_value <- function(x, name) {
    if (is.null(names(x)) || !name %in% names(x)) return(NA_real_)
    as.numeric(unname(x[[name]]))
  }

  .fit_pair <- function(df1, df2, level1, level2) {

    n_input_1 <- nrow(df1)
    n_input_2 <- nrow(df2)

    df1c <- .complete_cloud(df1)
    df2c <- .complete_cloud(df2)

    n1 <- nrow(df1c)
    n2 <- nrow(df2c)
    n_total <- n1 + n2
    p <- length(measure_nm)

    fail <- function(reason) {
      list(
        ok = FALSE,
        reason = reason,
        n_input_1 = n_input_1,
        n_input_2 = n_input_2,
        n_1 = n1,
        n_2 = n2,
        n_total = n_total,
        n_excluded = (n_input_1 + n_input_2) - n_total,
        data_1 = df1c,
        data_2 = df2c
      )
    }

    if (n1 == 0L) return(fail("empty_level_1"))
    if (n2 == 0L) return(fail("empty_level_2"))
    if (n1 < 2L || n2 < 2L) return(fail("class_with_lt_2_observations"))
    if (n_total < (p + 2L)) return(fail("insufficient_observations"))

    Y <- rbind(
      as.matrix(df1c[, measure_nm, drop = FALSE]),
      as.matrix(df2c[, measure_nm, drop = FALSE])
    )

    cls <- factor(
      c(rep(level1, n1), rep(level2, n2)),
      levels = c(level1, level2)
    )

    if (any(!is.finite(Y))) {
      return(fail("non_finite_values"))
    }

    direct_pillai <- .pillai_direct(Y, cls)
    if (!is.finite(direct_pillai)) {
      return(fail("pillai_not_finite"))
    }

    man <- tryCatch(
      stats::manova(Y ~ cls),
      error = function(e) e
    )

    if (inherits(man, "error")) {
      return(fail(paste0("manova_error: ", conditionMessage(man))))
    }

    sm <- tryCatch(
      suppressWarnings(summary(man, test = "Pillai")),
      error = function(e) e
    )

    if (inherits(sm, "error") || is.null(sm$stats) || nrow(sm$stats) < 1L) {
      rr <- if (inherits(sm, "error")) {
        conditionMessage(sm)
      } else {
        "MANOVA summary did not return statistics"
      }
      return(fail(paste0("manova_summary_error: ", rr)))
    }

    stat_row <- sm$stats[1L, , drop = TRUE]

    pillai_manova <- .stat_value(stat_row, "Pillai")
    approx_F <- .stat_value(stat_row, "approx F")
    num_df <- .stat_value(stat_row, "num Df")
    den_df <- .stat_value(stat_row, "den Df")
    manova_p <- .stat_value(stat_row, "Pr(>F)")

    # MANOVA's Pillai is the reported statistic. The direct calculation is
    # retained internally as a consistency check and as the permutation engine.
    pillai <- pillai_manova
    if (!is.finite(pillai)) pillai <- direct_pillai

    if (!is.finite(pillai)) {
      return(fail("pillai_not_finite_after_manova"))
    }

    list(
      ok = TRUE,
      reason = NA_character_,
      n_input_1 = n_input_1,
      n_input_2 = n_input_2,
      n_1 = n1,
      n_2 = n2,
      n_total = n_total,
      n_excluded = (n_input_1 + n_input_2) - n_total,
      data_1 = df1c,
      data_2 = df2c,
      Y = Y,
      cls = cls,
      pillai = as.numeric(pillai),
      pillai_direct = as.numeric(direct_pillai),
      approx_F = approx_F,
      num_df = num_df,
      den_df = den_df,
      manova_p = manova_p
    )
  }

  # ---------------------------------------------------------------------------
  # Calibration helpers
  # ---------------------------------------------------------------------------

  .ss_calibration <- function(fit) {

    out <- list(
      ss_null95 = NA_real_,
      ss_above_null95 = NA,
      ss_relation = NA_character_,
      ss_scope = "not_requested"
    )

    if (!calibration %in% c("stanley_sneller", "both")) {
      return(out)
    }

    if (length(measure_nm) != 2L) {
      out$ss_scope <- "not_applicable_non_bivariate"
      return(out)
    }

    n1 <- fit$n_1
    n2 <- fit$n_2
    n_total <- fit$n_total

    if (!is.finite(n_total) || n_total <= 0L) {
      out$ss_scope <- "not_applicable_no_sample"
      return(out)
    }

    threshold <- exp(1) / (n_total / 2)

    if (n1 >= 5L && n1 <= 100L && n2 >= 5L && n2 <= 100L) {
      scope <- "within_simulated_5_100_per_class"
    } else if (n1 < 5L || n2 < 5L) {
      if (n1 > 100L || n2 > 100L) {
        scope <- "outside_simulated_range"
      } else {
        scope <- "extrapolated_below_simulated_range"
      }
    } else if (n1 > 100L || n2 > 100L) {
      scope <- "extrapolated_above_simulated_range"
    } else {
      scope <- "outside_simulated_range"
    }

    above <- fit$pillai > threshold

    list(
      ss_null95 = as.numeric(threshold),
      ss_above_null95 = as.logical(above),
      ss_relation = if (isTRUE(above)) "above_null95" else "at_or_below_null95",
      ss_scope = scope
    )
  }

  .perm_calibration <- function(fit) {

    out <- list(
      perm_null95 = NA_real_,
      perm_p = NA_real_,
      perm_n_valid = NA_integer_,
      perm_values = list(NULL)
    )

    if (!calibration %in% c("permutation", "both")) {
      return(out)
    }

    labels <- as.character(fit$cls)
    vals <- rep(NA_real_, n_perm)

    for (i in seq_len(n_perm)) {
      vals[i] <- .pillai_direct(
        fit$Y,
        sample(labels, size = length(labels), replace = FALSE)
      )
    }

    valid <- vals[is.finite(vals)]
    n_valid <- length(valid)

    if (n_valid == 0L) {
      out$perm_n_valid <- 0L
      if (isTRUE(keep_resamples)) out$perm_values <- list(vals)
      return(out)
    }

    null95 <- as.numeric(stats::quantile(
      valid,
      probs = 0.95,
      names = FALSE,
      type = 8
    ))

    # Add-one correction.
    perm_p <- (1 + sum(valid >= fit$pillai)) / (1 + n_valid)

    out$perm_null95 <- null95
    out$perm_p <- as.numeric(perm_p)
    out$perm_n_valid <- as.integer(n_valid)
    if (isTRUE(keep_resamples)) out$perm_values <- list(vals)
    out
  }

  .pack_success <- function(
      fit,
      mode,
      factor_tested,
      level_1,
      level_2,
      requested_cloud,
      cloud,
      focal_cloud = NA_character_,
      reference_cloud = NA_character_,
      fallback_log = NA_character_
  ) {

    ss <- .ss_calibration(fit)
    pm <- .perm_calibration(fit)

    n_balance <- if (
      fit$n_1 > 0L &&
      fit$n_2 > 0L
    ) {
      min(fit$n_1, fit$n_2) / max(fit$n_1, fit$n_2)
    } else {
      NA_real_
    }

    out <- tibble::tibble(
      mode = mode,
      factor_tested = factor_tested,
      level_1 = as.character(level_1),
      level_2 = as.character(level_2),
      requested_cloud = requested_cloud,
      cloud = cloud,
      focal_cloud = focal_cloud,
      reference_cloud = reference_cloud,
      n_input = fit$n_input_1 + fit$n_input_2,
      n_excluded = fit$n_excluded,
      n_1 = fit$n_1,
      n_2 = fit$n_2,
      n_total = fit$n_total,
      mean_n_group = fit$n_total / 2,
      n_balance = n_balance,
      n_dimensions = length(measure_nm),
      pillai = fit$pillai,
      approx_F = fit$approx_F,
      num_df = fit$num_df,
      den_df = fit$den_df,
      manova_p = fit$manova_p,
      ss_null95 = ss$ss_null95,
      ss_above_null95 = ss$ss_above_null95,
      ss_relation = ss$ss_relation,
      ss_scope = ss$ss_scope,
      perm_null95 = pm$perm_null95,
      perm_p = pm$perm_p,
      perm_n_valid = pm$perm_n_valid,
      status = "ok",
      failure_reason = NA_character_,
      fallback_log = fallback_log
    )

    if (isTRUE(keep_resamples)) {
      out$perm_values <- pm$perm_values
    }

    out
  }

  .pack_failure <- function(
      mode,
      factor_tested,
      level_1,
      level_2,
      requested_cloud,
      failure_reason,
      fallback_log = NA_character_
  ) {

    out <- tibble::tibble(
      mode = mode,
      factor_tested = factor_tested,
      level_1 = as.character(level_1),
      level_2 = as.character(level_2),
      requested_cloud = requested_cloud,
      cloud = NA_character_,
      focal_cloud = NA_character_,
      reference_cloud = NA_character_,
      n_input = NA_integer_,
      n_excluded = NA_integer_,
      n_1 = NA_integer_,
      n_2 = NA_integer_,
      n_total = NA_integer_,
      mean_n_group = NA_real_,
      n_balance = NA_real_,
      n_dimensions = length(measure_nm),
      pillai = NA_real_,
      approx_F = NA_real_,
      num_df = NA_real_,
      den_df = NA_real_,
      manova_p = NA_real_,
      ss_null95 = NA_real_,
      ss_above_null95 = NA,
      ss_relation = NA_character_,
      ss_scope = if (
        calibration %in% c("stanley_sneller", "both") &&
        length(measure_nm) != 2L
      ) {
        "not_applicable_non_bivariate"
      } else if (calibration %in% c("stanley_sneller", "both")) {
        "not_calculated_failed_fit"
      } else {
        "not_requested"
      },
      perm_null95 = NA_real_,
      perm_p = NA_real_,
      perm_n_valid = if (
        calibration %in% c("permutation", "both")
      ) 0L else NA_integer_,
      status = "no_viable_cloud",
      failure_reason = failure_reason,
      fallback_log = fallback_log
    )

    if (isTRUE(keep_resamples)) {
      out$perm_values <- list(NULL)
    }

    out
  }

  # ---------------------------------------------------------------------------
  # WITHIN resolver
  # ---------------------------------------------------------------------------

  .resolve_within_one <- function(row_meta) {

    row_group <- if (group_user_supplied && group_nm %in% names(row_meta)) {
      row_meta[[group_nm]][1]
    } else {
      NULL
    }

    pool <- xdata

    if (!is.null(row_group)) {
      if (is.na(row_group)) {
        pool <- pool[is.na(pool[[group_nm]]), , drop = FALSE]
      } else {
        pool <- pool[
          !is.na(pool[[group_nm]]) &
            pool[[group_nm]] == row_group,
          ,
          drop = FALSE
        ]
      }
    }

    include_speaker <- (
      speaker_user_supplied &&
        isTRUE(initial_speaker_stratification) &&
        !isTRUE(exclude_focal_speaker_from_ref)
    )

    specs <- .make_specs(
      condition_cols = condition_nm,
      speaker_col = speaker_nm,
      include_speaker = include_speaker
    )

    requested_cloud <- .requested_cloud_label(
      row_meta,
      include_speaker = include_speaker
    )

    failures <- character(0)

    for (spec_name in names(specs)) {

      spec <- specs[[spec_name]]
      cand <- .subset_by_spec(pool, row_meta, spec)

      if (isTRUE(exclude_focal_speaker_from_ref)) {
        focal_sp <- row_meta[[speaker_nm]][1]

        if (is.na(focal_sp)) {
          cand <- cand[!is.na(cand[[speaker_nm]]), , drop = FALSE]
        } else {
          cand <- cand[
            is.na(cand[[speaker_nm]]) |
              cand[[speaker_nm]] != focal_sp,
            ,
            drop = FALSE
          ]
        }
      }

      d1 <- cand[
        !is.na(cand[[contrast_nm]]) &
          cand[[contrast_nm]] == contrast_levels_use[1],
        ,
        drop = FALSE
      ]
      d2 <- cand[
        !is.na(cand[[contrast_nm]]) &
          cand[[contrast_nm]] == contrast_levels_use[2],
        ,
        drop = FALSE
      ]

      fit <- .fit_pair(
        d1,
        d2,
        level1 = contrast_levels_use[1],
        level2 = contrast_levels_use[2]
      )

      effective <- .effective_spec_label(spec, row_meta)

      if (isTRUE(fit$ok)) {
        return(.pack_success(
          fit = fit,
          mode = "within",
          factor_tested = contrast_nm,
          level_1 = contrast_levels_use[1],
          level_2 = contrast_levels_use[2],
          requested_cloud = requested_cloud,
          cloud = effective,
          fallback_log = if (length(failures) == 0L) {
            NA_character_
          } else {
            paste(failures, collapse = " || ")
          }
        ))
      }

      failures <- c(
        failures,
        paste0(effective, ": ", fit$reason)
      )
    }

    .pack_failure(
      mode = "within",
      factor_tested = contrast_nm,
      level_1 = contrast_levels_use[1],
      level_2 = contrast_levels_use[2],
      requested_cloud = requested_cloud,
      failure_reason = "no_viable_within_cloud",
      fallback_log = paste(failures, collapse = " || ")
    )
  }

  # ---------------------------------------------------------------------------
  # BETWEEN resolver
  # ---------------------------------------------------------------------------

  .resolve_between_one <- function(row_meta, fg, ref_g) {

    contrast_value <- row_meta[[contrast_nm]][1]

    focal_pool <- xdata[
      !is.na(xdata[[group_nm]]) &
        xdata[[group_nm]] == fg &
        !is.na(xdata[[contrast_nm]]) &
        xdata[[contrast_nm]] == contrast_value,
      ,
      drop = FALSE
    ]

    ref_pool <- xdata[
      !is.na(xdata[[group_nm]]) &
        xdata[[group_nm]] == ref_g &
        !is.na(xdata[[contrast_nm]]) &
        xdata[[contrast_nm]] == contrast_value,
      ,
      drop = FALSE
    ]

    focal_specs <- .make_specs(
      condition_cols = condition_nm,
      speaker_col = speaker_nm,
      include_speaker = speaker_user_supplied
    )

    ref_specs <- .make_specs(
      condition_cols = condition_nm,
      speaker_col = NULL,
      include_speaker = FALSE
    )

    requested_cloud <- .requested_cloud_label(
      row_meta,
      include_speaker = speaker_user_supplied
    )

    failures <- character(0)

    for (focal_spec_name in names(focal_specs)) {

      focal_spec <- focal_specs[[focal_spec_name]]
      focal_cand <- .subset_by_spec(focal_pool, row_meta, focal_spec)
      focal_effective <- .effective_spec_label(focal_spec, row_meta)

      for (ref_spec_name in names(ref_specs)) {

        ref_spec <- ref_specs[[ref_spec_name]]
        ref_cand <- .subset_by_spec(ref_pool, row_meta, ref_spec)
        ref_effective <- .effective_spec_label(ref_spec, row_meta)

        fit <- .fit_pair(
          focal_cand,
          ref_cand,
          level1 = fg,
          level2 = ref_g
        )

        if (isTRUE(fit$ok)) {

          combined_cloud <- paste0(
            "focal: ", focal_effective,
            " || reference: ", ref_effective
          )

          return(.pack_success(
            fit = fit,
            mode = "between",
            factor_tested = group_nm,
            level_1 = fg,
            level_2 = ref_g,
            requested_cloud = requested_cloud,
            cloud = combined_cloud,
            focal_cloud = focal_effective,
            reference_cloud = paste0(ref_g, ": ", ref_effective),
            fallback_log = if (length(failures) == 0L) {
              NA_character_
            } else {
              paste(failures, collapse = " || ")
            }
          ))
        }

        failures <- c(
          failures,
          paste0(
            "focal=", focal_effective,
            "; reference=", ref_effective,
            ": ", fit$reason
          )
        )
      }
    }

    .pack_failure(
      mode = "between",
      factor_tested = group_nm,
      level_1 = fg,
      level_2 = ref_g,
      requested_cloud = requested_cloud,
      failure_reason = paste0(
        "no_viable_between_cloud_for_",
        contrast_nm, "=", contrast_value
      ),
      fallback_log = paste(failures, collapse = " || ")
    )
  }

  # ---------------------------------------------------------------------------
  # Build analysis requests and run
  # ---------------------------------------------------------------------------

  if (!isTRUE(inter_group)) {

    request_needs_speaker <- (
      speaker_user_supplied &&
        (
          isTRUE(initial_speaker_stratification) ||
            isTRUE(exclude_focal_speaker_from_ref)
        )
    )

    request_cols <- unique(c(
      group_nm,
      condition_nm,
      if (request_needs_speaker) speaker_nm else NULL
    ))

    if (length(request_cols) == 0L) {
      requests <- tibble::tibble(.vp_request_dummy = 1L)
    } else {
      requests <- dplyr::distinct(
        xdata,
        dplyr::across(dplyr::all_of(request_cols))
      )
    }

    results <- purrr::map_dfr(
      seq_len(nrow(requests)),
      function(i) {
        row_meta <- requests[i, , drop = FALSE]

        ans <- .resolve_within_one(row_meta)

        meta <- row_meta[
          ,
          setdiff(names(row_meta), ".vp_request_dummy"),
          drop = FALSE
        ]

        dplyr::bind_cols(meta, ans)
      }
    )

  } else {

    direction_results <- vector("list", nrow(directions))

    for (di in seq_len(nrow(directions))) {

      fg <- directions$focal_group[di]
      ref_g <- directions$reference_group[di]

      focal_data <- xdata[
        !is.na(xdata[[group_nm]]) &
          xdata[[group_nm]] == fg,
        ,
        drop = FALSE
      ]

      request_cols <- unique(c(
        group_nm,
        condition_nm,
        if (speaker_user_supplied) speaker_nm else NULL,
        contrast_nm
      ))

      requests <- dplyr::distinct(
        focal_data,
        dplyr::across(dplyr::all_of(request_cols))
      )

      direction_results[[di]] <- purrr::map_dfr(
        seq_len(nrow(requests)),
        function(i) {

          row_meta <- requests[i, , drop = FALSE]

          ans <- .resolve_between_one(
            row_meta = row_meta,
            fg = fg,
            ref_g = ref_g
          )

          ans$focal_group <- fg
          ans$reference_group <- ref_g

          dplyr::bind_cols(row_meta, ans)
        }
      )
    }

    results <- dplyr::bind_rows(direction_results)
  }

  # ---------------------------------------------------------------------------
  # Final cleanup / ordering
  # ---------------------------------------------------------------------------

  if (!isTRUE(diagnostics)) {
    results <- dplyr::select(
      results,
      -dplyr::any_of(c("failure_reason", "fallback_log"))
    )
  }

  if (!isTRUE(keep_resamples)) {
    results <- dplyr::select(
      results,
      -dplyr::any_of("perm_values")
    )
  }

  metadata_order <- unique(c(
    group_nm,
    condition_nm,
    speaker_nm,
    contrast_nm,
    "focal_group",
    "reference_group"
  ))
  metadata_order <- metadata_order[
    metadata_order %in% names(results)
  ]

  core_order <- c(
    "mode",
    "factor_tested",
    "level_1",
    "level_2",
    "requested_cloud",
    "cloud",
    "focal_cloud",
    "reference_cloud",
    "n_input",
    "n_excluded",
    "n_1",
    "n_2",
    "n_total",
    "mean_n_group",
    "n_balance",
    "n_dimensions",
    "pillai",
    "approx_F",
    "num_df",
    "den_df",
    "manova_p",
    "ss_null95",
    "ss_above_null95",
    "ss_relation",
    "ss_scope",
    "perm_null95",
    "perm_p",
    "perm_n_valid",
    "status",
    "failure_reason",
    "fallback_log",
    "perm_values"
  )
  core_order <- core_order[core_order %in% names(results)]

  other_cols <- setdiff(
    names(results),
    c(metadata_order, core_order)
  )

  results <- dplyr::select(
    results,
    dplyr::all_of(c(metadata_order, other_cols, core_order))
  )

  if (nrow(results) > 0L) {
    used_clouds <- unique(results$cloud[!is.na(results$cloud)])
    if (length(used_clouds) > 0L) {
      pretty_clouds <- gsub(" \\* ", " × ", used_clouds)
      message("clouds used: ", paste(pretty_clouds, collapse = " || "))
    }
  }

  attr(results, "vowelpillai_calibration") <- calibration
  attr(results, "dependent_vars") <- measure_nm
  attr(results, "contrast_var") <- contrast_nm
  attr(results, "stanley_sneller_reference") <-
    "Stanley & Sneller (2023), JASA 153(1):54-67, doi:10.1121/10.0016757"

  results
}
