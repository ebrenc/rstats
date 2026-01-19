# df1 = df %>% mahalanobis_outliers(vars = c(B0, B1), group_by = c(Subject, Vowel), chisq.prob = 0.99, verbose = T, act = T)

mahalanobis_outliers = function(data,
                                vars,
                                group_by = NULL,
                                chisq.prob = 0.99,
                                act = FALSE,
                                verbose = TRUE) {
  stopifnot(!missing(vars), length(vars) >= 1)
  
  # dependències mínimes
  require(dplyr)
  require(tidyr)
  require(purrr)
  
  vars = as.character(vars)
  group_by = if (is.null(group_by)) character(0) else as.character(group_by)
  
  missing_vars = setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("These vars are not in `data`: ", paste(missing_vars, collapse = ", "))
  }
  
  missing_groups = setdiff(group_by, names(data))
  if (length(missing_groups) > 0) {
    stop("These group_by vars are not in `data`: ", paste(missing_groups, collapse = ", "))
  }
  
  threshold = stats::qchisq(chisq.prob, df = length(vars))
  
  # funció segura per mahalanobis
  safe_mah = function(df_dep) {
    # df_dep: data frame només amb vars
    x = as.matrix(df_dep)
    
    # si la fila té NA, la mahalanobis() ja pot retornar NA (ok)
    # però si no hi ha prou dades o cov singular -> NA tot el grup
    out = tryCatch({
      center = colMeans(x, na.rm = TRUE)
      covmat = stats::cov(x, use = "pairwise.complete.obs")
      
      # comprovar singularitat / no finits
      if (any(!is.finite(center)) || any(!is.finite(covmat))) {
        rep(NA_real_, nrow(x))
      } else if (isTRUE(all.equal(det(covmat), 0))) {
        rep(NA_real_, nrow(x))
      } else {
        stats::mahalanobis(x, center = center, cov = covmat)
      }
    }, error = function(e) rep(NA_real_, nrow(x)))
    
    out
  }
  
  if (length(group_by) == 0) {
    # sense grups: calcular de cop
    data = data %>%
      mutate(
        mah = safe_mah(select(., all_of(vars))),
        mah_chisq = mah > threshold
      )
  } else {
    data = data %>%
      group_by(across(all_of(group_by))) %>%
      group_modify(~{
        .x %>%
          mutate(
            mah = safe_mah(select(., all_of(vars))),
            mah_chisq = mah > threshold
          )
      }) %>%
      ungroup()
  }
  
  if (verbose) {
    n_outliers = sum(data$mah_chisq, na.rm = TRUE)
    n_rows = sum(!is.na(data$mah))
    cat(n_outliers, "/", n_rows, "outliers", "\n")
    cat("Proportion:", ifelse(n_rows == 0, NA_real_, n_outliers / n_rows), "\n")
    cat("Distance threshold:", threshold, "\n")
  }
  
  if (isTRUE(act) && length(vars) == 1) {
    v = vars[[1]]
    data = data %>%
      mutate("{v}" := if_else(mah_chisq, NA_real_, .data[[v]])) %>%
      select(-mah, -mah_chisq)
  } else if (isTRUE(act)) {
    # si vols que act=TRUE funcioni també amb múltiples vars:
    data = data %>%
      mutate(across(all_of(vars), ~ if_else(mah_chisq, NA_real_, .x))) %>%
      select(-mah, -mah_chisq)
  }
  
  data
}
