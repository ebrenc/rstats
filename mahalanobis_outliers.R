# df1 = df %>% mahalanobis_outliers(vars = c(B0, B1), group_by = c(Subject, Vowel), chisq.prob = 0.99, verbose = T)

mahalanobis_outliers = function(data, vars = NULL, group_by = NULL, chisq.prob = 0.99, act = FALSE, verbose = FALSE) {
  require(tidyverse)
  
  vars_colnames = data %>% select({{vars}}) %>% colnames()
  data = data %>% 
    nest(groups = -c({{group_by}}, {{vars_colnames}}), dependent = c({{vars_colnames}})) %>%
    mutate(mah = map(dependent, ~mahalanobis(.x, center = colMeans(.x, na.rm = TRUE), cov = cov(.x, use = "pairwise.complete.obs")))) %>%
    unnest(c(groups, dependent, mah))
  
  threshold <- qchisq(chisq.prob, length(vars_colnames))
  data$mah_chisq <- data$mah > threshold
  
  if (verbose) {
    n_outliers = sum(data$mah > threshold, na.rm = TRUE)
    n_rows = length(!is.na(data$mah))
    cat(n_outliers, "/", n_rows, "outliers", "\n")
    cat("Proportion:", n_outliers / n_rows, "\n")
    cat("Distance threshold:", threshold, "\n")
  }
  
  if (act & length(vars_colnames) == 1) {
    data = data %>% 
      mutate(across({{vars_colnames}}, ~ if_else(mah_chisq, NA_real_, .))) %>%
      select(-c(mah, mah_chisq))
  }
  
  return(data)
}
