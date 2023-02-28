summarySE = function(data = NULL, measurevar, groupvars = NULL, na.rm = TRUE, conf.interval = .95, .drop = TRUE) {
  library(dplyr)
  
  n_custom = function(x) {sum(!is.na(x))}
  mode_custom = function(x) {x = x[!is.na(x)]; ux = unique(x); ux[which.max(tabulate(match(x, ux)))]}
  
  aggdata = data %>%
    group_by(across(all_of(groupvars))) %>%
    summarise(
      n = n_custom(!!sym(measurevar)),
      mode = mode_custom(!!sym(measurevar)),
      median = median(!!sym(measurevar), na.rm = na.rm),
      min = min(!!sym(measurevar), na.rm = na.rm),
      max = max(!!sym(measurevar), na.rm = na.rm),
      mean = mean(!!sym(measurevar), na.rm = na.rm),
      sd = sd(!!sym(measurevar), na.rm = na.rm),
      q25 = quantile(!!sym(measurevar), 0.25, na.rm = na.rm),
      q75 = quantile(!!sym(measurevar), 0.75, na.rm = na.rm),
      q1_3 = quantile(!!sym(measurevar), 1/3, na.rm = na.rm),
      q2_3 = quantile(!!sym(measurevar), 2/3, na.rm = na.rm)
    ) %>%
    mutate(se = sd / sqrt(n)) %>%
    mutate(ci = se * (qt(conf.interval/2 + .5, n-1))) %>%
    ungroup()
  return(aggdata)
}
