summarySE = function(data = NULL, measurevar, groupvars = NULL, na.rm = T, conf.interval = .95, .drop = TRUE) {
  library(dplyr)
  
  n_custom = function(x) {sum(!is.na(x))}
  mode_custom = function(x) {x = x[!is.na(x)]; ux = unique(x); ux[which.max(tabulate(match(x, ux)))]}
  
  aggdata = group_by_at(.tbl = data, .vars = groupvars) %>% 
    summarise(across(.cols = all_of(measurevar), .names = "{fn}",
                     .fns = list(
                       n = ~n_custom(.), 
                       mode = ~mode_custom(.), 
                       median = ~median(., na.rm = TRUE), 
                       min = ~min(., na.rm = TRUE), 
                       max = ~max(., na.rm = TRUE), 
                       mean = ~mean(., na.rm = TRUE), 
                       sd = ~sd(., na.rm = TRUE), 
                       q25 = ~quantile(., .25, na.rm = TRUE), 
                       q75 = ~quantile(., .75, na.rm = TRUE), 
                       q1_3 = ~quantile(., (1/3), na.rm = TRUE), 
                       q2_3 = ~quantile(., (2/3), na.rm = TRUE)
                     ))) %>% 
    mutate(se = sd / sqrt(n)) %>%
    mutate(ci = se * (qt(conf.interval/2 + .5, n-1))) %>%
    ungroup()
  return(aggdata)
}
