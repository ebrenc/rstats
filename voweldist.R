# inter_group = FALSE  # If true, calculates how different two groups produce a single segment
# segments_var  # Which segmental contras do we have to measure? I.e., /æ/ vs. /ʌ/. This variable must have 2 levels only
# tests_var  # This variable must have 2 levels only
# subjects_var
# items_var
# unique_items_var  # Combined with subjects_var, these variable(s) should allow us to uniquely identify each row of the data frame
# values_var  # You can put 2 or more variables here
# groups_var = NULL  # This variable must have 2 levels only
# learners_level = NULL
# natives_level = NULL

voweldist = function(
         x,
         inter_group = FALSE,
         segments_var,
         tests_var,
         subjects_var,
         items_var,
         unique_items_var,
         values_var,
         groups_var = NULL,
         learners_level = NULL,
         natives_level = NULL) {
  
  require(tidyverse)
  options(scipen=999)
  
  xdata = x
  
  # Preliminary tasks
  
  vars_to_char = c(segments_var, tests_var, groups_var, subjects_var, items_var, unique_items_var) %>% unique() %>% na.omit()
  xdata = xdata %>% 
    arrange(across(all_of(vars_to_char))) %>% 
    mutate(across(all_of(vars_to_char), as.character))
  rm(vars_to_char)
  
  if (!is.na(segments_var)) {segments_var_levels = xdata %>% select(all_of(segments_var)) %>% unique() %>% unlist() %>% as.character()} else {segments_var_levels = 1}
  if (!is.na(tests_var)) {tests_var_levels = xdata %>% select(all_of(tests_var)) %>% unique() %>% drop_na() %>% unlist() %>% as.character()} else {tests_var_levels = 1}
  n_values_var = values_var %>% length()
  if (!is.na(segments_var) & !is.na(tests_var)) {n_clusters_factors = length(segments_var_levels)*length(tests_var_levels)} else {n_clusters_factors = 1}
  
  # Filter out subjects with less than 3 cases in each combination of segments_var * tests_var
  
  df_target_subjects = xdata %>% select(all_of(c(groups_var, subjects_var, tests_var, segments_var) %>% na.omit()))
  if (inter_group == TRUE) {df_target_subjects = df_target_subjects %>% filter(!!as.name(groups_var) == learners_level)}
  subjects_out = c(
    df_target_subjects %>% distinct() %>% group_by_at(.vars = c(subjects_var)) %>% count() %>% filter(n < n_clusters_factors) %>% pluck(subjects_var),
    df_target_subjects %>% group_by_at(.vars = c(subjects_var, tests_var, segments_var) %>% na.omit()) %>% count() %>% filter(n <= n_values_var) %>% pluck(subjects_var)
  )
  if (length(subjects_out) > 0) {
    print(c("Subjects out: ", subjects_out))
    `%notin%` = Negate(`%in%`)
    # xdata = xdata %>% filter(eval(parse(text = subjects_var)) %notin% subjects_out)
    xdata = xdata %>% filter(!!as.name(subjects_var) %notin% subjects_out)
    rm(`%notin%`)
  }
  rm(df_target_subjects, subjects_out)
  
  # Euclidean within
  
  if (inter_group == FALSE) {
    
    # Euclidean distances 1 within 2
    
    d1 = xdata %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var))
    
    d2 = xdata %>% select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var) %>% na.omit()) %>% summarise_if(is.numeric, mean, na.rm = T)
    
    euc_segment_1 = left_join(d1, d2, by = c(subjects_var, tests_var) %>% na.omit())
    
    var_matrix = colnames(euc_segment_1) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment_1[var_matrix[i,1]] - euc_segment_1[var_matrix[i,2]])^2)
    }
    euc_segment_1 = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var)) %>%
      mutate(!!segments_var := segments_var_levels[1])
    rm(var_matrix, i, vector_to_sum, d1, d2)
    
    # Euclidean distances 2 within 1
    
    d1 = xdata %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var))
    
    d2 = xdata %>% select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var) %>% na.omit()) %>% summarise_if(is.numeric, mean, na.rm = T)
    
    euc_segment_2 = left_join(d1, d2, by = c(subjects_var, tests_var) %>% na.omit())
    
    var_matrix = colnames(euc_segment_2) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment_2[var_matrix[i,1]] - euc_segment_2[var_matrix[i,2]])^2)
    }
    euc_segment_2 = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var)) %>%
      mutate(!!segments_var := segments_var_levels[2])
    rm(var_matrix, i, vector_to_sum, d1, d2)
    
    # Incorporate Euclidean distances into the original data frame
    
    euc_distances = bind_rows(euc_segment_1, euc_segment_2) %>% rename(euc_dist = distance); rm(euc_segment_1, euc_segment_2)
    xdata = xdata %>% left_join(euc_distances, by = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit() %>% unique()); rm(euc_distances)
    
  }
  
  # Euclidean between
  
  if (inter_group == TRUE) {
    
    # Euclidean distances between
    
    d1 = xdata %>%
      select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(!!as.name(groups_var) == learners_level) %>% select(!all_of(groups_var))
    
    if (!is.na(segments_var)) {
      d2 = xdata %>% select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
        group_by_at(.vars = c(segments_var)) %>% 
        summarise_if(is.numeric, mean, na.rm = T)
      euc_segment = left_join(d1, d2, by = c(segments_var))
    } else {
      d2 = xdata %>% select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
        summarise_if(is.numeric, mean, na.rm = T)
      euc_segment = bind_cols(d1, d2)
    }
    
    var_matrix = colnames(euc_segment) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment[var_matrix[i,1]] - euc_segment[var_matrix[i,2]])^2)
    }
    euc_segment = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var))
    rm(var_matrix, i, vector_to_sum, d1, d2)
    
    # Incorporate Euclidean distances into the original data frame
    
    euc_segment = euc_segment %>% rename(euc_dist = distance)
    xdata = xdata %>% left_join(euc_segment, by = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit() %>% unique()); rm(euc_segment)
    
  }
  
  # Mahalanobis within
  
  if (inter_group == FALSE) {
    
    # Mahalanobis distances 1 within 2
    
    d1 = xdata %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, unique_items_var, tests_var) %>% na.omit() %>% unique()) %>% nest()
    
    d2 = xdata %>% select(all_of(c(subjects_var, segments_var, items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var) %>% na.omit()) %>% nest() %>%
      mutate(data = map(data, ~ select(.x, -all_of(items_var)) %>% distinct())) %>%
      filter(map_lgl(data, ~ ncol(.) < nrow(.))) %>%
      group_by_at(.vars = c(subjects_var)) %>% mutate(n_tests = n()) %>% ungroup() %>% filter(n_tests == length(tests_var_levels)) %>% select(-n_tests)
    
    mah_segment_1 = left_join(d1, d2, by = c(subjects_var, tests_var) %>% na.omit()) %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data.")) %>%
      mutate(!!segments_var := segments_var_levels[1])
    rm(d1, d2)
    
    # Mahalanobis distances 2 within 1
    
    d1 = xdata %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, unique_items_var, tests_var) %>% na.omit() %>% unique()) %>% nest()
    
    d2 = xdata %>% select(all_of(c(subjects_var, segments_var, items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var) %>% na.omit()) %>% nest() %>%
      mutate(data = map(data, ~ select(.x, -all_of(items_var)) %>% distinct())) %>%
      filter(map_lgl(data, ~ ncol(.) < nrow(.))) %>%
      group_by_at(.vars = c(subjects_var)) %>% mutate(n_tests = n()) %>% ungroup() %>% filter(n_tests == length(tests_var_levels)) %>% select(-n_tests)
    
    mah_segment_2 = left_join(d1, d2, by = c(subjects_var, tests_var) %>% na.omit()) %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data.")) %>%
      mutate(!!segments_var := segments_var_levels[2])
    rm(d1, d2)
    
    # Incorporate Mahalanobis distances into the original data frame
    
    mah_distances = bind_rows(mah_segment_1, mah_segment_2) %>% rename(mah_dist = distance); rm(mah_segment_1, mah_segment_2)
    xdata = xdata %>% left_join(mah_distances, by = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit() %>% unique()); rm(mah_distances)
    
  }
  
  # Mahalanobis between
  
  if (inter_group == TRUE) {
    
    # Mahalanobis distances between
    
    if (!is.na(segments_var)) {
      d1 = xdata %>%
        select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == learners_level) %>% select(!all_of(groups_var)) %>%
        group_by_at(.vars = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit() %>% unique()) %>% nest()
      d2 = xdata %>% select(all_of(c(groups_var, segments_var, items_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
        group_by_at(.vars = c(segments_var)) %>% nest() %>%
        mutate(data = map(data, ~ select(.x, -all_of(items_var))))
      mah_segment = left_join(d1, d2, by = c(segments_var))
    } else {
      d1 = xdata %>%
        select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == learners_level) %>% select(!all_of(groups_var)) %>%
        group_by_at(.vars = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit()) %>% nest() %>% rename(data.x = data)
      d2 = xdata %>% select(all_of(c(groups_var, segments_var, items_var, values_var) %>% na.omit() %>% unique())) %>%
        filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
        nest(data.y = everything()) %>%
        mutate(data.y = map(data.y, ~ select(.x, -all_of(items_var))))
      mah_segment = bind_cols(d1, d2)
    }
    
    mah_segment = mah_segment %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data."))
    rm(d1, d2)
    
    # Incorporate Mahalanobis distance into the original data frame
    
    mah_segment = mah_segment %>% rename(mah_dist = distance)
    xdata = xdata %>% left_join(mah_segment, by = c(subjects_var, segments_var, unique_items_var, tests_var) %>% na.omit() %>% unique()); rm(mah_segment)
    
  }
  
  # Pillai within
  
  if (inter_group == FALSE) {
    
    pillai = xdata %>%
      rename("segments_var" = all_of(segments_var)) %>%
      select(all_of(c(subjects_var, tests_var, values_var) %>% na.omit()), segments_var) %>% group_by_at(.vars = c(tests_var, subjects_var) %>% na.omit()) %>% nest() %>% 
      mutate(pillai = map(data, ~manova(formula = (.x %>% select(all_of(values_var)) %>% as.matrix()) ~ segments_var, data = .x) %>% 
                            summary() %>% pluck("stats") %>% as.data.frame() %>% rownames_to_column("factor") %>% 
                            filter(factor == "segments_var") %>% select(Pillai) %>% as.numeric()) %>% as.numeric()) %>%
      ungroup() %>% select(-c(data)) %>% rename(pil_dist = pillai)
    xdata = xdata %>% left_join(pillai, by = c(subjects_var, tests_var) %>% na.omit()); rm(pillai)
    
  }
  
  # Pillai between
  
  if (inter_group == TRUE) {
    
    vars_to_group = c(subjects_var, tests_var, segments_var) %>% na.omit()
    if (!is.na(segments_var)) {
      d1 = xdata %>%
        filter(!!as.name(groups_var) == learners_level) %>%
        rename("groups_var" = all_of(groups_var)) %>%
        select(all_of(vars_to_group), all_of(values_var)) %>% group_by_at(.vars = vars_to_group) %>% nest()
      d2 = xdata %>%
        filter(!!as.name(groups_var) == natives_level) %>%
        rename("groups_var" = all_of(groups_var)) %>%
        select(all_of(c(values_var, segments_var))) %>% group_by_at(.vars = c(segments_var)) %>% nest()
      pillai = left_join(d1, d2, by = c(segments_var))
    } else {
      d1 = xdata %>%
        filter(!!as.name(groups_var) == learners_level) %>%
        rename("groups_var" = all_of(groups_var)) %>%
        select(all_of(vars_to_group))
      if (length(vars_to_group) < ncol(d1)) {
        d1 = d1 %>% group_by_at(.vars = vars_to_group) %>% nest() %>% rename(data.x = data)
      } else {
        d1 = d1 %>% group_by_at(.vars = subjects_var) %>% nest() %>% rename(data.x = data)
      }
      d2 = xdata %>%
        filter(!!as.name(groups_var) == natives_level) %>%
        rename("groups_var" = all_of(groups_var)) %>%
        select(all_of(c(values_var))) %>% 
        nest(data.y = everything())
      pillai = bind_cols(d1, d2)
    }
    rm(vars_to_group)
    
    pillai = pillai %>% 
      rename(learners = data.x, natives = data.y) %>% pivot_longer(cols = c(learners, natives), names_to = "groups_var", values_to = "data") %>% unnest(data) %>% 
      group_by_at(.vars = c(subjects_var, tests_var, segments_var) %>% na.omit()) %>% nest() %>%
      mutate(pillai = map(data, ~manova(formula = (.x %>% select(all_of(values_var)) %>% as.matrix()) ~ groups_var, data = .x) %>% 
                            summary() %>% pluck("stats") %>% as.data.frame() %>% rownames_to_column("factor") %>% 
                            filter(factor == "groups_var") %>% select(Pillai) %>% as.numeric()) %>% as.numeric()) %>%
      ungroup() %>% select(-c(data))
    
    xdata = xdata %>% left_join(pillai, by = c(subjects_var, segments_var, tests_var) %>% na.omit()) %>% rename(pil_dist = pillai); rm(d1, d2, pillai)
    
  }
  
  # End
  
  if (inter_group == FALSE) {xdata = xdata %>% rename(dist_wit_euc = euc_dist, dist_wit_mah = mah_dist, dist_wit_pil = pil_dist)}
  if (inter_group == TRUE) {xdata = xdata %>% rename(dist_bet_euc = euc_dist, dist_bet_mah = mah_dist, dist_bet_pil = pil_dist)}
  
  return(xdata)
  
}
