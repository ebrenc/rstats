
voweldist = function(x, inter_group = FALSE, 
                     segments_var,  # By now, the script works if this variable has only 2 levels
                     tests_var,  # By now, the script works if this variable has only 2 levels
                     subjects_var,
                     items_var,
                     unique_items_var, # Combined with subjects_var, these variable(s) should allow us to uniquely identify each linguistic production present in the data frame
                     values_var, # You can put 2 or more variables here
                     groups_var = NULL,
                     learners_level = NULL,
                     natives_level = NULL) {
  
  require(tidyverse)
  options(scipen=999)
  
  df = x
  
  # Preliminary tasks
  
  vars_to_char = c(segments_var, tests_var, groups_var, subjects_var, items_var, unique_items_var) %>% unique()
  df = df %>% 
    arrange(across(all_of(vars_to_char))) %>% 
    mutate(across(all_of(vars_to_char), as.character))
  rm(vars_to_char)
  
  segments_var_levels = df %>% select(all_of(segments_var)) %>% unique() %>% unlist() %>% as.character()
  tests_var_levels = df %>% select(all_of(tests_var)) %>% unique() %>% drop_na() %>% unlist() %>% as.character()
  n_values_var = values_var %>% length()
  n_clusters_factors = length(segments_var_levels)*length(tests_var_levels)
  
  # Filter out subjects with less than 3 cases in each combination of segments_var * tests_var
  
  df_target_subjects = df %>% select(all_of(c(groups_var, subjects_var, tests_var, segments_var)))
  if (inter_group == TRUE) {df_target_subjects = df_target_subjects %>% filter(!!as.name(groups_var) == learners_level)}
  subjects_out = c(
    df_target_subjects %>% distinct() %>% group_by_at(.vars = c(subjects_var)) %>% count() %>% filter(n < n_clusters_factors) %>% pluck(subjects_var),
    df_target_subjects %>% group_by_at(.vars = c(subjects_var, tests_var, segments_var)) %>% count() %>% filter(n <= n_values_var) %>% pluck(subjects_var)
  )
  if (length(subjects_out) > 0) {
    print(c("Subjects out: ", subjects_out))
    `%notin%` = Negate(`%in%`)
    # df = df %>% filter(eval(parse(text = subjects_var)) %notin% subjects_out)
    df = df %>% filter(!!as.name(subjects_var) %notin% subjects_out)
    rm(`%notin%`)
  }
  rm(df_target_subjects, subjects_out)
  
  # Euclidean within

  if (inter_group == FALSE) {

    # Euclidean distances 1 within 2

    d1 = df %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var))

    d2 = df %>% select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var)) %>% summarise_if(is.numeric, mean, na.rm = T)

    euc_segment_1 = left_join(d1, d2, by = c(subjects_var, tests_var))

    var_matrix = colnames(euc_segment_1) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment_1[all_of(var_matrix[i,1])] - euc_segment_1[all_of(var_matrix[i,2])])^2)
    }
    euc_segment_1 = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var)) %>%
      mutate(!!segments_var := segments_var_levels[1])
    rm(var_matrix, i, vector_to_sum, d1, d2)

    # Euclidean distances 2 within 1

    d1 = df %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var))

    d2 = df %>% select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var)) %>% summarise_if(is.numeric, mean, na.rm = T)

    euc_segment_2 = left_join(d1, d2, by = c(subjects_var, tests_var))

    var_matrix = colnames(euc_segment_2) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment_2[all_of(var_matrix[i,1])] - euc_segment_2[all_of(var_matrix[i,2])])^2)
    }
    euc_segment_2 = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var)) %>%
      mutate(!!segments_var := segments_var_levels[2])
    rm(var_matrix, i, vector_to_sum, d1, d2)

    # Incorporate Euclidean distances into the original data frame

    euc_distances = bind_rows(euc_segment_1, euc_segment_2) %>% rename(euc_dist = distance); rm(euc_segment_1, euc_segment_2)
    df = df %>% left_join(euc_distances, by = c(subjects_var, segments_var, unique_items_var, tests_var)); rm(euc_distances)

  }

  # Euclidean between

  if (inter_group == TRUE) {

    # Euclidean distances between

    d1 = df %>%
      select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(!!as.name(groups_var) == learners_level) %>% select(!all_of(groups_var))

    d2 = df %>% select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
      group_by_at(.vars = c(segments_var)) %>% summarise_if(is.numeric, mean, na.rm = T)

    euc_segment = left_join(d1, d2, by = c(segments_var))

    var_matrix = colnames(euc_segment) %>% tail(n_values_var*2) %>% matrix(nrow = n_values_var)
    vector_to_sum = c()
    for (i in 1:n_values_var) {
      vector_to_sum = c(vector_to_sum, (euc_segment[all_of(var_matrix[i,1])] - euc_segment[all_of(var_matrix[i,2])])^2)
    }
    euc_segment = vector_to_sum %>% as.data.frame() %>% rowSums() %>% sqrt() %>% as.data.frame() %>% rename("distance" = ".") %>%
      bind_cols(d1, .) %>% select(!all_of(values_var))
    rm(var_matrix, i, vector_to_sum, d1, d2)

    # Incorporate Euclidean distances into the original data frame

    euc_segment = euc_segment %>% rename(euc_dist = distance)
    df = df %>% left_join(euc_segment, by = c(subjects_var, segments_var, unique_items_var, tests_var)); rm(euc_segment)

  }

  # Mahalanobis within

  if (inter_group == FALSE) {

    # Mahalanobis distances 1 within 2

    d1 = df %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, unique_items_var, tests_var)) %>% nest()

    d2 = df %>% select(all_of(c(subjects_var, segments_var, items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var)) %>% nest() %>%
      mutate(data = map(data, ~ select(.x, -all_of(items_var))))

    mah_segment_1 = left_join(d1, d2, by = c(subjects_var, tests_var)) %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data.")) %>%
      mutate(!!segments_var := segments_var_levels[1])
    rm(d1, d2)

    # Mahalanobis distances 2 within 1

    d1 = df %>%
      select(all_of(c(subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[2]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, unique_items_var, tests_var)) %>% nest()

    d2 = df %>% select(all_of(c(subjects_var, segments_var, items_var, tests_var, values_var))) %>%
      filter(eval(parse(text = segments_var)) == segments_var_levels[1]) %>% select(!all_of(segments_var)) %>%
      group_by_at(.vars = c(subjects_var, tests_var)) %>% nest() %>%
      mutate(data = map(data, ~ select(.x, -all_of(items_var))))

    mah_segment_2 = left_join(d1, d2, by = c(subjects_var, tests_var)) %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data.")) %>%
      mutate(!!segments_var := segments_var_levels[2])
    rm(d1, d2)

    # Incorporate Mahalanobis distances into the original data frame

    mah_distances = bind_rows(mah_segment_1, mah_segment_2) %>% rename(mah_dist = distance); rm(mah_segment_1, mah_segment_2)
    df = df %>% left_join(mah_distances, by = c(subjects_var, segments_var, unique_items_var, tests_var)); rm(mah_distances)

  }

  # Mahalanobis between

  if (inter_group == TRUE) {

    # Mahalanobis distances between

    d1 = df %>%
      select(all_of(c(groups_var, subjects_var, segments_var, unique_items_var, tests_var, values_var))) %>%
      filter(!!as.name(groups_var) == learners_level) %>% select(!all_of(groups_var)) %>%
      group_by_at(.vars = c(subjects_var, segments_var, unique_items_var, tests_var)) %>% nest()

    d2 = df %>% select(all_of(c(groups_var, segments_var, items_var, values_var))) %>%
      filter(!!as.name(groups_var) == natives_level) %>% select(!all_of(groups_var)) %>%
      group_by_at(.vars = c(segments_var)) %>% nest() %>%
      mutate(data = map(data, ~ select(.x, -all_of(items_var))))

    mah_segment = left_join(d1, d2, by = c(segments_var)) %>% rowwise() %>% filter(!is.null(data.y)) %>%
      mutate(distance = mahalanobis(data.x, MASS::cov.trob(data.y) %>% pluck("center"), MASS::cov.trob(data.y) %>% pluck("cov"))) %>%
      ungroup() %>% select(-contains("data."))
    rm(d1, d2)

    # Incorporate Mahalanobis distance into the original data frame

    mah_segment = mah_segment %>% rename(mah_dist = distance)
    df = df %>% left_join(mah_segment, by = c(subjects_var, segments_var, unique_items_var, tests_var)); rm(mah_segment)

  }
  
  # Pillai within
  
  if (inter_group == FALSE) {
    
    pillai = df %>%
      rename("segments_var" = all_of(segments_var)) %>%
      select(all_of(c(subjects_var, tests_var, values_var)), segments_var) %>% group_by_at(.vars = c(tests_var, subjects_var)) %>% nest() %>% 
      mutate(pillai = map(data, ~manova(formula = (.x %>% select(all_of(values_var)) %>% as.matrix()) ~ segments_var, data = .x) %>% 
                            summary() %>% pluck("stats") %>% as.data.frame() %>% rownames_to_column("factor") %>% 
                            filter(factor == "segments_var") %>% select(Pillai) %>% as.numeric()) %>% as.numeric()) %>%
      ungroup() %>% select(-c(data)) %>% rename(pil_dist = pillai)
    df = df %>% left_join(pillai, by = c(subjects_var, tests_var)); rm(pillai)
    
  }
  
  # Pillai between
  
  if (inter_group == TRUE) {
    
    d1 = df %>%
      filter(!!as.name(groups_var) == learners_level) %>%
      rename("groups_var" = all_of(groups_var)) %>%
      select(all_of(c(subjects_var, tests_var, values_var, segments_var))) %>% group_by_at(.vars = c(tests_var, subjects_var, segments_var)) %>% nest()
    
    d2 = df %>%
      filter(!!as.name(groups_var) == natives_level) %>%
      rename("groups_var" = all_of(groups_var)) %>%
      select(all_of(c(values_var, segments_var))) %>% group_by_at(.vars = c(segments_var)) %>% nest()
    
    pillai = left_join(d1, d2, by = c(segments_var)) %>% 
      rename(learners = data.x, natives = data.y) %>% pivot_longer(cols = c(learners, natives), names_to = "groups_var", values_to = "data") %>% unnest(data) %>% 
      group_by_at(.vars = c(subjects_var, tests_var, segments_var)) %>% nest() %>%
      mutate(pillai = map(data, ~manova(formula = (.x %>% select(all_of(values_var)) %>% as.matrix()) ~ groups_var, data = .x) %>% 
                            summary() %>% pluck("stats") %>% as.data.frame() %>% rownames_to_column("factor") %>% 
                            filter(factor == "groups_var") %>% select(Pillai) %>% as.numeric()) %>% as.numeric()) %>%
      ungroup() %>% select(-c(data))
    
    df = df %>% left_join(pillai, by = c(subjects_var, segments_var, tests_var)) %>% rename(pil_dist = pillai); rm(d1, d2, pillai)
    
  }
  
  # End
  
  if (inter_group == FALSE) {df = df %>% rename(dist_wit_euc = euc_dist, dist_wit_mah = mah_dist, dist_wit_pil = pil_dist)}
  if (inter_group == TRUE) {df = df %>% rename(dist_bet_euc = euc_dist, dist_bet_mah = mah_dist, dist_bet_pil = pil_dist)}
  
  return(df)
  
}
