
glmmTable = function(model, path = NA, title = "Model", extract = FALSE) {
  
  require(tidyverse); require(car)
  
  model.anova = car::Anova(model, type = c("II", "III", 2, 3), test.statistic = c("Chisq", "F")) %>%
    tibble::rownames_to_column("term") %>% 
    dplyr::rename("p.Chisq" = "Pr(>Chisq)") %>% 
    select(term, Chisq, Df, p.Chisq)
  
  require(emmeans)
  emm_options(lmerTest.limit = 10000, disable.pbkrtest = T, lmer.df = "satterthwaite", msg.interaction = F)
  
  if (class(model) == "glmmTMB") {
    terms = attr(model$modelInfo$reTrms$cond$terms$fixed, "dataClasses") %>% enframe() %>% tail(-1)
  } else {
    terms = model %>% terms() %>% attr("term.labels") %>% str_subset("^[^:]+$") %>%
      sapply(., function(effect) {
        class(model@frame[[effect]])
      }) %>% enframe()
  }
  
  factor_names = terms %>% pluck("name")
  categorical_factor_names = terms %>% filter(value %in% c("factor", "character")) %>% pluck("name")
  non_categorical_factor_names = factor_names[!factor_names %in% categorical_factor_names]
  rm(terms)
  
  n = length(factor_names)
  for (f in 1:n) {assign(str_c("f", f), factor_names[f])}; rm(f)
  
  combinations = c()
  for (r in 1:n) {
    newset = gtools::permutations(n=n,r=r,v=as.character(1:n)) %>% data.frame() %>% rename_with(~str_replace_all(., "X", "V"), .cols = everything())
    if(ncol(newset) > 2){
      for (v in 3:ncol(newset)) {
        maintain = newset[,v] > newset[,v-1]
        newset = newset %>% filter(maintain)
        rm(maintain)
      }
      rm(v)
    }
    newset = newset %>% unite(col="col",sep="") %>% pluck("col")
    combinations = c(combinations, newset)
    rm(newset)
  }
  rm(r, n)
  
  for (c in combinations) {
    
    for (cnumber in 1:str_count(c)) {
      assign(x = str_c("number", cnumber), value = str_c("f", str_sub(c, cnumber, cnumber)) %>% parse(text = .) %>% eval())
    }; rm(cnumber)
    
    if (str_c("emmeans(model, ~ ", number1, ")") %>% parse(text = .) %>% eval() %>% data.frame() %>% as_tibble() %>% select(all_of(number1)) %>% dplyr::count() > 1) {
      
      lineofresults = str_c("emmeans(model, ~ ", number1)
      if (str_count(c) == 1) {lineofresults = str_c(lineofresults, ")")}
      if (str_count(c) == 2) {lineofresults = str_c(lineofresults, " | ", number2, ")")}
      if (str_count(c) > 2) {
        numbers = c()
        for (number in 3:str_count(c)) { numbers = c(numbers, " * ", str_c("number", number) %>% parse(text = .) %>% eval()) }  
        numbers = numbers %>% paste(collapse="")
        lineofresults = str_c(lineofresults, " | ", number2, numbers, ")")
        rm(numbers)
      }
      
      term = lineofresults %>% str_replace_all(fixed(" | "), ":") %>% str_replace_all(fixed(" * "), ":") %>% str_remove_all(fixed("emmeans(model, ~ ")) %>% str_remove_all(fixed(")"))
      result = lineofresults %>% parse(text = .) %>% eval() %>% 
        pairs(adjust = 'bonferroni') %>% summary() %>% as.data.frame() %>% mutate(contrastfield = number1) %>% 
        mutate(cohensd = reduce(select(summary(eff_size(eval(parse(text = lineofresults)), sigma = sigma(model), edf = df.residual(model))), effect.size), .f=full_join)) %>% 
        mutate(combinations = c) %>% mutate(term = term)
      rm(lineofresults, term)
      
    } else { result = tibble(
      contrast = character(),
      estimate = numeric(),
      SE = numeric(),
      df = numeric(),
      t.ratio = numeric(),
      p.value = numeric(),
      contrastfield = character(),
      cohensd = numeric(),
      combinations = character(),
      term = character())}
    
    assign(str_c("pair", c), result)
    
  }
  
  for (number in 1:str_count(c)) {rm(list = str_c("number", str_sub(c, number, number)))}
  rm(c, result, number)
  
  model.emmeans = tibble()
  for (pair in 1:length(combinations)) {
    model.emmeans = bind_rows(model.emmeans, str_c("pair", eval(parse(text = combinations[pair]))) %>% parse(text = .) %>% eval())
  }
  
  str_arrange = function(x){
    x %>%
      stringr::str_split("") %>% # Split string into letters
      purrr::map(~sort(.) %>% paste(collapse = "")) %>% # Sort and re-combine
      as_vector() # Convert list into vector
  }
  
  model.emmeans = model.emmeans %>% 
    mutate(order1 = str_count(combinations)) %>% 
    mutate(order2 = str_arrange(combinations)) %>% 
    arrange(order1, order2)
  
  rm(list = str_c("pair", eval(parse(text = "combinations"))))
  
  term = model.emmeans %>% select(combinations, order2, term) %>% distinct()
  term2 = term %>% filter(combinations == order2) %>% select(-combinations)
  term = term %>% select(-term) %>% left_join(term2, by = "order2") %>% select(-order2)
  rm(term2)
  
  model.emmeans = model.emmeans %>% select(-term) %>% left_join(term, by = "combinations") %>% select(-combinations) %>% relocate(term, .before = contrast)
  model.emmeans = model.emmeans %>% select(-c(order1, order2))
  
  if (length(non_categorical_factor_names) != 0) {
    model.emmeans = model.emmeans %>% filter(!str_detect(term, non_categorical_factor_names)) %>% janitor::remove_empty("cols")
  }
  rm(non_categorical_factor_names)
  
  if (nrow(model.emmeans) != 0) {
    
    model.emmeans = model.emmeans %>% 
      mutate(df = round(df, 2)) %>%
      separate(contrast, into = c("cfield.1","cfield.2"), sep=" - ") %>%
      mutate(contrast = case_when(p.value < .05 & estimate > 0 ~ paste(cfield.1, " > ", cfield.2),
                                  p.value < .05 & estimate < 0 ~ paste(cfield.2, " > ", cfield.1),
                                  p.value >= .05 & p.value < .1 & estimate > 0 ~ paste(cfield.1, " *> ", cfield.2),
                                  p.value >= .05 & p.value < .1 & estimate < 0 ~ paste(cfield.2, " *> ", cfield.1),
                                  p.value >= .1 & estimate > 0 ~ paste(cfield.1, " = ", cfield.2),
                                  p.value >= .1 & estimate < 0 ~ paste(cfield.2, " = ", cfield.1),
                                  p.value >= .1 & estimate < 0 ~ paste(cfield.1, " == ", cfield.2),
                                  TRUE ~ 'not computable'))
    
    if (length(unique(categorical_factor_names)) == 1) {
      model.emmeans = model.emmeans %>%select(term, contrastfield, contrast, cohensd, p.value)
    } else {
      model.emmeans = model.emmeans %>%
        unite(col = "insidelevel", all_of(categorical_factor_names), sep = " ", remove = T, na.rm = T) %>%
        select(term, contrastfield, insidelevel, contrast, cohensd, p.value)
    }
    
    results.table = model.emmeans %>% 
      tidyr::nest(tests = -c(term, contrastfield)) %>% 
      pivot_wider(names_from = contrastfield, values_from = tests)
    
    for (cat_colname in colnames(results.table)[-1]) {
      results.table = results.table %>% rename_with(~str_c("f", which(categorical_factor_names == cat_colname)), all_of(cat_colname))
    }
    rm(cat_colname)
    
    for (catfactor in 1:length(categorical_factor_names)) {
      assign(str_c("results.table", catfactor), 
             results.table %>% select(term, str_c("f", all_of(catfactor))) %>% 
               janitor::remove_empty("cols") %>% unnest(cols = str_c("f", catfactor), names_sep = "."))
    }
    rm(catfactor)
    
    results.table = tibble()
    for (pair in 1:length(categorical_factor_names)) {
      results.table = bind_rows(results.table, str_c("results.table", pair) %>% parse(text = .) %>% eval())
    }
    rm(pair)
    rm(list = str_c("results.table", 1:length(categorical_factor_names)))
    
    results.table = results.table %>% 
      full_join(model.anova, by = "term") %>% 
      relocate(Chisq, Df, p.Chisq, .after = term)
    
  } else if (nrow(model.emmeans) == 0) {results.table = model.anova}
  
  rm(model.anova, model.emmeans)
  
  term = term %>% mutate(order1 = str_count(combinations))
  term = term %>% mutate(order2 = str_arrange(combinations))
  results.table = results.table %>% left_join(term %>% select(-combinations), by = "term", relationship = "many-to-many") %>% arrange(order1, order2) %>% select(-c(order1, order2))
  results.table = results.table %>% distinct()
  rm(term, str_arrange, combinations)
  results.table = results.table %>% filter(!is.na(Chisq)) 
  
  # Repair rounding and effect sizes
  
  results.table = results.table %>%
    mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>%
    mutate(across(contains("cohensd"), abs))
  
  if (!is.na(path)) {
    # Do the flextable
    
    require(flextable)
    
    set_flextable_defaults(font.family = "Arial", font.size = 9, digits = 3)
    
    number_of_columns = ncol(results.table)
    number_of_fs = round((number_of_columns-4)/4)
    
    results = results.table %>% flextable() %>% theme_vanilla()
    
    results = results %>% 
      colformat_num(j = "Df", digits = 0) %>%
      merge_v(j = "term", target = c(1:4)) %>% 
      bg(~ (p.Chisq < .1 & p.Chisq >= .05), ~ term + p.Chisq, bg = "#d7eef3") %>%
      bg(~ p.Chisq < .05, ~ term + p.Chisq, bg = "#90ccde")
    if(number_of_columns > 4) {results = results %>% merge_v(j = c(5:number_of_columns))}
    
    if (number_of_columns >= 7) {
      
      for (f in 1:number_of_fs) {
        f.p.value = str_c("f", f, ".p.value")
        results = results %>% 
          bg(i = (results$body$dataset[[f.p.value]] < .1 & results$body$dataset[[f.p.value]] >= .05), j = which(results$body$col_keys == f.p.value), bg = "#d7eef3") %>% 
          bg(i = results$body$dataset[[f.p.value]] < .05, j = which(results$body$col_keys == f.p.value), bg = "#90ccde")
      }
      rm(f, f.p.value)
    }
    
    if (number_of_fs == 1) {
      values = c("Model results", "", "", "", categorical_factor_names[1], "", "")
    }
    
    if (number_of_fs > 1) {
      values = c("Model results", "", "", "")
      for (f in 1:length(categorical_factor_names)) {
        values = c(values, categorical_factor_names[f], "", "", "")
      }
      rm(f)
    }
    
    ###
    
    for (alias in 1:number_of_columns) {
      content = results$header$content$data[[alias]]$txt
      if (content == "term") {results$header$content$data[[alias]]$txt = "Term"}
      else if (content == "Chisq") {results$header$content$data[[alias]]$txt = "χ²"}
      else if (content == "Df") {results$header$content$data[[alias]]$txt = "df"}
      else if (content == "p.Chisq") {results$header$content$data[[alias]]$txt = "p"}
      else if (grepl("insidelevel", content, fixed = T)) {results$header$content$data[[alias]]$txt = "Levels"}
      else if (grepl("contrast", content, fixed = T)) {results$header$content$data[[alias]]$txt = "Contrast"}
      else if (grepl("cohensd", content, fixed = T)) {results$header$content$data[[alias]]$txt = "Coh.d"}
      else if (grepl("p.value", content, fixed = T)) {results$header$content$data[[alias]]$txt = "Sig"}
      rm(content)
    }; rm(alias, number_of_columns)
    
    results = results %>% add_header_row(values = values, top = TRUE); rm(values)
    results = results %>% merge_at(i = 1, j = 1:4, part = "header")
    if (number_of_fs == 1) {results = results %>% merge_at(i = 1, j = 5:7, part = "header")}
    if (number_of_fs > 1) {for (f in 1:length(categorical_factor_names)) {results = results %>% merge_at(i = 1, j = (1:4)+4*f, part = "header")}; rm(f)}
    
    results = results %>% set_formatter_type(fmt_double="%.01f") %>% padding(padding = 4, part = "all") %>% autofit()
    
    rm(list = str_c("f", 1:number_of_fs))
    rm(categorical_factor_names, factor_names, number_of_fs)
    
    results %>% save_as_html(path = path, title = title)
  }
  
  if (extract == TRUE) {return(results.table)}
  
}
