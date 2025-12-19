
glmmTable = function(model, path = NA, title = "Model", extract = FALSE) {
  
  library(tidyverse); library(car); library(janitor)
  library(emmeans)
  emm_options(lmerTest.limit = 10000, disable.pbkrtest = TRUE,
              lmer.df = "satterthwaite", msg.interaction = FALSE)
  
  # -------------------------
  # 0) Helpers (robustness)
  # -------------------------
  
  is_glmmTMB <- inherits(model, "glmmTMB")
  is_lm      <- inherits(model, "lm") && !inherits(model, "glm")
  is_glm     <- inherits(model, "glm") && !inherits(model, "glmerMod") && !is_glmmTMB
  is_lmer    <- inherits(model, "lmerMod")
  is_glmer   <- inherits(model, "glmerMod")
  
  # robust model frame (avoid any @ usage)
  model_frame <- tryCatch(
    model.frame(model),
    error = function(e) {
      if (is_glmmTMB && !is.null(model$frame)) return(model$frame)
      if (!is.null(model$model)) return(model$model)
      NULL
    }
  )
  
  # robust residual df (some models don't have it cleanly)
  safe_df_resid <- function(m) {
    out <- tryCatch(df.residual(m), error = function(e) NA_real_)
    as.numeric(out)
  }
  
  # robust sigma (for eff_size); may be undefined for some GLM/GLMM
  safe_sigma <- function(m) {
    out <- tryCatch(sigma(m), error = function(e) NA_real_)
    as.numeric(out)
  }
  
  safe_eff_size <- function(emm_obj, m) {
    out <- tryCatch({
      es <- eff_size(emm_obj, sigma = safe_sigma(m), edf = safe_df_resid(m))
      s  <- summary(es)
      if ("effect.size" %in% names(s)) return(as.numeric(s$effect.size[1]))
      if ("effect.size" %in% names(es)) return(as.numeric(es$effect.size[1]))
      NA_real_
    }, error = function(e) NA_real_)
    as.numeric(out)
  }
  
  # robust car::Anova call: try several settings
  safe_anova <- function(m) {
    # Choose sensible default:
    # - lm: F
    # - glm/glmer/glmmTMB: Chisq often
    # - lmer: F (or Chisq sometimes supported depending on method)
    test_try <- list(
      if (is_lm) "F" else NULL,
      if (is_glm || is_glmer || is_glmmTMB) "Chisq" else NULL,
      c("Chisq", "F"),
      "F",
      "Chisq"
    ) %>% compact()
    
    for (ts in test_try) {
      out <- tryCatch({
        car::Anova(m, type = c("II", "III", 2, 3), test.statistic = ts) %>%
          data.frame() %>%
          janitor::clean_names() %>%
          tibble::rownames_to_column("term")
      }, error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    
    # If everything fails, return empty tibble with required cols
    tibble(term = character(), chisq = numeric(), df = numeric(), p_chisq = numeric())
  }
  
  model.anova <- safe_anova(model)
  
  # pick whatever p-value column exists and rename to p_chisq
  p_candidates <- c("pr_chisq", "pr_f", "pr_gt_chisq", "pr_gt_f", "p")
  p_col <- p_candidates[p_candidates %in% names(model.anova)][1]
  if (!is.na(p_col)) {
    model.anova <- model.anova %>% dplyr::rename(p_chisq = !!rlang::sym(p_col))
  } else {
    model.anova <- model.anova %>% mutate(p_chisq = NA_real_)
  }
  
  # also standardize chisq/df column names if needed
  # (car::Anova sometimes gives "chisq" or "f" depending on test)
  if (!"chisq" %in% names(model.anova)) {
    f_candidates <- c("f", "f_value", "statistic")
    f_col <- f_candidates[f_candidates %in% names(model.anova)][1]
    if (!is.na(f_col)) model.anova <- model.anova %>% rename(chisq = !!rlang::sym(f_col))
  }
  if (!"df" %in% names(model.anova)) {
    df_candidates <- c("df", "num_df", "den_df")
    df_col <- df_candidates[df_candidates %in% names(model.anova)][1]
    if (!is.na(df_col)) model.anova <- model.anova %>% rename(df = !!rlang::sym(df_col))
    if (!"df" %in% names(model.anova)) model.anova <- model.anova %>% mutate(df = NA_real_)
  }
  if (!"chisq" %in% names(model.anova)) model.anova <- model.anova %>% mutate(chisq = NA_real_)
  
  # -------------------------
  # 1) Terms + data classes (compatible across model types)
  # -------------------------
  
  # Preferred: use dataClasses from terms()
  terms_tbl <- tryCatch({
    dc <- attr(terms(model), "dataClasses")
    dc <- dc[names(dc) != "(Intercept)"]
    tibble(name = names(dc), value = unname(dc))
  }, error = function(e) NULL)
  
  # glmmTMB special-case fallback (your original approach)
  if (is.null(terms_tbl) && is_glmmTMB) {
    terms_tbl <- tryCatch({
      attr(model$modelInfo$reTrms$cond$terms$fixed, "dataClasses") %>%
        enframe(name = "name", value = "value") %>%
        tail(-1)
    }, error = function(e) NULL)
  }
  
  # Final fallback: derive from model_frame
  if (is.null(terms_tbl)) {
    if (!is.null(model_frame)) {
      terms_labels <- tryCatch(attr(terms(model), "term.labels"), error = function(e) character())
      main_terms   <- stringr::str_subset(terms_labels, "^[^:]+$")
      
      terms_tbl <- tibble(name = main_terms) %>%
        mutate(value = purrr::map_chr(name, \(eff) {
          if (!is.null(model_frame) && eff %in% names(model_frame)) class(model_frame[[eff]])[1] else NA_character_
        }))
    } else {
      terms_tbl <- tibble(name = character(), value = character())
    }
  }
  
  terms <- terms_tbl
  rm(terms_tbl)
  
  factor_names <- terms %>% pluck("name")
  
  if ("value" %in% names(terms)) {
    categorical_factor_names <-
      terms %>%
      filter(value %in% c("factor", "character")) %>%
      pluck("name")
  } else {
    # fallback: assume all factors are categorical
    categorical_factor_names <- factor_names
  }
  
  non_categorical_factor_names <-
    factor_names[!factor_names %in% categorical_factor_names]
  rm(terms)
  
  n = length(factor_names)
  for (f in 1:n) {assign(str_c("f", f), factor_names[f])}; rm(f)
  
  # -------------------------
  # 2) Build combinations (unchanged)
  # -------------------------
  
  combinations = c()
  for (r in 1:n) {
    newset = gtools::permutations(n=n,r=r,v=as.character(1:n)) %>% data.frame() %>%
      rename_with(~str_replace_all(., "X", "V"), .cols = everything())
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
  
  # -------------------------
  # 3) Emmeans + pairs (robust)
  # -------------------------
  
  for (c in combinations) {
    
    for (cnumber in 1:str_count(c)) {
      assign(x = str_c("number", cnumber),
             value = str_c("f", str_sub(c, cnumber, cnumber)) %>% parse(text = .) %>% eval())
    }; rm(cnumber)
    
    # robust check: only proceed if emmeans exists and has >1 level
    ok_to_compute <- tryCatch({
      emm1 <- emmeans(model, as.formula(paste0("~ ", number1)))
      df1  <- as.data.frame(emm1) %>% as_tibble()
      if (!(number1 %in% names(df1))) return(FALSE)
      dplyr::count(df1, .data[[number1]]) %>% nrow() > 1
    }, error = function(e) FALSE)
    
    if (ok_to_compute) {
      
      lineofresults = str_c("emmeans(model, ~ ", number1)
      if (str_count(c) == 1) {lineofresults = str_c(lineofresults, ")")}
      if (str_count(c) == 2) {lineofresults = str_c(lineofresults, " | ", number2, ")")}
      if (str_count(c) > 2) {
        numbers = c()
        for (number in 3:str_count(c)) {
          numbers = c(numbers, " * ", str_c("number", number) %>% parse(text = .) %>% eval())
        }
        numbers = numbers %>% paste(collapse="")
        lineofresults = str_c(lineofresults, " | ", number2, numbers, ")")
        rm(numbers)
      }
      
      term = lineofresults %>%
        str_replace_all(fixed(" | "), ":") %>%
        str_replace_all(fixed(" * "), ":") %>%
        str_remove_all(fixed("emmeans(model, ~ ")) %>%
        str_remove_all(fixed(")"))
      
      emm_obj <- tryCatch(lineofresults %>% parse(text = .) %>% eval(),
                          error = function(e) NULL)
      
      if (!is.null(emm_obj)) {
        pair_df <- tryCatch({
          pairs(emm_obj, adjust = 'bonferroni') %>%
            summary() %>% as.data.frame()
        }, error = function(e) data.frame())
        
        result <- pair_df %>%
          as_tibble() %>%
          mutate(contrastfield = number1) %>%
          mutate(cohensd = safe_eff_size(emm_obj, model)) %>%
          mutate(combinations = c) %>%
          mutate(term = term)
        
      } else {
        result = tibble(
          contrast = character(),
          estimate = numeric(),
          SE = numeric(),
          df = numeric(),
          t.ratio = numeric(),
          p.value = numeric(),
          contrastfield = character(),
          cohensd = numeric(),
          combinations = character(),
          term = character()
        )
      }
      
      rm(lineofresults, term)
      
    } else {
      result = tibble(
        contrast = character(),
        estimate = numeric(),
        SE = numeric(),
        df = numeric(),
        t.ratio = numeric(),
        p.value = numeric(),
        contrastfield = character(),
        cohensd = numeric(),
        combinations = character(),
        term = character()
      )
    }
    
    if ("contrast1" %in% colnames(result) && !"contrast" %in% colnames(result)) {
      result <- result %>% dplyr::rename(contrast = contrast1)
    }
    
    assign(str_c("pair", c), result)
  }
  
  for (number in 1:str_count(c)) {rm(list = str_c("number", str_sub(c, number, number)))}
  rm(c, result, number)
  
  model.emmeans = tibble()
  for (pair in 1:length(combinations)) {
    model.emmeans = bind_rows(model.emmeans,
                              str_c("pair", eval(parse(text = combinations[pair]))) %>%
                                parse(text = .) %>% eval())
  }
  
  str_arrange = function(x){
    x %>%
      stringr::str_split("") %>%
      purrr::map(~sort(.) %>% paste(collapse = "")) %>%
      as_vector()
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
  
  model.emmeans = model.emmeans %>%
    select(-term) %>%
    left_join(term, by = "combinations") %>%
    select(-combinations) %>%
    relocate(term, .before = contrast)
  
  model.emmeans = model.emmeans %>% select(-c(order1, order2))
  
  if (length(non_categorical_factor_names) != 0) {
    pattern_nc <- paste0("\\b(", paste(non_categorical_factor_names, collapse = "|"), ")\\b")
    model.emmeans <- model.emmeans %>%
      filter(!str_detect(term, pattern_nc)) %>%
      janitor::remove_empty("cols")
    rm(pattern_nc)
  }
  rm(non_categorical_factor_names)
  
  if (nrow(model.emmeans) != 0) {
    
    model.emmeans = model.emmeans %>%
      mutate(df = round(df, 2)) %>%
      separate(contrast, into = c("cfield.1","cfield.2"), sep=" - ") %>%
      mutate(contrast = case_when(
        p.value < .05 & estimate > 0 ~ paste(cfield.1, " > ", cfield.2),
        p.value < .05 & estimate < 0 ~ paste(cfield.2, " > ", cfield.1),
        p.value >= .05 & p.value < .1 & estimate > 0 ~ paste(cfield.1, " *> ", cfield.2),
        p.value >= .05 & p.value < .1 & estimate < 0 ~ paste(cfield.2, " *> ", cfield.1),
        p.value >= .1 & estimate > 0 ~ paste(cfield.1, " = ", cfield.2),
        p.value >= .1 & estimate < 0 ~ paste(cfield.2, " = ", cfield.1),
        TRUE ~ 'not computable'
      ))
    
    if (length(unique(categorical_factor_names)) == 1) {
      model.emmeans = model.emmeans %>% select(term, contrastfield, contrast, cohensd, p.value)
    } else {
      model.emmeans = model.emmeans %>%
        unite(col = "insidelevel", all_of(categorical_factor_names),
              sep = " ", remove = TRUE, na.rm = TRUE) %>%
        select(term, contrastfield, insidelevel, contrast, cohensd, p.value)
    }
    
    results.table = model.emmeans %>%
      tidyr::nest(tests = -c(term, contrastfield)) %>%
      pivot_wider(names_from = contrastfield, values_from = tests)
    
    for (cat_colname in colnames(results.table)[-1]) {
      results.table = results.table %>%
        rename_with(~str_c("f", which(categorical_factor_names == cat_colname)), all_of(cat_colname))
    }
    rm(cat_colname)
    
    for (catfactor in 1:length(categorical_factor_names)) {
      assign(str_c("results.table", catfactor),
             results.table %>% select(term, str_c("f", all_of(catfactor))) %>%
               janitor::remove_empty("cols") %>%
               unnest(cols = str_c("f", catfactor), names_sep = "."))
    }
    rm(catfactor)
    
    results.table = tibble()
    for (pair in 1:length(categorical_factor_names)) {
      results.table = bind_rows(results.table,
                                str_c("results.table", pair) %>% parse(text = .) %>% eval())
    }
    rm(pair)
    rm(list = str_c("results.table", 1:length(categorical_factor_names)))
    
    results.table = results.table %>%
      full_join(model.anova, by = "term") %>%
      relocate(chisq, df, p_chisq, .after = term)
    
  } else {
    results.table = model.anova
  }
  
  rm(model.anova, model.emmeans)
  
  term = term %>% mutate(order1 = str_count(combinations))
  term = term %>% mutate(order2 = str_arrange(combinations))
  
  results.table = results.table %>%
    left_join(term %>% select(-combinations), by = "term", relationship = "many-to-many") %>%
    arrange(order1, order2) %>%
    select(-c(order1, order2)) %>%
    distinct()
  
  rm(term, str_arrange, combinations)
  
  results.table = results.table %>% filter(!is.na(chisq))
  
  # Repair rounding and effect sizes
  results.table = results.table %>%
    mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>%
    mutate(across(contains("cohensd"), abs))
  
  # -------------------------
  # 4) Output formatting (unchanged)
  # -------------------------
  if (!is.na(path)) {
    
    library(flextable)
    set_flextable_defaults(font.family = "Arial", font.size = 9, digits = 3)
    
    number_of_columns = ncol(results.table)
    number_of_fs = round((number_of_columns-4)/4)
    
    results = results.table %>% flextable() %>% theme_vanilla()
    
    results = results %>%
      colformat_num(j = "df", digits = 0) %>%
      merge_v(j = "term", target = c(1:4)) %>%
      bg(~ (p_chisq < .1 & p_chisq >= .05), ~ term + p_chisq, bg = "#d7eef3") %>%
      bg(~ p_chisq < .05, ~ term + p_chisq, bg = "#90ccde")
    if(number_of_columns > 4) {results = results %>% merge_v(j = c(5:number_of_columns))}
    
    if (number_of_columns >= 7) {
      for (f in 1:number_of_fs) {
        f.p.value = str_c("f", f, ".p.value")
        if (f.p.value %in% results$body$col_keys) {
          results = results %>%
            bg(i = (results$body$dataset[[f.p.value]] < .1 & results$body$dataset[[f.p.value]] >= .05),
               j = which(results$body$col_keys == f.p.value), bg = "#d7eef3") %>%
            bg(i = results$body$dataset[[f.p.value]] < .05,
               j = which(results$body$col_keys == f.p.value), bg = "#90ccde")
        }
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
    
    for (alias in 1:number_of_columns) {
      content = results$header$content$data[[alias]]$txt
      if (content == "term") {results$header$content$data[[alias]]$txt = "Term"}
      else if (content == "chisq") {results$header$content$data[[alias]]$txt = "χ²"}
      else if (content == "df") {results$header$content$data[[alias]]$txt = "df"}
      else if (content == "p_chisq") {results$header$content$data[[alias]]$txt = "p"}
      else if (grepl("insidelevel", content, fixed = TRUE)) {results$header$content$data[[alias]]$txt = "Levels"}
      else if (grepl("contrast", content, fixed = TRUE)) {results$header$content$data[[alias]]$txt = "Contrast"}
      else if (grepl("cohensd", content, fixed = TRUE)) {results$header$content$data[[alias]]$txt = "Coh.d"}
      else if (grepl("p.value", content, fixed = TRUE)) {results$header$content$data[[alias]]$txt = "Sig"}
      rm(content)
    }
    rm(alias, number_of_columns)
    
    results = results %>% add_header_row(values = values, top = TRUE); rm(values)
    results = results %>% merge_at(i = 1, j = 1:4, part = "header")
    if (number_of_fs == 1) {results = results %>% merge_at(i = 1, j = 5:7, part = "header")}
    if (number_of_fs > 1) {
      for (f in 1:length(categorical_factor_names)) {
        results = results %>% merge_at(i = 1, j = (1:4)+4*f, part = "header")
      }
      rm(f)
    }
    
    results = results %>% colformat_double(digits = 3) %>% padding(padding = 4, part = "all") %>% autofit()
    
    rm(list = str_c("f", 1:number_of_fs))
    rm(categorical_factor_names, factor_names, number_of_fs)
    
    results %>% save_as_html(path = path, title = title)
  }
  
  if (extract == TRUE) {return(results.table)}
}

