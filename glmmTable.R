glmmTable <- function(model, path = NA, title = "Model", extract = FALSE, adjust = "bonferroni") {
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(car)
    library(janitor)
    library(emmeans)
    library(lmerTest)
  })
  
  emm_options(
    lmerTest.limit = 10000,
    disable.pbkrtest = TRUE,
    lmer.df = "satterthwaite",
    msg.interaction = FALSE
  )
  
  # -------------------------
  # 0) Helpers (robustness)
  # -------------------------
  
  is_glmmTMB <- inherits(model, "glmmTMB")
  is_lm      <- inherits(model, "lm") && !inherits(model, "glm")
  is_glm     <- inherits(model, "glm") && !inherits(model, "glmerMod") && !is_glmmTMB
  is_lmer    <- inherits(model, "lmerMod")
  is_glmer   <- inherits(model, "glmerMod")
  
  model_frame <- tryCatch(
    model.frame(model),
    error = function(e) {
      if (is_glmmTMB && !is.null(model$frame)) return(model$frame)
      if (!is.null(model$model)) return(model$model)
      NULL
    }
  )
  
  safe_df_resid <- function(m) {
    out <- tryCatch(df.residual(m), error = function(e) NA_real_)
    as.numeric(out)
  }
  
  safe_sigma <- function(m) {
    out <- tryCatch(sigma(m), error = function(e) NA_real_)
    as.numeric(out)
  }
  
  safe_eff_size_tbl <- function(emm_obj, m) {
    tryCatch({
      es <- eff_size(emm_obj, sigma = safe_sigma(m), edf = safe_df_resid(m))
      as_tibble(summary(es)) %>%
        select(contrast, cohensd = effect.size)
    }, error = function(e) {
      # Si falla, retornem taula buida amb les columnes esperades
      tibble(contrast = character(), cohensd = numeric())
    })
  }
  
  safe_anova <- function(m) {
    
    test_try <- list(
      if (is_lm) "F" else NULL,
      if (is_glm || is_glmer || is_glmmTMB) "Chisq" else NULL,
      c("Chisq", "F"),
      "F",
      "Chisq"
    ) %>% purrr::compact()
    
    for (ts in test_try) {
      out <- tryCatch({
        car::Anova(m, type = c("II", "III", 2, 3), test.statistic = ts) %>%
          data.frame() %>%
          janitor::clean_names() %>%
          tibble::rownames_to_column("term")
      }, error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    
    tibble(term = character(), chisq = numeric(), df = numeric(), p_chisq = numeric())
  }
  
  model.anova <- safe_anova(model)
  
  # Standardize p-value column
  p_candidates <- c("pr_chisq", "pr_f", "pr_gt_chisq", "pr_gt_f", "p")
  p_col <- p_candidates[p_candidates %in% names(model.anova)][1]
  if (!is.na(p_col)) {
    model.anova <- model.anova %>% dplyr::rename(p_chisq = !!rlang::sym(p_col))
  } else {
    model.anova <- model.anova %>% mutate(p_chisq = NA_real_)
  }
  
  # Standardize chisq/df
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
  if (!"chisq" %in% names(model.anova)) {
    model.anova <- model.anova %>% mutate(chisq = NA_real_)
  }
  
  # -------------------------
  # 1) Identify categorical predictors robustly
  # -------------------------
  
  terms_tbl <- tryCatch({
    dc <- attr(terms(model), "dataClasses")
    dc <- dc[names(dc) != "(Intercept)"]
    tibble(name = names(dc), value = unname(dc))
  }, error = function(e) NULL)
  
  if (is.null(terms_tbl) && is_glmmTMB) {
    terms_tbl <- tryCatch({
      attr(model$modelInfo$reTrms$cond$terms$fixed, "dataClasses") %>%
        enframe(name = "name", value = "value") %>%
        tail(-1)
    }, error = function(e) NULL)
  }
  
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
  
  # Safety: if something weird came back, enforce column names
  if (!all(c("name", "value") %in% names(terms_tbl))) {
    if (ncol(terms_tbl) >= 2) {
      names(terms_tbl)[1:2] <- c("name", "value")
    } else {
      terms_tbl <- tibble(name = character(), value = character())
    }
  }
  
  factor_names <- terms_tbl %>% dplyr::pull(name)
  categorical_factor_names <- terms_tbl %>%
    filter(value %in% c("factor", "character", "ordered")) %>%
    pull(name)
  
  # If no categorical predictors, just return ANOVA
  if (length(categorical_factor_names) < 1) {
    results.table <- model.anova
    if (extract) return(results.table)
    return(invisible(NULL))
  }
  
  # -------------------------
  # 2) Extract ANOVA terms and map them to "simple-effects" specs
  # -------------------------
  
  anova_terms <- model.anova %>%
    filter(!is.na(term)) %>%
    pull(term) %>%
    unique()
  
  anova_terms <- anova_terms[anova_terms != "(Intercept)"]
  anova_terms <- anova_terms[!tolower(anova_terms) %in% c("residuals", "residual")]
  
  split_term_vars <- function(term) {
    term <- stringr::str_replace_all(term, fixed("*"), ":")
    unlist(stringr::str_split(term, fixed(":")))
  }
  
  # Build a list of emmeans specs to compute, for EACH ANOVA term:
  # - main: focal
  # - interaction: focal | other*other...
  specs_tbl <- tibble()
  
  for (tt in anova_terms) {
    vars <- split_term_vars(tt)
    
    # keep only terms where ALL vars are categorical
    if (!all(vars %in% categorical_factor_names)) next
    
    # for each variable in the term, compute simple effects of that focal variable
    for (focal in vars) {
      others <- setdiff(vars, focal)
      
      if (length(others) == 0) {
        spec <- paste0("~ ", focal)
        inside <- NA_character_
      } else {
        spec <- paste0("~ ", focal, " | ", paste(others, collapse = "*"))
        inside <- paste(others, collapse = "*")
      }
      
      specs_tbl <- bind_rows(specs_tbl, tibble(
        term = tt,               # ANOVA term
        contrastfield = focal,   # which factor we're contrasting
        inside_by = inside,      # conditioning factors (if any)
        spec = spec
      ))
    }
  }
  
  # If ANOVA had no categorical terms we can contrast
  if (nrow(specs_tbl) == 0) {
    results.table <- model.anova
    if (extract) return(results.table)
    return(invisible(NULL))
  }
  
  # -------------------------
  # 3) Compute emmeans + pairs for EVERY spec (ALWAYS, if computable)
  # -------------------------
  
  model.emmeans <- tibble()
  
  for (i in seq_len(nrow(specs_tbl))) {
    
    tt   <- specs_tbl$term[i]
    cf   <- specs_tbl$contrastfield[i]
    iby  <- specs_tbl$inside_by[i]
    spec <- specs_tbl$spec[i]
    
    emm_obj <- tryCatch(
      emmeans::emmeans(model, specs = as.formula(spec)),
      error = function(e) NULL
    )
    if (is.null(emm_obj)) next
    
    # If only 1 EMM row, no pairs possible
    emm_df <- tryCatch(as.data.frame(emm_obj), error = function(e) NULL)
    if (is.null(emm_df) || nrow(emm_df) < 2) next
    
    pair_df <- tryCatch({
      pairs(emm_obj, adjust = adjust) %>%
        summary() %>%
        as.data.frame()
    }, error = function(e) NULL)
    
    if (is.null(pair_df) || nrow(pair_df) == 0) next
    
    res <- pair_df %>% as_tibble()
    
    # ensure 'contrast' column name
    if ("contrast1" %in% names(res) && !"contrast" %in% names(res)) {
      res <- res %>% rename(contrast = contrast1)
    }
    
    res <- res %>%
      left_join(safe_eff_size_tbl(emm_obj, model), by = "contrast") %>%
      mutate(
        term = tt,
        contrastfield = cf,
        inside_by = iby
      )
    
    model.emmeans <- bind_rows(model.emmeans, res)
  }
  
  # If no contrasts computed, return ANOVA
  if (nrow(model.emmeans) == 0) {
    results.table <- model.anova
    if (extract) return(results.table)
    return(invisible(NULL))
  }
  
  # -------------------------
  # 3b) Keep a readable directional label BUT also keep numeric columns
  # -------------------------
  
  model.emmeans <- model.emmeans %>%
    mutate(df = round(df, 2)) %>%
    tidyr::separate(contrast, into = c("cfield.1", "cfield.2"), sep = " - ", remove = FALSE) %>%
    mutate(
      contrast_dir = case_when(
        p.value < .05 & estimate > 0 ~ paste(cfield.1, " > ", cfield.2),
        p.value < .05 & estimate < 0 ~ paste(cfield.2, " > ", cfield.1),
        p.value >= .05 & p.value < .1 & estimate > 0 ~ paste(cfield.1, " *> ", cfield.2),
        p.value >= .05 & p.value < .1 & estimate < 0 ~ paste(cfield.2, " *> ", cfield.1),
        p.value >= .1 & estimate > 0 ~ paste(cfield.1, " = ", cfield.2),
        p.value >= .1 & estimate < 0 ~ paste(cfield.2, " = ", cfield.1),
        TRUE ~ "not computable"
      ),
      estimate = abs(estimate),
      t.ratio  = abs(t.ratio),
      cohensd  = abs(cohensd)
    ) %>%
    select(
      term, contrastfield, inside_by,
      contrast_dir,
      estimate, SE, df, t.ratio, p.value,
      cohensd
    )
  
  # -------------------------
  # 3c) Wide table per ANOVA term
  #     We'll keep one nested column per contrastfield, inside each term.
  # -------------------------
  
  results.table <- model.emmeans %>%
    nest(tests = -c(term, contrastfield)) %>%
    pivot_wider(names_from = contrastfield, values_from = tests) %>%
    full_join(model.anova, by = "term") %>%
    relocate(chisq, df, p_chisq, .after = term) %>%
    filter(!is.na(chisq)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(across(contains("cohensd"), abs))
  
  # -------------------------
  # 4) Output formatting (your original style, kept minimal)
  # -------------------------
  
  results.table_clean <- results.table %>%
    mutate(
      across(
        where(is.list),
        ~ map(.x, \(df) df %>%
                janitor::remove_empty("cols") %>%
                transmute(
                  insidelevel = dplyr::if_else("insidelevel" %in% names(df),
                                               as.character(df$insidelevel),
                                               NA_character_) |> (\(x) NULL)(),
                  contrast_dir,
                  d = abs(cohensd),
                  t = abs(t.ratio),
                  p = p.value
                ))
      )
    )
  
  list_cols <- names(results.table_clean)[sapply(results.table_clean, is.list)]
  
  results.table_long <- reduce(
    list_cols,
    .init = results.table_clean,
    .f = function(acc, col) {
      tidyr::unnest(acc, cols = all_of(col), names_sep = ".")
    }
  )
  
  if (!is.na(path)) {
    
    library(flextable)
    
    set_flextable_defaults(font.family = "Arial", font.size = 9, digits = 3)
    
    global_cols <- c("term", "chisq", "df", "p_chisq")
    all_cols    <- names(results.table_long)
    
    # detecta prefixos (abans del punt) que NO són globals
    detail_cols <- setdiff(all_cols, global_cols)
    prefixes <- detail_cols %>%
      str_extract("^[^.]+") %>%
      unique()
    
    # quantes columnes té cada prefix (en l'ordre en què apareixen)
    prefix_ncols <- vapply(prefixes, function(p) sum(str_detect(detail_cols, paste0("^", p, "\\."))), numeric(1))
    
    # Header superior (agrupació)
    top_values   <- c("Model results", prefixes)
    top_widths   <- c(length(global_cols), as.integer(prefix_ncols))
    
    pretty_names <- all_cols
    
    # globals (bonics)
    pretty_names[pretty_names == "term"]    <- "Term"
    pretty_names[pretty_names == "chisq"]   <- "χ²"
    pretty_names[pretty_names == "df"]      <- "df"
    pretty_names[pretty_names == "p_chisq"] <- "p"
    
    # subheaders curts (sense prefix)
    pretty_names <- ifelse(
      all_cols %in% global_cols,
      pretty_names,
      sub("^[^.]+\\.", "", all_cols)  # treu "group." / "phone." / etc.
    )
    
    # i ara canvia noms interns a labels APA
    pretty_names <- pretty_names %>%
      str_replace("^insidelevel$", "Levels") %>%
      str_replace("^contrast_dir$", "Contrast") %>%
      str_replace("^d$", "d") %>%
      str_replace("^t$", "t") %>%
      str_replace("^p$", "p")
    
    # Ara crea la flextable
    ft <- flextable(results.table_long, col_keys = all_cols) %>%
      theme_vanilla() %>%
      set_header_labels(values = setNames(pretty_names, all_cols)) %>%
      add_header_row(values = top_values, colwidths = top_widths, top = TRUE) %>%
      merge_h(part = "header") %>%
      colformat_double(digits = 3) %>%
      colformat_num(j = "df", digits = 0) %>%
      merge_v(j = global_cols) %>%
      bg(~ (p_chisq < .1 & p_chisq >= .05), j = c("term","p_chisq"), bg = "#d7eef3") %>%
      bg(~ p_chisq < .05,                  j = c("term","p_chisq"), bg = "#90ccde") %>%
      padding(padding = 4, part = "all") %>%
      autofit()
    
    p_cols <- grep("\\.p$", names(results.table_long), value = TRUE)
    if (length(p_cols) > 0) { ft <- ft %>% merge_v(j = p_cols) }
    
    ft %>% save_as_html(path = path, title = title)
  }
  
  if (extract) return(results.table_long)
  invisible(NULL)
}
