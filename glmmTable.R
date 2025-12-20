glmmTable = function(model, path = NA, title = "Model", extract = FALSE, adjust = "bonferroni") {
  
  library(tidyverse); library(car); library(janitor)
  library(emmeans)
  
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
  
  # robust model frame (avoid @)
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
  
  # pick whatever p-value column exists and rename to p_chisq
  p_candidates <- c("pr_chisq", "pr_f", "pr_gt_chisq", "pr_gt_f", "p")
  p_col <- p_candidates[p_candidates %in% names(model.anova)][1]
  if (!is.na(p_col)) {
    model.anova <- model.anova %>% dplyr::rename(p_chisq = !!rlang::sym(p_col))
  } else {
    model.anova <- model.anova %>% mutate(p_chisq = NA_real_)
  }
  
  # standardize statistic column name to chisq if needed
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
  
  # Fallback: derive from model_frame
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
  
  # --- FIX: ensure columns exist ---
  if (!"name" %in% names(terms_tbl))  terms_tbl$name  <- character()
  if (!"value" %in% names(terms_tbl)) terms_tbl$value <- NA_character_
  
  factor_names <- terms_tbl %>% dplyr::pull(name) %>% unique()
  
  categorical_factor_names <- terms_tbl %>%
    dplyr::filter(.data$value %in% c("factor", "character")) %>%
    dplyr::pull(name) %>%
    unique()
  
  non_categorical_factor_names <- setdiff(factor_names, categorical_factor_names)
  
  # If no categorical factors: return ANOVA only
  if (length(categorical_factor_names) < 1) {
    results.table <- model.anova
    if (extract) return(results.table)
    return(invisible(NULL))
  }
  
  # -------------------------
  # 1b) Drop unused levels in model_frame (critical for pairwise)
  # -------------------------
  if (!is.null(model_frame)) {
    for (v in categorical_factor_names) {
      if (v %in% names(model_frame)) {
        if (is.factor(model_frame[[v]])) model_frame[[v]] <- droplevels(model_frame[[v]])
        if (is.character(model_frame[[v]])) model_frame[[v]] <- factor(model_frame[[v]]) %>% droplevels()
      }
    }
  }
  
  # -------------------------
  # 2) Terms to contrast = ANOVA terms (categorical only)
  # -------------------------
  
  anova_terms <- model.anova %>%
    dplyr::filter(!is.na(term)) %>%
    dplyr::pull(term) %>%
    unique()
  
  anova_terms <- anova_terms[anova_terms != "(Intercept)"]
  anova_terms <- anova_terms[!tolower(anova_terms) %in% c("residuals", "residual")]
  
  split_term_vars <- function(term) {
    term <- stringr::str_replace_all(term, fixed("*"), ":")
    unlist(stringr::str_split(term, fixed(":")))
  }
  
  # -------------------------
  # 3) Emmeans + pairs for EACH categorical ANOVA term
  # -------------------------
  
  model.emmeans <- tibble::tibble()
  
  for (tt in anova_terms) {
    
    vars <- split_term_vars(tt)
    
    # only if ALL variables in this ANOVA term are categorical
    if (!all(vars %in% categorical_factor_names)) next
    
    # emmeans formula: ~ A  or ~ A:B etc.
    fml <- stats::as.formula(paste0("~ ", paste(vars, collapse = ":")))
    
    # compute emmeans (try to force use of the same data)
    emm_obj <- tryCatch(
      emmeans::emmeans(model, fml),
      error = function(e) NULL
    )
    if (is.null(emm_obj)) next
    
    emm_df <- tryCatch(as.data.frame(emm_obj), error = function(e) NULL)
    if (is.null(emm_df) || nrow(emm_df) < 2) next
    
    # pairs
    pair_df <- tryCatch({
      emmeans::pairs(emm_obj, adjust = adjust) %>%
        summary() %>%
        as.data.frame()
    }, error = function(e) NULL)
    
    if (is.null(pair_df) || nrow(pair_df) == 0) next
    
    res <- pair_df %>% tibble::as_tibble()
    
    # standardize column name
    if ("contrast1" %in% names(res) && !"contrast" %in% names(res)) {
      res <- res %>% dplyr::rename(contrast = contrast1)
    }
    
    # Attach term metadata + effect size (one per term)
    res <- res %>%
      dplyr::mutate(
        term = tt,
        contrastfield = tt,
        cohensd = safe_eff_size(emm_obj, model)
      )
    
    model.emmeans <- dplyr::bind_rows(model.emmeans, res)
  }
  
  # If no contrasts: return ANOVA only
  if (nrow(model.emmeans) == 0) {
    results.table <- model.anova
    if (extract) return(results.table)
    return(invisible(NULL))
  }
  
  # -------------------------
  # 3b) Format contrasts (keep both raw + directional label)
  # -------------------------
  
  model.emmeans <- model.emmeans %>%
    dplyr::mutate(df = round(df, 2)) %>%
    tidyr::separate(contrast, into = c("cfield.1","cfield.2"), sep = " - ", remove = FALSE) %>%
    dplyr::mutate(
      contrast_dir = dplyr::case_when(
        p.value < .05 & estimate > 0 ~ paste(cfield.1, " > ", cfield.2),
        p.value < .05 & estimate < 0 ~ paste(cfield.2, " > ", cfield.1),
        p.value >= .05 & p.value < .1 & estimate > 0 ~ paste(cfield.1, " *> ", cfield.2),
        p.value >= .05 & p.value < .1 & estimate < 0 ~ paste(cfield.2, " *> ", cfield.1),
        p.value >= .1 & estimate > 0 ~ paste(cfield.1, " = ", cfield.2),
        p.value >= .1 & estimate < 0 ~ paste(cfield.2, " = ", cfield.1),
        TRUE ~ "not computable"
      )
    ) %>%
    dplyr::select(
      term, contrastfield,
      contrast_dir,
      estimate, SE, df, t.ratio, p.value,
      cohensd
    )
  
  # -------------------------
  # 3c) Wide table: one block per ANOVA term
  # -------------------------
  
  results.table <- model.emmeans %>%
    tidyr::nest(tests = -c(term, contrastfield)) %>%
    tidyr::pivot_wider(names_from = contrastfield, values_from = tests) %>%
    dplyr::full_join(model.anova, by = "term") %>%
    dplyr::relocate(chisq, df, p_chisq, .after = term)
  
  results.table <- results.table %>%
    dplyr::filter(!is.na(chisq)) %>%
    dplyr::mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>%
    dplyr::mutate(across(contains("cohensd"), abs))
  
  # -------------------------
  # 4) Output formatting (kept compatible with your current)
  # -------------------------
  if (!is.na(path)) {
    
    library(flextable)
    set_flextable_defaults(font.family = "Arial", font.size = 9, digits = 3)
    
    number_of_columns <- ncol(results.table)
    number_of_fs <- round((number_of_columns - 4) / 4)
    
    results <- results.table %>% flextable() %>% theme_vanilla()
    
    results <- results %>%
      colformat_num(j = "df", digits = 0) %>%
      merge_v(j = "term", target = c(1:4)) %>%
      bg(~ (p_chisq < .1 & p_chisq >= .05), ~ term + p_chisq, bg = "#d7eef3") %>%
      bg(~ p_chisq < .05, ~ term + p_chisq, bg = "#90ccde")
    
    if (number_of_columns > 4) {
      results <- results %>% merge_v(j = c(5:number_of_columns))
    }
    
    if (number_of_columns >= 7) {
      for (f in 1:number_of_fs) {
        f.p.value <- str_c("f", f, ".p.value")
        if (f.p.value %in% results$body$col_keys) {
          results <- results %>%
            bg(i = (results$body$dataset[[f.p.value]] < .1 & results$body$dataset[[f.p.value]] >= .05),
               j = which(results$body$col_keys == f.p.value), bg = "#d7eef3") %>%
            bg(i = results$body$dataset[[f.p.value]] < .05,
               j = which(results$body$col_keys == f.p.value), bg = "#90ccde")
        }
      }
      rm(f, f.p.value)
    }
    
    values <- c("Model results", "", "", "")
    results <- results %>% add_header_row(values = values, top = TRUE)
    results <- results %>% merge_at(i = 1, j = 1:4, part = "header")
    
    results <- results %>% colformat_double(digits = 3) %>%
      padding(padding = 4, part = "all") %>%
      autofit()
    
    results %>% save_as_html(path = path, title = title)
  }
  
  if (extract) return(results.table)
}
