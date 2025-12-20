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
      es  <- eff_size(emm_obj, sigma = safe_sigma(m), edf = safe_df_resid(m))
      out <- as_tibble(summary(es))
      
      if ("contrast1" %in% names(out) && !"contrast" %in% names(out)) {
        out <- out %>% rename(contrast = contrast1)
      }
      
      out %>%
        rename(cohensd = effect.size) %>%
        select(-any_of(c("SE","df","lower.CL","upper.CL","t.ratio","z.ratio","p.value")))
      # IMPORTANT: NO eliminis les columnes by
    }, error = function(e) tibble(contrast = character(), cohensd = numeric()))
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
  
  if (is_empty(terms_tbl) && is_glmmTMB) {
    terms_tbl <- tryCatch({
      attr(model$modelInfo$reTrms$cond$terms$fixed, "dataClasses") %>%
        enframe(name = "name", value = "value") %>%
        tail(-1)
    }, error = function(e) NULL)
  }
  
  if (is_empty(terms_tbl)) {
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
    
    eff_tbl <- safe_eff_size_tbl(emm_obj, model)
    
    if ("contrast1" %in% names(eff_tbl) && !"contrast" %in% names(eff_tbl)) {
      eff_tbl <- eff_tbl %>% rename(contrast = contrast1)
    }
    
    # Columnes que apareixen tant a res com a eff_tbl (inclourà contrast + by vars)
    join_by <- intersect(names(res), names(eff_tbl))
    
    # evita que s'hi coli estimate/SE/etc al join:
    stats <- c("estimate","SE","df","t.ratio","z.ratio","p.value",
               "lower.CL","upper.CL","asymp.LCL","asymp.UCL")
    join_by <- setdiff(join_by, stats)
    if (!"contrast" %in% join_by) join_by <- c("contrast", join_by)
    
    by_vars <- if (!is.na(iby)) strsplit(iby, "\\*")[[1]] else character()
    
    res <- res %>%
      mutate(
        insidelevel = if (length(by_vars) > 0 && all(by_vars %in% names(.))) {
          pmap_chr(across(all_of(by_vars)), \(...) {
            vals <- c(...)
            paste(paste0(by_vars, "=", vals), collapse = ", ")
          })
        } else {
          NA_character_
        }
      )
    
    n_before <- nrow(res)
    
    # Si l'eff_tbl no és usable per estrats, no fem join (evita matches incorrectes)
    if (nrow(eff_tbl) == 0 || setequal(sort(names(eff_tbl)), sort(c("contrast","cohensd")))) {
      res$cohensd <- NA_real_
    } else {
      res <- res %>% left_join(eff_tbl, by = join_by)
    }
    
    # AQUEST mutate ha d'anar SEMPRE
    res <- res %>%
      mutate(
        term = tt,
        contrastfield = cf,
        inside_by = iby
      )
    
    if (nrow(res) != n_before) {
      warning("Join duplicated rows for spec: ", spec,
              " | join_by = ", paste(join_by, collapse = ", "),
              " | eff_tbl cols = ", paste(names(eff_tbl), collapse = ", "))
    }
    
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
  
  has_t <- "t.ratio" %in% names(model.emmeans)
  has_z <- "z.ratio" %in% names(model.emmeans)
  
  model.emmeans <- model.emmeans %>%
    mutate(df = round(df, 2)) %>%
    tidyr::separate(contrast, into = c("cfield.1", "cfield.2"), sep = " - ", remove = FALSE) %>%
    mutate(
      contrast_dir = case_when(
        p.value < .05 & estimate > 0 ~ paste(cfield.1, " > ",  cfield.2),
        p.value < .05 & estimate < 0 ~ paste(cfield.2, " > ",  cfield.1),
        p.value >= .05 & p.value < .1 & estimate > 0 ~ paste(cfield.1, " *> ", cfield.2),
        p.value >= .05 & p.value < .1 & estimate < 0 ~ paste(cfield.2, " *> ", cfield.1),
        p.value >= .1 & estimate > 0 ~ paste(cfield.1, " ≈ ",  cfield.2),
        p.value >= .1 & estimate < 0 ~ paste(cfield.2, " ≈ ",  cfield.1),
        TRUE ~ "not computable"
      ),
      estimate = abs(estimate),
      cohensd  = abs(cohensd)
    ) %>%
    { 
      df <- .
      if (has_t) {
        df <- df %>%
          mutate(
            t.ratio = abs(.data[["t.ratio"]]),
            t = abs(.data[["t.ratio"]])
          )
      } else {
        df <- df %>% mutate(t = NA_real_)
      }
      if (has_z) {
        df <- df %>%
          mutate(
            z.ratio = abs(.data[["z.ratio"]]),
            z = abs(.data[["z.ratio"]])
          )
      } else {
        df <- df %>% mutate(z = NA_real_)
      }
      df
    }
  
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
  
  # --- 1) normalitza qualsevol minitaula (o NULL) a un format estàndard ---
  clean_nested_tests <- function(df) {
    if (is.null(df) || !inherits(df, "data.frame") || nrow(df) == 0) {
      return(tibble(
        insidelevel  = character(),
        contrast_dir = character(),
        d            = numeric(),
        t            = numeric(),
        z            = numeric(),
        p            = numeric()
      ))
    }
    
    df <- janitor::remove_empty(df, "cols")
    
    tibble(
      insidelevel  = if ("insidelevel" %in% names(df)) as.character(df$insidelevel) else NA_character_,
      contrast_dir = if ("contrast_dir" %in% names(df)) df$contrast_dir else NA_character_,
      d            = if ("cohensd" %in% names(df)) abs(df$cohensd) else NA_real_,
      t            = if ("t.ratio" %in% names(df)) abs(df$t.ratio) else if ("t" %in% names(df)) df$t else NA_real_,
      z            = if ("z.ratio" %in% names(df)) abs(df$z.ratio) else if ("z" %in% names(df)) df$z else NA_real_,
      p            = if ("p.value" %in% names(df)) df$p.value else NA_real_
    )
  }
  
  # --- 2) results.table (nested) -> excel_long (tidy) ---
  excel_long <- results.table %>%
    mutate(across(where(is.list), ~ purrr::map(.x, clean_nested_tests))) %>%
    tidyr::pivot_longer(where(is.list), names_to = "contrastfield", values_to = "tests") %>%
    tidyr::unnest(tests, keep_empty = TRUE) %>%
    filter(!is.na(contrast_dir) & !is.na(p)) %>%
    mutate(.stat = dplyr::coalesce(t, z)) %>%
    select(term, chisq, df, p_chisq, contrastfield,
           insidelevel, contrast_dir, d, t, z, p, .stat)
  
  # decideix quin estadístic tens (mirant les columnes disponibles al teu excel_long actual)
  stat_name <- if ("t" %in% names(excel_long) && any(!is.na(excel_long$t))) "t" else "z"
  
  excel_long <- excel_long %>%
    mutate(.stat = dplyr::coalesce(t, z)) %>%
    select(-any_of(c("t","z"))) %>%
    rename(!!stat_name := .stat)
  
  # --- 3) excel_long -> excel_wide (real wide per Excel) + ordre columnes ---
  
  global_cols <- c("term", "chisq", "df", "p_chisq")
  
  term_levels <- excel_long %>% distinct(term) %>% pull(term)
  
  factor_levels <- excel_long %>%
    mutate(term = factor(term, levels = term_levels)) %>%
    arrange(term) %>%
    distinct(contrastfield) %>%
    pull(contrastfield)
  
  measures_pref <- c("insidelevel", "contrast_dir", "d", "t", "z", "p")
  measures <- intersect(measures_pref, names(excel_long))
  
  excel_wide <- excel_long %>%
    mutate(term = factor(term, levels = term_levels)) %>%
    group_by(term, contrastfield) %>%
    arrange(insidelevel, contrast_dir, .by_group = TRUE) %>%
    mutate(row_id = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols    = c(term, chisq, df, p_chisq, row_id),
      names_from = contrastfield,
      values_from = all_of(measures),
      names_glue = "{contrastfield}.{.value}"   # <<< CLAU: factor.mesura
    ) %>%
    arrange(term, row_id) %>%
    select(-row_id)
  
  wide_cols <- purrr::map(factor_levels, \(f) paste(f, measures, sep = ".")) %>%
    unlist(use.names = FALSE) %>%
    intersect(names(excel_wide))
  
  excel_wide <- excel_wide %>%
    select(all_of(global_cols), all_of(wide_cols))
  
  # -------------------------
  # 5) Flextable
  # -------------------------
  
  if (!is.na(path)) {
    
    library(flextable)
    library(officer)
    
    set_flextable_defaults(font.family = "Arial", font.size = 9, digits = 3)
    
    global_cols   <- c("term", "chisq", "df", "p_chisq")
    all_cols      <- names(excel_wide)
    detail_cols   <- setdiff(all_cols, global_cols)
    
    # prefixos (Factor) per fer el header superior
    prefixes <- str_extract(detail_cols, "^[^.]+") %>% unique()
    prefix_ncols <- vapply(
      prefixes,
      function(p) sum(str_detect(detail_cols, paste0("^", p, "\\."))),
      numeric(1)
    )
    
    top_values <- c("Model results", prefixes)
    top_widths <- c(length(global_cols), as.integer(prefix_ncols))
    
    # noms "bonics"
    pretty_names <- all_cols
    pretty_names[pretty_names == "term"]    <- "Term"
    pretty_names[pretty_names == "chisq"]   <- "χ²"
    pretty_names[pretty_names == "df"]      <- "df"
    pretty_names[pretty_names == "p_chisq"] <- "p"
    
    pretty_names <- ifelse(all_cols %in% global_cols, pretty_names, sub("^[^.]+\\.", "", all_cols))
    
    pretty_names <- pretty_names %>%
      str_replace("^insidelevel$", "Levels") %>%
      str_replace("^contrast_dir$", "Contrast") %>%
      str_replace("^d$", "d") %>%
      str_replace("^t$", "t") %>%
      str_replace("^z$", "z") %>%
      str_replace("^p$", "p")
    
    # 1) crea flextable DIRECTAMENT des de excel_wide (manté numèrics!)
    ft <- flextable(excel_wide, col_keys = all_cols) %>%
      theme_vanilla() %>%
      set_header_labels(values = setNames(pretty_names, all_cols)) %>%
      add_header_row(values = top_values, colwidths = top_widths, top = TRUE) %>%
      merge_h(part = "header")
    
    # 2) MERGE GLOBALS "SEGONS term" (com el teu script original)
    #    - això fa que chisq/df/p_chisq s’estenguin cap avall dins de cada Term
    #    - i MAI no es barregin entre termes diferents, encara que tinguin el mateix número
    ft <- ft %>%
      merge_v(j = "term", target = c("chisq", "df", "p_chisq"), part = "body") %>%
      merge_v(j = "term", part = "body")
    
    # 3) format numèric (sense convertir a text)
    #    - globals:
    ft <- ft %>%
      colformat_double(j = "chisq", digits = 3, na_str = "") %>%
      colformat_num(j = "df", digits = 0, na_str = "") %>%
      colformat_double(j = "p_chisq", digits = 3, na_str = "")
    
    #    - detalls (numèrics) amb 3 decimals i NA buit
    num_detail <- detail_cols[vapply(excel_wide[detail_cols], is.numeric, logical(1))]
    if (length(num_detail) > 0) {
      ft <- ft %>% colformat_double(j = num_detail, digits = 3, na_str = "")
    }
    
    ft <- ft %>%
      padding(padding = 4, part = "all") %>%
      autofit()
    
    # 4) colors globals (com abans)
    ft <- ft %>%
      bg(~ (p_chisq < .1 & p_chisq >= .05), j = c("term", "p_chisq"), bg = "#d7eef3") %>%
      bg(~ p_chisq < .05,                  j = c("term", "p_chisq"), bg = "#90ccde")
    
    # 5) colors també per a les p de les minitaules (totes les columnes que acaben en ".p")
    p_detail_cols <- grep("\\.p$", detail_cols, value = TRUE)
    if (length(p_detail_cols) > 0) {
      for (pc in p_detail_cols) {
        # important: funciona encara que hi hagi NA
        ft <- ft %>%
          bg(i = which(!is.na(excel_wide[[pc]]) & excel_wide[[pc]] < .1 & excel_wide[[pc]] >= .05),
             j = pc, bg = "#d7eef3") %>%
          bg(i = which(!is.na(excel_wide[[pc]]) & excel_wide[[pc]] < .05),
             j = pc, bg = "#90ccde")
      }
    }
    
    # 6) treu les ratlles internes a cel·les buides (sense mergejar-les)
    transparent <- fp_border(color = "transparent", width = 0)
    
    for (j in detail_cols) {
      x <- excel_wide[[j]]
      i_empty <- if (is.numeric(x)) which(is.na(x)) else which(is.na(x) | x == "")
      if (length(i_empty) > 0) {
        ft <- border(ft, i = i_empty, j = j, border.top = transparent, part = "body")
        ft <- border(ft, i = i_empty, j = j, border.bottom = transparent, part = "body")
      }
    }
    
    term_breaks <- which(
      excel_wide$term != dplyr::lag(excel_wide$term)
    )
    
    thick_border <- fp_border(color = "black", width = 2)
    
    ft <- ft %>%
      border(
        i = term_breaks,
        j = seq_len(ncol(excel_wide)),
        border.top = thick_border,
        part = "body"
      )
    
    # --- separadors gruixuts verticals (blocs) ---
    thick_v <- fp_border(color = "black", width = 2)
    
    # columnes globals i de detall (ja les tens abans)
    global_cols <- c("term", "chisq", "df", "p_chisq")
    all_cols    <- names(excel_wide)
    detail_cols <- setdiff(all_cols, global_cols)
    
    # 1) separador gruixut després de p_chisq
    j_after_p <- match("p_chisq", all_cols)
    
    ft <- ft %>%
      border(j = j_after_p, border.right = thick_v, part = "all")
    
    # 2) separadors gruixuts entre minitaules (entre blocs de factors)
    #    - detecta prefixos (factor) en l'ordre en què apareixen a excel_wide
    prefixes <- stringr::str_extract(detail_cols, "^[^.]+") %>% unique()
    
    # per a cada prefix, troba l'última columna del seu bloc
    last_col_each_block <- purrr::map_chr(prefixes, \(p) {
      cols_p <- detail_cols[stringr::str_detect(detail_cols, paste0("^", p, "\\."))]
      dplyr::last(cols_p)
    })
    
    # índexs d'aquestes columnes
    j_block_ends <- match(last_col_each_block, all_cols)
    
    # aplica el border dret gruixut a finals de bloc
    ft <- ft %>%
      border(j = j_block_ends, border.right = thick_v, part = "all")
    
    ft %>% save_as_html(path = path, title = title)
  }
  
  if (extract) return(excel_wide)
  invisible(NULL)
}
