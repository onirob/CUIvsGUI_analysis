suppressPackageStartupMessages({
  library(ordinal)
  library(emmeans)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================================================
# SETTINGS
# ============================================================
CSV_PATH <- "reliance_clmm.csv"
OUT_DIR  <- "results_H5_reliance_clmm"

PARTICIPANT_COL <- "participant_id"
COND_COL        <- "interface_cond"
COMPLEXITY_COL  <- "task_complexity"
ITEM_COL        <- "item_key"
Y_RAW_COL       <- "reliance_value"

REF_COND          <- "chatbot"
OTHER_COND        <- "dashboard"
ITEM_LEVELS       <- c("r1", "r2", "r3")
COMPLEXITY_LEVELS <- c("easy", "mid", "hard")
REF_COMPLEXITY    <- "easy"

LINK_FUN       <- "probit"
CONF_LEVEL     <- 0.95
EMM_WEIGHTS    <- "equal"
NAGQ_SINGLE_RE <- 10

OUT_DIR_H5A <- file.path(OUT_DIR, "H5A_interface")
OUT_DIR_H5B <- file.path(OUT_DIR, "H5B_interaction")

# ============================================================
# HELPERS
# ============================================================
iqr_num <- function(x) {
  stats::IQR(x, na.rm = TRUE, type = 7)
}

save_csv <- function(x, file) {
  readr::write_csv(as_tibble(x), file)
}

make_coef_table <- function(model, conf.level = 0.95) {
  sm <- summary(model)
  tab <- as.data.frame(coef(sm))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  zcrit <- qnorm(1 - (1 - conf.level) / 2)

  tab <- tab %>%
    mutate(
      conf.low  = Estimate - zcrit * `Std. Error`,
      conf.high = Estimate + zcrit * `Std. Error`,
      component = case_when(
        str_detect(term, "\\|") ~ "threshold",
        TRUE ~ "location"
      )
    ) %>%
    relocate(term, component)

  as_tibble(tab)
}

dashboard_minus_chatbot_contrast <- function() {
  list(dashboard_minus_chatbot = c(-1, 1))
}

complexity_simple_contrasts <- function() {
  list(
    mid_minus_easy  = c(-1, 1, 0),
    hard_minus_easy = c(-1, 0, 1),
    hard_minus_mid  = c(0, -1, 1)
  )
}

did_like_contrasts <- function() {
  list(
    did_mid_vs_easy  = c(-1, 1, 0),
    did_hard_vs_easy = c(-1, 0, 1),
    did_hard_vs_mid  = c(0, -1, 1)
  )
}

standardize_emm_ci_cols <- function(x) {
  df <- as.data.frame(x)

  if (!("lower.CL" %in% names(df))) {
    if ("asymp.LCL" %in% names(df)) df$lower.CL <- df$asymp.LCL
    if ("LCL" %in% names(df)) df$lower.CL <- df$LCL
  }

  if (!("upper.CL" %in% names(df))) {
    if ("asymp.UCL" %in% names(df)) df$upper.CL <- df$asymp.UCL
    if ("UCL" %in% names(df)) df$upper.CL <- df$UCL
  }

  if (!("z.ratio" %in% names(df)) && "t.ratio" %in% names(df)) {
    df$z.ratio <- df$t.ratio
  }

  if (!("lower.CL" %in% names(df))) df$lower.CL <- NA_real_
  if (!("upper.CL" %in% names(df))) df$upper.CL <- NA_real_

  as_tibble(df)
}

add_interface_direction_flags <- function(df) {
  df %>%
    mutate(
      dashboard_higher_reliance = estimate > 0,
      chatbot_higher_reliance   = estimate < 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

add_did_direction_flags <- function(df) {
  df %>%
    mutate(
      dashboard_chatbot_gap_larger_at_higher_complexity = estimate > 0,
      dashboard_chatbot_gap_smaller_at_higher_complexity = estimate < 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

empty_random_tbl <- function() {
  tibble(
    item_key = character(),
    model = character(),
    group = character(),
    coef_name = character(),
    variance_random_effect = numeric(),
    sd_random_effect = numeric(),
    source = character()
  )
}

extract_from_stdev_vector <- function(stdev, item_key, model_label, source) {
  if (is.null(stdev) || length(stdev) == 0) {
    return(empty_random_tbl())
  }

  groups <- names(stdev)
  if (is.null(groups) || length(groups) != length(stdev) || any(groups == "")) {
    groups <- paste0("random_term_", seq_along(stdev))
  }

  tibble(
    item_key = item_key,
    model = model_label,
    group = groups,
    coef_name = "(Intercept)",
    variance_random_effect = unname(stdev)^2,
    sd_random_effect = unname(stdev),
    source = source
  )
}

extract_from_st_object <- function(ST_obj, item_key, model_label, source) {
  if (is.null(ST_obj) || length(ST_obj) == 0) {
    return(empty_random_tbl())
  }

  if (!is.list(ST_obj)) {
    ST_obj <- list(ST_obj)
  }

  groups <- names(ST_obj)
  if (is.null(groups) || length(groups) != length(ST_obj) || any(groups == "")) {
    groups <- paste0("random_term_", seq_along(ST_obj))
  }

  out <- vector("list", length(ST_obj))

  for (i in seq_along(ST_obj)) {
    block <- ST_obj[[i]]
    grp   <- groups[[i]]

    mat <- as.matrix(block)
    if (length(mat) == 0) next

    cov_mat <- mat %*% t(mat)
    sds <- sqrt(diag(cov_mat))

    coef_names <- rownames(cov_mat)
    if (is.null(coef_names) || any(coef_names == "")) {
      coef_names <- paste0("coef_", seq_along(sds))
    }

    out[[i]] <- tibble(
      item_key = item_key,
      model = model_label,
      group = grp,
      coef_name = coef_names,
      variance_random_effect = as.numeric(sds)^2,
      sd_random_effect = as.numeric(sds),
      source = source
    )
  }

  bind_rows(out)
}

extract_from_tau <- function(model, item_key, model_label, source = "exp(model$tau)") {
  if (is.null(model$tau) || length(model$tau) == 0) {
    return(empty_random_tbl())
  }

  sd_vals <- exp(model$tau)
  groups <- names(sd_vals)

  if (is.null(groups) || length(groups) != length(sd_vals) || any(groups == "")) {
    if (!is.null(model$ST) && length(model$ST) == length(sd_vals)) {
      groups <- names(model$ST)
    }
  }

  if (is.null(groups) || length(groups) != length(sd_vals) || any(groups == "")) {
    groups <- paste0("random_term_", seq_along(sd_vals))
  }

  tibble(
    item_key = item_key,
    model = model_label,
    group = groups,
    coef_name = "(Intercept)",
    variance_random_effect = unname(sd_vals)^2,
    sd_random_effect = unname(sd_vals),
    source = source
  )
}

extract_from_summary_text <- function(model, item_key, model_label, source = "summary_text") {
  txt <- capture.output(summary(model))
  if (length(txt) == 0) return(empty_random_tbl())

  start_idx <- which(str_detect(txt, "^Random effects:"))
  end_idx   <- which(str_detect(txt, "^Number of groups:"))

  if (length(start_idx) == 0 || length(end_idx) == 0) {
    return(empty_random_tbl())
  }

  start_idx <- start_idx[1]
  end_idx   <- end_idx[end_idx > start_idx][1]

  if (is.na(end_idx) || (end_idx - start_idx) < 2) {
    return(empty_random_tbl())
  }

  block <- txt[(start_idx + 2):(end_idx - 1)]
  block <- block[str_squish(block) != ""]

  if (length(block) == 0) {
    return(empty_random_tbl())
  }

  out <- lapply(block, function(line) {
    parts <- str_split(str_squish(line), "\\s+")[[1]]
    if (length(parts) < 4) return(NULL)

    var_val <- suppressWarnings(as.numeric(parts[3]))
    sd_val  <- suppressWarnings(as.numeric(parts[4]))

    if (is.na(var_val) || is.na(sd_val)) return(NULL)

    tibble(
      item_key = item_key,
      model = model_label,
      group = parts[1],
      coef_name = parts[2],
      variance_random_effect = var_val,
      sd_random_effect = sd_val,
      source = source
    )
  })

  bind_rows(out)
}

extract_random_effect_sd <- function(model, item_key, model_label) {
  out <- extract_from_stdev_vector(model$stDev, item_key, model_label, "model$stDev")
  if (nrow(out) > 0) return(out)

  sm <- summary(model)

  out <- extract_from_stdev_vector(sm$stDev, item_key, model_label, "summary(model)$stDev")
  if (nrow(out) > 0) return(out)

  out <- extract_from_st_object(model$ST, item_key, model_label, "model$ST")
  if (nrow(out) > 0) return(out)

  out <- extract_from_st_object(sm$ST, item_key, model_label, "summary(model)$ST")
  if (nrow(out) > 0) return(out)

  out <- extract_from_tau(model, item_key, model_label, "exp(model$tau)")
  if (nrow(out) > 0) return(out)

  out <- extract_from_summary_text(model, item_key, model_label, "summary_text")
  if (nrow(out) > 0) return(out)

  empty_random_tbl()
}

build_item_dataset <- function(df, item) {
  dat <- df %>%
    filter(.data[[ITEM_COL]] == item) %>%
    mutate(
      y_num = as.numeric(.data[[Y_RAW_COL]])
    )

  y_levels <- sort(unique(dat$y_num))
  if (length(y_levels) < 2) {
 stop(paste0("For item ", item, " ci are less of 2 categories observed."))
  }

  dat$y_ord <- ordered(dat$y_num, levels = y_levels)
  dat$y_ord_num <- as.integer(dat$y_ord)

  dat
}

# ============================================================
# OUTPUT DIRS
# ============================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_H5A, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_H5B, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================
if (!file.exists(CSV_PATH)) {
 stop(paste0("CSV not found: ", CSV_PATH))
}

df <- read_csv(CSV_PATH, show_col_types = FALSE)

needed_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL, ITEM_COL, Y_RAW_COL)
missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0) {
 stop("Missing these columns in the CSV: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# CLEANING
# ============================================================
df <- df %>%
  mutate(
    !!PARTICIPANT_COL := as.character(.data[[PARTICIPANT_COL]]),
    !!COND_COL        := str_to_lower(str_trim(as.character(.data[[COND_COL]]))),
    !!COMPLEXITY_COL  := str_to_lower(str_trim(as.character(.data[[COMPLEXITY_COL]]))),
    !!ITEM_COL        := str_to_lower(str_trim(as.character(.data[[ITEM_COL]]))),
    !!Y_RAW_COL       := as.numeric(.data[[Y_RAW_COL]])
  ) %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[COMPLEXITY_COL]]),
    !is.na(.data[[ITEM_COL]]),
    !is.na(.data[[Y_RAW_COL]])
  )

# participant assigned to only one interface
n_conds <- df %>%
  group_by(.data[[PARTICIPANT_COL]]) %>%
  summarise(n_cond = n_distinct(.data[[COND_COL]]), .groups = "drop")

bad_ids <- n_conds %>%
  filter(n_cond > 1) %>%
  pull(.data[[PARTICIPANT_COL]])

if (length(bad_ids) > 0) {
 message("ATTENZIONE: elimino participants assegnati to more of a interface: ",
          paste(head(bad_ids, 10), collapse = ", "))
  df <- df %>%
    filter(!(.data[[PARTICIPANT_COL]] %in% bad_ids))
}

df <- df %>%
  filter(
    .data[[COND_COL]] %in% c(REF_COND, OTHER_COND),
    .data[[COMPLEXITY_COL]] %in% COMPLEXITY_LEVELS,
    .data[[ITEM_COL]] %in% ITEM_LEVELS
  )

conds_present <- sort(unique(df[[COND_COL]]))
if (!all(c(REF_COND, OTHER_COND) %in% conds_present)) {
 stop("Levels expected in interface_cond not present. Found: ",
       paste(conds_present, collapse = ", "))
}
if (length(unique(df[[COND_COL]])) != 2) {
 stop("Mi aspetto esattamente 2 levels of interface_cond.")
}

present_complexities <- sort(unique(df[[COMPLEXITY_COL]]))
if (!all(COMPLEXITY_LEVELS %in% present_complexities)) {
 stop("Are used all and 3 the levels of task_complexity: ",
       paste(COMPLEXITY_LEVELS, collapse = ", "),
       ". Trovati: ", paste(present_complexities, collapse = ", "))
}

present_items <- sort(unique(df[[ITEM_COL]]))
if (!all(ITEM_LEVELS %in% present_items)) {
 stop("Are used all and 3 the item: ",
       paste(ITEM_LEVELS, collapse = ", "),
       ". Trovati: ", paste(present_items, collapse = ", "))
}

# enforce integer ordinal outcome
if (any(abs(df[[Y_RAW_COL]] - round(df[[Y_RAW_COL]])) > 1e-8, na.rm = TRUE)) {
 stop("reliance_value contiene values not interi. The CLMM here richiede ordinal categories discrete.")
}
df[[Y_RAW_COL]] <- as.integer(round(df[[Y_RAW_COL]]))

# duplicates
dup_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL, ITEM_COL)
dups <- duplicated(df[dup_cols]) | duplicated(df[dup_cols], fromLast = TRUE)
if (any(dups)) {
 stop("Ci are rows duplicate in the same cell logic (participant x interface x complexity x item).")
}

df[[COND_COL]]       <- factor(df[[COND_COL]], levels = c(REF_COND, OTHER_COND))
df[[COMPLEXITY_COL]] <- factor(df[[COMPLEXITY_COL]], levels = COMPLEXITY_LEVELS)
df[[ITEM_COL]]       <- factor(df[[ITEM_COL]], levels = ITEM_LEVELS)
df[[PARTICIPANT_COL]] <- factor(df[[PARTICIPANT_COL]])

# global ordinal coding for descriptives
global_y_levels <- sort(unique(df[[Y_RAW_COL]]))
df$y_ord_global <- ordered(df[[Y_RAW_COL]], levels = global_y_levels)
df$y_ord_num_global <- as.integer(df$y_ord_global)

save_csv(df, file.path(OUT_DIR, "h5_reliance_cleaned_analysis_data.csv"))

# ============================================================
# DESCRIPTIVES
# ============================================================
h5a_desc_raw <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
    .groups = "drop"
  )

h5a_desc_ord <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]]) %>%
  summarise(
    n            = n(),
    mean_class   = mean(y_ord_num_global, na.rm = TRUE),
    sd_class     = sd(y_ord_num_global, na.rm = TRUE),
    median_class = median(y_ord_num_global, na.rm = TRUE),
    iqr_class    = iqr_num(y_ord_num_global),
    .groups = "drop"
  )

h5a_n_pp <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]]) %>%
  summarise(
    n_participants = n_distinct(.data[[PARTICIPANT_COL]]),
    .groups = "drop"
  )

save_csv(h5a_desc_raw, file.path(OUT_DIR_H5A, "h5a_descriptives_raw_by_item_interface.csv"))
save_csv(h5a_desc_ord, file.path(OUT_DIR_H5A, "h5a_descriptives_ordinal_by_item_interface.csv"))
save_csv(h5a_n_pp,     file.path(OUT_DIR_H5A, "h5a_n_participants_by_item_interface.csv"))

h5b_desc_raw <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]], .data[[COMPLEXITY_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
    .groups = "drop"
  )

h5b_desc_ord <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]], .data[[COMPLEXITY_COL]]) %>%
  summarise(
    n            = n(),
    mean_class   = mean(y_ord_num_global, na.rm = TRUE),
    sd_class     = sd(y_ord_num_global, na.rm = TRUE),
    median_class = median(y_ord_num_global, na.rm = TRUE),
    iqr_class    = iqr_num(y_ord_num_global),
    .groups = "drop"
  )

h5b_n_pp <- df %>%
  group_by(.data[[ITEM_COL]], .data[[COND_COL]], .data[[COMPLEXITY_COL]]) %>%
  summarise(
    n_participants = n_distinct(.data[[PARTICIPANT_COL]]),
    .groups = "drop"
  )

save_csv(h5b_desc_raw, file.path(OUT_DIR_H5B, "h5b_descriptives_raw_by_cell.csv"))
save_csv(h5b_desc_ord, file.path(OUT_DIR_H5B, "h5b_descriptives_ordinal_by_cell.csv"))
save_csv(h5b_n_pp,     file.path(OUT_DIR_H5B, "h5b_n_participants_by_cell.csv"))

# ============================================================
# CONTRAST CODING
# ============================================================
options(contrasts = c("contr.treatment", "contr.treatment"))

# ============================================================
# H5A
# Model: y ~ interface + task_complexity + (1|participant)
# ============================================================
h5a_fit_tbl <- list()
h5a_lrt_tbl <- list()
h5a_coef_all_tbl <- list()
h5a_coef_location_tbl <- list()
h5a_random_tbl <- list()
h5a_emm_interface_latent_tbl <- list()
h5a_emm_interface_meanclass_tbl <- list()
h5a_overall_contrast_latent_tbl <- list()
h5a_overall_contrast_meanclass_tbl <- list()
h5a_prob_interface_tbl <- list()

for (item in ITEM_LEVELS) {
  dat <- build_item_dataset(df, item)

  form_add <- as.formula(
    paste("y_ord ~", COND_COL, "+", COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")", sep = "")
  )

  form_no_interface <- as.formula(
    paste("y_ord ~", COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")", sep = "")
  )

  form_no_complexity <- as.formula(
    paste("y_ord ~", COND_COL, "+ (1|", PARTICIPANT_COL, ")", sep = "")
  )

  m_add <- clmm(
    formula = form_add,
    data    = dat,
    link    = LINK_FUN,
    Hess    = TRUE,
    nAGQ    = NAGQ_SINGLE_RE
  )

  m_no_interface <- clmm(
    formula = form_no_interface,
    data    = dat,
    link    = LINK_FUN,
    Hess    = TRUE,
    nAGQ    = NAGQ_SINGLE_RE
  )

  m_no_complexity <- clmm(
    formula = form_no_complexity,
    data    = dat,
    link    = LINK_FUN,
    Hess    = TRUE,
    nAGQ    = NAGQ_SINGLE_RE
  )

  capture.output(
    summary(m_add),
    file = file.path(OUT_DIR_H5A, paste0("h5a_", item, "_additive_model_summary.txt"))
  )

  fit_tbl_item <- tibble(
    item_key = item,
    model  = c("no_interface", "no_complexity", "additive"),
    nobs   = c(nobs(m_no_interface), nobs(m_no_complexity), nobs(m_add)),
    logLik = c(as.numeric(logLik(m_no_interface)), as.numeric(logLik(m_no_complexity)), as.numeric(logLik(m_add))),
    AIC    = c(AIC(m_no_interface), AIC(m_no_complexity), AIC(m_add)),
    BIC    = c(BIC(m_no_interface), BIC(m_no_complexity), BIC(m_add))
  )
  h5a_fit_tbl[[item]] <- fit_tbl_item

  lrt_interface <- as.data.frame(anova(m_no_interface, m_add))
  lrt_interface$item_key <- item
  lrt_interface$effect <- "interface_main_effect_adjusted_for_complexity"

  lrt_complexity <- as.data.frame(anova(m_no_complexity, m_add))
  lrt_complexity$item_key <- item
  lrt_complexity$effect <- "complexity_main_effect_adjusted_for_interface"

  h5a_lrt_tbl[[item]] <- bind_rows(
    as_tibble(lrt_interface),
    as_tibble(lrt_complexity)
  )

  coef_add_all <- make_coef_table(m_add, conf.level = CONF_LEVEL) %>%
    mutate(item_key = item, .before = 1)

  coef_add_location <- coef_add_all %>%
    filter(component == "location")

  h5a_coef_all_tbl[[item]] <- coef_add_all
  h5a_coef_location_tbl[[item]] <- coef_add_location

  rand_add <- extract_random_effect_sd(m_add, item, "additive")
  h5a_random_tbl[[item]] <- rand_add

  emm_overall_interface_latent <- emmeans(
    m_add,
    specs   = as.formula(paste("~", COND_COL)),
    mode    = "latent",
    weights = EMM_WEIGHTS
  )

  emm_overall_interface_meanclass <- emmeans(
    m_add,
    specs   = as.formula(paste("~", COND_COL)),
    mode    = "mean.class",
    weights = EMM_WEIGHTS
  )

  emm_interface_latent_tbl <- summary(
    emm_overall_interface_latent,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  emm_interface_meanclass_tbl <- summary(
    emm_overall_interface_meanclass,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5a_emm_interface_latent_tbl[[item]] <- emm_interface_latent_tbl
  h5a_emm_interface_meanclass_tbl[[item]] <- emm_interface_meanclass_tbl

  overall_interface_contrast_latent <- summary(
    contrast(
      emm_overall_interface_latent,
      method = dashboard_minus_chatbot_contrast()
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_interface_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  overall_interface_contrast_meanclass <- summary(
    contrast(
      emm_overall_interface_meanclass,
      method = dashboard_minus_chatbot_contrast()
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_interface_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  h5a_overall_contrast_latent_tbl[[item]] <- overall_interface_contrast_latent
  h5a_overall_contrast_meanclass_tbl[[item]] <- overall_interface_contrast_meanclass

  emm_prob_interface <- emmeans(
    m_add,
    specs   = as.formula(paste("~ y_ord |", COND_COL)),
    mode    = "prob",
    weights = EMM_WEIGHTS
  )

  prob_interface_tbl <- confint(
    emm_prob_interface,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5a_prob_interface_tbl[[item]] <- prob_interface_tbl
}

h5a_model_fit <- bind_rows(h5a_fit_tbl)
h5a_lrt       <- bind_rows(h5a_lrt_tbl)
h5a_coef_all  <- bind_rows(h5a_coef_all_tbl)
h5a_coef_loc  <- bind_rows(h5a_coef_location_tbl)
h5a_random    <- bind_rows(h5a_random_tbl)
h5a_emm_latent <- bind_rows(h5a_emm_interface_latent_tbl)
h5a_emm_meanclass <- bind_rows(h5a_emm_interface_meanclass_tbl)
h5a_contrast_latent <- bind_rows(h5a_overall_contrast_latent_tbl)
h5a_contrast_meanclass <- bind_rows(h5a_overall_contrast_meanclass_tbl)
h5a_prob_interface <- bind_rows(h5a_prob_interface_tbl)

save_csv(h5a_model_fit, file.path(OUT_DIR_H5A, "h5a_model_fit.csv"))
save_csv(h5a_lrt, file.path(OUT_DIR_H5A, "h5a_omnibus_lrt.csv"))
save_csv(
  h5a_lrt %>% filter(effect == "interface_main_effect_adjusted_for_complexity"),
  file.path(OUT_DIR_H5A, "h5a_interface_lrt.csv")
)
save_csv(
  h5a_lrt %>% filter(effect == "complexity_main_effect_adjusted_for_interface"),
  file.path(OUT_DIR_H5A, "h5a_complexity_lrt.csv")
)

capture.output(h5a_lrt, file = file.path(OUT_DIR_H5A, "h5a_omnibus_lrt.txt"))

save_csv(h5a_coef_all, file.path(OUT_DIR_H5A, "h5a_additive_coef_all.csv"))
save_csv(h5a_coef_loc, file.path(OUT_DIR_H5A, "h5a_additive_coef_location.csv"))
save_csv(h5a_random, file.path(OUT_DIR_H5A, "h5a_random_effect_sd.csv"))
save_csv(h5a_emm_latent, file.path(OUT_DIR_H5A, "h5a_emm_interface_latent.csv"))
save_csv(h5a_contrast_latent, file.path(OUT_DIR_H5A, "h5a_overall_interface_contrast_latent.csv"))
save_csv(h5a_emm_meanclass, file.path(OUT_DIR_H5A, "h5a_emm_interface_meanclass.csv"))
save_csv(h5a_contrast_meanclass, file.path(OUT_DIR_H5A, "h5a_overall_interface_contrast_meanclass.csv"))
save_csv(h5a_prob_interface, file.path(OUT_DIR_H5A, "h5a_category_probabilities_by_interface.csv"))

writeLines(
  c(
    "Random-effect extraction diagnostics - H5A",
    paste0("rows: ", nrow(h5a_random)),
    if (nrow(h5a_random) > 0) paste0("source(s): ", paste(unique(h5a_random$source), collapse = ", ")) else "source(s): none"
  ),
  con = file.path(OUT_DIR_H5A, "h5a_random_effect_sd_diagnostics.txt")
)

# ============================================================
# H5B
# Model: y ~ interface * task_complexity + (1|participant)
# ============================================================
h5b_fit_tbl <- list()
h5b_lrt_tbl <- list()
h5b_coef_all_tbl <- list()
h5b_coef_location_tbl <- list()
h5b_coef_main_tbl <- list()
h5b_coef_interaction_tbl <- list()
h5b_random_tbl <- list()

h5b_emm_interface_by_comp_latent_tbl <- list()
h5b_simple_interface_by_comp_latent_tbl <- list()
h5b_emm_interface_by_comp_meanclass_tbl <- list()
h5b_simple_interface_by_comp_meanclass_tbl <- list()

h5b_emm_comp_by_interface_latent_tbl <- list()
h5b_simple_comp_by_interface_latent_tbl <- list()
h5b_emm_comp_by_interface_meanclass_tbl <- list()
h5b_simple_comp_by_interface_meanclass_tbl <- list()

h5b_did_like_latent_tbl <- list()
h5b_did_like_meanclass_tbl <- list()
h5b_all_cells_latent_tbl <- list()
h5b_all_cells_meanclass_tbl <- list()
h5b_prob_cells_tbl <- list()

for (item in ITEM_LEVELS) {
  dat <- build_item_dataset(df, item)

  form_add <- as.formula(
    paste("y_ord ~", COND_COL, "+", COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")", sep = "")
  )

  form_int <- as.formula(
    paste("y_ord ~", COND_COL, "*", COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")", sep = "")
  )

  m_add <- clmm(
    formula = form_add,
    data    = dat,
    link    = LINK_FUN,
    Hess    = TRUE,
    nAGQ    = NAGQ_SINGLE_RE
  )

  m_int <- clmm(
    formula = form_int,
    data    = dat,
    link    = LINK_FUN,
    Hess    = TRUE,
    nAGQ    = NAGQ_SINGLE_RE
  )

  capture.output(
    summary(m_add),
    file = file.path(OUT_DIR_H5B, paste0("h5b_", item, "_additive_model_summary.txt"))
  )
  capture.output(
    summary(m_int),
    file = file.path(OUT_DIR_H5B, paste0("h5b_", item, "_interaction_model_summary.txt"))
  )

  fit_tbl_item <- tibble(
    item_key = item,
    model  = c("additive", "interaction"),
    nobs   = c(nobs(m_add), nobs(m_int)),
    logLik = c(as.numeric(logLik(m_add)), as.numeric(logLik(m_int))),
    AIC    = c(AIC(m_add), AIC(m_int)),
    BIC    = c(BIC(m_add), BIC(m_int))
  )
  h5b_fit_tbl[[item]] <- fit_tbl_item

  lrt_interaction <- as.data.frame(anova(m_add, m_int))
  lrt_interaction$item_key <- item
  lrt_interaction$effect <- "interface_x_complexity_interaction"

  h5b_lrt_tbl[[item]] <- as_tibble(lrt_interaction)

  coef_int_all <- make_coef_table(m_int, conf.level = CONF_LEVEL) %>%
    mutate(item_key = item, .before = 1)

  coef_int_location <- coef_int_all %>%
    filter(component == "location")

  coef_int_main_only <- coef_int_location %>%
    filter(!str_detect(term, ":"))

  coef_int_interaction_only <- coef_int_location %>%
    filter(str_detect(term, ":"))

  h5b_coef_all_tbl[[item]] <- coef_int_all
  h5b_coef_location_tbl[[item]] <- coef_int_location
  h5b_coef_main_tbl[[item]] <- coef_int_main_only
  h5b_coef_interaction_tbl[[item]] <- coef_int_interaction_only

  rand_int <- extract_random_effect_sd(m_int, item, "interaction")
  h5b_random_tbl[[item]] <- rand_int

  emm_interface_by_comp_latent <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COND_COL, "|", COMPLEXITY_COL)),
    mode    = "latent",
    weights = EMM_WEIGHTS
  )

  gap_by_comp_latent <- contrast(
    emm_interface_by_comp_latent,
    method = dashboard_minus_chatbot_contrast()
  )

  simple_effects_interface_latent <- summary(
    gap_by_comp_latent,
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_interface_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  h5b_emm_interface_by_comp_latent_tbl[[item]] <- summary(
    emm_interface_by_comp_latent,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_simple_interface_by_comp_latent_tbl[[item]] <- simple_effects_interface_latent

  emm_interface_by_comp_meanclass <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COND_COL, "|", COMPLEXITY_COL)),
    mode    = "mean.class",
    weights = EMM_WEIGHTS
  )

  gap_by_comp_meanclass <- contrast(
    emm_interface_by_comp_meanclass,
    method = dashboard_minus_chatbot_contrast()
  )

  simple_effects_interface_meanclass <- summary(
    gap_by_comp_meanclass,
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_interface_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  h5b_emm_interface_by_comp_meanclass_tbl[[item]] <- summary(
    emm_interface_by_comp_meanclass,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_simple_interface_by_comp_meanclass_tbl[[item]] <- simple_effects_interface_meanclass

  emm_complexity_by_interface_latent <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COMPLEXITY_COL, "|", COND_COL)),
    mode    = "latent",
    weights = EMM_WEIGHTS
  )

  simple_effects_complexity_latent <- summary(
    contrast(
      emm_complexity_by_interface_latent,
      method = complexity_simple_contrasts(),
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_emm_comp_by_interface_latent_tbl[[item]] <- summary(
    emm_complexity_by_interface_latent,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_simple_comp_by_interface_latent_tbl[[item]] <- simple_effects_complexity_latent

  emm_complexity_by_interface_meanclass <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COMPLEXITY_COL, "|", COND_COL)),
    mode    = "mean.class",
    weights = EMM_WEIGHTS
  )

  simple_effects_complexity_meanclass <- summary(
    contrast(
      emm_complexity_by_interface_meanclass,
      method = complexity_simple_contrasts(),
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_emm_comp_by_interface_meanclass_tbl[[item]] <- summary(
    emm_complexity_by_interface_meanclass,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_simple_comp_by_interface_meanclass_tbl[[item]] <- simple_effects_complexity_meanclass

  did_like_latent <- summary(
    contrast(
      gap_by_comp_latent,
      method = did_like_contrasts(),
      by = NULL,
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_did_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  did_like_meanclass <- summary(
    contrast(
      gap_by_comp_meanclass,
      method = did_like_contrasts(),
      by = NULL,
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE),
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    add_did_direction_flags() %>%
    mutate(item_key = item, .before = 1)

  h5b_did_like_latent_tbl[[item]] <- did_like_latent
  h5b_did_like_meanclass_tbl[[item]] <- did_like_meanclass

  emm_all_cells_latent <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COND_COL, "*", COMPLEXITY_COL)),
    mode    = "latent",
    weights = EMM_WEIGHTS
  )

  emm_all_cells_meanclass <- emmeans(
    m_int,
    specs   = as.formula(paste("~", COND_COL, "*", COMPLEXITY_COL)),
    mode    = "mean.class",
    weights = EMM_WEIGHTS
  )

  h5b_all_cells_latent_tbl[[item]] <- summary(
    emm_all_cells_latent,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_all_cells_meanclass_tbl[[item]] <- summary(
    emm_all_cells_meanclass,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  emm_prob_cells <- emmeans(
    m_int,
    specs   = as.formula(paste("~ y_ord |", COND_COL, "*", COMPLEXITY_COL)),
    mode    = "prob",
    weights = EMM_WEIGHTS
  )

  prob_cells_tbl <- confint(
    emm_prob_cells,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5b_prob_cells_tbl[[item]] <- prob_cells_tbl
}

h5b_model_fit <- bind_rows(h5b_fit_tbl)
h5b_lrt <- bind_rows(h5b_lrt_tbl)
h5b_coef_all <- bind_rows(h5b_coef_all_tbl)
h5b_coef_loc <- bind_rows(h5b_coef_location_tbl)
h5b_coef_main <- bind_rows(h5b_coef_main_tbl)
h5b_coef_interaction <- bind_rows(h5b_coef_interaction_tbl)
h5b_random <- bind_rows(h5b_random_tbl)

h5b_emm_interface_by_comp_latent <- bind_rows(h5b_emm_interface_by_comp_latent_tbl)
h5b_simple_interface_by_comp_latent <- bind_rows(h5b_simple_interface_by_comp_latent_tbl)
h5b_emm_interface_by_comp_meanclass <- bind_rows(h5b_emm_interface_by_comp_meanclass_tbl)
h5b_simple_interface_by_comp_meanclass <- bind_rows(h5b_simple_interface_by_comp_meanclass_tbl)

h5b_emm_comp_by_interface_latent <- bind_rows(h5b_emm_comp_by_interface_latent_tbl)
h5b_simple_comp_by_interface_latent <- bind_rows(h5b_simple_comp_by_interface_latent_tbl)
h5b_emm_comp_by_interface_meanclass <- bind_rows(h5b_emm_comp_by_interface_meanclass_tbl)
h5b_simple_comp_by_interface_meanclass <- bind_rows(h5b_simple_comp_by_interface_meanclass_tbl)

h5b_did_like_latent <- bind_rows(h5b_did_like_latent_tbl)
h5b_did_like_meanclass <- bind_rows(h5b_did_like_meanclass_tbl)
h5b_all_cells_latent <- bind_rows(h5b_all_cells_latent_tbl)
h5b_all_cells_meanclass <- bind_rows(h5b_all_cells_meanclass_tbl)
h5b_prob_cells <- bind_rows(h5b_prob_cells_tbl)

save_csv(h5b_model_fit, file.path(OUT_DIR_H5B, "h5b_model_fit.csv"))
save_csv(h5b_lrt, file.path(OUT_DIR_H5B, "h5b_interaction_lrt.csv"))
capture.output(h5b_lrt, file = file.path(OUT_DIR_H5B, "h5b_interaction_lrt.txt"))

save_csv(h5b_coef_all, file.path(OUT_DIR_H5B, "h5b_interaction_coef_all.csv"))
save_csv(h5b_coef_loc, file.path(OUT_DIR_H5B, "h5b_interaction_coef_location.csv"))
save_csv(h5b_coef_main, file.path(OUT_DIR_H5B, "h5b_interaction_main_terms_only.csv"))
save_csv(h5b_coef_interaction, file.path(OUT_DIR_H5B, "h5b_interaction_terms_only.csv"))

save_csv(h5b_random, file.path(OUT_DIR_H5B, "h5b_random_effect_sd.csv"))
writeLines(
  c(
    "Random-effect extraction diagnostics - H5B",
    paste0("rows: ", nrow(h5b_random)),
    if (nrow(h5b_random) > 0) paste0("source(s): ", paste(unique(h5b_random$source), collapse = ", ")) else "source(s): none"
  ),
  con = file.path(OUT_DIR_H5B, "h5b_random_effect_sd_diagnostics.txt")
)

save_csv(
  h5b_emm_interface_by_comp_latent,
  file.path(OUT_DIR_H5B, "h5b_emm_interface_by_complexity_latent.csv")
)
save_csv(
  h5b_simple_interface_by_comp_latent,
  file.path(OUT_DIR_H5B, "h5b_simple_effects_interface_by_complexity_latent.csv")
)
save_csv(
  h5b_emm_interface_by_comp_meanclass,
  file.path(OUT_DIR_H5B, "h5b_emm_interface_by_complexity_meanclass.csv")
)
save_csv(
  h5b_simple_interface_by_comp_meanclass,
  file.path(OUT_DIR_H5B, "h5b_simple_effects_interface_by_complexity_meanclass.csv")
)

save_csv(
  h5b_emm_comp_by_interface_latent,
  file.path(OUT_DIR_H5B, "h5b_emm_complexity_by_interface_latent.csv")
)
save_csv(
  h5b_simple_comp_by_interface_latent,
  file.path(OUT_DIR_H5B, "h5b_simple_effects_complexity_by_interface_latent.csv")
)
save_csv(
  h5b_emm_comp_by_interface_meanclass,
  file.path(OUT_DIR_H5B, "h5b_emm_complexity_by_interface_meanclass.csv")
)
save_csv(
  h5b_simple_comp_by_interface_meanclass,
  file.path(OUT_DIR_H5B, "h5b_simple_effects_complexity_by_interface_meanclass.csv")
)

save_csv(
  h5b_did_like_latent,
  file.path(OUT_DIR_H5B, "h5b_did_like_latent.csv")
)
save_csv(
  h5b_did_like_meanclass,
  file.path(OUT_DIR_H5B, "h5b_did_like_meanclass.csv")
)

save_csv(
  h5b_all_cells_latent,
  file.path(OUT_DIR_H5B, "h5b_all_cells_latent.csv")
)
save_csv(
  h5b_all_cells_meanclass,
  file.path(OUT_DIR_H5B, "h5b_all_cells_meanclass.csv")
)

save_csv(
  h5b_prob_cells,
  file.path(OUT_DIR_H5B, "h5b_category_probabilities_by_cell.csv")
)

cat("\n============================================================\n")
cat("H5 reliance CLMM completata.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("============================================================\n\n")

cat("Riferimenti dei fattori:\n")
cat("- interface_cond reference =", REF_COND, "\n")
cat("- task_complexity reference =", REF_COMPLEXITY, "\n\n")

cat("File principali H5A:\n")
cat("- H5A_interface/h5a_interface_lrt.csv\n")
cat("- H5A_interface/h5a_additive_coef_location.csv\n")
cat("- H5A_interface/h5a_overall_interface_contrast_latent.csv\n")
cat("- H5A_interface/h5a_category_probabilities_by_interface.csv\n\n")

cat("File principali H5B:\n")
cat("- H5B_interaction/h5b_interaction_lrt.csv\n")
cat("- H5B_interaction/h5b_interaction_terms_only.csv\n")
cat("- H5B_interaction/h5b_simple_effects_interface_by_complexity_latent.csv\n")
cat("- H5B_interaction/h5b_simple_effects_complexity_by_interface_latent.csv\n")
cat("- H5B_interaction/h5b_did_like_latent.csv\n")
cat("- H5B_interaction/h5b_category_probabilities_by_cell.csv\n\n")

cat("Direction of the contrasts of interface: dashboard_minus_chatbot\n")
cat("estimate > 0 => dashboard higher reliance => chatbot lower reliance\n\n")

cat("Direction of the DiD-like:\n")
cat("estimate > 0 => the dashboard–chatbot gap is larger at the level of complexity higher of the contrast\n")