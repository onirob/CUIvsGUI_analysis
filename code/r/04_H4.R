suppressPackageStartupMessages({
  library(ordinal)
  library(emmeans)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================================================
# H4a - Task complexity -> Mental Workload (ordinal)
#
# HYPOTHESIS 4(a):
#   Higher task complexity is associated with higher MWL.
#
# MODEL STRATEGY
#   - cumulative link mixed model (CLMM)
#   - ordinalized MWL outcome
#   - task_complexity as PRIMARY predictor
#   - interface_cond as nuisance additive covariate
#   - participant random intercept
#
# PRIMARY TEST
#   - omnibus LRT: model with task_complexity vs model without it
# ============================================================

# ============================================================
# PATHS RELATIVE TO SCRIPT DIRECTORY
# ============================================================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  return(normalizePath(getwd()))
}

SCRIPT_DIR <- get_script_dir()

CSV_PATH <- file.path(SCRIPT_DIR, "nasa_tlx_clmm.csv")
OUT_DIR  <- file.path(SCRIPT_DIR, "results_H4")

# ============================================================
# SETTINGS
# ============================================================
PARTICIPANT_COL <- "participant_id"
COND_COL        <- "interface_cond"
COMPLEXITY_COL  <- "task_complexity"
Y_RAW_COL       <- "tlx_mean"

REF_COND          <- "chatbot"
OTHER_COND        <- "dashboard"
COMPLEXITY_LEVELS <- c("easy", "mid", "hard")
REF_COMPLEXITY    <- "easy"

ROUND_TO       <- 5
LINK_FUN       <- "probit"
CONF_LEVEL     <- 0.95
EMM_WEIGHTS    <- "equal"
NAGQ_SINGLE_RE <- 10

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

complexity_simple_contrasts <- function() {
  list(
    mid_minus_easy  = c(-1, 1, 0),
    hard_minus_easy = c(-1, 0, 1),
    hard_minus_mid  = c(0, -1, 1)
  )
}

chatbot_minus_dashboard_contrast <- function() {
  list(chatbot_minus_dashboard = c(1, -1))
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

add_complexity_direction_flags <- function(df) {
  df %>%
    mutate(
      supports_h4_higher_mwl = estimate > 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

add_interface_direction_flags <- function(df) {
  df %>%
    mutate(
      supports_chatbot_lower_mwl = estimate < 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

empty_random_tbl <- function() {
  tibble(
    model = character(),
    group = character(),
    coef_name = character(),
    variance_random_effect = numeric(),
    sd_random_effect = numeric(),
    source = character()
  )
}

extract_from_stdev_vector <- function(stdev, model_label, source) {
  if (is.null(stdev) || length(stdev) == 0) {
    return(empty_random_tbl())
  }

  groups <- names(stdev)
  if (is.null(groups) || length(groups) != length(stdev) || any(groups == "")) {
    groups <- paste0("random_term_", seq_along(stdev))
  }

  tibble(
    model = model_label,
    group = groups,
    coef_name = "(Intercept)",
    variance_random_effect = unname(stdev)^2,
    sd_random_effect = unname(stdev),
    source = source
  )
}

extract_from_st_object <- function(ST_obj, model_label, source) {
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

extract_from_tau <- function(model, model_label, source = "exp(model$tau)") {
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
    model = model_label,
    group = groups,
    coef_name = "(Intercept)",
    variance_random_effect = unname(sd_vals)^2,
    sd_random_effect = unname(sd_vals),
    source = source
  )
}

extract_from_summary_text <- function(model, model_label, source = "summary_text") {
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

extract_random_effect_sd <- function(model, model_label) {
  out <- extract_from_stdev_vector(model$stDev, model_label, "model$stDev")
  if (nrow(out) > 0) return(out)

  sm <- summary(model)

  out <- extract_from_stdev_vector(sm$stDev, model_label, "summary(model)$stDev")
  if (nrow(out) > 0) return(out)

  out <- extract_from_st_object(model$ST, model_label, "model$ST")
  if (nrow(out) > 0) return(out)

  out <- extract_from_st_object(sm$ST, model_label, "summary(model)$ST")
  if (nrow(out) > 0) return(out)

  out <- extract_from_tau(model, model_label, "exp(model$tau)")
  if (nrow(out) > 0) return(out)

  out <- extract_from_summary_text(model, model_label, "summary_text")
  if (nrow(out) > 0) return(out)

  empty_random_tbl()
}

# ============================================================
# OUTPUT DIR
# ============================================================
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================
if (!file.exists(CSV_PATH)) {
 stop(paste0("CSV not found: ", CSV_PATH))
}

df <- read_csv(CSV_PATH, show_col_types = FALSE)

needed_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL, Y_RAW_COL)
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
    !!Y_RAW_COL       := as.numeric(.data[[Y_RAW_COL]])
  ) %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[COMPLEXITY_COL]]),
    !is.na(.data[[Y_RAW_COL]])
  )

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

conds_present <- sort(unique(df[[COND_COL]]))
if (!all(c(REF_COND, OTHER_COND) %in% conds_present)) {
 stop("Levels expected in interface_cond not present. Found: ",
       paste(conds_present, collapse = ", "))
}
if (length(unique(df[[COND_COL]])) != 2) {
 stop("Mi aspetto esattamente 2 levels of interface_cond.")
}

df[[COND_COL]] <- factor(df[[COND_COL]], levels = c(REF_COND, OTHER_COND))

present_complexities <- sort(unique(df[[COMPLEXITY_COL]]))
if (!all(COMPLEXITY_LEVELS %in% present_complexities)) {
 stop("Are used all the levels of task_complexity: ",
       paste(COMPLEXITY_LEVELS, collapse = ", "),
       ". Trovati: ", paste(present_complexities, collapse = ", "))
}

df[[COMPLEXITY_COL]] <- factor(df[[COMPLEXITY_COL]], levels = COMPLEXITY_LEVELS)

dup_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL)
dups <- duplicated(df[dup_cols]) | duplicated(df[dup_cols], fromLast = TRUE)
if (any(dups)) {
 stop("Ci are rows duplicate in the same cell logic. Controlla the CSV.")
}

# ============================================================
# ORDINALIZATION OUTCOME
# ============================================================
df <- df %>%
  mutate(
    tlx_round = round(.data[[Y_RAW_COL]] / ROUND_TO) * ROUND_TO
  )

y_levels <- sort(unique(df$tlx_round))
if (length(y_levels) < 3) {
 stop("After l'ordinalization at least are required 3 categories. Changes ROUND_TO.")
}

df$y_ord <- ordered(df$tlx_round, levels = y_levels)
df$y_ord_num <- as.integer(df$y_ord)
df[[PARTICIPANT_COL]] <- factor(df[[PARTICIPANT_COL]])

save_csv(df, file.path(OUT_DIR, "h4a_cleaned_analysis_data.csv"))

# ============================================================
# DESCRIPTIVES
# ============================================================
desc_raw_by_complexity <- df %>%
  group_by(.data[[COMPLEXITY_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
    .groups = "drop"
  )

desc_raw_by_cell <- df %>%
  group_by(.data[[COND_COL]], .data[[COMPLEXITY_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
    .groups = "drop"
  )

desc_ord_by_complexity <- df %>%
  group_by(.data[[COMPLEXITY_COL]]) %>%
  summarise(
    n            = n(),
    mean_class   = mean(y_ord_num, na.rm = TRUE),
    sd_class     = sd(y_ord_num, na.rm = TRUE),
    median_class = median(y_ord_num, na.rm = TRUE),
    iqr_class    = iqr_num(y_ord_num),
    .groups = "drop"
  )

desc_ord_by_cell <- df %>%
  group_by(.data[[COND_COL]], .data[[COMPLEXITY_COL]]) %>%
  summarise(
    n            = n(),
    mean_class   = mean(y_ord_num, na.rm = TRUE),
    sd_class     = sd(y_ord_num, na.rm = TRUE),
    median_class = median(y_ord_num, na.rm = TRUE),
    iqr_class    = iqr_num(y_ord_num),
    .groups = "drop"
  )

n_participants_by_interface <- df %>%
  distinct(.data[[PARTICIPANT_COL]], .data[[COND_COL]]) %>%
  count(.data[[COND_COL]], name = "n_participants")

save_csv(desc_raw_by_complexity, file.path(OUT_DIR, "h4a_descriptives_raw_by_complexity.csv"))
save_csv(desc_raw_by_cell,       file.path(OUT_DIR, "h4a_descriptives_raw_by_cell.csv"))
save_csv(desc_ord_by_complexity, file.path(OUT_DIR, "h4a_descriptives_ordinal_by_complexity.csv"))
save_csv(desc_ord_by_cell,       file.path(OUT_DIR, "h4a_descriptives_ordinal_by_cell.csv"))
save_csv(n_participants_by_interface, file.path(OUT_DIR, "h4a_n_participants_by_interface.csv"))

# ============================================================
# MODELS
# ============================================================
form_null <- as.formula(
  paste("y_ord ~", COND_COL, "+ (1|", PARTICIPANT_COL, ")")
)

form_add <- as.formula(
  paste("y_ord ~", COND_COL, "+", COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")")
)

options(contrasts = c("contr.treatment", "contr.treatment"))

m_null <- clmm(
  formula = form_null,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = NAGQ_SINGLE_RE
)

m_add <- clmm(
  formula = form_add,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = NAGQ_SINGLE_RE
)

capture.output(summary(m_null), file = file.path(OUT_DIR, "h4a_null_model_summary.txt"))
capture.output(summary(m_add),  file = file.path(OUT_DIR, "h4a_additive_model_summary.txt"))

# ============================================================
# MODEL FIT / LRT
# ============================================================
fit_tbl <- tibble(
  model  = c("null_interface_only", "additive_complexity"),
  nobs   = c(nobs(m_null), nobs(m_add)),
  logLik = c(as.numeric(logLik(m_null)), as.numeric(logLik(m_add))),
  AIC    = c(AIC(m_null), AIC(m_add)),
  BIC    = c(BIC(m_null), BIC(m_add))
)
save_csv(fit_tbl, file.path(OUT_DIR, "h4a_model_fit.csv"))

lrt_tbl <- as.data.frame(anova(m_null, m_add))
save_csv(lrt_tbl, file.path(OUT_DIR, "h4a_complexity_lrt.csv"))
capture.output(anova(m_null, m_add), file = file.path(OUT_DIR, "h4a_complexity_lrt.txt"))

# ============================================================
# COEFFICIENTS
# ============================================================
coef_null_all <- make_coef_table(m_null, conf.level = CONF_LEVEL)
coef_add_all  <- make_coef_table(m_add,  conf.level = CONF_LEVEL)

coef_null_location <- coef_null_all %>% filter(component == "location")
coef_add_location  <- coef_add_all  %>% filter(component == "location")

save_csv(coef_null_all,      file.path(OUT_DIR, "h4a_null_coef_all.csv"))
save_csv(coef_add_all,       file.path(OUT_DIR, "h4a_additive_coef_all.csv"))
save_csv(coef_null_location, file.path(OUT_DIR, "h4a_null_coef_location.csv"))
save_csv(coef_add_location,  file.path(OUT_DIR, "h4a_additive_coef_location.csv"))

# ============================================================
# RANDOM EFFECT SD
# ============================================================
rand_null <- extract_random_effect_sd(m_null, "null_interface_only")
rand_add  <- extract_random_effect_sd(m_add,  "additive_complexity")
rand_tbl  <- bind_rows(rand_null, rand_add)

save_csv(rand_tbl, file.path(OUT_DIR, "h4a_random_effect_sd.csv"))

# ============================================================
# EMM: COMPLEXITY (PRIMARY FOLLOW-UP)
# ============================================================
emm_complexity_latent <- emmeans(
  m_add,
  specs   = ~ task_complexity,
  mode    = "latent",
  weights = EMM_WEIGHTS
)

contrasts_complexity_latent <- summary(
  contrast(
    emm_complexity_latent,
    method = complexity_simple_contrasts(),
    adjust = "holm"
  ),
  infer = c(TRUE, TRUE),
  level = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_complexity_direction_flags()

emm_complexity_meanclass <- emmeans(
  m_add,
  specs   = ~ task_complexity,
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

contrasts_complexity_meanclass <- summary(
  contrast(
    emm_complexity_meanclass,
    method = complexity_simple_contrasts(),
    adjust = "holm"
  ),
  infer = c(TRUE, TRUE),
  level = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_complexity_direction_flags()

save_csv(
  summary(emm_complexity_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h4a_emm_complexity_latent.csv")
)

save_csv(
  contrasts_complexity_latent,
  file.path(OUT_DIR, "h4a_contrasts_complexity_latent.csv")
)

save_csv(
  summary(emm_complexity_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h4a_emm_complexity_meanclass.csv")
)

save_csv(
  contrasts_complexity_meanclass,
  file.path(OUT_DIR, "h4a_contrasts_complexity_meanclass.csv")
)

# ============================================================
# EMM: INTERFACE (EXPLORATORY)
# ============================================================
emm_interface_latent <- emmeans(
  m_add,
  specs   = ~ interface_cond,
  mode    = "latent",
  weights = EMM_WEIGHTS
)

overall_interface_contrast_latent <- summary(
  contrast(
    emm_interface_latent,
    method = chatbot_minus_dashboard_contrast()
  ),
  infer = c(TRUE, TRUE),
  level = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_interface_direction_flags()

save_csv(
  summary(emm_interface_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h4a_emm_interface_latent.csv")
)

save_csv(
  overall_interface_contrast_latent,
  file.path(OUT_DIR, "h4a_overall_interface_contrast_latent_from_additive_model.csv")
)

# ============================================================
# ALL CELL MEANS
# ============================================================
emm_all_cells_latent <- emmeans(
  m_add,
  specs   = ~ interface_cond * task_complexity,
  mode    = "latent",
  weights = EMM_WEIGHTS
)

emm_all_cells_meanclass <- emmeans(
  m_add,
  specs   = ~ interface_cond * task_complexity,
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

save_csv(
  summary(emm_all_cells_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h4a_all_cells_latent.csv")
)

save_csv(
  summary(emm_all_cells_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h4a_all_cells_meanclass.csv")
)

# ============================================================
# CATEGORY PROBABILITIES BY CELL
# ============================================================
emm_prob_cells <- emmeans(
  m_add,
  specs   = ~ y_ord | interface_cond * task_complexity,
  mode    = "prob",
  weights = EMM_WEIGHTS
)

prob_cells_tbl <- confint(emm_prob_cells, level = CONF_LEVEL) %>%
  standardize_emm_ci_cols()

save_csv(prob_cells_tbl, file.path(OUT_DIR, "h4a_category_probabilities_by_cell.csv"))

# ============================================================
# FINAL CONSOLE NOTES
# ============================================================
cat("\n============================================================\n")
cat("H4a additive complexity model completato.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("============================================================\n\n")

cat("Riferimenti dei fattori:\n")
cat("- interface_cond reference =", REF_COND, "\n")
cat("- task_complexity reference =", REF_COMPLEXITY, "\n\n")

cat("Interpretation key of the coefficients of the additive model:\n")
cat("- interface_conddashboard = dashboard vs chatbot, adjusted for complexity\n")
cat("- task_complexitymid = mid vs easy, adjusted for interface\n")
cat("- task_complexityhard = hard vs easy, adjusted for interface\n\n")