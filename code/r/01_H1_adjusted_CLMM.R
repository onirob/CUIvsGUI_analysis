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
CSV_PATH <- "nasa_tlx_clmm.csv"
OUT_DIR  <- "results_H1_adjusted"

PARTICIPANT_COL <- "participant_id"
COND_COL        <- "interface_cond"
COMPLEXITY_COL  <- "task_complexity"
Y_RAW_COL       <- "tlx_mean"

# Opzionale: attiva only if you have a column task_id ben fatta
USE_TASK_RE <- FALSE
TASK_COL    <- "task_id"

REF_COND   <- "chatbot"
OTHER_COND <- "dashboard"

COMPLEXITY_LEVELS <- c("easy", "mid", "hard")

ROUND_TO   <- 5
LINK_FUN   <- "probit"
CONF_LEVEL <- 0.95
EMM_WEIGHTS <- "equal"

# with a only random effect you can keep >1; if you have problems of stability try 5 or 1
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

  est_col <- "Estimate"
  se_col  <- "Std. Error"

  if (!all(c(est_col, se_col) %in% names(tab))) {
 stop("Not trovo the columns Estimate / Std. Error in the summary of the model.")
  }

  tab <- tab %>%
    mutate(
      conf.low  = .data[[est_col]] - zcrit * .data[[se_col]],
      conf.high = .data[[est_col]] + zcrit * .data[[se_col]],
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

add_direction_flags <- function(df) {
  df %>%
    mutate(
      supports_chatbot_lower = estimate > 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
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
if (USE_TASK_RE) needed_cols <- c(needed_cols, TASK_COL)

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
  )

if (USE_TASK_RE) {
  df <- df %>%
    mutate(!!TASK_COL := as.character(.data[[TASK_COL]]))
}

df <- df %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[COMPLEXITY_COL]]),
    !is.na(.data[[Y_RAW_COL]])
  )

# ogni participant must stare in a only interface
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

present_complexities <- COMPLEXITY_LEVELS[COMPLEXITY_LEVELS %in% unique(df[[COMPLEXITY_COL]])]
if (length(present_complexities) < 2) {
 stop("Are used almeno 2 levels of task_complexity present.")
}
df[[COMPLEXITY_COL]] <- ordered(df[[COMPLEXITY_COL]], levels = present_complexities)

dup_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL)
if (USE_TASK_RE) dup_cols <- c(dup_cols, TASK_COL)

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
if (USE_TASK_RE) df[[TASK_COL]] <- factor(df[[TASK_COL]])

save_csv(df, file.path(OUT_DIR, "h1_cleaned_analysis_data.csv"))

# ============================================================
# DESCRIPTIVES
# ============================================================
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

desc_raw_by_interface <- df %>%
  group_by(.data[[COND_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
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

save_csv(desc_raw_by_cell,      file.path(OUT_DIR, "h1_descriptives_raw_by_cell.csv"))
save_csv(desc_raw_by_interface, file.path(OUT_DIR, "h1_descriptives_raw_by_interface.csv"))
save_csv(desc_ord_by_cell,      file.path(OUT_DIR, "h1_descriptives_ordinal_by_cell.csv"))

# ============================================================
# MODEL
# ============================================================
random_terms <- c(paste0("(1|", PARTICIPANT_COL, ")"))
if (USE_TASK_RE) random_terms <- c(random_terms, paste0("(1|", TASK_COL, ")"))

rhs_h1 <- paste(c(COND_COL, COMPLEXITY_COL, random_terms), collapse = " + ")
form_h1 <- as.formula(paste("y_ord ~", rhs_h1))

nAGQ_to_use <- if (length(random_terms) == 1) NAGQ_SINGLE_RE else 1

options(contrasts = c("contr.treatment", "contr.poly"))

m_h1 <- clmm(
  formula = form_h1,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = nAGQ_to_use
)

capture.output(summary(m_h1), file = file.path(OUT_DIR, "h1_model_summary.txt"))

# ============================================================
# FIXED EFFECTS + RANDOM EFFECTS
# ============================================================
coef_all <- make_coef_table(m_h1, conf.level = CONF_LEVEL)
coef_location  <- coef_all %>% filter(component == "location")
coef_threshold <- coef_all %>% filter(component == "threshold")

save_csv(coef_all,       file.path(OUT_DIR, "h1_coef_all.csv"))
save_csv(coef_location,  file.path(OUT_DIR, "h1_coef_location.csv"))
save_csv(coef_threshold, file.path(OUT_DIR, "h1_coef_thresholds.csv"))

# ---- random effects: estrazione robusta ----
vc <- tryCatch(VarCorr(m_h1), error = function(e) NULL)

if (!is.null(vc) && length(vc) > 0) {
  rand_tbl <- tibble(
    group = names(vc),
    variance_random_effect = vapply(vc, function(M) M[1, 1], numeric(1)),
    sd_random_effect = vapply(
      vc,
      function(M) {
        sd_attr <- attr(M, "stddev")
        if (!is.null(sd_attr) && length(sd_attr) >= 1) {
          as.numeric(sd_attr[1])
        } else {
          sqrt(as.numeric(M[1, 1]))
        }
      },
      numeric(1)
    )
  )
} else if (!is.null(m_h1$stDev) && length(m_h1$stDev) > 0) {
  grp <- names(m_h1$stDev)
  if (is.null(grp) || length(grp) != length(m_h1$stDev)) {
    grp <- paste0("RE_", seq_along(m_h1$stDev))
  }
  rand_tbl <- tibble(
    group = grp,
    sd_random_effect = as.numeric(unname(m_h1$stDev))
  )
} else if (!is.null(m_h1$tau) && length(m_h1$tau) > 0) {
  grp <- names(m_h1$tau)
  if (is.null(grp) || length(grp) != length(m_h1$tau)) {
    grp <- paste0("RE_", seq_along(m_h1$tau))
  }
  rand_tbl <- tibble(
    group = grp,
    sd_random_effect = as.numeric(exp(unname(m_h1$tau)))
  )
} else {
  rand_tbl <- tibble(
    note = "Random-effect SD non disponibile: controllare VarCorr(m_h1), m_h1$stDev, m_h1$tau e summary(m_h1)."
  )
}

save_csv(rand_tbl, file.path(OUT_DIR, "h1_random_effect_sd.csv"))

fit_tbl <- tibble(
  model  = "H1_adjusted",
  nobs   = nobs(m_h1),
  logLik = as.numeric(logLik(m_h1)),
  AIC    = AIC(m_h1),
  BIC    = BIC(m_h1)
)
save_csv(fit_tbl, file.path(OUT_DIR, "h1_model_fit.csv"))

# ============================================================
# EFFECT SIZES / CONTRASTS
# dashboard_minus_chatbot > 0 ==> dashboard higher MWL ==> chatbot lower MWL
# ============================================================
emm_interface_latent <- emmeans(
  m_h1,
  specs   = as.formula(paste("~", COND_COL)),
  mode    = "latent",
  weights = EMM_WEIGHTS
)

contrast_interface_latent <- summary(
  contrast(
    emm_interface_latent,
    method = dashboard_minus_chatbot_contrast()
  ),
  infer = c(TRUE, TRUE),
  level = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_direction_flags()

emm_interface_meanclass <- emmeans(
  m_h1,
  specs   = as.formula(paste("~", COND_COL)),
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

contrast_interface_meanclass <- summary(
  contrast(
    emm_interface_meanclass,
    method = dashboard_minus_chatbot_contrast()
  ),
  infer = c(TRUE, TRUE),
  level = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_direction_flags()

# predicted probabilities for category ordinal, aggregate on the complexity
emm_probs <- emmeans(
  m_h1,
  specs   = as.formula(paste("~ y_ord |", COND_COL)),
  mode    = "prob",
  weights = EMM_WEIGHTS
)

prob_tbl <- confint(emm_probs, level = CONF_LEVEL) %>%
  standardize_emm_ci_cols()

save_csv(
  summary(emm_interface_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h1_emm_interface_latent.csv")
)

save_csv(
  contrast_interface_latent,
  file.path(OUT_DIR, "h1_contrast_interface_latent.csv")
)

save_csv(
  summary(emm_interface_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h1_emm_interface_meanclass.csv")
)

save_csv(
  contrast_interface_meanclass,
  file.path(OUT_DIR, "h1_contrast_interface_meanclass.csv")
)

save_csv(
  prob_tbl,
  file.path(OUT_DIR, "h1_category_probabilities.csv")
)

cat("\n============================================================\n")
cat("H1 adjusted completata.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("============================================================\n\n")

cat("File principali:\n")
cat("- h1_coef_location.csv\n")
cat("- h1_contrast_interface_latent.csv\n")
cat("- h1_contrast_interface_meanclass.csv\n")
cat("- h1_category_probabilities.csv\n\n")

cat("Direction of the contrast: dashboard_minus_chatbot\n")
cat("estimate > 0 => dashboard higher => chatbot lower (supports H1)\n")