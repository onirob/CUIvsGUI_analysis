suppressPackageStartupMessages({
  library(ordinal)
  library(emmeans)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================================================
# H3a - Data Literacy x Interface -> Mental Workload (ordinal)
#
# PRIMARY HYPOTHESIS
#   Data literacy moderates the effect of interface type on MWL,
#   such that LOWER data literacy disadvantages the graphical
#   interface (dashboard) more than the conversational interface
#   (chatbot).
#
# MODEL
#   - cumulative link mixed model (CLMM)
#   - ordinal MWL outcome
#   - random intercept for participant
#   - task complexity included as nuisance fixed effect
#
# DL HANDLING (FIXED CATEGORIES FROM BDLI MEAN)
#   - use bdli_mean instead of bdli_total
#   - dl_ord categories:
#       low    = 1 <= mean <= 5
#       medium = 5 < mean <= 6
#       high   = 6 < mean <= 7
#
# NOTE ON BOUNDARIES
#   - exact 5.00 -> low
#   - exact 6.00 -> medium
#
# PRIMARY TEST
#   - additive model vs interaction model
#   - focus on interface_cond * dl_ord
#
# PARAMETER LABELS (updated)
#   beta1 = interface_conddashboard
#   beta2 = dl_ordmedium
#   beta3 = dl_ordhigh
#   beta4 = interface_conddashboard:dl_ordmedium
#   beta5 = interface_conddashboard:dl_ordhigh
#   beta6 = task_complexitymid
#   beta7 = task_complexityhard
#
# OUTPUT
#   - saved in a folder "results_H3" in the same directory as this script
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

# expects the csv to be in the same folder as this script
CSV_PATH <- file.path(SCRIPT_DIR, "h3a_mwl_dl_clmm.csv")
OUT_DIR  <- file.path(SCRIPT_DIR, "results_H3")

# ============================================================
# SETTINGS
# ============================================================
PARTICIPANT_COL <- "participant_id"
COND_COL        <- "interface_cond"
TASK_COMPLEXITY_COL <- "task_complexity"
TASK_ID_COL     <- "task_id"
Y_RAW_COL       <- "tlx_mean"

DL_MEAN_COL     <- "bdli_mean"
DL_N_COL        <- "bdli_n_answered"

REF_COND        <- "chatbot"
OTHER_COND      <- "dashboard"

TASK_COMPLEXITY_LEVELS <- c("easy", "mid", "hard")
REF_COMPLEXITY <- "easy"

DL_LEVELS <- c("low", "medium", "high")
REF_DL    <- "low"

# BDLI validity
REQUIRE_ALL_BDLI_ITEMS <- TRUE
BDLI_ITEMS_REQUIRED    <- 12

# MWL ordinalization
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

add_holm_correction <- function(df,
                                group_cols = NULL,
                                raw_p_col = "p.value.raw",
                                new_col = "p.value.holm",
                                alpha = 0.05) {
  out <- as_tibble(df)

  if (!(raw_p_col %in% names(out))) {
    out[[raw_p_col]] <- NA_real_
  }

  valid_groups <- character(0)
  if (!is.null(group_cols) && length(group_cols) > 0) {
    valid_groups <- group_cols[group_cols %in% names(out)]
  }

  if (length(valid_groups) == 0) {
    out[[new_col]] <- p.adjust(out[[raw_p_col]], method = "holm")
  } else {
    out <- out %>%
      group_by(across(all_of(valid_groups))) %>%
      mutate(!!new_col := p.adjust(.data[[raw_p_col]], method = "holm")) %>%
      ungroup()
  }

  out %>%
    mutate(
      holm_reject_0.05 = !is.na(.data[[new_col]]) & (.data[[new_col]] < alpha)
    )
}

make_parameter_label_lookup <- function() {
  tibble(
    term = c(
      paste0(COND_COL, OTHER_COND),
      "dl_ordmedium",
      "dl_ordhigh",
      paste0(COND_COL, OTHER_COND, ":dl_ordmedium"),
      paste0("dl_ordmedium:", COND_COL, OTHER_COND),
      paste0(COND_COL, OTHER_COND, ":dl_ordhigh"),
      paste0("dl_ordhigh:", COND_COL, OTHER_COND),
      paste0(TASK_COMPLEXITY_COL, "mid"),
      paste0(TASK_COMPLEXITY_COL, "hard")
    ),
    symbolic_label = c(
      "beta1",
      "beta2",
      "beta3",
      "beta4",
      "beta4",
      "beta5",
      "beta5",
      "beta6",
      "beta7"
    ),
    parameter_meaning = c(
      "dashboard vs chatbot at low DL",
      "medium vs low DL within chatbot",
      "high vs low DL within chatbot",
      "change in dashboard-chatbot gap at medium DL vs low DL",
      "change in dashboard-chatbot gap at medium DL vs low DL",
      "change in dashboard-chatbot gap at high DL vs low DL",
      "change in dashboard-chatbot gap at high DL vs low DL",
      "task complexity mid vs easy",
      "task complexity hard vs easy"
    )
  ) %>%
    distinct(term, .keep_all = TRUE)
}

make_parameter_label_export <- function() {
  tibble(
    symbolic_label = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7"),
    model_term = c(
      paste0(COND_COL, OTHER_COND),
      "dl_ordmedium",
      "dl_ordhigh",
      paste0(COND_COL, OTHER_COND, ":dl_ordmedium"),
      paste0(COND_COL, OTHER_COND, ":dl_ordhigh"),
      paste0(TASK_COMPLEXITY_COL, "mid"),
      paste0(TASK_COMPLEXITY_COL, "hard")
    ),
    meaning = c(
      "dashboard vs chatbot at low DL",
      "medium vs low DL within chatbot",
      "high vs low DL within chatbot",
      "change in dashboard-chatbot gap at medium DL vs low DL",
      "change in dashboard-chatbot gap at high DL vs low DL",
      "task complexity mid vs easy",
      "task complexity hard vs easy"
    )
  )
}

add_symbolic_labels <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !("term" %in% names(df))) {
    return(as_tibble(df))
  }

  lbl <- make_parameter_label_lookup()

  df %>%
    left_join(lbl, by = "term") %>%
    relocate(symbolic_label, .after = term) %>%
    relocate(parameter_meaning, .after = symbolic_label)
}

make_coef_table <- function(model, conf.level = 0.95) {
  sm <- summary(model)
  tab <- as.data.frame(coef(sm))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL

  zcrit <- qnorm(1 - (1 - conf.level) / 2)

  est_val <- if ("Estimate" %in% names(tab)) tab$Estimate else rep(NA_real_, nrow(tab))
  se_val  <- if ("Std. Error" %in% names(tab)) tab$`Std. Error` else rep(NA_real_, nrow(tab))
  z_val   <- if ("z value" %in% names(tab)) tab$`z value` else rep(NA_real_, nrow(tab))
  p_raw   <- if ("Pr(>|z|)" %in% names(tab)) tab$`Pr(>|z|)` else rep(NA_real_, nrow(tab))

  tab <- tab %>%
    mutate(
      estimate    = est_val,
      SE          = se_val,
      z.ratio     = z_val,
      p.value     = p_raw,
      p.value.raw = p_raw,
      conf.low    = estimate - zcrit * SE,
      conf.high   = estimate + zcrit * SE,
      lower.CL    = conf.low,
      upper.CL    = conf.high,
      component = case_when(
        str_detect(term, "\\|") ~ "threshold",
        TRUE ~ "location"
      )
    ) %>%
    relocate(term, component)

  add_symbolic_labels(as_tibble(tab))
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

  if (!("estimate" %in% names(df)) && "Estimate" %in% names(df)) {
    df$estimate <- df$Estimate
  }

  if (!("SE" %in% names(df)) && "Std. Error" %in% names(df)) {
    df$SE <- df[["Std. Error"]]
  }

  if (!("p.value.raw" %in% names(df))) {
    if ("p.value" %in% names(df)) {
      df$p.value.raw <- df$p.value
    } else if ("Pr(>|z|)" %in% names(df)) {
      df$p.value.raw <- df[["Pr(>|z|)"]]
    } else if ("Pr(>|t|)" %in% names(df)) {
      df$p.value.raw <- df[["Pr(>|t|)"]]
    } else {
      df$p.value.raw <- NA_real_
    }
  }

  if (!("p.value" %in% names(df))) df$p.value <- df$p.value.raw
  if (!("lower.CL" %in% names(df))) df$lower.CL <- NA_real_
  if (!("upper.CL" %in% names(df))) df$upper.CL <- NA_real_

  as_tibble(df)
}

summarize_contrast_with_holm <- function(x,
                                         conf.level = 0.95,
                                         holm_by = NULL) {
  summary(
    x,
    infer  = c(TRUE, TRUE),
    adjust = "none",
    level  = conf.level
  ) %>%
    standardize_emm_ci_cols() %>%
    add_holm_correction(group_cols = holm_by, raw_p_col = "p.value.raw")
}

add_interface_direction_flags <- function(df) {
  df %>%
    mutate(
      dashboard_higher_than_chatbot = estimate > 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

add_dl_direction_flags <- function(df) {
  df %>%
    mutate(
      higher_dl_lower_mwl = estimate < 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

add_did_direction_flags <- function(df) {
  df %>%
    mutate(
      supports_h3a_if_higher_dl_lowers_mwl = estimate < 0,
      ci_excludes_zero = (!is.na(lower.CL) & !is.na(upper.CL)) &
        ((lower.CL > 0) | (upper.CL < 0))
    )
}

dashboard_minus_chatbot_contrast <- function() {
  list(dashboard_minus_chatbot = c(-1, 1))
}

# for 3 DL levels: low, medium, high
# high_minus_low = [-1, 0, 1]
low_high_dl_contrast <- function() {
  list(high_minus_low = c(-1, 0, 1))
}

dashboard_minus_chatbot_two_estimates <- function() {
  list(dashboard_minus_chatbot = c(-1, 1))
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

make_dl_groups_from_mean <- function(x) {
  case_when(
    !is.na(x) & x >= 1 & x <= 5 ~ "low",
    !is.na(x) & x > 5 & x <= 6 ~ "medium",
    !is.na(x) & x > 6 & x <= 7 ~ "high",
    TRUE ~ NA_character_
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

needed_cols <- c(
  PARTICIPANT_COL,
  COND_COL,
  TASK_COMPLEXITY_COL,
  TASK_ID_COL,
  Y_RAW_COL,
  DL_MEAN_COL,
  DL_N_COL
)

missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0) {
 stop("Missing these columns in the CSV: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# CLEANING
# ============================================================
df <- df %>%
  mutate(
    !!PARTICIPANT_COL     := as.character(.data[[PARTICIPANT_COL]]),
    !!COND_COL            := str_to_lower(str_trim(as.character(.data[[COND_COL]]))),
    !!TASK_COMPLEXITY_COL := str_to_lower(str_trim(as.character(.data[[TASK_COMPLEXITY_COL]]))),
    !!TASK_ID_COL         := as.character(.data[[TASK_ID_COL]]),
    !!Y_RAW_COL           := as.numeric(.data[[Y_RAW_COL]]),
    !!DL_MEAN_COL         := as.numeric(.data[[DL_MEAN_COL]]),
    !!DL_N_COL            := as.numeric(.data[[DL_N_COL]])
  ) %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[TASK_COMPLEXITY_COL]]),
    !is.na(.data[[TASK_ID_COL]]),
    !is.na(.data[[Y_RAW_COL]])
  )

if (REQUIRE_ALL_BDLI_ITEMS) {
  df <- df %>% filter(.data[[DL_N_COL]] >= BDLI_ITEMS_REQUIRED)
}

df <- df %>% filter(!is.na(.data[[DL_MEAN_COL]]))

# remove participants assigned to more than one interface
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

present_complexities <- sort(unique(df[[TASK_COMPLEXITY_COL]]))
if (!all(TASK_COMPLEXITY_LEVELS %in% present_complexities)) {
 stop("Are used all the levels of task_complexity: ",
       paste(TASK_COMPLEXITY_LEVELS, collapse = ", "),
       ". Trovati: ", paste(present_complexities, collapse = ", "))
}
df[[TASK_COMPLEXITY_COL]] <- factor(df[[TASK_COMPLEXITY_COL]], levels = TASK_COMPLEXITY_LEVELS)

# keep only participants with all 3 task-complexity levels
counts_task <- df %>%
  group_by(.data[[PARTICIPANT_COL]]) %>%
  summarise(n_task = n_distinct(.data[[TASK_COMPLEXITY_COL]]), .groups = "drop")

keep_ids <- counts_task %>%
  filter(n_task == length(TASK_COMPLEXITY_LEVELS)) %>%
  pull(.data[[PARTICIPANT_COL]])

df <- df %>% filter(.data[[PARTICIPANT_COL]] %in% keep_ids)

dup_cols <- c(PARTICIPANT_COL, COND_COL, TASK_COMPLEXITY_COL, TASK_ID_COL)
dups <- duplicated(df[dup_cols]) | duplicated(df[dup_cols], fromLast = TRUE)
if (any(dups)) {
 stop("Ci are rows duplicate in the same cell logic. Controlla the CSV.")
}

# ============================================================
# BUILD FIXED ORDINAL DL MODERATOR FROM BDLI MEAN
# ============================================================
df <- df %>%
  mutate(
    dl_ord = make_dl_groups_from_mean(.data[[DL_MEAN_COL]])
  ) %>%
  filter(!is.na(dl_ord))

df$dl_ord <- factor(df$dl_ord, levels = DL_LEVELS, ordered = TRUE)

present_dl <- levels(droplevels(df$dl_ord))
if (!all(DL_LEVELS %in% present_dl)) {
 stop("After the categorizzazione of the BDLI mean are used all and 3 the levels of dl_ord: ",
       paste(DL_LEVELS, collapse = ", "),
       ". Trovati: ", paste(present_dl, collapse = ", "))
}

dl_map <- tibble(
  dl_ord = factor(DL_LEVELS, levels = DL_LEVELS, ordered = TRUE),
  lower_inclusive = c(1, NA_real_, NA_real_),
  lower_exclusive = c(NA_real_, 5, 6),
  upper_inclusive = c(5, 6, 7),
  upper_exclusive = c(NA_real_, NA_real_, NA_real_),
  rule = c(
    "1 <= bdli_mean <= 5",
    "5 < bdli_mean <= 6",
    "6 < bdli_mean <= 7"
  )
)

save_csv(dl_map, file.path(OUT_DIR, "h3a_dl_group_mapping.csv"))
save_csv(make_parameter_label_export(), file.path(OUT_DIR, "h3a_parameter_label_mapping.csv"))

# ============================================================
# ORDINALIZE MWL OUTCOME
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

save_csv(df, file.path(OUT_DIR, "h3a_cleaned_analysis_data.csv"))

# ============================================================
# DESCRIPTIVES
# ============================================================
desc_raw_by_cell <- df %>%
  group_by(.data[[COND_COL]], dl_ord, .data[[TASK_COMPLEXITY_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    sd     = sd(.data[[Y_RAW_COL]], na.rm = TRUE),
    median = median(.data[[Y_RAW_COL]], na.rm = TRUE),
    iqr    = iqr_num(.data[[Y_RAW_COL]]),
    .groups = "drop"
  )

desc_ord_by_cell <- df %>%
  group_by(.data[[COND_COL]], dl_ord, .data[[TASK_COMPLEXITY_COL]]) %>%
  summarise(
    n            = n(),
    mean_class   = mean(y_ord_num, na.rm = TRUE),
    sd_class     = sd(y_ord_num, na.rm = TRUE),
    median_class = median(y_ord_num, na.rm = TRUE),
    iqr_class    = iqr_num(y_ord_num),
    .groups = "drop"
  )

desc_n_by_interface_dl <- df %>%
  distinct(.data[[PARTICIPANT_COL]], .data[[COND_COL]], dl_ord) %>%
  group_by(.data[[COND_COL]], dl_ord) %>%
  summarise(n_participants = n(), .groups = "drop")

save_csv(desc_raw_by_cell, file.path(OUT_DIR, "h3a_descriptives_raw_by_cell.csv"))
save_csv(desc_ord_by_cell, file.path(OUT_DIR, "h3a_descriptives_ordinal_by_cell.csv"))
save_csv(desc_n_by_interface_dl, file.path(OUT_DIR, "h3a_n_participants_by_interface_dl.csv"))

# ============================================================
# MODELS
# ============================================================
options(contrasts = c("contr.treatment", "contr.treatment"))

form_add <- as.formula(
  paste("y_ord ~", COND_COL, "+ dl_ord +", TASK_COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")")
)

form_int <- as.formula(
  paste("y_ord ~", COND_COL, "* dl_ord +", TASK_COMPLEXITY_COL, "+ (1|", PARTICIPANT_COL, ")")
)

nAGQ_to_use <- NAGQ_SINGLE_RE

m_add <- clmm(
  formula = form_add,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = nAGQ_to_use
)

m_int <- clmm(
  formula = form_int,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = nAGQ_to_use
)

capture.output(summary(m_add), file = file.path(OUT_DIR, "h3a_additive_model_summary.txt"))
capture.output(summary(m_int), file = file.path(OUT_DIR, "h3a_interaction_model_summary.txt"))

# ============================================================
# MODEL FIT / LRT
# ============================================================
fit_tbl <- tibble(
  model  = c("additive", "interaction"),
  nobs   = c(nobs(m_add), nobs(m_int)),
  logLik = c(as.numeric(logLik(m_add)), as.numeric(logLik(m_int))),
  AIC    = c(AIC(m_add), AIC(m_int)),
  BIC    = c(BIC(m_add), BIC(m_int))
)
save_csv(fit_tbl, file.path(OUT_DIR, "h3a_model_fit.csv"))

lrt_tbl <- as.data.frame(anova(m_add, m_int))
save_csv(lrt_tbl, file.path(OUT_DIR, "h3a_interaction_lrt.csv"))
capture.output(anova(m_add, m_int), file = file.path(OUT_DIR, "h3a_interaction_lrt.txt"))

# ============================================================
# COEFFICIENTS
# ============================================================
coef_add_all <- make_coef_table(m_add, conf.level = CONF_LEVEL)
coef_int_all <- make_coef_table(m_int, conf.level = CONF_LEVEL)

coef_add_location <- coef_add_all %>% filter(component == "location")
coef_int_location <- coef_int_all %>% filter(component == "location")
coef_int_main_terms_only <- coef_int_location %>% filter(!str_detect(term, ":"))

# internal DiDs = beta4 e beta5
coef_int_interaction_only <- coef_int_location %>%
  filter(str_detect(term, ":")) %>%
  add_holm_correction(group_cols = NULL, raw_p_col = "p.value.raw") %>%
  add_did_direction_flags()

save_csv(coef_add_all,  file.path(OUT_DIR, "h3a_additive_coef_all.csv"))
save_csv(coef_int_all,  file.path(OUT_DIR, "h3a_interaction_coef_all.csv"))
save_csv(coef_add_location, file.path(OUT_DIR, "h3a_additive_coef_location.csv"))
save_csv(coef_int_location, file.path(OUT_DIR, "h3a_interaction_coef_location.csv"))
save_csv(coef_int_main_terms_only, file.path(OUT_DIR, "h3a_interaction_main_terms_only.csv"))
save_csv(coef_int_interaction_only, file.path(OUT_DIR, "h3a_interaction_terms_only.csv"))
save_csv(coef_int_interaction_only, file.path(OUT_DIR, "h3a_internal_did_terms_only.csv"))

# ============================================================
# RANDOM EFFECT SD
# ============================================================
rand_add <- extract_random_effect_sd(m_add, "additive")
rand_int <- extract_random_effect_sd(m_int, "interaction")
rand_tbl <- bind_rows(rand_add, rand_int)

if (nrow(rand_tbl) == 0) {
 warning("These are not riuscito to estrarre the random effects; controlla the summary testuali.")
}

save_csv(rand_tbl, file.path(OUT_DIR, "h3a_random_effect_sd.csv"))

writeLines(
  c(
    "Random-effect extraction diagnostics",
    paste0("additive rows: ", nrow(rand_add)),
    paste0("interaction rows: ", nrow(rand_int)),
    if (nrow(rand_add) > 0) paste0("additive source(s): ", paste(unique(rand_add$source), collapse = ", ")) else "additive source(s): none",
    if (nrow(rand_int) > 0) paste0("interaction source(s): ", paste(unique(rand_int$source), collapse = ", ")) else "interaction source(s): none"
  ),
  con = file.path(OUT_DIR, "h3a_random_effect_sd_diagnostics.txt")
)

# ============================================================
# SIMPLE EFFECTS OF INTERFACE WITHIN EACH DL LEVEL
# dashboard_minus_chatbot > 0 ==> dashboard higher MWL
# H3a pattern: this gap should be largest at LOW DL
# Holm on the family of the 3 test
# ============================================================
emm_interface_by_dl_latent <- emmeans(
  m_int,
  specs   = as.formula(paste("~", COND_COL, "| dl_ord")),
  mode    = "latent",
  weights = EMM_WEIGHTS
)

gap_by_dl_latent <- contrast(
  emm_interface_by_dl_latent,
  method = dashboard_minus_chatbot_contrast()
)

simple_effects_interface_latent <- summarize_contrast_with_holm(
  gap_by_dl_latent,
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_interface_direction_flags()

emm_interface_by_dl_meanclass <- emmeans(
  m_int,
  specs   = as.formula(paste("~", COND_COL, "| dl_ord")),
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

gap_by_dl_meanclass <- contrast(
  emm_interface_by_dl_meanclass,
  method = dashboard_minus_chatbot_contrast()
)

simple_effects_interface_meanclass <- summarize_contrast_with_holm(
  gap_by_dl_meanclass,
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_interface_direction_flags()

save_csv(
  summary(emm_interface_by_dl_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_emm_interface_by_dl_latent.csv")
)

save_csv(
  simple_effects_interface_latent,
  file.path(OUT_DIR, "h3a_simple_effects_interface_by_dl_latent.csv")
)

save_csv(
  summary(emm_interface_by_dl_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_emm_interface_by_dl_meanclass.csv")
)

save_csv(
  simple_effects_interface_meanclass,
  file.path(OUT_DIR, "h3a_simple_effects_interface_by_dl_meanclass.csv")
)

# ============================================================
# SIMPLE EFFECTS OF DL (HIGH - LOW) WITHIN EACH INTERFACE
# high_minus_low < 0 ==> higher DL associated with lower MWL
# equivalently: LOW DL is associated with higher MWL
# H3a pattern: this contrast should be more negative in dashboard
# Holm on the family of the 2 test
# ============================================================
emm_dl_by_interface_latent <- emmeans(
  m_int,
  specs   = as.formula(paste("~ dl_ord |", COND_COL)),
  mode    = "latent",
  weights = EMM_WEIGHTS
)

dl_lowhigh_by_interface_latent <- contrast(
  emm_dl_by_interface_latent,
  method = low_high_dl_contrast()
)

simple_effects_dl_lowhigh_latent <- summarize_contrast_with_holm(
  dl_lowhigh_by_interface_latent,
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_dl_direction_flags()

emm_dl_by_interface_meanclass <- emmeans(
  m_int,
  specs   = as.formula(paste("~ dl_ord |", COND_COL)),
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

dl_lowhigh_by_interface_meanclass <- contrast(
  emm_dl_by_interface_meanclass,
  method = low_high_dl_contrast()
)

simple_effects_dl_lowhigh_meanclass <- summarize_contrast_with_holm(
  dl_lowhigh_by_interface_meanclass,
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_dl_direction_flags()

save_csv(
  summary(emm_dl_by_interface_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_emm_dl_by_interface_latent.csv")
)

save_csv(
  simple_effects_dl_lowhigh_latent,
  file.path(OUT_DIR, "h3a_simple_effects_dl_lowhigh_by_interface_latent.csv")
)

save_csv(
  summary(emm_dl_by_interface_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_emm_dl_by_interface_meanclass.csv")
)

save_csv(
  simple_effects_dl_lowhigh_meanclass,
  file.path(OUT_DIR, "h3a_simple_effects_dl_lowhigh_by_interface_meanclass.csv")
)

# ============================================================
# DID-LIKE CONTRAST ON THE LOW-HIGH DL EFFECT
# dashboard_minus_chatbot < 0 ==> the high-low improvement is more negative
# in dashboard than in chatbot
# equivalently: LOW DL penalizes dashboard more than chatbot
# a only test -> Holm coincide with raw
# ============================================================
did_like_latent <- summarize_contrast_with_holm(
  contrast(
    dl_lowhigh_by_interface_latent,
    method = dashboard_minus_chatbot_two_estimates(),
    by = NULL
  ),
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_did_direction_flags()

did_like_meanclass <- summarize_contrast_with_holm(
  contrast(
    dl_lowhigh_by_interface_meanclass,
    method = dashboard_minus_chatbot_two_estimates(),
    by = NULL
  ),
  conf.level = CONF_LEVEL,
  holm_by    = NULL
) %>%
  add_did_direction_flags()

save_csv(
  did_like_latent,
  file.path(OUT_DIR, "h3a_did_like_lowhigh_latent.csv")
)

save_csv(
  did_like_meanclass,
  file.path(OUT_DIR, "h3a_did_like_lowhigh_meanclass.csv")
)

# ============================================================
# ALL CELL MEANS (AVERAGED OVER TASK)
# ============================================================
emm_all_cells_latent <- emmeans(
  m_int,
  specs   = as.formula(paste("~", COND_COL, "* dl_ord")),
  mode    = "latent",
  weights = EMM_WEIGHTS
)

emm_all_cells_meanclass <- emmeans(
  m_int,
  specs   = as.formula(paste("~", COND_COL, "* dl_ord")),
  mode    = "mean.class",
  weights = EMM_WEIGHTS
)

save_csv(
  summary(emm_all_cells_latent, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_all_cells_latent.csv")
)

save_csv(
  summary(emm_all_cells_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  file.path(OUT_DIR, "h3a_all_cells_meanclass.csv")
)

# ============================================================
# CATEGORY PROBABILITIES BY CELL (AVERAGED OVER TASK)
# ============================================================
emm_prob_cells <- emmeans(
  m_int,
  specs   = as.formula(paste("~ y_ord |", COND_COL, "* dl_ord")),
  mode    = "prob",
  weights = EMM_WEIGHTS
)

prob_cells_tbl <- confint(emm_prob_cells, level = CONF_LEVEL) %>%
  standardize_emm_ci_cols()

save_csv(prob_cells_tbl, file.path(OUT_DIR, "h3a_category_probabilities_by_cell.csv"))

# ============================================================
# OPTIONAL: OVERALL INTERFACE EFFECT FROM FULL MODEL
# averaged over DL and task
# ============================================================
emm_overall_interface_latent <- emmeans(
  m_int,
  specs   = as.formula(paste("~", COND_COL)),
  mode    = "latent",
  weights = EMM_WEIGHTS
)

overall_interface_contrast_latent <- summary(
  contrast(
    emm_overall_interface_latent,
    method = dashboard_minus_chatbot_contrast()
  ),
  infer  = c(TRUE, TRUE),
  adjust = "none",
  level  = CONF_LEVEL
) %>%
  standardize_emm_ci_cols() %>%
  add_interface_direction_flags()

save_csv(
  overall_interface_contrast_latent,
  file.path(OUT_DIR, "h3a_overall_interface_contrast_latent_from_full_model.csv")
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("H3a interaction (Data Literacy mean categories x Interface) completata.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("============================================================\n\n")

cat("Riferimenti dei fattori:\n")
cat("- interface_cond reference =", REF_COND, "\n")
cat("- task_complexity reference =", REF_COMPLEXITY, "\n")
cat("- dl_ord reference =", REF_DL, "\n\n")

cat("Regole of categorizzazione of the Data Literacy (bdli_mean):\n")
cat("- low    : 1 <= bdli_mean <= 5\n")
cat("- medium : 5 < bdli_mean <= 6\n")
cat("- high   : 6 < bdli_mean <= 7\n\n")

cat("Labels parameter of the fixed effects:\n")
cat("- beta1 = interface_conddashboard = dashboard vs chatbot at the level 'low' of DL\n")
cat("- beta2 = dl_ordmedium = difference medium vs low in the chatbot\n")
cat("- beta3 = dl_ordhigh = difference high vs low in the chatbot\n")
cat("- beta4 = interface_conddashboard:dl_ordmedium = change of the dashboard–chatbot gap to DL medium vs low\n")
cat("- beta5 = interface_conddashboard:dl_ordhigh = change of the dashboard–chatbot gap to DL high vs low\n")
cat("- beta6 = task_complexitymid = effect of control mid vs easy\n")
cat("- beta7 = task_complexityhard = effect of control hard vs easy\n\n")

cat("Correzioni Holm implementate:\n")
cat("- internal DiDs of the model (beta4, beta5)\n")
cat("- simple effects of the interface within DL -> family of 3 test\n")
cat("- simple effects high-minus-low DL within interface -> family of 2 test\n")
cat("- DID-like exact contrast -> a only test (Holm coincide col raw)\n\n")

cat("Columns aggiunte in the file corretti:\n")
cat("- p.value.raw\n")
cat("- p.value.holm\n")
cat("- holm_reject_0.05\n\n")

cat("Pattern expected for H3a:\n")
cat("- the dashboard–chatbot gap should be largest to DL low\n")
cat("- the contrast high_minus_low should be more negative in the dashboard\n")
cat("- the DID-like on the contrast high_minus_low should be < 0\n")
cat(" => this equivale to indicate that low DL penalizza of more the dashboard that the chatbot\n\n")

cat("File principali:\n")
cat("- h3a_interaction_lrt.csv\n")
cat("- h3a_interaction_terms_only.csv\n")
cat("- h3a_internal_did_terms_only.csv\n")
cat("- h3a_simple_effects_interface_by_dl_latent.csv\n")
cat("- h3a_simple_effects_dl_lowhigh_by_interface_latent.csv\n")
cat("- h3a_did_like_lowhigh_latent.csv\n")
cat("- h3a_category_probabilities_by_cell.csv\n")
cat("- h3a_random_effect_sd.csv\n")
cat("- h3a_dl_group_mapping.csv\n")
cat("- h3a_parameter_label_mapping.csv\n\n")

cat("Direction of the contrasts of interface: dashboard_minus_chatbot\n")
cat("estimate > 0 => dashboard higher in MWL => chatbot lower in MWL\n")
cat("For H3a, ci si aspetta that this estimate sia larger to DL low.\n\n")

cat("Direction of the contrasts DL low-high: high_minus_low\n")
cat("estimate < 0 => higher DL associated to MWL lower\n")
cat(" => equivalentemente, low DL associated to MWL higher\n\n")

cat("Direction of the DID-like on the contrast high-low:\n")
cat("estimate < 0 => the effect high-low DL is more negative in the dashboard that in the chatbot\n")
cat(" => equivalentemente, low DL penalizza more the dashboard of the chatbot\n")