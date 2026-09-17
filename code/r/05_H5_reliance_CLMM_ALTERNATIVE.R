suppressPackageStartupMessages({
  library(ordinal)
  library(emmeans)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(rlang)
})

# ============================================================
# H5 ALTERNATIVE - Reliance CLMM
#
# Alternative/confounder model:
#   industry_role is added as a covariate/confounder.
#
# H5A adjusted:
#   Interface type is associated with differences in reliance,
#   after controlling for task_complexity and industry_role.
#
# H5B adjusted:
#   The interface effect on reliance depends on task_complexity,
#   after controlling for industry_role.
#
# industry_role is NOT a moderator here.
# It is included only as a main-effect control variable.
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

  normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()

CSV_PATH <- file.path(SCRIPT_DIR, "reliance_clmm.csv")
OUT_DIR  <- file.path(SCRIPT_DIR, "results_H5_reliance_clmm_ALTERNATIVE")

QUESTIONNAIRE_CSV_CANDIDATES <- c(
  file.path(SCRIPT_DIR, "questionnaire_responses.csv"),
  file.path(SCRIPT_DIR, "..", "questionnaire_responses_20260203094443.csv"),
  file.path(SCRIPT_DIR, "..", "..", "analysis_data", "questionnaire_responses_20260203094443.csv")
)

resolve_existing_path <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) return(normalizePath(p))
  }
  NA_character_
}

QUESTIONNAIRE_CSV_PATH <- resolve_existing_path(QUESTIONNAIRE_CSV_CANDIDATES)

# ============================================================
# SETTINGS
# ============================================================
PARTICIPANT_COL <- "participant_id"
COND_COL        <- "interface_cond"
COMPLEXITY_COL  <- "task_complexity"
ITEM_COL        <- "item_key"
Y_RAW_COL       <- "reliance_value"

INDUSTRY_ROLE_COL      <- "industry_role"
INDUSTRY_ROLE_CODE_COL <- "industry_role_code"

REF_COND          <- "chatbot"
OTHER_COND        <- "dashboard"
ITEM_LEVELS       <- c("r1", "r2", "r3")
COMPLEXITY_LEVELS <- c("easy", "mid", "hard")
REF_COMPLEXITY    <- "easy"

ROLE_ORDER <- c(
  "Junior Management",
  "Middle Management",
  "Upper Management"
)
REF_ROLE <- "Junior Management"

LINK_FUN       <- "probit"
CONF_LEVEL     <- 0.95
EMM_WEIGHTS    <- "equal"
NAGQ_SINGLE_RE <- 10

OUT_DIR_H5A <- file.path(OUT_DIR, "H5A_interface_ALTERNATIVE")
OUT_DIR_H5B <- file.path(OUT_DIR, "H5B_interaction_ALTERNATIVE")

# ============================================================
# HELPERS
# ============================================================
alt_file <- function(dir, filename) {
  file.path(
    dir,
    sub("(\\.[^.]+)$", "_ALTERNATIVE\\1", filename)
  )
}

iqr_num <- function(x) {
  stats::IQR(x, na.rm = TRUE, type = 7)
}

save_csv <- function(x, file) {
  readr::write_csv(as_tibble(x), file)
}

epsilon_sq_kw <- function(H, n, k) {
  denom <- n - k
  if (denom <= 0) return(NA_real_)
  max(0, min(1, as.numeric((H - k + 1) / denom)))
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

add_industry_role_if_needed <- function(df) {
  if (INDUSTRY_ROLE_COL %in% names(df)) {
    return(df)
  }

  if (is.na(QUESTIONNAIRE_CSV_PATH) || !file.exists(QUESTIONNAIRE_CSV_PATH)) {
    stop(
 "The industry_role column is not in the primary file and no questionnaire CSV was found. ",
      "Candidati controllati: ",
      paste(QUESTIONNAIRE_CSV_CANDIDATES, collapse = " | ")
    )
  }

  message("Carico industry_role da: ", QUESTIONNAIRE_CSV_PATH)

  qr <- read_csv(QUESTIONNAIRE_CSV_PATH, show_col_types = FALSE)

  needed_qr_cols <- c(PARTICIPANT_COL, "questionnaire_name", "item_key", "value_text")
  missing_qr_cols <- setdiff(needed_qr_cols, names(qr))

  if (length(missing_qr_cols) > 0) {
 stop("Missing columns in the questionnaire CSV: ", paste(missing_qr_cols, collapse = ", "))
  }

  role_df <- qr %>%
    filter(
      .data[["questionnaire_name"]] == "PrescreenStart",
      .data[["item_key"]] == "industry_role"
    ) %>%
    transmute(
      !!PARTICIPANT_COL := as.character(.data[[PARTICIPANT_COL]]),
      !!INDUSTRY_ROLE_COL := as.character(.data[["value_text"]])
    ) %>%
    distinct(.data[[PARTICIPANT_COL]], .keep_all = TRUE)

  df %>%
    mutate(!!PARTICIPANT_COL := as.character(.data[[PARTICIPANT_COL]])) %>%
    left_join(role_df, by = PARTICIPANT_COL)
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
df <- add_industry_role_if_needed(df)

needed_cols <- c(
  PARTICIPANT_COL,
  COND_COL,
  COMPLEXITY_COL,
  ITEM_COL,
  Y_RAW_COL,
  INDUSTRY_ROLE_COL
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
    !!PARTICIPANT_COL   := as.character(.data[[PARTICIPANT_COL]]),
    !!COND_COL          := str_to_lower(str_trim(as.character(.data[[COND_COL]]))),
    !!COMPLEXITY_COL    := str_to_lower(str_trim(as.character(.data[[COMPLEXITY_COL]]))),
    !!ITEM_COL          := str_to_lower(str_trim(as.character(.data[[ITEM_COL]]))),
    !!Y_RAW_COL         := as.numeric(.data[[Y_RAW_COL]]),
    !!INDUSTRY_ROLE_COL := str_squish(as.character(.data[[INDUSTRY_ROLE_COL]]))
  ) %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[COMPLEXITY_COL]]),
    !is.na(.data[[ITEM_COL]]),
    !is.na(.data[[Y_RAW_COL]]),
    !is.na(.data[[INDUSTRY_ROLE_COL]])
  )

unexpected_roles <- sort(setdiff(
  unique(df[[INDUSTRY_ROLE_COL]][!is.na(df[[INDUSTRY_ROLE_COL]])]),
  ROLE_ORDER
))

if (length(unexpected_roles) > 0) {
  message("ATTENZIONE: industry_role inattesi esclusi: ",
          paste(unexpected_roles, collapse = ", "))
}

df <- df %>%
  filter(.data[[INDUSTRY_ROLE_COL]] %in% ROLE_ORDER)

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

if (any(abs(df[[Y_RAW_COL]] - round(df[[Y_RAW_COL]])) > 1e-8, na.rm = TRUE)) {
 stop("reliance_value contiene values not interi. The CLMM here richiede ordinal categories discrete.")
}

df[[Y_RAW_COL]] <- as.integer(round(df[[Y_RAW_COL]]))

dup_cols <- c(PARTICIPANT_COL, COND_COL, COMPLEXITY_COL, ITEM_COL)
dups <- duplicated(df[dup_cols]) | duplicated(df[dup_cols], fromLast = TRUE)

if (any(dups)) {
 stop("Ci are rows duplicate in the same cell logic: participant x interface x complexity x item.")
}

df[[COND_COL]]        <- factor(df[[COND_COL]], levels = c(REF_COND, OTHER_COND))
df[[COMPLEXITY_COL]]  <- factor(df[[COMPLEXITY_COL]], levels = COMPLEXITY_LEVELS)
df[[ITEM_COL]]        <- factor(df[[ITEM_COL]], levels = ITEM_LEVELS)
df[[INDUSTRY_ROLE_COL]] <- factor(df[[INDUSTRY_ROLE_COL]], levels = ROLE_ORDER)
df[[INDUSTRY_ROLE_CODE_COL]] <- as.integer(df[[INDUSTRY_ROLE_COL]]) - 1
df[[PARTICIPANT_COL]] <- factor(df[[PARTICIPANT_COL]])

global_y_levels <- sort(unique(df[[Y_RAW_COL]]))
df$y_ord_global <- ordered(df[[Y_RAW_COL]], levels = global_y_levels)
df$y_ord_num_global <- as.integer(df$y_ord_global)

save_csv(df, alt_file(OUT_DIR, "h5_reliance_cleaned_analysis_data.csv"))

# ============================================================
# INDUSTRY_ROLE DIAGNOSTICS
# ============================================================
role_participant <- df %>%
  group_by(
    .data[[PARTICIPANT_COL]],
    .data[[INDUSTRY_ROLE_COL]],
    .data[[INDUSTRY_ROLE_CODE_COL]],
    .data[[COND_COL]]
  ) %>%
  summarise(
    reliance_mean_participant = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    .groups = "drop"
  )

role_counts <- role_participant %>%
  count(.data[[INDUSTRY_ROLE_COL]], name = "n") %>%
  mutate(prop = n / sum(n))

role_by_interface <- role_participant %>%
  count(.data[[COND_COL]], .data[[INDUSTRY_ROLE_COL]], name = "n") %>%
  group_by(.data[[COND_COL]]) %>%
  mutate(prop_within_interface = n / sum(n)) %>%
  ungroup()

desc_raw_by_role <- role_participant %>%
  group_by(.data[[INDUSTRY_ROLE_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(reliance_mean_participant, na.rm = TRUE),
    sd     = sd(reliance_mean_participant, na.rm = TRUE),
    median = median(reliance_mean_participant, na.rm = TRUE),
    iqr    = iqr_num(reliance_mean_participant),
    .groups = "drop"
  )

role_kw <- kruskal.test(
  reliance_mean_participant ~ industry_role,
  data = role_participant
)

role_spearman <- suppressWarnings(
  cor.test(
    role_participant$industry_role_code,
    role_participant$reliance_mean_participant,
    method = "spearman",
    exact = FALSE
  )
)

role_interface_chisq <- suppressWarnings(
  chisq.test(table(role_participant[[INDUSTRY_ROLE_COL]], role_participant[[COND_COL]]))
)

industry_role_omnibus <- tibble(
  analysis = c(
    "Industry role vs participant-level reliance",
    "Ordinal trend: industry role code vs participant-level reliance",
    "Industry role balance across interface condition"
  ),
  test = c(
    "Kruskal-Wallis",
    "Spearman rank correlation",
    "Chi-square"
  ),
  statistic_name = c("H", "rho", "X-squared"),
  statistic = c(
    as.numeric(role_kw$statistic),
    as.numeric(role_spearman$estimate),
    as.numeric(role_interface_chisq$statistic)
  ),
  df = c(
    as.numeric(role_kw$parameter),
    NA_real_,
    as.numeric(role_interface_chisq$parameter)
  ),
  p.value = c(
    as.numeric(role_kw$p.value),
    as.numeric(role_spearman$p.value),
    as.numeric(role_interface_chisq$p.value)
  ),
  effect_name = c(
    "epsilon_squared",
    "rho",
    NA_character_
  ),
  effect_size = c(
    epsilon_sq_kw(
      role_kw$statistic,
      n = nrow(role_participant),
      k = n_distinct(role_participant[[INDUSTRY_ROLE_COL]])
    ),
    as.numeric(role_spearman$estimate),
    NA_real_
  ),
  n = c(
    nrow(role_participant),
    nrow(role_participant),
    nrow(role_participant)
  )
)

save_csv(role_counts, alt_file(OUT_DIR, "h5_industry_role_counts.csv"))
save_csv(role_by_interface, alt_file(OUT_DIR, "h5_industry_role_by_interface.csv"))
save_csv(desc_raw_by_role, alt_file(OUT_DIR, "h5_descriptives_raw_by_industry_role.csv"))
save_csv(industry_role_omnibus, alt_file(OUT_DIR, "h5_industry_role_omnibus.csv"))

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

save_csv(h5a_desc_raw, alt_file(OUT_DIR_H5A, "h5a_descriptives_raw_by_item_interface.csv"))
save_csv(h5a_desc_ord, alt_file(OUT_DIR_H5A, "h5a_descriptives_ordinal_by_item_interface.csv"))
save_csv(h5a_n_pp, alt_file(OUT_DIR_H5A, "h5a_n_participants_by_item_interface.csv"))

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

save_csv(h5b_desc_raw, alt_file(OUT_DIR_H5B, "h5b_descriptives_raw_by_cell.csv"))
save_csv(h5b_desc_ord, alt_file(OUT_DIR_H5B, "h5b_descriptives_ordinal_by_cell.csv"))
save_csv(h5b_n_pp, alt_file(OUT_DIR_H5B, "h5b_n_participants_by_cell.csv"))

# ============================================================
# CONTRAST CODING
# ============================================================
options(contrasts = c("contr.treatment", "contr.treatment"))

# ============================================================
# H5A
# Model: y ~ industry_role + interface + task_complexity + (1|participant)
# ============================================================
h5a_fit_tbl <- list()
h5a_lrt_tbl <- list()
h5a_coef_all_tbl <- list()
h5a_coef_location_tbl <- list()
h5a_coef_industry_role_tbl <- list()
h5a_random_tbl <- list()
h5a_emm_interface_latent_tbl <- list()
h5a_emm_interface_meanclass_tbl <- list()
h5a_overall_contrast_latent_tbl <- list()
h5a_overall_contrast_meanclass_tbl <- list()
h5a_prob_interface_tbl <- list()

for (item in ITEM_LEVELS) {
  dat <- build_item_dataset(df, item)

  form_add <- as.formula(
    paste(
      "y_ord ~",
      INDUSTRY_ROLE_COL,
      "+",
      COND_COL,
      "+",
      COMPLEXITY_COL,
      "+ (1|",
      PARTICIPANT_COL,
      ")"
    )
  )

  form_no_interface <- as.formula(
    paste(
      "y_ord ~",
      INDUSTRY_ROLE_COL,
      "+",
      COMPLEXITY_COL,
      "+ (1|",
      PARTICIPANT_COL,
      ")"
    )
  )

  form_no_complexity <- as.formula(
    paste(
      "y_ord ~",
      INDUSTRY_ROLE_COL,
      "+",
      COND_COL,
      "+ (1|",
      PARTICIPANT_COL,
      ")"
    )
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
    file = alt_file(OUT_DIR_H5A, paste0("h5a_", item, "_additive_model_summary.txt"))
  )

  fit_tbl_item <- tibble(
    item_key = item,
    model  = c("no_interface_adjusted_for_role", "no_complexity_adjusted_for_role", "additive_adjusted_for_role"),
    nobs   = c(nobs(m_no_interface), nobs(m_no_complexity), nobs(m_add)),
    logLik = c(as.numeric(logLik(m_no_interface)), as.numeric(logLik(m_no_complexity)), as.numeric(logLik(m_add))),
    AIC    = c(AIC(m_no_interface), AIC(m_no_complexity), AIC(m_add)),
    BIC    = c(BIC(m_no_interface), BIC(m_no_complexity), BIC(m_add))
  )

  h5a_fit_tbl[[item]] <- fit_tbl_item

  lrt_interface <- as.data.frame(anova(m_no_interface, m_add))
  lrt_interface$item_key <- item
  lrt_interface$effect <- "interface_main_effect_adjusted_for_complexity_and_industry_role"

  lrt_complexity <- as.data.frame(anova(m_no_complexity, m_add))
  lrt_complexity$item_key <- item
  lrt_complexity$effect <- "complexity_main_effect_adjusted_for_interface_and_industry_role"

  h5a_lrt_tbl[[item]] <- bind_rows(
    as_tibble(lrt_interface),
    as_tibble(lrt_complexity)
  )

  coef_add_all <- make_coef_table(m_add, conf.level = CONF_LEVEL) %>%
    mutate(item_key = item, .before = 1)

  coef_add_location <- coef_add_all %>%
    filter(component == "location")

  coef_add_industry_role <- coef_add_location %>%
    filter(str_detect(term, INDUSTRY_ROLE_COL)) %>%
    add_holm_correction(raw_p_col = "p.value.raw")

  h5a_coef_all_tbl[[item]] <- coef_add_all
  h5a_coef_location_tbl[[item]] <- coef_add_location
  h5a_coef_industry_role_tbl[[item]] <- coef_add_industry_role

  rand_add <- extract_random_effect_sd(m_add, item, "additive_adjusted_for_role")
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

  h5a_emm_interface_latent_tbl[[item]] <- summary(
    emm_overall_interface_latent,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5a_emm_interface_meanclass_tbl[[item]] <- summary(
    emm_overall_interface_meanclass,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)

  h5a_overall_contrast_latent_tbl[[item]] <- summary(
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

  h5a_overall_contrast_meanclass_tbl[[item]] <- summary(
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

  emm_prob_interface <- emmeans(
    m_add,
    specs   = as.formula(paste("~ y_ord |", COND_COL)),
    mode    = "prob",
    weights = EMM_WEIGHTS
  )

  h5a_prob_interface_tbl[[item]] <- confint(
    emm_prob_interface,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)
}

h5a_model_fit <- bind_rows(h5a_fit_tbl)
h5a_lrt       <- bind_rows(h5a_lrt_tbl)
h5a_coef_all  <- bind_rows(h5a_coef_all_tbl)
h5a_coef_loc  <- bind_rows(h5a_coef_location_tbl)
h5a_coef_industry_role <- bind_rows(h5a_coef_industry_role_tbl)
h5a_random    <- bind_rows(h5a_random_tbl)
h5a_emm_latent <- bind_rows(h5a_emm_interface_latent_tbl)
h5a_emm_meanclass <- bind_rows(h5a_emm_interface_meanclass_tbl)
h5a_contrast_latent <- bind_rows(h5a_overall_contrast_latent_tbl)
h5a_contrast_meanclass <- bind_rows(h5a_overall_contrast_meanclass_tbl)
h5a_prob_interface <- bind_rows(h5a_prob_interface_tbl)

save_csv(h5a_model_fit, alt_file(OUT_DIR_H5A, "h5a_model_fit.csv"))
save_csv(h5a_lrt, alt_file(OUT_DIR_H5A, "h5a_omnibus_lrt.csv"))

save_csv(
  h5a_lrt %>% filter(effect == "interface_main_effect_adjusted_for_complexity_and_industry_role"),
  alt_file(OUT_DIR_H5A, "h5a_interface_lrt.csv")
)

save_csv(
  h5a_lrt %>% filter(effect == "complexity_main_effect_adjusted_for_interface_and_industry_role"),
  alt_file(OUT_DIR_H5A, "h5a_complexity_lrt.csv")
)

capture.output(h5a_lrt, file = alt_file(OUT_DIR_H5A, "h5a_omnibus_lrt.txt"))

save_csv(h5a_coef_all, alt_file(OUT_DIR_H5A, "h5a_additive_coef_all.csv"))
save_csv(h5a_coef_loc, alt_file(OUT_DIR_H5A, "h5a_additive_coef_location.csv"))
save_csv(h5a_coef_industry_role, alt_file(OUT_DIR_H5A, "h5a_additive_industry_role_terms_only.csv"))
save_csv(h5a_random, alt_file(OUT_DIR_H5A, "h5a_random_effect_sd.csv"))
save_csv(h5a_emm_latent, alt_file(OUT_DIR_H5A, "h5a_emm_interface_latent.csv"))
save_csv(h5a_contrast_latent, alt_file(OUT_DIR_H5A, "h5a_overall_interface_contrast_latent.csv"))
save_csv(h5a_emm_meanclass, alt_file(OUT_DIR_H5A, "h5a_emm_interface_meanclass.csv"))
save_csv(h5a_contrast_meanclass, alt_file(OUT_DIR_H5A, "h5a_overall_interface_contrast_meanclass.csv"))
save_csv(h5a_prob_interface, alt_file(OUT_DIR_H5A, "h5a_category_probabilities_by_interface.csv"))

writeLines(
  c(
    "Random-effect extraction diagnostics - H5A ALTERNATIVE",
    paste0("rows: ", nrow(h5a_random)),
    if (nrow(h5a_random) > 0) paste0("source(s): ", paste(unique(h5a_random$source), collapse = ", ")) else "source(s): none"
  ),
  con = alt_file(OUT_DIR_H5A, "h5a_random_effect_sd_diagnostics.txt")
)

# ============================================================
# H5B
# Model: y ~ industry_role + interface * task_complexity + (1|participant)
# ============================================================
h5b_fit_tbl <- list()
h5b_lrt_tbl <- list()
h5b_coef_all_tbl <- list()
h5b_coef_location_tbl <- list()
h5b_coef_main_tbl <- list()
h5b_coef_industry_role_tbl <- list()
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
    paste(
      "y_ord ~",
      INDUSTRY_ROLE_COL,
      "+",
      COND_COL,
      "+",
      COMPLEXITY_COL,
      "+ (1|",
      PARTICIPANT_COL,
      ")"
    )
  )

  form_int <- as.formula(
    paste(
      "y_ord ~",
      INDUSTRY_ROLE_COL,
      "+",
      COND_COL,
      "*",
      COMPLEXITY_COL,
      "+ (1|",
      PARTICIPANT_COL,
      ")"
    )
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
    file = alt_file(OUT_DIR_H5B, paste0("h5b_", item, "_additive_model_summary.txt"))
  )

  capture.output(
    summary(m_int),
    file = alt_file(OUT_DIR_H5B, paste0("h5b_", item, "_interaction_model_summary.txt"))
  )

  fit_tbl_item <- tibble(
    item_key = item,
    model  = c("additive_adjusted_for_role", "interaction_adjusted_for_role"),
    nobs   = c(nobs(m_add), nobs(m_int)),
    logLik = c(as.numeric(logLik(m_add)), as.numeric(logLik(m_int))),
    AIC    = c(AIC(m_add), AIC(m_int)),
    BIC    = c(BIC(m_add), BIC(m_int))
  )

  h5b_fit_tbl[[item]] <- fit_tbl_item

  lrt_interaction <- as.data.frame(anova(m_add, m_int))
  lrt_interaction$item_key <- item
  lrt_interaction$effect <- "interface_x_complexity_interaction_adjusted_for_industry_role"

  h5b_lrt_tbl[[item]] <- as_tibble(lrt_interaction)

  coef_int_all <- make_coef_table(m_int, conf.level = CONF_LEVEL) %>%
    mutate(item_key = item, .before = 1)

  coef_int_location <- coef_int_all %>%
    filter(component == "location")

  coef_int_main_only <- coef_int_location %>%
    filter(!str_detect(term, ":"))

  coef_int_industry_role <- coef_int_location %>%
    filter(str_detect(term, INDUSTRY_ROLE_COL)) %>%
    add_holm_correction(raw_p_col = "p.value.raw")

  coef_int_interaction_only <- coef_int_location %>%
    filter(str_detect(term, ":")) %>%
    add_holm_correction(raw_p_col = "p.value.raw") %>%
    add_did_direction_flags()

  h5b_coef_all_tbl[[item]] <- coef_int_all
  h5b_coef_location_tbl[[item]] <- coef_int_location
  h5b_coef_main_tbl[[item]] <- coef_int_main_only
  h5b_coef_industry_role_tbl[[item]] <- coef_int_industry_role
  h5b_coef_interaction_tbl[[item]] <- coef_int_interaction_only

  rand_int <- extract_random_effect_sd(m_int, item, "interaction_adjusted_for_role")
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

  h5b_prob_cells_tbl[[item]] <- confint(
    emm_prob_cells,
    level = CONF_LEVEL
  ) %>%
    standardize_emm_ci_cols() %>%
    mutate(item_key = item, .before = 1)
}

h5b_model_fit <- bind_rows(h5b_fit_tbl)
h5b_lrt <- bind_rows(h5b_lrt_tbl)
h5b_coef_all <- bind_rows(h5b_coef_all_tbl)
h5b_coef_loc <- bind_rows(h5b_coef_location_tbl)
h5b_coef_main <- bind_rows(h5b_coef_main_tbl)
h5b_coef_industry_role <- bind_rows(h5b_coef_industry_role_tbl)
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

save_csv(h5b_model_fit, alt_file(OUT_DIR_H5B, "h5b_model_fit.csv"))
save_csv(h5b_lrt, alt_file(OUT_DIR_H5B, "h5b_interaction_lrt.csv"))
capture.output(h5b_lrt, file = alt_file(OUT_DIR_H5B, "h5b_interaction_lrt.txt"))

save_csv(h5b_coef_all, alt_file(OUT_DIR_H5B, "h5b_interaction_coef_all.csv"))
save_csv(h5b_coef_loc, alt_file(OUT_DIR_H5B, "h5b_interaction_coef_location.csv"))
save_csv(h5b_coef_main, alt_file(OUT_DIR_H5B, "h5b_interaction_main_terms_only.csv"))
save_csv(h5b_coef_industry_role, alt_file(OUT_DIR_H5B, "h5b_interaction_industry_role_terms_only.csv"))
save_csv(h5b_coef_interaction, alt_file(OUT_DIR_H5B, "h5b_interaction_terms_only.csv"))
save_csv(h5b_coef_interaction, alt_file(OUT_DIR_H5B, "h5b_internal_did_terms_only.csv"))

save_csv(h5b_random, alt_file(OUT_DIR_H5B, "h5b_random_effect_sd.csv"))

writeLines(
  c(
    "Random-effect extraction diagnostics - H5B ALTERNATIVE",
    paste0("rows: ", nrow(h5b_random)),
    if (nrow(h5b_random) > 0) paste0("source(s): ", paste(unique(h5b_random$source), collapse = ", ")) else "source(s): none"
  ),
  con = alt_file(OUT_DIR_H5B, "h5b_random_effect_sd_diagnostics.txt")
)

save_csv(
  h5b_emm_interface_by_comp_latent,
  alt_file(OUT_DIR_H5B, "h5b_emm_interface_by_complexity_latent.csv")
)

save_csv(
  h5b_simple_interface_by_comp_latent,
  alt_file(OUT_DIR_H5B, "h5b_simple_effects_interface_by_complexity_latent.csv")
)

save_csv(
  h5b_emm_interface_by_comp_meanclass,
  alt_file(OUT_DIR_H5B, "h5b_emm_interface_by_complexity_meanclass.csv")
)

save_csv(
  h5b_simple_interface_by_comp_meanclass,
  alt_file(OUT_DIR_H5B, "h5b_simple_effects_interface_by_complexity_meanclass.csv")
)

save_csv(
  h5b_emm_comp_by_interface_latent,
  alt_file(OUT_DIR_H5B, "h5b_emm_complexity_by_interface_latent.csv")
)

save_csv(
  h5b_simple_comp_by_interface_latent,
  alt_file(OUT_DIR_H5B, "h5b_simple_effects_complexity_by_interface_latent.csv")
)

save_csv(
  h5b_emm_comp_by_interface_meanclass,
  alt_file(OUT_DIR_H5B, "h5b_emm_complexity_by_interface_meanclass.csv")
)

save_csv(
  h5b_simple_comp_by_interface_meanclass,
  alt_file(OUT_DIR_H5B, "h5b_simple_effects_complexity_by_interface_meanclass.csv")
)

save_csv(h5b_did_like_latent, alt_file(OUT_DIR_H5B, "h5b_did_like_latent.csv"))
save_csv(h5b_did_like_meanclass, alt_file(OUT_DIR_H5B, "h5b_did_like_meanclass.csv"))
save_csv(h5b_all_cells_latent, alt_file(OUT_DIR_H5B, "h5b_all_cells_latent.csv"))
save_csv(h5b_all_cells_meanclass, alt_file(OUT_DIR_H5B, "h5b_all_cells_meanclass.csv"))
save_csv(h5b_prob_cells, alt_file(OUT_DIR_H5B, "h5b_category_probabilities_by_cell.csv"))

cat("\n============================================================\n")
cat("H5 reliance CLMM ALTERNATIVE completata: adjusted for industry_role.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("All the file have suffisso _ALTERNATIVE.\n")
cat("============================================================\n\n")

cat("H5A adjusted model:\n")
cat("y_ord ~ industry_role + interface_cond + task_complexity + (1|participant_id)\n\n")

cat("H5B adjusted model:\n")
cat("y_ord ~ industry_role + interface_cond * task_complexity + (1|participant_id)\n\n")

cat("Riferimenti dei fattori:\n")
cat("- interface_cond reference =", REF_COND, "\n")
cat("- task_complexity reference =", REF_COMPLEXITY, "\n")
cat("- industry_role reference =", REF_ROLE, "\n\n")

cat("File principali H5A ALTERNATIVE:\n")
cat("- H5A_interface_ALTERNATIVE/h5a_interface_lrt_ALTERNATIVE.csv\n")
cat("- H5A_interface_ALTERNATIVE/h5a_additive_coef_location_ALTERNATIVE.csv\n")
cat("- H5A_interface_ALTERNATIVE/h5a_additive_industry_role_terms_only_ALTERNATIVE.csv\n")
cat("- H5A_interface_ALTERNATIVE/h5a_overall_interface_contrast_latent_ALTERNATIVE.csv\n")
cat("- H5A_interface_ALTERNATIVE/h5a_category_probabilities_by_interface_ALTERNATIVE.csv\n\n")

cat("File principali H5B ALTERNATIVE:\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_interaction_lrt_ALTERNATIVE.csv\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_interaction_terms_only_ALTERNATIVE.csv\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_interaction_industry_role_terms_only_ALTERNATIVE.csv\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_simple_effects_interface_by_complexity_latent_ALTERNATIVE.csv\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_did_like_latent_ALTERNATIVE.csv\n")
cat("- H5B_interaction_ALTERNATIVE/h5b_category_probabilities_by_cell_ALTERNATIVE.csv\n\n")

cat("Direction of the contrasts of interface: dashboard_minus_chatbot\n")
cat("estimate > 0 => dashboard higher reliance => chatbot lower reliance\n\n")

cat("Direction of the DiD-like:\n")
cat("estimate > 0 => the dashboard–chatbot gap is larger at the level of complexity higher of the contrast\n")