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
# H3a ALTERNATIVE
# Data Literacy x Interface -> Mental Workload
#
# Alternative/confounder model:
#   industry_role is added as a covariate/confounder.
#
# Original H3a:
#   Data literacy moderates the interface effect on MWL.
#
# Alternative H3a:
#   Data literacy moderates the interface effect on MWL
#   AFTER controlling for industry_role.
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

CSV_PATH <- file.path(SCRIPT_DIR, "h3a_mwl_dl_clmm.csv")
OUT_DIR  <- file.path(SCRIPT_DIR, "results_H3_ALTERNATIVE")

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
TASK_COMPLEXITY_COL <- "task_complexity"
TASK_ID_COL     <- "task_id"
Y_RAW_COL       <- "tlx_mean"

DL_MEAN_COL     <- "bdli_mean"
DL_N_COL        <- "bdli_n_answered"

INDUSTRY_ROLE_COL      <- "industry_role"
INDUSTRY_ROLE_CODE_COL <- "industry_role_code"

REF_COND        <- "chatbot"
OTHER_COND      <- "dashboard"

TASK_COMPLEXITY_LEVELS <- c("easy", "mid", "hard")
REF_COMPLEXITY <- "easy"

DL_LEVELS <- c("low", "medium", "high")
REF_DL    <- "low"

ROLE_ORDER <- c(
  "Junior Management",
  "Middle Management",
  "Upper Management"
)
REF_ROLE <- "Junior Management"

REQUIRE_ALL_BDLI_ITEMS <- TRUE
BDLI_ITEMS_REQUIRED    <- 12

ROUND_TO       <- 5
LINK_FUN       <- "probit"
CONF_LEVEL     <- 0.95
EMM_WEIGHTS    <- "equal"
NAGQ_SINGLE_RE <- 10
ALPHA          <- 0.05

# ============================================================
# HELPERS
# ============================================================
alt_file <- function(filename) {
  file.path(
    OUT_DIR,
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
      paste0(TASK_COMPLEXITY_COL, "hard"),
      paste0(INDUSTRY_ROLE_COL, "Middle Management"),
      paste0(INDUSTRY_ROLE_COL, "Upper Management")
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
      "beta7",
      "control_role_middle",
      "control_role_upper"
    ),
    parameter_meaning = c(
      "dashboard vs chatbot at low DL, adjusted for industry_role",
      "medium vs low DL within chatbot, adjusted for industry_role",
      "high vs low DL within chatbot, adjusted for industry_role",
      "change in dashboard-chatbot gap at medium DL vs low DL, adjusted for industry_role",
      "change in dashboard-chatbot gap at medium DL vs low DL, adjusted for industry_role",
      "change in dashboard-chatbot gap at high DL vs low DL, adjusted for industry_role",
      "change in dashboard-chatbot gap at high DL vs low DL, adjusted for industry_role",
      "task complexity mid vs easy, adjusted for interface, DL, and industry_role",
      "task complexity hard vs easy, adjusted for interface, DL, and industry_role",
      "Middle Management vs Junior Management, control effect",
      "Upper Management vs Junior Management, control effect"
    )
  ) %>%
    distinct(term, .keep_all = TRUE)
}

make_parameter_label_export <- function() {
  tibble(
    symbolic_label = c(
      "beta1",
      "beta2",
      "beta3",
      "beta4",
      "beta5",
      "beta6",
      "beta7",
      "control_role_middle",
      "control_role_upper"
    ),
    model_term = c(
      paste0(COND_COL, OTHER_COND),
      "dl_ordmedium",
      "dl_ordhigh",
      paste0(COND_COL, OTHER_COND, ":dl_ordmedium"),
      paste0(COND_COL, OTHER_COND, ":dl_ordhigh"),
      paste0(TASK_COMPLEXITY_COL, "mid"),
      paste0(TASK_COMPLEXITY_COL, "hard"),
      paste0(INDUSTRY_ROLE_COL, "Middle Management"),
      paste0(INDUSTRY_ROLE_COL, "Upper Management")
    ),
    meaning = c(
      "dashboard vs chatbot at low DL, adjusted for industry_role",
      "medium vs low DL within chatbot, adjusted for industry_role",
      "high vs low DL within chatbot, adjusted for industry_role",
      "change in dashboard-chatbot gap at medium DL vs low DL, adjusted for industry_role",
      "change in dashboard-chatbot gap at high DL vs low DL, adjusted for industry_role",
      "task complexity mid vs easy, adjusted for interface, DL, and industry_role",
      "task complexity hard vs easy, adjusted for interface, DL, and industry_role",
      "Middle Management vs Junior Management, control effect",
      "Upper Management vs Junior Management, control effect"
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
df <- add_industry_role_if_needed(df)

needed_cols <- c(
  PARTICIPANT_COL,
  COND_COL,
  TASK_COMPLEXITY_COL,
  TASK_ID_COL,
  Y_RAW_COL,
  DL_MEAN_COL,
  DL_N_COL,
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
    !!PARTICIPANT_COL     := as.character(.data[[PARTICIPANT_COL]]),
    !!COND_COL            := str_to_lower(str_trim(as.character(.data[[COND_COL]]))),
    !!TASK_COMPLEXITY_COL := str_to_lower(str_trim(as.character(.data[[TASK_COMPLEXITY_COL]]))),
    !!TASK_ID_COL         := as.character(.data[[TASK_ID_COL]]),
    !!Y_RAW_COL           := as.numeric(.data[[Y_RAW_COL]]),
    !!DL_MEAN_COL         := as.numeric(.data[[DL_MEAN_COL]]),
    !!DL_N_COL            := as.numeric(.data[[DL_N_COL]]),
    !!INDUSTRY_ROLE_COL   := str_squish(as.character(.data[[INDUSTRY_ROLE_COL]]))
  ) %>%
  filter(
    !is.na(.data[[PARTICIPANT_COL]]),
    !is.na(.data[[COND_COL]]),
    !is.na(.data[[TASK_COMPLEXITY_COL]]),
    !is.na(.data[[TASK_ID_COL]]),
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

if (REQUIRE_ALL_BDLI_ITEMS) {
  df <- df %>% filter(.data[[DL_N_COL]] >= BDLI_ITEMS_REQUIRED)
}

df <- df %>% filter(!is.na(.data[[DL_MEAN_COL]]))

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

df[[TASK_COMPLEXITY_COL]] <- factor(
  df[[TASK_COMPLEXITY_COL]],
  levels = TASK_COMPLEXITY_LEVELS
)

df[[INDUSTRY_ROLE_COL]] <- factor(
  df[[INDUSTRY_ROLE_COL]],
  levels = ROLE_ORDER
)

df[[INDUSTRY_ROLE_CODE_COL]] <- as.integer(df[[INDUSTRY_ROLE_COL]]) - 1

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

save_csv(dl_map, alt_file("h3a_dl_group_mapping.csv"))
save_csv(make_parameter_label_export(), alt_file("h3a_parameter_label_mapping.csv"))

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

save_csv(df, alt_file("h3a_cleaned_analysis_data.csv"))

# ============================================================
# INDUSTRY ROLE DIAGNOSTICS
# ============================================================
role_participant <- df %>%
  group_by(
    .data[[PARTICIPANT_COL]],
    .data[[INDUSTRY_ROLE_COL]],
    .data[[INDUSTRY_ROLE_CODE_COL]],
    .data[[COND_COL]],
    dl_ord
  ) %>%
  summarise(
    tlx_mean_participant = mean(.data[[Y_RAW_COL]], na.rm = TRUE),
    bdli_mean_participant = mean(.data[[DL_MEAN_COL]], na.rm = TRUE),
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

role_by_dl <- role_participant %>%
  count(dl_ord, .data[[INDUSTRY_ROLE_COL]], name = "n") %>%
  group_by(dl_ord) %>%
  mutate(prop_within_dl = n / sum(n)) %>%
  ungroup()

role_by_interface_dl <- role_participant %>%
  count(.data[[COND_COL]], dl_ord, .data[[INDUSTRY_ROLE_COL]], name = "n") %>%
  group_by(.data[[COND_COL]], dl_ord) %>%
  mutate(prop_within_interface_dl = n / sum(n)) %>%
  ungroup()

desc_raw_by_role <- role_participant %>%
  group_by(.data[[INDUSTRY_ROLE_COL]]) %>%
  summarise(
    n      = n(),
    mean   = mean(tlx_mean_participant, na.rm = TRUE),
    sd     = sd(tlx_mean_participant, na.rm = TRUE),
    median = median(tlx_mean_participant, na.rm = TRUE),
    iqr    = iqr_num(tlx_mean_participant),
    .groups = "drop"
  )

role_kw <- kruskal.test(
  tlx_mean_participant ~ industry_role,
  data = role_participant
)

role_spearman <- suppressWarnings(
  cor.test(
    role_participant$industry_role_code,
    role_participant$tlx_mean_participant,
    method = "spearman",
    exact = FALSE
  )
)

role_interface_chisq <- suppressWarnings(
  chisq.test(table(role_participant[[INDUSTRY_ROLE_COL]], role_participant[[COND_COL]]))
)

role_dl_chisq <- suppressWarnings(
  chisq.test(table(role_participant[[INDUSTRY_ROLE_COL]], role_participant$dl_ord))
)

industry_role_omnibus <- tibble(
  analysis = c(
    "Industry role vs participant-level NASA-TLX",
    "Ordinal trend: industry role code vs participant-level NASA-TLX",
    "Industry role balance across interface condition",
    "Industry role balance across DL groups"
  ),
  test = c(
    "Kruskal-Wallis",
    "Spearman rank correlation",
    "Chi-square",
    "Chi-square"
  ),
  statistic_name = c("H", "rho", "X-squared", "X-squared"),
  statistic = c(
    as.numeric(role_kw$statistic),
    as.numeric(role_spearman$estimate),
    as.numeric(role_interface_chisq$statistic),
    as.numeric(role_dl_chisq$statistic)
  ),
  df = c(
    as.numeric(role_kw$parameter),
    NA_real_,
    as.numeric(role_interface_chisq$parameter),
    as.numeric(role_dl_chisq$parameter)
  ),
  p.value = c(
    as.numeric(role_kw$p.value),
    as.numeric(role_spearman$p.value),
    as.numeric(role_interface_chisq$p.value),
    as.numeric(role_dl_chisq$p.value)
  ),
  effect_name = c(
    "epsilon_squared",
    "rho",
    NA_character_,
    NA_character_
  ),
  effect_size = c(
    epsilon_sq_kw(
      role_kw$statistic,
      n = nrow(role_participant),
      k = n_distinct(role_participant[[INDUSTRY_ROLE_COL]])
    ),
    as.numeric(role_spearman$estimate),
    NA_real_,
    NA_real_
  ),
  n = c(
    nrow(role_participant),
    nrow(role_participant),
    nrow(role_participant),
    nrow(role_participant)
  )
)

save_csv(role_counts, alt_file("h3a_industry_role_counts.csv"))
save_csv(role_by_interface, alt_file("h3a_industry_role_by_interface.csv"))
save_csv(role_by_dl, alt_file("h3a_industry_role_by_dl.csv"))
save_csv(role_by_interface_dl, alt_file("h3a_industry_role_by_interface_dl.csv"))
save_csv(desc_raw_by_role, alt_file("h3a_descriptives_raw_by_industry_role.csv"))
save_csv(industry_role_omnibus, alt_file("h3a_industry_role_omnibus.csv"))

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

save_csv(desc_raw_by_cell, alt_file("h3a_descriptives_raw_by_cell.csv"))
save_csv(desc_ord_by_cell, alt_file("h3a_descriptives_ordinal_by_cell.csv"))
save_csv(desc_n_by_interface_dl, alt_file("h3a_n_participants_by_interface_dl.csv"))

# ============================================================
# MODELS
# ============================================================
options(contrasts = c("contr.treatment", "contr.treatment"))

form_add <- as.formula(
  paste(
    "y_ord ~",
    INDUSTRY_ROLE_COL,
    "+",
    COND_COL,
    "+ dl_ord +",
    TASK_COMPLEXITY_COL,
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
    "* dl_ord +",
    TASK_COMPLEXITY_COL,
    "+ (1|",
    PARTICIPANT_COL,
    ")"
  )
)

m_add <- clmm(
  formula = form_add,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = NAGQ_SINGLE_RE
)

m_int <- clmm(
  formula = form_int,
  data    = df,
  link    = LINK_FUN,
  Hess    = TRUE,
  nAGQ    = NAGQ_SINGLE_RE
)

capture.output(summary(m_add), file = alt_file("h3a_additive_model_summary.txt"))
capture.output(summary(m_int), file = alt_file("h3a_interaction_model_summary.txt"))

# ============================================================
# MODEL FIT / LRT
# ============================================================
fit_tbl <- tibble(
  model  = c("additive_adjusted_for_industry_role", "interaction_adjusted_for_industry_role"),
  nobs   = c(nobs(m_add), nobs(m_int)),
  logLik = c(as.numeric(logLik(m_add)), as.numeric(logLik(m_int))),
  AIC    = c(AIC(m_add), AIC(m_int)),
  BIC    = c(BIC(m_add), BIC(m_int))
)

save_csv(fit_tbl, alt_file("h3a_model_fit.csv"))

lrt_tbl <- as.data.frame(anova(m_add, m_int))
save_csv(lrt_tbl, alt_file("h3a_interaction_lrt.csv"))
capture.output(anova(m_add, m_int), file = alt_file("h3a_interaction_lrt.txt"))

# ============================================================
# COEFFICIENTS
# ============================================================
coef_add_all <- make_coef_table(m_add, conf.level = CONF_LEVEL)
coef_int_all <- make_coef_table(m_int, conf.level = CONF_LEVEL)

coef_add_location <- coef_add_all %>% filter(component == "location")
coef_int_location <- coef_int_all %>% filter(component == "location")

coef_int_main_terms_only <- coef_int_location %>%
  filter(!str_detect(term, ":"))

coef_int_industry_role_terms_only <- coef_int_location %>%
  filter(str_detect(term, INDUSTRY_ROLE_COL)) %>%
  add_holm_correction(group_cols = NULL, raw_p_col = "p.value.raw")

coef_int_interaction_only <- coef_int_location %>%
  filter(str_detect(term, ":")) %>%
  add_holm_correction(group_cols = NULL, raw_p_col = "p.value.raw") %>%
  add_did_direction_flags()

save_csv(coef_add_all, alt_file("h3a_additive_coef_all.csv"))
save_csv(coef_int_all, alt_file("h3a_interaction_coef_all.csv"))
save_csv(coef_add_location, alt_file("h3a_additive_coef_location.csv"))
save_csv(coef_int_location, alt_file("h3a_interaction_coef_location.csv"))
save_csv(coef_int_main_terms_only, alt_file("h3a_interaction_main_terms_only.csv"))
save_csv(coef_int_industry_role_terms_only, alt_file("h3a_interaction_industry_role_terms_only.csv"))
save_csv(coef_int_interaction_only, alt_file("h3a_interaction_terms_only.csv"))
save_csv(coef_int_interaction_only, alt_file("h3a_internal_did_terms_only.csv"))

# ============================================================
# RANDOM EFFECT SD
# ============================================================
rand_add <- extract_random_effect_sd(m_add, "additive_adjusted")
rand_int <- extract_random_effect_sd(m_int, "interaction_adjusted")
rand_tbl <- bind_rows(rand_add, rand_int)

if (nrow(rand_tbl) == 0) {
 warning("These are not riuscito to estrarre the random effects; controlla the summary testuali.")
}

save_csv(rand_tbl, alt_file("h3a_random_effect_sd.csv"))

writeLines(
  c(
    "Random-effect extraction diagnostics - H3a ALTERNATIVE adjusted model",
    paste0("additive rows: ", nrow(rand_add)),
    paste0("interaction rows: ", nrow(rand_int)),
    if (nrow(rand_add) > 0) paste0("additive source(s): ", paste(unique(rand_add$source), collapse = ", ")) else "additive source(s): none",
    if (nrow(rand_int) > 0) paste0("interaction source(s): ", paste(unique(rand_int$source), collapse = ", ")) else "interaction source(s): none"
  ),
  con = alt_file("h3a_random_effect_sd_diagnostics.txt")
)

# ============================================================
# SIMPLE EFFECTS OF INTERFACE WITHIN EACH DL LEVEL
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
  alt_file("h3a_emm_interface_by_dl_latent.csv")
)

save_csv(
  simple_effects_interface_latent,
  alt_file("h3a_simple_effects_interface_by_dl_latent.csv")
)

save_csv(
  summary(emm_interface_by_dl_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  alt_file("h3a_emm_interface_by_dl_meanclass.csv")
)

save_csv(
  simple_effects_interface_meanclass,
  alt_file("h3a_simple_effects_interface_by_dl_meanclass.csv")
)

# ============================================================
# SIMPLE EFFECTS OF DL HIGH - LOW WITHIN EACH INTERFACE
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
  alt_file("h3a_emm_dl_by_interface_latent.csv")
)

save_csv(
  simple_effects_dl_lowhigh_latent,
  alt_file("h3a_simple_effects_dl_lowhigh_by_interface_latent.csv")
)

save_csv(
  summary(emm_dl_by_interface_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  alt_file("h3a_emm_dl_by_interface_meanclass.csv")
)

save_csv(
  simple_effects_dl_lowhigh_meanclass,
  alt_file("h3a_simple_effects_dl_lowhigh_by_interface_meanclass.csv")
)

# ============================================================
# DID-LIKE CONTRAST ON HIGH-LOW DL EFFECT
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

save_csv(did_like_latent, alt_file("h3a_did_like_lowhigh_latent.csv"))
save_csv(did_like_meanclass, alt_file("h3a_did_like_lowhigh_meanclass.csv"))

# ============================================================
# ALL CELL MEANS
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
  alt_file("h3a_all_cells_latent.csv")
)

save_csv(
  summary(emm_all_cells_meanclass, level = CONF_LEVEL) %>% standardize_emm_ci_cols(),
  alt_file("h3a_all_cells_meanclass.csv")
)

# ============================================================
# CATEGORY PROBABILITIES BY CELL
# ============================================================
emm_prob_cells <- emmeans(
  m_int,
  specs   = as.formula(paste("~ y_ord |", COND_COL, "* dl_ord")),
  mode    = "prob",
  weights = EMM_WEIGHTS
)

prob_cells_tbl <- confint(emm_prob_cells, level = CONF_LEVEL) %>%
  standardize_emm_ci_cols()

save_csv(prob_cells_tbl, alt_file("h3a_category_probabilities_by_cell.csv"))

# ============================================================
# OVERALL INTERFACE EFFECT FROM FULL MODEL
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
  alt_file("h3a_overall_interface_contrast_latent_from_full_model.csv")
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("H3a ALTERNATIVE completata: adjusted model for industry_role.\n")
cat("Output salvati in:", OUT_DIR, "\n")
cat("All the file have suffisso _ALTERNATIVE.\n")
cat("============================================================\n\n")

cat("Additive model adjusted:\n")
cat("y_ord ~ industry_role + interface_cond + dl_ord + task_complexity + (1|participant_id)\n\n")

cat("Interaction model adjusted:\n")
cat("y_ord ~ industry_role + interface_cond * dl_ord + task_complexity + (1|participant_id)\n\n")

cat("Riferimenti dei fattori:\n")
cat("- interface_cond reference =", REF_COND, "\n")
cat("- task_complexity reference =", REF_COMPLEXITY, "\n")
cat("- dl_ord reference =", REF_DL, "\n")
cat("- industry_role reference =", REF_ROLE, "\n\n")

cat("Come cambia H3a:\n")
cat("- Before: the moderation interface x data literacy was estimated without control for industry_role.\n")
cat("- Ora: the moderation interface x data literacy must rimanere after controlling for industry_role.\n")
cat("- industry_role is to covariates/confounder, not to moderator.\n\n")

cat("Pattern expected for H3a adjusted:\n")
cat("- the dashboard–chatbot gap should be largest to DL low\n")
cat("- the contrast high_minus_low should be more negative in the dashboard\n")
cat("- the DID-like on the contrast high_minus_low should be < 0\n")
cat(" => low DL penalizza more the dashboard that the chatbot, to the netto of industry_role\n\n")

cat("File principali ALTERNATIVE:\n")
cat("- h3a_interaction_lrt_ALTERNATIVE.csv\n")
cat("- h3a_interaction_terms_only_ALTERNATIVE.csv\n")
cat("- h3a_interaction_industry_role_terms_only_ALTERNATIVE.csv\n")
cat("- h3a_simple_effects_interface_by_dl_latent_ALTERNATIVE.csv\n")
cat("- h3a_simple_effects_dl_lowhigh_by_interface_latent_ALTERNATIVE.csv\n")
cat("- h3a_did_like_lowhigh_latent_ALTERNATIVE.csv\n")
cat("- h3a_industry_role_omnibus_ALTERNATIVE.csv\n")
cat("- h3a_category_probabilities_by_cell_ALTERNATIVE.csv\n")