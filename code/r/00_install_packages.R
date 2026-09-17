repos <- "https://cloud.r-project.org"

pkgs <- c(
  "ordinal",
  "emmeans",
  "readr",
  "dplyr",
  "stringr",
  "tibble"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))

if (length(to_install) > 0) {
  install.packages(to_install, repos = repos, dependencies = TRUE)
}

invisible(lapply(pkgs, library, character.only = TRUE))

cat("\nPacchetti installati e caricati correttamente.\n")
sessionInfo()