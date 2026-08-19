# ============================================================================
# COMPASS shared utility functions
# Sourced by all analysis scripts to avoid duplication.
# ============================================================================

normalize_gene <- function(x) toupper(trimws(as.character(x)))

normalize_tissue <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("_", "-", z, fixed = TRUE)
  z <- gsub("[[:space:]]+", "-", z)
  z <- gsub("-+", "-", z)
  z[z == "testis"] <- "testicle"
  z[z == "head-and neck"] <- "head-and-neck"
  z
}

read_gmt_genes <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nchar(ln) > 0][1]
  normalize_gene(strsplit(ln, "\t")[[1]][-(1:2)])
}

load_rda_list <- function(path) {
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  mget(nm, envir = e, inherits = FALSE)
}

average_precision_grouped <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- score[ok]
  npos <- sum(y == 1)
  if (npos == 0) return(NA_real_)

  o <- order(score, decreasing = TRUE)
  y <- y[o]
  score <- score[o]

  grp <- cumsum(c(TRUE, diff(score) != 0))
  pos_by <- as.numeric(tapply(y, grp, sum))
  n_by <- as.numeric(tapply(y, grp, length))

  cum_pos <- cumsum(pos_by)
  cum_n <- cumsum(n_by)
  precision <- cum_pos / cum_n
  recall <- cum_pos / npos
  delta_recall <- c(recall[1], diff(recall))

  sum(delta_recall * precision)
}

require_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "))
  }
}

rna_conf_score <- function(x) {
  dplyr::case_when(
    tolower(as.character(x)) == "high"   ~ 1.0,
    tolower(as.character(x)) == "medium" ~ 0.4,
    tolower(as.character(x)) == "low"    ~ 0.2,
    TRUE                                 ~ 0
  )
}

chip_bind_score <- function(x) {
  dplyr::case_when(
    as.character(x) == "High confidence binding"   ~ 1.0,
    as.character(x) == "Medium-high confidence"     ~ 0.8,
    as.character(x) == "Medium confidence"          ~ 0.4,
    as.character(x) == "Low-medium confidence"      ~ 0.3,
    as.character(x) == "Low confidence"             ~ 0.15,
    TRUE                                            ~ 0
  )
}
