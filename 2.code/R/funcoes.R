# ============================================================================
# P06 - Funções auxiliares do pipeline
# ============================================================================

# Leitura/gravação. Por defeito usa arrow; options(p06.io = "csv") troca para CSV
# (usado só para testar o pipeline num ambiente sem arrow).
ler_parquet <- function(path, cols = NULL) {
  if (identical(getOption("p06.io", "arrow"), "csv")) {
    d <- utils::read.csv(sub("\\.parquet$", ".csv", path), stringsAsFactors = FALSE,
                         encoding = "UTF-8", na.strings = c("", "NA"))
    d <- tibble::as_tibble(d)
    if (!is.null(cols)) d <- dplyr::select(d, dplyr::any_of(cols))
    return(d)
  }
  if (is.null(cols)) arrow::read_parquet(path)
  else arrow::read_parquet(path, col_select = dplyr::any_of(cols))
}

gravar_parquet <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (identical(getOption("p06.io", "arrow"), "csv")) {
    utils::write.csv(x, sub("\\.parquet$", ".csv", path), na = "", row.names = FALSE, fileEncoding = "UTF-8")
    # guarda a ordem dos níveis dos factores para reconstruir o Parquet
    niveis <- lapply(Filter(is.factor, x), levels)
    jsonlite::write_json(niveis, sub("\\.parquet$", "_niveis.json", path), auto_unbox = FALSE)
  } else {
    arrow::write_parquet(x, path)
  }
  invisible(path)
}

# Códigos DHS de "não sabe / em falta" (8, 9, 98, 99, 998, 999) -> NA
na_dhs <- function(x, codigos = c(8, 9)) {
  x <- as.numeric(x)
  dplyr::if_else(x %in% codigos, NA_real_, x)
}

# Proporção ponderada simples (para controlo de qualidade; não substitui survey)
prop_ponderada <- function(y, w) sum(y * w, na.rm = TRUE) / sum(w[!is.na(y)])

# Pára o pipeline com mensagem clara
verificar <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
