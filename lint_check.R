x <- lintr::lint_package()
df <- as.data.frame(x)
exclude_pat <- "IWTimage|ITPaov|ITPlm|ITP1|ITP2|ITPimage|data\\.R"
src <- df[grepl("^R/", df$filename) & !grepl(exclude_pat, df$filename), ]
src <- src[order(src$filename, src$line_number), ]
if (nrow(src) == 0) {
  cat("No lint in main R/ source files\n")
} else {
  out <- paste0(
    src$filename,
    ":",
    src$line_number,
    " [",
    src$linter,
    "] ",
    src$message
  )
  cat(out, sep = "\n")
}
