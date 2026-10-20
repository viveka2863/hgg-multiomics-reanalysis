# Supplement round 1: one-time conversion of the needed sheets of the paper's supplementary workbooks to CSV.
# Run from the repo root: Rscript analysis/C2_differential/improved/c2b_convert_xlsx.R
# Input (gitignored, not committed): data/external/*.xlsx. Output (gitignored): data/external/csv/*.csv.
# Everything is read as text so full numeric precision is kept; header rows are kept (parsed in the c2b/c5b/c7b scripts).
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(readxl))
stopifnot(file.exists("CLAUDE.md"))
dir.create("data/external/csv", recursive = TRUE, showWarnings = FALSE)
wb <- function(n) sprintf("data/external/41467_2019_11661_%s_ESM.xlsx", n)
conv <- function(file, sheet, out) {
  stopifnot(file.exists(wb(file)))
  x <- as.data.frame(suppressMessages(read_excel(wb(file), sheet = sheet, col_names = FALSE, col_types = "text")))
  write.table(x, file.path("data/external/csv", out), sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE, na = "")
  cat(sprintf("%-6s sheet %-14s -> %-10s %6d x %2d\n", file, sheet, out, nrow(x), ncol(x)))
}
conv("MOESM4", 1, "sd1.csv")       # Supp Data 1 (sheet name has a trailing space, so read by index)
conv("MOESM5", 1, "sd2a.csv")      # Supp Data 2a (sheet 2b not needed)
conv("MOESM7", "4e", "t4e.csv")    # kinase activity ANOVA, 90 kinases
conv("MOESM7", "4h", "t4h.csv")    # a(P) inputs, 22 sites
conv("MOESM9", "4a", "src4a.csv")  # Source Data Fig 4a, 41 kinases
