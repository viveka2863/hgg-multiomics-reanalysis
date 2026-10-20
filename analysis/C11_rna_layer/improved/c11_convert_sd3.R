# C11: one-time conversion of Supp Data 3 (MOESM6) to CSV, same approach as C2/c2b_convert_xlsx.R (which is left unmodified).
# Run from the repo root: Rscript analysis/C11_rna_layer/improved/c11_convert_sd3.R. Output (gitignored): data/external/csv/sd3.csv
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(readxl))
stopifnot(file.exists("CLAUDE.md"))
f <- "data/external/41467_2019_11661_MOESM6_ESM.xlsx"; stopifnot(file.exists(f))
print(excel_sheets(f))
x <- as.data.frame(suppressMessages(read_excel(f, sheet = 1, col_names = FALSE, col_types = "text")))
dir.create("data/external/csv", recursive = TRUE, showWarnings = FALSE)
write.table(x, "data/external/csv/sd3.csv", sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE, na = "")
cat(sprintf("sd3.csv %d x %d\n", nrow(x), ncol(x)))
