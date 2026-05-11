library(openxlsx)

projdir <- "/Users/gracekenney/Documents/UNC/Research/AGAL"
dat <- read.csv(file.path(projdir, "results/RNA_all_res.csv"))

dat <- dat[, c(colnames(dat)[1:10], sort(colnames(dat)[11:22]))]

# get comp names
comp_names <- colnames(dat)[grepl(pattern = "log2", x = colnames(dat))] |>
  sub(pattern = "log2FC_", replacement = "") |>
  sub(pattern = "_vs_", replacement = " vs ") |>
  gsub(pattern = "C", replacement = "WT")

# get sample names
samp_names <- colnames(dat)[11:22] |>
  sub(pattern = "AGAL_", replacement = "") |>
  sub(pattern = "control_", replacement = "WT_") |>
  sub(pattern = "knockout_", replacement = "KO_")

## Build report
wb <- createWorkbook()
addWorksheet(wb, 'RNA_Analysis')

col_border <- createStyle(
  border = "right", numFmt = "@"
)

bottom_border <- createStyle(
  border = "bottom"
)

right_border <- createStyle(
  border = "right"
)

writeData(wb, 1, dat, startRow = 3, startCol = 1, colNames = FALSE)

# write headers in row 1
writeData(wb, 1, "Normalized Counts", 
          startRow = 1, startCol = 11)
writeData(wb, 1, "ENSEMBL ID", 
          startRow = 2, startCol = 1)
writeData(wb, 1, "Gene Name", 
          startRow = 2, startCol = 2)

for (i in seq_along(11:22)) {
  writeData(wb, 1, samp_names[i], 
            startRow = 2, startCol = (11:22)[i])
}

col_nums <- seq(3, 10, 2)
for (i in seq_along(comp_names)) {
  writeData(wb, 1, comp_names[i], 
            startRow = 1, startCol = col_nums[i])
  writeData(wb, 1, "log2FC", 
            startRow = 2, startCol = col_nums[i])
  writeData(wb, 1, "padj", 
            startRow = 2, startCol = col_nums[i]+1)
  addStyle(wb, 1, col_border, rows = 1:(nrow(dat)+2), cols = col_nums[i]+1)
}


# Create header style
header_style <- createStyle(
  fontSize = 12,
  fontColour = "#000000",
  halign = "center",
  textDecoration = "bold",
  border = "bottomright"
)

mergeCells(wb, sheet=1, cols=3:4, rows=1)
mergeCells(wb, sheet=1, cols=5:6, rows=1)
mergeCells(wb, sheet=1, cols=7:8, rows=1)
mergeCells(wb, sheet=1, cols=9:10, rows=1)
mergeCells(wb, sheet=1, cols=11:22, rows=1)

addStyle(wb, 1, col_border, rows = 1:(nrow(dat)+2), cols = 2)
addStyle(wb, 1, header_style, rows = 1, cols = 1:22)
addStyle(wb, 1, bottom_border, rows = 1, cols = 1)
addStyle(wb, 1, right_border, rows = 2:(nrow(dat)+2), cols = 22)

saveWorkbook(wb, file.path(projdir, "results/RNA_res.xlsx"),
             overwrite = TRUE)
