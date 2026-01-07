library(gprofiler2)

getSigRes <- function(res, padj = 0.01, lfc = 1) {
  sig_res <- res[!is.na(res$padj), ]
  sig_res <- sig_res[(sig_res$padj < padj), ]
  sig_res <- sig_res[(abs(sig_res$log2FoldChange) > lfc), ] |>
    as.data.frame()
  return(sig_res)
}
