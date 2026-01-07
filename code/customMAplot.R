# create an MA plot for DESeq results
# color the points by significance over some adjusted pvalue (thresh)
# only label these significant points according to thresh  cutoff
# mo = max overlaps to provide to geom_text_repel for gene labeling
customMAplot <- function(res, mo = 10, pval_thresh = 0.01, lfc_thresh = 1) {
  res_df <- as.data.frame(res)
  
  iv.sig <- res_df$padj < pval_thresh
  res_df$point_color <- (abs(res_df$log2FoldChange) > lfc_thresh) & iv.sig
  
  res_df$labels <- apply(res_df, 1, 
                         function(x) ifelse(x["point_color"], x["gene_name"], ""))
  
  ma_plt <- ggplot(res_df, aes(x = log2(baseMean + 1),
                               y = log2FoldChange,
                               label = labels)) +
    geom_point(aes(color = point_color), size = 0.5) +
    scale_colour_manual(values = c("black", "dodgerblue"),
                        name="padj < 0.01 & \n|log2FC| > 1") +
    geom_text_repel(col = "dodgerblue", max.overlaps = mo) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "grey", linetype="dashed")
  return(ma_plt)
}
