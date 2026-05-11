library(ggplot2)
library(ggpubr)

datadir <- "/work/users/g/k/gkenney/AGAL/data"
outdir <- "/work/users/g/k/gkenney/AGAL/results"

df <- read.csv(file.path(datadir, "Fig_3B_data.csv"))
colnames(df) <- c("sample", "perc")

df$group <- lapply(df$sample, 
                   function(x) strsplit(x, split = " ")[[1]][1]) |>
  unlist()

pdf(file = file.path(outdir, "figures/fig3b.pdf"),
    height = 5, width = 3)
ggplot(df, aes(x = group, y = perc, col = group)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  xlab("") +
  ylab("CD123+ CD63+ cells (%)") +
  scale_color_manual(values = c("dodgerblue2", "#bc3754")) +
  theme(legend.position = "") +
  stat_compare_means(method = "t.test",
                     method.args = list(var.equal = FALSE), # Force Welch's T-test
                     label = "p.signif",                   # Shows stars (*, **, ns)
                     comparisons = list(c("Control", "KO"))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_x_discrete(labels = c("Control", "GGTA1 KO")) +
  ylim(c(0, 21))
dev.off()

