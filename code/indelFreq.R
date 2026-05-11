library(ggplot2)

datadir <- "/work/users/g/k/gkenney/AGAL/data"
outdir <- "/work/users/g/k/gkenney/AGAL/results"

df <- read.csv(file.path(datadir, "Fig_1A_data.csv"))
colnames(df) <- c("indel_size", "indel_freq")

pdf(file = file.path(outdir, "figures/indelFreq.pdf"))
ggplot(df, aes(x = indel_size, y = indel_freq)) +
  geom_col(fill = "dodgerblue2") +
  scale_x_continuous(breaks = seq(-12, 2, by = 2)) +
  theme_classic() +
  xlab("Indel Size (bp)") +
  ylab("Indel Frequency (%)")
dev.off()
  
