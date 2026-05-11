library(DESeq2)
library(ggplot2)

outdir <- "/work/users/g/k/gkenney/AGAL/results"

dds <- readRDS(file.path(outdir, "dds.rda"))

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Genotype", "Treatment"), returnData=TRUE)
# mirror flip PC1
pcaData$PC1f <- pcaData$PC1 * (-1)

percentVar <- round(100 * attr(pcaData, "percentVar"))
group_pca <- ggplot(pcaData, aes(PC1f, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#bc3754", "dodgerblue2"),
                     guide = guide_legend(override.aes = list(size = 4)),
                     labels = c("d0", "d3")) +
  scale_shape_manual(values = scales::shape_pal()(2), 
                     labels = c("Control", "Knockout"))
group_pca |> 
  ggsave(filename = file.path(outdir, "figures/PCA.pdf"),
         height = 5, width = 5)

