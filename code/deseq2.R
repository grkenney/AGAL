library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(org.Bt.eg.db)
library(clusterProfiler)


# ---- Setup Directories ---- 

# base directory of pipeline outputs and output directory for 
# figures and objects
basedir <- "/proj/phanstiel_lab/Data/processed/AGAL/rna"
outdir <- "/work/users/g/k/gkenney/AGAL/results"
codedir <- "/work/users/g/k/gkenney/AGAL/code"

source(file.path(codedir, "customMAplot.R"))
source(file.path(codedir, "getSigRes.R"))

# make an output directory for figures
figdir <- file.path(outdir, "figures")
if (!dir.exists(figdir)) {
  dir.create(figdir)
}


# ---- Load Data ----

# load gtf to convert gene ids to gene names
gtf <- rtracklayer::import('Bos_taurus.ARS-UCD1.3.113.gtf')

gns <- mcols(gtf)[, c("gene_id", "gene_name")] |> unique()

# if there isn't a gene name, just use the id
gns$gene_name <- apply(gns, 1, function(x) ifelse(is.na(x[2]), x[1], x[2]))
rownames(gns) <- gns$gene_id

# import tximport object from bagpipes
txi <- readRDS(file.path(basedir, "bagPipes/output/quant/AGAL_tximport.rds"))

# load the sample sheet
samplesheet <- read.delim(file.path(basedir, "bagPipes/RNAsamplesheet.txt")) |>
  as.data.frame()

# reorder factors for treatment and genotype
samplesheet$Treatment <- factor(samplesheet$Treatment, 
                                levels = c("undifferentiated", "differentiated"))
samplesheet$Genotype <- factor(samplesheet$Genotype, 
                               levels = c("control", "knockout"))
# make a new group combining both genotype and treatment
samplesheet$group <- factor(paste0(samplesheet$Genotype, "_", 
                                   samplesheet$Treatment))


# ---- DESeq2 ----

# make a DESeqDataSet object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samplesheet,
                                   design = ~ 0 + group)

# keep only rows that have 10 or more counts in at least 6 samples
paste("Rows before filter: ", nrow(ddsTxi))
keep <- rowSums(counts(ddsTxi) >= 10) >= 6
ddsTxi <- ddsTxi[keep,]
paste("Rows after filter: ", nrow(ddsTxi))

# run DESeq2
dds <- DESeq(ddsTxi)
# resultsNames(dds)

# output normalized counts
counts(dds)
write.csv(counts(dds, normalized=T), 
          file = file.path(outdir, "norm_counts.csv"))
saveRDS(dds, file.path(outdir, "dds.rda"))


## ---- undiff KOvC ----

### ---- MA Plot ----
res_undiff <- results(dds,
                        contrast = c("group",
                                     "knockout_undifferentiated",
                                     "control_undifferentiated"),
                      alpha = 0.01)
res_undiff$gene_name <- gns[rownames(res_undiff), "gene_name"]

customMAplot(res_undiff, pval_thresh = 0.01, lfc_thresh = 1) |> 
  ggsave(filename = file.path(outdir, "figures/MAplot_undiff_KOvC.pdf"),
         height = 7, width = 10)


resLFC_undiff <- lfcShrink(dds,
                           contrast = c("group",
                                        "knockout_undifferentiated",
                                        "control_undifferentiated"),
                           alpha = 0.01,
                           type = "ashr")
resLFC_undiff$gene_name <- gns[rownames(resLFC_undiff), "gene_name"]

customMAplot(resLFC_undiff, mo=Inf, pval_thresh = 0.01, lfc_thresh = 1) |> 
  ggsave(filename = file.path(outdir, "figures/MAplotShrink_undiff_KOvC.pdf"),
         height = 7, width = 10)


### ---- DE ----
res_undiff_sig <- getSigRes(res_undiff)
write.csv(res_undiff_sig, 
          file = file.path(outdir, "DE_undiff_KOvC.csv"),
          row.names = T)


## ---- diff KOvC ----

### ---- MA Plot ----
res_diff <- results(dds,
                    contrast = c("group",
                                 "knockout_differentiated",
                                 "control_differentiated"),
                    alpha = 0.01)
res_diff$gene_name <- gns[rownames(res_diff), "gene_name"]

customMAplot(res_diff, mo = Inf, pval_thresh = 0.01, lfc_thresh = 1) |> 
  ggsave(filename = file.path(outdir, "figures/MAplot_diff_KOvC.pdf"),
         height = 7, width = 10)

resLFC_diff <- lfcShrink(dds,
                         contrast = c("group",
                                      "knockout_differentiated",
                                      "control_differentiated"),
                         alpha = 0.01,
                         type = "ashr")
resLFC_diff$gene_name <- gns[rownames(resLFC_diff), "gene_name"]

customMAplot(resLFC_diff, mo = Inf, pval_thresh = 0.01, lfc_thresh = 1) |> 
  ggsave(filename = file.path(outdir, "figures/MAplotShrink_diff_KOvC.pdf"),
         height = 7, width = 10)


### ---- DE ----
res_diff_sig <- getSigRes(res_diff)
write.csv(res_diff_sig, 
          file = file.path(outdir, "DE_diff_KOvC.csv"),
          row.names = F)


## ---- control diff v undiff ----

res_cntrl <- results(dds,
                    contrast = c("group",
                                 "control_differentiated",
                                 "control_undifferentiated"),
                    alpha = 0.01)
res_cntrl$gene_name <- gns[rownames(res_cntrl), "gene_name"]

res_cntrl_sig <- getSigRes(res_cntrl)

customMAplot(res_cntrl, mo = 1, pval_thresh = 0.01, lfc_thresh = 1)

write.csv(as.data.frame(res_cntrl), 
          file = file.path(outdir, "DE_C_diff_v_undiff.csv"),
          row.names = T)


## ---- KO diff v undiff ----

res_ko <- results(dds,
                     contrast = c("group",
                                  "knockout_differentiated",
                                  "knockout_undifferentiated"),
                  alpha = 0.01)
res_ko$gene_name <- gns[rownames(res_ko), "gene_name"]

res_ko_sig <- getSigRes(res_ko)

write.csv(as.data.frame(res_ko), 
          file = file.path(outdir, "DE_KO_diff_v_undiff.csv"),
          row.names = T)


# ---- Results Table ----

all_res <- counts(dds, normalized=TRUE) |>
  as.data.frame()

# add gene id column
all_res$ensembl_id <- rownames(all_res)

# Add DEG results

# KO_undiff_vs_C_undiff
all_res <- merge(data.frame(ensembl_id = rownames(res_undiff),
                            log2FC_KO_undiff_vs_C_undiff = res_undiff$log2FoldChange,
                            padj_KO_undiff_vs_C_undiff = res_undiff$padj), 
                 all_res,
                 by = "ensembl_id")

# KO_diff_vs_C_diff
all_res <- merge(data.frame(ensembl_id = rownames(res_diff),
                            log2FC_KO_diff_vs_C_diff = res_diff$log2FoldChange,
                            padj_KO_diff_vs_C_diff = res_diff$padj), 
                 all_res,
                 by = "ensembl_id")

# C_diff_vs_C_undiff
all_res <- merge(data.frame(ensembl_id = rownames(res_cntrl),
                            log2FC_C_diff_vs_C_undiff = res_cntrl$log2FoldChange,
                            padj_C_diff_vs_C_undiff = res_cntrl$padj), 
                 all_res,
                 by = "ensembl_id")

# KO_diff_vs_KO_undiff
all_res <- merge(data.frame(ensembl_id = rownames(res_ko),
                            log2FC_KO_diff_vs_KO_undiff = res_ko$log2FoldChange,
                            padj_KO_diff_vs_KO_undiff = res_ko$padj), 
                 all_res,
                 by = "ensembl_id")

# add gene name column
all_res$gene_name <- gns[all_res$ensembl_id, "gene_name"]
# move ensembl_id and gene_name to first column
all_res <- dplyr::relocate(all_res, c(ensembl_id, gene_name))

# rename undiff to d0 and diff to d3
colnames(all_res) <- gsub(pattern = "undifferentiated", replacement = "d0", x = colnames(all_res))
colnames(all_res) <- gsub(pattern = "differentiated", replacement = "d3", x = colnames(all_res))
colnames(all_res) <- gsub(pattern = "undiff", replacement = "d0", x = colnames(all_res))
colnames(all_res) <- gsub(pattern = "diff", replacement = "d3", x = colnames(all_res))

write.csv(all_res, 
          file = file.path(outdir, "RNA_all_res.csv"),
          row.names = F)

