library(dplyr)
library(gprofiler2)
library(viridis)

outdir <- "/work/users/g/k/gkenney/AGAL/results"


# ---- read data ----

dds <- readRDS(file.path(outdir, "dds.rda"))
res_cntrl_df <- read.csv(file.path(outdir, "DE_C_diff_v_undiff.csv"))
res_ko_df <- read.csv(file.path(outdir, "DE_KO_diff_v_undiff.csv"))

colnames(res_cntrl_df)[1] <- "ENSEMBL"
colnames(res_ko_df)[1] <- "ENSEMBL"

# set filter cutoffs
padj_cutoff <- 0.01
logfc_cutoff <- 1


# ---- get sig results ----

res_cntrl_sig <- res_cntrl_df %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::filter(padj < padj_cutoff)

res_ko_sig <- res_ko_df %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::filter(padj < padj_cutoff)


# ---- get up and down genes ----

up_c_gns <- res_cntrl_sig %>% 
  filter(log2FoldChange > logfc_cutoff) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select(ENSEMBL) %>% 
  unlist()

up_ko_gns <- res_ko_sig %>% 
  filter(log2FoldChange > logfc_cutoff) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select(ENSEMBL) %>% 
  unlist()

down_c_gns <- res_cntrl_sig %>% 
  filter(log2FoldChange < (-logfc_cutoff)) %>%
  arrange(desc(log2FoldChange)) %>% 
  select(ENSEMBL) %>% 
  unlist()

down_ko_gns <- res_ko_sig %>% 
  filter(log2FoldChange < (-logfc_cutoff)) %>%
  arrange(desc(log2FoldChange)) %>% 
  select(ENSEMBL) %>% 
  unlist()


# ---- go enrichment ----

# upregulated genes
go_c_up <- gost(query = up_c_gns, organism = "btaurus", ordered_query = FALSE,
                correction_method = "fdr", custom_bg = rownames(dds),
                sources = "GO:BP")
go_ko_up <- gost(query = up_ko_gns, organism = "btaurus", ordered_query = FALSE,
                 correction_method = "fdr", custom_bg = rownames(dds),
                 sources = "GO:BP")

# downregulated genes
go_c_down <- gost(query = down_c_gns, organism = "btaurus", ordered_query = FALSE,
                  correction_method = "fdr", custom_bg = rownames(dds),
                  sources = "GO:BP")
go_ko_down <- gost(query = down_ko_gns, organism = "btaurus", ordered_query = FALSE,
                   correction_method = "fdr", custom_bg = rownames(dds),
                   sources = "GO:BP")

# save results
go_out_dir <- file.path(outdir, "GO")
if(!file.exists(go_out_dir)) { dir.create(go_out_dir) }

apply(go_c_up$result, 2, as.character) |>
  write.csv(file = file.path(go_out_dir, "GO_cntrl_up.csv"),
            row.names = F)

apply(go_ko_up$result, 2, as.character) |>
  write.csv(file = file.path(go_out_dir, "GO_ko_up.csv"),
            row.names = F)

apply(go_c_down$result, 2, as.character) |>
  write.csv(file = file.path(go_out_dir, "GO_cntrl_down.csv"),
            row.names = F)

apply(go_ko_down$result, 2, as.character) |>
  write.csv(file = file.path(go_out_dir, "GO_ko_down.csv"),
            row.names = F)


# ---- select go terms ----
select_go_terms <- c("cell differentiation", # up terms
                     "anatomical structure development",
                     "muscle cell differentiation", 
                     "striated muscle cell differentiation",
                     "myotube differentiation",
                     "cell cycle", # down terms
                     "mitotic cell cycle",
                     "regulation of cell cycle",
                     "metaphase chromosome alignment",
                     "mitotic DNA replication")

# subset to selected go terms
go_c_up_select <- go_c_up$result %>% 
  filter(term_name %in% select_go_terms) %>% 
  mutate(genotype = "WT", direction = "up")
go_ko_up_select <- go_ko_up$result %>% 
  filter(term_name %in% select_go_terms) %>% 
  mutate(genotype = "KO", direction = "up")
go_c_down_select <- go_c_down$result %>% 
  filter(term_name %in% select_go_terms) %>% 
  mutate(genotype = "WT", direction = "down")
go_ko_down_select <- go_ko_down$result %>% 
  filter(term_name %in% select_go_terms) %>% 
  mutate(genotype = "KO", , direction = "down")

# combine into single data frame
go_select <- rbind(go_c_up_select, go_ko_up_select,
                   go_c_down_select, go_ko_down_select)


# ---- plot ----

# reorder fators for plotting
go_select$genotype <- factor(go_select$genotype, 
                             levels = c("WT", "KO"))
go_select$direction <- factor(go_select$direction, 
                              levels = c("up", "down"))
go_select$term_name <- factor(go_select$term_name, 
                              levels = select_go_terms)

# plot dotplot
pdf(file = file.path(outdir, "figures/go_dot_plot.pdf"),
    height = 5, width = 5)
ggplot(go_select, 
       aes(x = term_name, y = genotype,
           size = -log10(p_value), col = direction)) +
  geom_point() +
  coord_flip() +
  xlab("") + ylab("Genotype") +
  theme_classic() +
  scale_size(range = c(4, 14), breaks = c(5,10,15,20)) +
  scale_x_discrete(limits = rev) +
  labs(color = "Direction", size = "-log10(adj. pval)") +
  scale_color_manual(values = c("#bc3754", "dodgerblue2"),
                     guide = guide_legend(override.aes = list(size = 4)))
dev.off()
