
# Get Data --------------------------------------------------------------------

data_pathways <- read_csv(here::here("MitoData", "Fibroblasts_Mitopathways.csv")) %>%
  select(RNAseq_sampleID, Condition, Pathway, score)
all_pathways <- unique(data_pathways$Pathway)
meta <- read_csv(here::here("MitoData", "Fibroblasts_Mitopathways.csv")) %>%  
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Treatments, Clinical_condition, Condition) %>%
  unique()

# Create pathway pairs ----------------------------------------------------

pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 
colnames(pw_pairs) <- c("pw1", "pw2")

data_long  <- data_pathways

janitor::tabyl(data_long, Pathway)
data_added1 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw1) %>%
  dplyr::rename(Pathway = pw1) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value1 = score, Pathway1 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(Pathway1) %>%
  group_by(RNAseq_sampleID, Condition) %>%
  nest() %>%
  dplyr::rename(pw1 = data) 

data_added2 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw2) %>%
  dplyr::rename(Pathway = pw2) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value2 = score, Pathway2 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(desc(Pathway2)) %>%
  unique() %>%
  group_by(RNAseq_sampleID, Condition) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 


# Calculate ratios --------------------------------------------------------

data_ratios =full_join(data_added1, data_added2, by = c("RNAseq_sampleID", "Condition")) %>%
  dplyr::mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  dplyr::select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  dplyr::mutate(ratio = value1 / value2)
sapply(data_ratios, class)
ctrl_ratios = data_ratios %>%
  group_by(Pathway1, Pathway2) %>%
  mutate(ctrl_ratio_av = mean(ratio)) %>%
  select(Pathway1, Pathway2, ctrl_ratio_av) %>%
  unique() 

data_ratios_corrected <- data_ratios %>%
  full_join(ctrl_ratios) %>%
  select(-c(value1, value2)) %>%
  dplyr::mutate(corrected_ratio = ratio / ctrl_ratio_av)


# Calculate mtPPS ---------------------------------------------------------

mtPPS <-data_ratios_corrected %>%
  group_by(RNAseq_sampleID, Pathway1) %>%
  mutate(mtPPS = mean(corrected_ratio)) %>%
  select(-c(Pathway2, corrected_ratio, ratio, ctrl_ratio_av)) %>%
  unique() %>%
  full_join(meta) %>%
  pivot_wider(names_from = "Pathway1", values_from = "mtPPS")  


# Figure 6C PCA -------------------------------------------------------------------

## hFB12 / HC1 -------------------------------------------------------------------

theme_pca <- theme(axis.text = element_text(size = 6),
                   axis.title = element_text(size = 7),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none"
)
data_pca <- mtPPS %>%
  filter(Line %in% "hFB12") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation

p <- autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca
print(p)

## hFB13 / HC2 -------------------------------------------------------------------

data_pca <- mtPPS %>%
  filter(Line %in% "hFB13") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation

p <- autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca
print(p)

## hFB14 / HC3 -------------------------------------------------------------------

data_pca <- mtPPS %>%
  filter(Line %in% "hFB14") %>%
  filter(Study_part %in% 2)
pca <- prcomp(data_pca[,-c(1:8)] , scale. = T)
loadings <- pca$rotation

p <- autoplot(pca, data = data_pca, colour = 'Passage', alpha = 0.6,
         fill = 'Passage', size = 2.3, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca
print(p)

# Figure 6D Heatmap -------------------------------------------------------
all_pathways <- colnames(mtPPS[,9:ncol(mtPPS)])
cor_res <- mtPPS %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway") %>%
  group_by(Pathway) %>%
  nest(-Pathway) %>%
  mutate(cor=map(data,~cor.test(.x$value, .x$Passage, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied, .drop = T) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

## Top and bottom 10 pathways
top_10 <- cor_res %>%
  dplyr::filter( (padj < 0.05) | (padj = 0.05) ) %>%
  ungroup() %>%
  arrange(desc(estimate)) %>%
  dplyr::slice(1:10) %>%
  pull(Pathway)

bottom_10 <- cor_res %>%
  dplyr::filter( (padj < 0.05) | (padj = 0.05) ) %>%
  ungroup() %>%
  arrange(estimate) %>%
  dplyr::slice(1:10) %>%
  pull(Pathway)

test <- mtPPS  %>%
  ungroup() %>%
  dplyr::select(Passage, Days_grown_Udays, Line, all_of(c(top_10, bottom_10))) %>%
  arrange(Line, Passage) %>%
  dplyr::rename(`Mt dynamics and surveillance` = `Mitochondrial dynamics and surveillance`,
                `NAD biosynthesis and metab.` = `NAD biosynthesis and metabolism`)

# Meta and expression data
meta  <- test[,1:3] #%>% column_to_rownames("Unique_variable_name")

exprs <- test[4:ncol(test)]
exprs = log10(exprs)
exprs = t(exprs)
rownames(exprs)
colnames(exprs)

# Col annotation
col_anno = columnAnnotation(
  `Cell_line`=meta$Line,
  Passage = anno_barplot(meta$Passage, axis_param=list(gp=gpar(fontsize = 4))),
  col=list(`Cell_line`  = c("hFB12" = "#636363", "hFB13" = "#969696", "hFB14" = "#D9D9D9")),
  show_annotation_name = F,
  show_legend =  T
)

col_anno <- re_size(col_anno, annotation_height = unit(c(0.1,0.3), "in"))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue","white", "#ff4242")) #"#ffa76b",
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")
HM2 <- Heatmap(exprs,
               name = "Pathway weight",
               cluster_rows = rowdend,
               show_column_names =F,
               cluster_columns =F,
               row_title="",
               row_dend_side = "right",
               row_names_side = "left",
               show_row_dend=F,
               show_row_names=T,
               top_annotation=col_anno,
               row_names_gp = grid::gpar(fontsize = 5),
               column_names_gp = grid::gpar(fontsize = 6),
               heatmap_legend_param = list(
                 title = "mtPPS (log10)",
                 title_position = "lefttop-rot",
                 legend_height = unit(0.9, "in"),
                 labels_gp = gpar(fontsize = 8),
                 title_gp = gpar(fontsize = 7, fontface= "bold")
               )
)

draw(HM2)

# Figure 6E Scatter plots ---------------------------------------------------------------

p <- mtPPS %>%
  ungroup() %>%
  dplyr::select(Line, Passage, `Fission`, `Mitophagy`, `mtDNA repair`, `mtDNA maintenance`) %>%
  pivot_longer(cols = -c(Line, Passage), names_to = "Pathway", values_to = "mtPPS") %>%
  ggplot(aes(x = Passage, y = mtPPS, color = Pathway)) +
  geom_point(size =0.8, alpha = 0.7) +
  scale_color_manual(values = c("Fission" ="#CE9332", "Mitophagy" = "#55BC82",
                                "mtDNA repair" ="#56BCC2", "mtDNA maintenance" = "#4FADF0")) +
  theme_bw() +
  facet_wrap(~Pathway, scales = "free", ncol = 4) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position =  "none",
        strip.text = element_text(size = 6))
print(p)
cor_res %>%
  filter(Pathway %in% c("Fission", "Mitophagy", "mtDNA repair", "mtDNA maintenance"))

