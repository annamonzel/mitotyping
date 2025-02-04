rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)
library(ggfortify)
# Setup -------------------------------------------------------------------

color_conditions = c("Control" ="#929292", 
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")

# Figure 6H Heatmap -------------------------------------------------------

mtPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mtPPS_All.csv"))
data_heatmap <- mtPPS %>%
  column_to_rownames("RNAseq_sampleID") %>%
  arrange(fct_relevel(Condition, c("Control",
                                   
                                   "Hypoxia","Contact_Inhibition", "Contact_Inhibition_Hypoxia",
                                   
                                   "2DG", "Galactose", "BHB",
                                   "SURF1","Oligomycin","MitoNUITs", 
                                   "DEX", "MitoNUITs_DEX", "Oligomycin_DEX", "SURF1_DEX")), Passage)

exprs <- log10(data_heatmap[,8:ncol(data_heatmap)])
anno_df <- data_heatmap[,1:7]

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(exprs, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


col_anno <- columnAnnotation(Passage = anno_barplot(anno_df$Passage),
                             Condition = anno_df$Condition,
                             col=list(#Treatment  = color_conditions,
                               Condition = color_conditions
                             ))

set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              row_names_gp = gpar(fontsize = 6), row_km = 4)

hm <- draw(hm)


row_order <- data.table::melt(row_order(hm)) %>%
  dplyr::rename(Position = value, Cluster = L1) %>%
  mutate(row_order_in_heatmap = seq(1:149))
clusters <- data.frame(Pathway = colnames(exprs), Position = seq(1:149)) %>%
  inner_join(row_order) %>%
  arrange()

hm

pdf(here::here("Figures","Figure6","Fig6G_Heatmap_mtPPS_noLabels_small.pdf"), height = 6, width = 10, "in")
set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              show_row_names = F, row_km = 4)
draw(hm)
dev.off()



# Suppl FigS8 ------------------------------------------------------------


## B Heatmap ---------------------------------------------------------------

pdf(here::here("Figures","Figure6","Suppl_FigS8B_Heatmap_mtPPS_withLabels_large.pdf"), height = 12, width = 10, "in")
set.seed(123)
hm <- ComplexHeatmap::Heatmap(t(exprs),
                              column_title = "",
                              column_title_side = c("top"),
                              cluster_columns = F,
                              cluster_rows = T,
                              top_annotation = col_anno,
                              show_column_names = F,
                              show_row_names = T, row_km = 4,
                              row_names_gp = gpar(fontsize = 6))

draw(hm)
dev.off()

## C PCA - combine datasets ---------------------------------------------------

mitocarta_sheet4 <- readxl::read_xls(here::here("Data","MitoCarta", "OriginalData","HumanMitoCarta3_0.xls"), sheet = 4) %>%
  dplyr::select(MitoPathway, Genes) %>%
  na.omit()
gene_to_pathway <- splitstackshape::cSplit(mitocarta_sheet4, 'Genes', ',') %>%
  column_to_rownames("MitoPathway") %>%
  t() %>%
  as.data.frame()
gene_to_pathway <- gene_to_pathway %>%
  pivot_longer(cols = colnames(gene_to_pathway), names_to = "Pathway", values_to = "Gene") %>%
  na.omit()

all_pathways <- unique(gene_to_pathway$Pathway)
mtPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mtPPS_All.csv"))
tissue_mtPPS <- read_csv(here::here("Data", "HumanProteinAtlas", "ProcessedData", "mtPPS_all_tissues.csv"))

combined <- mtPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  #  inner_join(meta) %>%
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Tissue, Group, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, Tissue, Group), names_to = "Pathway") %>%
  full_join(tissue_mtPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway")%>%
              mutate(Passage =NA) ) %>%
  pivot_wider(names_from = "Pathway", values_from = "value") %>%
  arrange(desc(Group))

theme_pca <- theme(axis.text = element_text(size = 6),
                   axis.title = element_text(size = 7),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none"
)


pca <- prcomp(combined[, -c(1:7)], scale = T)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Suppl_Fig8C_loadings.csv"))

summary(pca)
autoplot(pca, data = combined , colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 1, shape = 21, lwd = 0.5, )+
  viridis::scale_color_viridis(direction = -1, na.value = "grey") +
  viridis::scale_fill_viridis(direction = -1, na.value = "grey") +
  theme_minimal() +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_pca## Note the switch in sign in PCA, axis reversed for visualization purposes. 

ggsave(here::here("Figures", "Figure6", "Suppl_FigS8C_PCA_Tissue_Fibroblasts_AllTreatments_AllPassages.png"), height = 1.35, width = 1.5,unit = "in", dpi = 1200)


pca$x <- pca$x[which(combined$Group == "Fibroblasts"), ]
autoplot(pca, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 2, shape = 21, lwd = 0.5)+
  viridis::scale_color_viridis(direction = -1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_minimal() +
  theme_pca + 
  scale_x_reverse()+
  scale_y_reverse() 
ggsave(here::here("Figures", "Figure6", "Suppl_FigS8C_PCA_Tissue_Fibroblasts_AllTreatments_AllPassages_subset.png"), 
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)

color_conditions = c("Control" ="#929292", 
                     "Tissue" = "gray",
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")



## D PCA All passages averaged ------------------------------------------------

combined <- mtPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
  select(RNAseq_sampleID, Line, Study_part, Passage, Days_grown_Udays, 
         Tissue, Group,Condition, Treatments, Days_grown_Udays, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line, Study_part, 
                         Passage, Days_grown_Udays, Tissue, 
                         Group, Condition, Treatments, Days_grown_Udays), names_to = "Pathway") %>%
  unique() %>%
  group_by(Line, Condition, Treatments, Pathway ) %>%
  mutate(mean = mean(value)) %>% 
  select(-c(value, Passage, Days_grown_Udays, RNAseq_sampleID, Study_part)) %>%
  unique() %>%
  full_join(tissue_mtPPS %>%pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mean")%>%
              mutate(Passage =NA) %>%
              mutate(Condition = "Tissue")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mean") %>%
  arrange(desc(Group))

pca <- prcomp(combined[, -c(1:7)], scale = T)
summary(pca)
loadings <- pca$rotation
write_csv(as.data.frame(loadings) %>% rownames_to_column("Pathway"), 
          here::here("Supplemental/Suppl_Fig8D_loadings.csv"))

autoplot(pca, data = combined , colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 1, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) +
  scale_fill_manual(values = color_conditions  ) +
  theme_minimal() +
  theme_pca +
  scale_y_reverse()
ggsave(here::here("Figures", "Figure6", "Suppl_FigS8D_PCA_Tissue_Fibroblasts_AllTreatments_AllPassagesAverage.png"),
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)


pca$x <- pca$x[which(combined$Group == "Fibroblasts"), ]
autoplot(pca, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 
ggsave(here::here("Figures", "Figure6", "Suppl_Fig8D_PCA_Tissue_Fibroblasts_AllTreatments_AllPassagesAveraged_subset.png"), 
       height = 1.35, width = 1.5,unit = "in", dpi = 1200)


## PC3 vs 4
autoplot(pca, x = 3, y = 4, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 

## PC3 vs 4
autoplot(pca, x = 5, y = 6, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 


## E 3D plot -----------------------------------------------------------------
#pca <- prcomp(combined[, -c(1:7)], scale = T)
summary_pca <- summary(pca)
color_conditions = c("Control" ="#929292", 
                     #"Tissue" = "gray",
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5")
var_1 <- round(summary_pca$importance[2,1]*100,2)
var_2 <- round(summary_pca$importance[2,2]*100,2)
var_3 <- round(summary_pca$importance[2,3]*100,2)
group_color_df <- data.frame(Color = color_conditions) %>%
  rownames_to_column("Condition")
df <- pca$x
conditions <- combined %>% filter(!Condition %in% "Tissue") %>%
  pull(Condition)
df <- data.frame(PC1=df[,1], PC2=df[,2], PC3=df[,3], Condition = as.factor(conditions)) %>%
  full_join(group_color_df, by = "Condition") %>%
  na.omit() 

rgl::par3d(cex=0.001)
rgl::bg3d(color = "white")
with(df, rgl::plot3d(PC1,PC2,PC3, col= Color, alpha = 0.6, size = 15, type = "p", 
                     xlab = "", ylab = "", zlab = "",
                     ann = FALSE, axes = FALSE
))
rgl::box3d()


# Suppl FigS9 Radar charts ---------------------------------------------------

## Fibroblast data
fibs <- mtPPS %>%
  filter(!Condition %in% "Control") %>%
  select(Condition, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(Condition), names_to = "Pathway") %>%
  unique() %>%
  group_by(Condition, Pathway ) %>%
  mutate(mean = mean(value)) %>% 
  select(-value) %>%
  unique() 
## Axis label with percent changes
label_fibs <- fibs %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mean)) %>%
  mutate(pct_change = round((mean - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Condition, label, sep = " ", remove= FALSE) %>%
  select(Condition, label, Pathway)

## Tissue data
tissues <- tissue_mtPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mtPPS")%>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 
## Axis label with percent changes
label_tissues <- tissue_mtPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mtPPS)) %>%
  mutate(pct_change = round((mtPPS - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Tissue, label, sep = " ", remove= FALSE) %>%
  select(Tissue, label, Pathway) %>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 


## Complex I ---------------------------------------------------------------
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "Complex I") 
fibs_sub <- fibs %>%
  filter(Pathway %in% "Complex I") #%>%
fibs_average <- mean(fibs_sub$mean)
tissues_average <- tissue_mtPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Pathway) %>%
  mutate(Grand_mean = mean(mtPPS)) %>%
  filter(Pathway %in% "Complex I") %>%
  pull(Grand_mean) %>%
  unique()
min <- min(tissues_sub$mtPPS, fibs_sub$mean)
min <- min + min/100
max <- max(tissues_sub$mtPPS, fibs_sub$mean)
max <- max + max/100
data <- fibs_sub %>%
  pivot_wider(names_from = "Condition", values_from = "mean") %>%
  ungroup() %>%
  select(-Pathway) %>%
  select(Hypoxia, SURF1_DEX, Oligomycin_DEX,
         MitoNUITs_DEX, DEX, Oligomycin, MitoNUITs, 
         SURF1, BHB, Galactose, `2DG`, Contact_Inhibition_Hypoxia,
         Contact_Inhibition, Hypoxia) %>%
  as.data.frame()
new_colnames <- as.data.frame(colnames(data)) %>%
  dplyr::rename(Condition = `colnames(data)`) %>%
  full_join(label_fibs %>%
              filter(Pathway %in% "Complex I"))
colnames(data) <- new_colnames$label

data <- rbind(rep(max,13) , rep(min,13), rep(fibs_average, 13),data) 
rownames(data) <- c("Max", "Min", "Average", "Fibroblasts")
colors_border=c( 
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(73/255, 116/255, 159/255,0.9) 
  )
colors_in=c( 
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(73/255, 116/255, 159/255,0.4) 
) 
png("Figures/Figure6/Suppl_FigS9A_radar_Fibs_CI.png", width = 6.4, height = 6, units = "in", res = 1200)
radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Fibroblasts", calcex=0.8, seg = 3
) 
dev.off()

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mtPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure6/Suppl_FigS9A_radar_Tissues_CI.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Tissues", calcex=0.8, seg = 3
) 
dev.off()

## FAO ---------------------------------------------------------------
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "Fatty acid oxidation") 
fibs_sub <- fibs %>%
  filter(Pathway %in% "Fatty acid oxidation") #%>%
fibs_average <- mean(fibs_sub$mean)
tissues_average <- tissue_mtPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Pathway) %>%
  mutate(Grand_mean = mean(mtPPS)) %>%
  filter(Pathway %in% "Fatty acid oxidation") %>%
  pull(Grand_mean) %>%
  unique()
min <- min(tissues_sub$mtPPS, fibs_sub$mean)
min <- min + min/100
max <- max(tissues_sub$mtPPS, fibs_sub$mean)
max <- max + max/100
data <- fibs_sub %>%
  pivot_wider(names_from = "Condition", values_from = "mean") %>%
  ungroup() %>%
  select(-Pathway) %>%
  select(Hypoxia, SURF1_DEX, Oligomycin_DEX,
         MitoNUITs_DEX, DEX, Oligomycin, MitoNUITs, 
         SURF1, BHB, Galactose, `2DG`, Contact_Inhibition_Hypoxia,
         Contact_Inhibition, Hypoxia) %>%
  as.data.frame()
new_colnames <- as.data.frame(colnames(data)) %>%
  dplyr::rename(Condition = `colnames(data)`) %>%
  full_join(label_fibs %>%
              filter(Pathway %in% "Fatty acid oxidation"))
colnames(data) <- new_colnames$label

data <- rbind(rep(max,13) , rep(min,13), rep(fibs_average, 13),data) 
rownames(data) <- c("Max", "Min", "Average", "Fibroblasts")
colors_border=c( 
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(73/255, 116/255, 159/255,0.9) 
)
colors_in=c( 
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(73/255, 116/255, 159/255,0.4) 
) 
png("Figures/Figure6/Suppl_FigS9B_radar_Fibs_FAO.png", width = 6.4, height = 6, units = "in", res = 1200)
radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Fibroblasts", calcex=0.8, seg = 3
) 
dev.off()

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mtPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure6/Suppl_FigS9B_radar_Tissues_FAO.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Tissues", calcex=0.8, seg = 3
) 
dev.off()


## mtDNA repair ---------------------------------------------------------------
## Fibroblast data only controls!!

fibs <- mtPPS %>%
  filter(Condition %in% "Control") %>%
  filter(Study_part %in% 2) %>%
  select(Line, Passage, all_of(all_pathways)) %>%
  pivot_longer(cols = -c(Line, Passage), names_to = "Pathway") %>%
  unique() %>%
  filter(Pathway %in% "mtDNA repair") 

## scatterplot
fibs %>% 
  ggplot(aes(x = Passage, y = value)) +
  geom_point(color ="darkmagenta") +
  facet_wrap(~Line) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 8))
ggsave("Figures/Figure6/Suppl_FigS9C_mtDNA_repair_scatter.png", width = 3.7, height = 2, units = "in", dpi = 1200)


label_tissues <- tissue_mtPPS %>%
  pivot_longer(cols = -c(Tissue, Group) , names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Pathway) %>%
  mutate(mean_all = mean(mtPPS)) %>%
  mutate(pct_change = round((mtPPS - mean_all) / mean_all * 100)) %>% 
  mutate(pct = "%") %>%
  mutate(sign = case_when(
    pct_change <0 ~ "",
    TRUE ~ "+"
  )) %>%
  unite(label, sign, pct_change, pct, sep = "") %>%
  unite(label, Tissue, label, sep = " ", remove= FALSE) %>%
  select(Tissue, label, Pathway) %>%
  filter(Tissue %in% c("Liver", "Cerebral cortex", "Skeletal muscle", "Bone marrow", "Small intestine")) 


### Tissue radar
tissues_sub <- tissues %>%
  full_join(label_tissues) %>%
  ungroup() %>%
  select(-c(Group, Tissue)) %>%
  dplyr::rename(Tissue = label) %>%
  filter(Pathway %in% "mtDNA repair") 

# Tissues
data <- tissues_sub %>%
  pivot_wider(names_from = "Tissue", values_from = "mtPPS") %>%
  ungroup() %>%
  select(-Pathway) %>%
  as.data.frame() 
data <- rbind(rep(max,5) , rep(min,5), rep(tissues_average, 5), data ) 
rownames(data) <- c("Max", "Min", "Average", "Tissues")
colors_border=c(
  rgb(210/255, 210/255, 210/255,0.9), 
  rgb(0,0,0,0.6) 
)
colors_in=c(
  rgb(210/255, 210/255, 210/255,0.4), 
  rgb(0,0,0,0.2) 
) 
png("Figures/Figure6/Suppl_FigS9C_radar_Tissues_mtDNA.png", width = 6.4, height = 6, units = "in", res = 1200)

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="mtDNA repair in Tissues", calcex=0.8, seg = 3
) 
dev.off()





# Suppl FigS10 Dynamic range -----------------------------------------------------------
mtPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mtPPS_All.csv")) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line,  Study_part, Passage,
                         Days_grown_Udays, Treatments, Clinical_condition, Condition), 
               names_to = "Pathway", values_to = "mtPPS") 

mtPPS_Tissues <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mtPPS_all_tissues.csv" )) %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "mtPPS")

delta_controls <- mtPPS %>%
  filter(Condition %in% "Control") %>%
  filter(Study_part %in% 2) %>%
  group_by(Pathway) %>%
  mutate(delta_controls = max(mtPPS) - min(mtPPS)) %>%
  select(Pathway, delta_controls) %>%
  unique()

delta_treatments <- mtPPS %>%
  filter(!Condition %in% "Control") %>%
  group_by(Condition, Pathway, Line) %>%
  mutate(mtPPS = mean(mtPPS)) %>%
  select(Condition, Pathway, Line, mtPPS) %>%
  unique() %>%
  group_by(Pathway) %>%
  mutate(delta_treatments = max(mtPPS) - min(mtPPS)) %>%
  select(Pathway, delta_treatments) %>%
  unique()

delta_tissues <- mtPPS_Tissues %>%
  group_by(Pathway) %>%
  mutate(delta_tissues = max(mtPPS) - min(mtPPS)) %>%
  select(Pathway, delta_tissues) %>%
  unique()

combined <- delta_controls %>%
  full_join(delta_treatments) %>%
  full_join(delta_tissues) %>%
  dplyr::rename(`Treated Fibroblasts` = delta_treatments, Tissues = delta_tissues, `Control Fibroblasts` =delta_controls ) %>%
  pivot_longer(cols = -Pathway, names_to = "dataset", values_to = "delta") %>%
  unite(Group, Pathway, dataset, remove = F) %>%
  mutate(dataset = as.factor(dataset))
combined$dataset <- factor(combined$dataset, levels = c("Tissues", "Control Fibroblasts", "Treated Fibroblasts"))
order <- combined %>%
  ungroup() %>%
  filter(dataset %in% "Tissues") %>%
  arrange(desc(delta)) %>%
  mutate(order = seq(1:149)) %>%
  select(Pathway, order)

data1 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Tissues")
data2 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Treated Fibroblasts")
data3 <- combined %>%
  inner_join(order) %>%
  arrange(desc(order)) %>%
  mutate(Pathway = as.factor(Pathway)) %>%
  arrange(desc(order)) %>% 
  filter(dataset %in% "Control Fibroblasts")

ggplot() +
  geom_bar(data = data1, aes(x = Pathway, y = delta), stat = "identity", fill = "grey", alpha =0.7) + 
  geom_bar(data = data2, aes(x = Pathway, y = delta), stat = "identity", fill = "dodgerblue4", alpha =0.7) + # must include argument label "data"
  geom_bar(data = data3, aes(x = Pathway, y = delta), stat = "identity", fill = "darkmagenta", alpha =0.7) + 
  theme_classic() +
  scale_x_discrete(limits=data1$Pathway) +
  coord_flip() +
  labs(x = "MitoPathways",y = "Dynamic range of mtPPS") + 
  theme(axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 7))
ggsave(here::here("Figures", "Figure6", "Fig6H_dynamic_range.png"), 
       height = 5, width = 1.81,unit = "in", dpi = 1200)

ggplot() +
  geom_bar(data = data1, aes(x = Pathway, y = delta), stat = "identity", fill = "grey", alpha =0.7) + 
  geom_bar(data = data2, aes(x = Pathway, y = delta), stat = "identity", fill = "dodgerblue4", alpha =0.7) + # must include argument label "data"
  geom_bar(data = data3, aes(x = Pathway, y = delta), stat = "identity", fill = "darkmagenta", alpha =0.7) + 
  theme_bw() +
  scale_x_discrete(limits=data1$Pathway) +
  coord_flip() +
  labs(x = "MitoPathways",y = "Dynamic range of mtPPS") + 
  theme(axis.text.y= element_text(size =6),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 7))
ggsave(here::here("Figures", "Figure6", "Suppl_Fig10A_dynamic_range_large.png"), 
       height = 10, width = 5,unit = "in", dpi = 1200)


## Bivariate mtPPS ---------------------------------------------------------
color_conditions = c("Control" ="#929292", 
                     "DEX"= "#E7298A", 
                     "Oligomycin"= "#1E78B4", 
                     "Oligomycin_DEX"= "#A5CEE3", 
                     "SURF1"= "#E3211D", 
                     "SURF1_DEX"= "#FC9A99", 
                     "MitoNUITs"="#33A12C", 
                     "MitoNUITs_DEX"= "#B2DF8A", 
                     "Galactose" = "#D96004", 
                     "2DG" = "#FFDA2E", 
                     "BHB"= "#B05928",
                     "Hypoxia"= "#F681C0",
                     "Contact_Inhibition" = "#6A3E9A",  
                     "Contact_Inhibition_Hypoxia"= "#CBB2D5",
                     "Tissue" = "lightgrey")

mtPPS <- read_csv(here::here("Data", "Fibroblasts", 
                                        "ProcessedData", "mtPPS_All.csv")) %>%
  pivot_longer(cols = -c(RNAseq_sampleID, Line,  Study_part, Passage,
                         Days_grown_Udays, Treatments, Clinical_condition, Condition), 
               names_to = "Pathway", values_to = "mtPPS") %>%
  select(RNAseq_sampleID, Condition, Pathway, Line, mtPPS) %>%
  unique() 

mtPPS_Tissues <- read_csv(here::here("Data", "HumanProteinAtlas","ProcessedData", "mtPPS_all_tissues.csv" )) %>%
  pivot_longer(cols = -c(Tissue, Group), names_to = "Pathway", values_to = "mtPPS") %>%
  mutate(Condition = "Tissue")

combined <- mtPPS %>%
  full_join(mtPPS_Tissues)

test <- combined %>%
  filter(Pathway %in% c("Phospholipid metabolism", "Catechol metabolism")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mtPPS") 

test %>%
  ggplot(aes(x = `Phospholipid metabolism`, y= `Catechol metabolism`, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "mtPPS") + 
  scale_color_manual(values = color_conditions)  +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 8),
        legend.position = "none")
ggsave(here::here("Figures", "Figure6", "Suppl_FigS10B_Tissues_vs_fibs_mtPPS_example1.png"), 
       height = 2.9, width = 2.61,unit = "in", dpi = 1200)

test <- combined %>%
  filter(Pathway %in% c("Cholesterol-associated", "Iron homeostasis")) %>%
  pivot_wider(names_from = "Pathway", values_from = "mtPPS") 
test %>%
  ggplot(aes(x = `Cholesterol-associated`, y= `Iron homeostasis`, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "mtPPS") + 
  scale_color_manual(values = color_conditions)  +
  theme_bw() +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 8),
        legend.position = "none")
ggsave(here::here("Figures", "Figure6", "Suppl_FigS10C_Tissues_vs_fibs_mtPPS_example2.png"), 
       height = 2.9, width = 2.61,unit = "in", dpi = 1200)
