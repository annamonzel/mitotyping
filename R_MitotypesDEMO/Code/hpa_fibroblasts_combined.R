
# Get Fibroblast mtPPS --------------------------------------------------------

mtPPS <- read_csv(here::here("MitoData", "Fibroblasts_mtPPS.csv"))
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



## Figure 6G Heatmap -------------------------------------------------------

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

draw(hm)


# Calculate Tissue mtPPS ------------------------------------------------------------
processed_data <- read_csv(here::here("MitoData", "HPA_Mitogenes.csv"))
gene_to_pathway <- read_csv(here::here("MitoData", "gene_to_pathway_human.csv"))

data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique()%>%
  group_by(Tissue, Pathway) %>%
  mutate(score = mean(value, na.rm = T)) %>%
  select(-c(Gene, value)) %>%
  unique() 

## Create pathway pairs ----------------------------------------------------------

all_pathways <- unique(data_pathways$Pathway)
data_long <- data_pathways %>%
  select(Tissue, Group, Pathway, score) 

pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 
colnames(pw_pairs) <- c("pw1", "pw2")

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
  group_by(Tissue, Group) %>%
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
  group_by(Tissue, Group) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 

## Calculate pathway ratios ----------------------------------------------------------

data_ratios =full_join(data_added1, data_added2, by = c("Tissue", "Group")) %>%
  mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  mutate(ratio = value1 / value2) %>%
  group_by(Pathway1, Pathway2) %>%
  mutate(average_ratio = mean(ratio)) %>%
  mutate(corrected_ratio = ratio / average_ratio)

## Calculate mtPPS ----------------------------------------------------------

tissue_mtPPS <- data_ratios %>%
  filter(!is.na(corrected_ratio)) %>%
  group_by(Tissue, Pathway1) %>%
  mutate(mtPPS = mean(corrected_ratio)) %>%
  select(Tissue, Group, Pathway1, mtPPS) %>%
  unique() %>%
  pivot_wider(names_from = "Pathway1", values_from = "mtPPS")

# Combine both datasets ----------------------------------------------------

combined <- mtPPS %>%
  mutate(Group = "Fibroblasts") %>%
  mutate(Tissue = "Control") %>%
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


# Combined PCA ---------------------------------------------------------------------


## Colored by passage ------------------------------------------------------


pca <- prcomp(combined[, -c(1:7)], scale = T)
loadings <- pca$rotation

summary(pca)
p <- autoplot(pca, data = combined , colour = 'Passage', alpha = 0.6, scale = 1,
         fill = 'Passage', size = 1, shape = 21, lwd = 0.5, )+
  viridis::scale_color_viridis(direction = -1, na.value = "grey") +
  viridis::scale_fill_viridis(direction = -1, na.value = "grey") +
  theme_minimal() +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_pca## Note the switch in sign in PCA, axis reversed for visualization purposes 
print(p)

## Colored by treatment (average across lifespan) -------------------------------

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

p <- autoplot(pca, data = combined , colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 1, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) +
  scale_fill_manual(values = color_conditions  ) +
  theme_minimal() +
  theme_pca +
  scale_y_reverse()
print(p)

## Fibroblast subset -------------------------------
pca$x <- pca$x[which(combined$Group == "Fibroblasts"), ]
p <- autoplot(pca, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 
print(p)


## PC2 vs 3
p <- autoplot(pca, x = 3, y = 2, data = combined %>%
           filter(Group %in% "Fibroblasts"), colour = 'Condition', alpha = 0.6, scale = 1,
         fill = 'Condition', size = 2, shape = 21, lwd = 0.5, )+
  scale_color_manual(values = color_conditions  ) + 
  scale_fill_manual(values = color_conditions  ) + 
  theme_minimal() +
  theme_pca + 
  scale_y_reverse() 
print(p)

# Combined Radar charts ---------------------------------------------------

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


radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Fibroblasts", calcex=0.8, seg = 3
) 

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


radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Complex I in Tissues", calcex=0.8, seg = 3
) 


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

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Fibroblasts", calcex=0.8, seg = 3
) 


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


radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="Fatty acid oxidation in Tissues", calcex=0.8, seg = 3
) 



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
p <- fibs %>% 
  ggplot(aes(x = Passage, y = value)) +
  geom_point(color ="darkmagenta") +
  facet_wrap(~Line) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 8))
print(p)


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

radarchart(data  , axistype=0, pcol = colors_border, pfcol = colors_in,
           plwd=c(0.5,1.5),  pty=16,plty=c(1,3),
           cglcol="grey",cglwd = 1, cglty=3, axislabcol="grey",
           vlcex=0.8, title="mtDNA repair in Tissues", calcex=0.8, seg = 3
) 



# Suppl FigS10 Dynamic range -----------------------------------------------------------


delta_controls <- mtPPS %>%
  filter(Condition %in% "Control") %>%
  filter(Study_part %in% 2) %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Pathway) %>%
  mutate(delta_controls = max(mtPPS) - min(mtPPS)) %>%
  select(Pathway, delta_controls) %>%
  unique()

delta_treatments <- mtPPS %>%
  filter(!Condition %in% "Control") %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway", values_to = "mtPPS") %>%
  group_by(Condition, Pathway, Line) %>%
  mutate(mtPPS = mean(mtPPS)) %>%
  select(Condition, Pathway, Line, mtPPS) %>%
  unique() %>%
  group_by(Pathway) %>%
  mutate(delta_treatments = max(mtPPS) - min(mtPPS)) %>%
  select(Pathway, delta_treatments) %>%
  unique()

delta_tissues <- tissue_mtPPS %>%
  pivot_longer(cols = all_of(all_pathways), names_to = "Pathway", values_to = "mtPPS") %>%
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




p <- ggplot() +
  geom_bar(data = data1, aes(x = Pathway, y = delta), stat = "identity", fill = "grey", alpha =0.7) + 
  geom_bar(data = data2, aes(x = Pathway, y = delta), stat = "identity", fill = "dodgerblue4", alpha =0.7) + # must include argument label "data"
  geom_bar(data = data3, aes(x = Pathway, y = delta), stat = "identity", fill = "darkmagenta", alpha =0.7) + 
  theme_bw() +
  scale_x_discrete(limits=data1$Pathway) +
  coord_flip() +
  labs(x = "MitoPathways",y = "Dynamic range of mtPPS", title = "Tissues vs Fibroblasts - dynamic range") + 
  theme(axis.text.y= element_text(size =6),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 7))
print(p)


