
color_groups <- c( "CNS"= "#1DB100",
                   "Contractile"="#F8BA00",
                   "Reproductive"="#EE220C",
                   "Digestive"="#265A8C",
                   "Anabolic"="#EF5FA7",
                   "Secretory"="#A8A8A8",
                   "Other"="khaki",
                   "Immune"="#00A89D"
)


tissue_to_group_ms <- function(x) {
  mutate(x, Group = case_when(
    Tissue == "Cerebrum" ~ "CNS",
    Tissue == "Cerebellum" ~ "CNS",
    Tissue == "Brainstem" ~ "CNS",
    Tissue == "Spinal cord" ~ "CNS",
    Tissue == "Kidney" ~ "Anabolic",
    Tissue == "Liver" ~ "Anabolic",
    Tissue == "Heart" ~ "Contractile",
    Tissue == "Skeletal muscle" ~ "Contractile",
    Tissue == "Adipose" ~ "Secretory",
    Tissue == "Small intestine" ~ "Digestive",
    Tissue == "Large intestine" ~ "Digestive",
    Tissue == "Stomach" ~ "Digestive",
    Tissue == "Placenta" ~ "Reproductive",
    Tissue == "Testis" ~ "Reproductive"
  ), .after = "Tissue")
}

## Read the data
data <- read_csv(here::here("MitoData/MitoCartaMouse_processed.csv"))
data <- data %>%
  pivot_longer(cols = colnames(data[2:ncol(data)])) %>%
  pivot_wider(names_from = "Gene", values_from = "value") %>%
  column_to_rownames("name")

## Heatmap data preparation
exprs <- log10(data)
col_fun = colorRamp2(c(4, 8,  12),c("#fff899", "#ffa76b", "#ff4242"))

### Clustering Columns
distance    = dist(t(as.matrix(exprs)), method="euclidean" ) 
coldend     = hclust(distance, method="ward.D2")  
### Clustering Rows
row_dist    = dist(as.matrix(exprs), method="euclidean" ) 
rowdend     = hclust(row_dist, method="ward.D2")  

## Column annotation
row_anno_df <- data %>%
  as.data.frame() %>%
  rownames_to_column("Tissue") %>%
  tissue_to_group_ms() %>%
  select(Tissue, Group)

row_anno = rowAnnotation(
  `Tissue group`=row_anno_df$Group,
  col=list(`Tissue group`  = color_groups),
  show_annotation_name = F,
  show_legend =  F,
  simple_anno_size = unit(0.1, "in"),
  annotation_legend_param = list(nrow=6))

## Heatmap
HM <- Heatmap(exprs, 
              name = "Log(protein exprs)", 
              col=col_fun,
              column_title="977 mitochondrial proteins",
              show_column_dend=F,
              show_column_names = F,
              cluster_columns =coldend,
              row_title="",
              row_dend_side = "right",
              row_names_side = "left",
              show_row_dend=TRUE,
              show_row_names=T,
              cluster_rows =rowdend,
              row_names_gp = grid::gpar(fontsize = 6),
              column_title_gp = grid::gpar(fontsize = 6),
              row_dend_width = unit(0.15, "in"),
              right_annotation=row_anno,
              width = unit(0.0013, "in")*977,
              heatmap_legend_param = list(
                at = c(4, 8, 12),
                title = "Log10 exprs",
                grid_height = unit(0.6, "in"),
                grid_width = unit(0.1, "in"),
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 6)
              )
)

draw(HM, heatmap_legend_side = "right", annotation_legend_side = "right")


rm(list = setdiff(ls(), c("processed_data", "color_groups", "color_tissues", "tissue_to_group_ms")))

# Figure 2B PCA Mouse -----------------------------------------------------

data <- read_csv(here::here("MitoData/MitoCartaMouse_processed.csv")) %>%
  pivot_longer(cols = -Gene, names_to = "Tissue") %>%
  tissue_to_group_ms() %>%
  mutate(value = log10(value)) %>%
  pivot_wider(names_from = "Gene", values_from = "value")
pca <- prcomp(data[,-c(1,2)], scale. = T)
summary(pca)
p <- autoplot(pca, data = data, colour = 'Group', size = 2.8, alpha = 0.7)+#,loadings = T, loadings.label = TRUE, loadings.label.size  = 1) +
  theme_bw() +
  scale_color_manual(values = c(c("CNS"= "#1DB100",
                                  "Contractile"="#F8BA00",
                                  "Reproductive"="#EE220C",
                                  "Digestive"="#265A8C",
                                  "Anabolic"="#EF5FA7",
                                  "Secretory"="#A8A8A8",
                                  "Other"="khaki",
                                  "Immune"="#00A89D"))) +
  theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

print(p)



