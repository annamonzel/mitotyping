color_groups <- c( "Brain"= "#1DB100",
                   "Contractile"="#F8BA00",
                   "Reproductive"="#EE220C",
                   "Digestive"="#265A8C",
                   "Anabolic"="#EF5FA7",
                   "Secretory"="#A8A8A8",
                   "Other"="khaki",
                   "Immune"="#00A89D"
)

color_tissue <- c(
  "Adipose tissue"= "#A8A8A8","Adrenal gland"= "#A8A8A8", "Amygdala"= "#1DB100",             
  "Appendix"="#265A8C", "Basal ganglia"="#1DB100","Bone marrow"="#00A89D", 
  "Breast"="khaki", "Cerebellum"= "#1DB100", "Cerebral cortex"= "#1DB100",        
  "Cervix"="#EE220C", "Choroid plexus"= "#1DB100", "Colon"="#265A8C", "Duodenum" ="#265A8C",           
  "Endometrium"="#EE220C", "Epididymis"="#EE220C", "Esophagus"="#265A8C", "Fallopian tube"="#EE220C",       
  "Gallbladder" ="#265A8C", "Heart muscle"="#F8BA00", "Hippocampal formation"= "#1DB100", 
  "Hypothalamus" = "#1DB100", "Kidney"="#EF5FA7", "Liver"="#EF5FA7","Lung"="khaki",                
  "Lymph node"="#00A89D", "Medulla oblongata"= "#1DB100", "Midbrain"= "#1DB100",
  "Olfactory bulb"= "#1DB100", "Ovary" ="#EE220C","Pancreas"="khaki", 
  "Parathyroid gland"="#A8A8A8", "Pituitary gland"="#A8A8A8","Placenta"="#EE220C",              
  "Pons"= "#1DB100","Prostate"="#EE220C","Rectum"= "#265A8C","Retina"= "khaki",             
  "Salivary gland"="#A8A8A8", "Seminal vesicle" ="#EE220C", "Skeletal muscle"="#F8BA00",       
  "Skin" ="khaki","Small intestine"="#265A8C","Smooth muscle" ="khaki","Spinal cord"= "#1DB100",            
  "Spleen"="#00A89D","Stomach"="#265A8C", "Testis"="#EE220C","Thalamus"= "#1DB100",         
  "Thymus"="#00A89D", "Thyroid gland"="#A8A8A8","Tongue"="khaki","Tonsil"="#00A89D",                
  "Urinary bladder"="khaki", "Vagina"="#EE220C", "White matter"= "#1DB100"         
)
theme_pca <-   theme(
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.text = element_blank(),
  legend.title = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())

## Data normalization ------------------------------------------------------

data <- read_csv(here::here("MitoData", "HPA_Mitogenes.csv"))
data_norm <- data %>%
  pivot_longer(cols= colnames(data[3:ncol(data)]), names_to = "Gene", values_to = "nTPM")  %>%
  group_by(Gene) %>%
  mutate(mean = mean(nTPM, na.rm= T)) %>%
  mutate(sd = sd(nTPM, na.rm= T)) %>%
  mutate(norm = (nTPM-mean)/sd) %>%
  select(-c(mean, sd, nTPM, Group)) %>%
  unique() %>%
  pivot_wider(names_from = "Tissue", values_from = "norm")



# Figure2D Heatmap --------------------------------------------------------

## Prepare the heatmap data
meta  <- data %>% select(Tissue, Group)
exprs <- data_norm[,2:ncol(data_norm)]
rownames(exprs) <- data_norm$Gene
range(exprs)
exprs = t(exprs)
## Column annotation
row_anno = rowAnnotation(
  `Tissue Group`=meta$Group,
  col=list(`Tissue Group`  = color_groups),
  show_annotation_name = F,
  show_legend =  T,
  simple_anno_size = unit(0.1, "in"),
  annotation_legend_param = list(nrow=8, labels_gp = gpar(fontsize = 6),
                                 title_gp = gpar(fontsize = 7, fontface = "bold"),
                                 grid_height = unit(0.02, "in"),
                                 grid_width = unit(0.1, "in")
  )
)

col_fun = colorRamp2(c(-3, 0, 3), c("#fff899", "#ffa76b", "#ff4242"))
col_fun(seq(-3, 3))

distance    = dist(t(as.matrix(exprs)), method="euclidean" ) 
coldend     = hclust(distance, method="ward.D2")  
row_dist    = dist(as.matrix(exprs), method="euclidean" )
rowdend     = hclust(row_dist, method="ward.D2")

HM2 <- Heatmap(exprs, 
               name = "Scaled mt transcripts", 
               col=col_fun,
               column_title="1134 mitochondrial genes",
               #row_title_side="left",
               #row_title_rot=90,
               show_column_dend=F,
               show_column_names = F,
               cluster_columns =coldend,
               row_title="",
               row_dend_side = "right",
               row_names_side = "left",
               show_row_dend=TRUE,
               show_row_names=T,
               cluster_rows=rowdend,
               right_annotation=row_anno,
               row_names_gp = grid::gpar(fontsize = 6),
               column_title_gp = grid::gpar(fontsize = 6),
               row_dend_width = unit(0.15, "in"),
               width = unit(0.95, "in"),
               heatmap_legend_param = list(
                 at = c(-3,  3),
                 labels = c("low",  "high"),
                 title = "Rel. expression",
                 title_position = "lefttop-rot",
                 grid_height = unit(0.6, "in"),
                 grid_width = unit(0.1, "in"),
                 title_gp = gpar(fontsize = 6, face = "bold"),
                 labels_gp = gpar(fontsize = 6),
                 legend_height = unit(0.57, "in")
               )
)

draw(HM2, heatmap_legend_side = "left", annotation_legend_side = "bottom")#, height = unit(0.05, "mm")*1135)

rm(list = setdiff(ls(), c("data", "tissue_to_group_hm", "color_groups", "theme_pca", "color_tissue")))



# Figure2E PCA ------------------------------------------------------------

## PCA 1 -------------------------------------------------------------------

pca <- prcomp(data[,-c(1,2)], scale. = T)
summary(pca)
p <- p <- autoplot(pca, data = data, colour = 'Group',
                   fill = 'Group', size = 2.8, alpha = 0.6, shape = 21)+
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = color_groups) +
  theme_bw() +
  theme_pca 
print(p)

## PCA 2 -------------------------------------------------------------------

data_sub <- data %>%
  filter(!Tissue %in% c("Tongue",  "Heart muscle" ,"Skeletal muscle","Cerebellum",
                        "White matter" ,"Cerebral cortex","Choroid plexus" ,
                        "Thalamus", "Hypothalamus", "Medulla oblongata", 
                        "Basal ganglia", 
                        "Pons", "Spinal cord", "Midbrain", "Amygdala",
                        "Hippocampal formation", 
                        "Olfactory bulb", "Liver", "Kidney", 
                        "Adrenal gland", "Parathyroid gland"))
pca <- prcomp(data_sub[,-c(1,2)], scale. = T)
summary(pca)
p <-  p <- autoplot(pca, data = data_sub, colour = 'Tissue',
                    fill = 'Group', size = 2.8, alpha = 0.6, shape = 21)+
  scale_color_manual(values = color_tissue) +
  scale_fill_manual(values = color_groups) +
  theme_pca +
  theme(panel.background=element_rect(color = "black", fill = "white", linetype = "dotted"))

print(p)

