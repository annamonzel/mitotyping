


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
  "Adipose tissue"= "#A8A8A8","Adrenal gland"= "#A8A8A8", "Amygdala"= "#A8A8A8",             
  "Appendix"="#A8A8A8", "Basal ganglia"="#A8A8A8","Bone marrow"="#A8A8A8", 
  "Breast"="#A8A8A8", "Cerebellum"= "#A8A8A8", "Cerebral cortex"= "#1DB100",        
  "Cervix"="#A8A8A8", "Choroid plexus"= "#A8A8A8", "Colon"="#A8A8A8", "Duodenum" ="#A8A8A8",           
  "Endometrium"="#A8A8A8", "Epididymis"="#A8A8A8", "Esophagus"="#A8A8A8", "Fallopian tube"="#A8A8A8",       
  "Gallbladder" ="#A8A8A8", "Heart muscle"="#A8A8A8", "Hippocampal formation"= "#A8A8A8", 
  "Hypothalamus" = "#A8A8A8", "Kidney"="#A8A8A8", "Liver"="#EF5FA7","Lung"="#A8A8A8",                
  "Lymph node"="#A8A8A8", "Medulla oblongata"= "#A8A8A8", "Midbrain"="#A8A8A8",
  "Olfactory bulb"= "#A8A8A8", "Ovary" ="#A8A8A8","Pancreas"="#A8A8A8", 
  "Parathyroid gland"="#A8A8A8", "Pituitary gland"="#A8A8A8","Placenta"="#A8A8A8",              
  "Pons"= "#A8A8A8","Prostate"="#A8A8A8","Rectum"= "#A8A8A8","Retina"="#A8A8A8",# "#1DB100",             
  "Salivary gland"="#A8A8A8", "Seminal vesicle" ="#A8A8A8", "Skeletal muscle"="#A8A8A8",       
  "Skin" ="#A8A8A8","Small intestine"="#A8A8A8","Smooth muscle" ="#A8A8A8","Spinal cord"= "#A8A8A8",            
  "Spleen"="#A8A8A8","Stomach"="#A8A8A8", "Testis"="#A8A8A8","Thalamus"= "#A8A8A8",         
  "Thymus"="#A8A8A8", "Thyroid gland"="#A8A8A8","Tongue"="#A8A8A8","Tonsil"="#A8A8A8",                
  "Urinary bladder"="#A8A8A8", "Vagina"="#A8A8A8", "White matter"= "#A8A8A8"       
)


processed_data <- read_csv(here::here("MitoData", "HPA_Mitogenes.csv"))
gene_to_pathway <- read_csv(here::here("MitoData", "gene_to_pathway_human.csv"))

## Combine data
data_pathways <- processed_data %>%
  pivot_longer(cols = colnames(processed_data[3:ncol(processed_data)]), names_to = "Gene") %>%
  full_join(gene_to_pathway, by = "Gene") %>%
  filter(!is.na(Pathway)) %>%
  unique() %>%
  group_by(Tissue, Pathway) %>%
  mutate(score= mean(value)) %>%
  ungroup() %>%
  select(Tissue, Group, Pathway, score) %>%
  unique() 


# Suppl FigS5C_Bivariate score CI FAO --------------------------------------------------

data_plot <- data_pathways %>%
  filter(Pathway %in% c("Complex I", "Fatty acid oxidation")) %>%
  pivot_wider(names_from = "Pathway", values_from = "score") %>%
  mutate(ratio = round(`Complex I` / `Fatty acid oxidation`), 1) %>%
  mutate(space = "             ") %>%
  unite(ratio, space, ratio, sep = "") %>%
  select(Tissue, Group, `Complex I`, `Fatty acid oxidation`, ratio) %>%
  unique() %>%
  mutate(average_ratio = mean(as.numeric(ratio))) %>%
  mutate(line_average_ratio_x = runif(55, 0,180) ) %>%
  mutate(line_average_ratio_y = line_average_ratio_x * 81)


p <- data_plot %>%
  ggplot(aes(x = `Fatty acid oxidation`, y = `Complex I`, color = Tissue)) +
  geom_point() +
  geom_line(aes(x = line_average_ratio_x, y= line_average_ratio_y), inherit.aes = F, linetype = "dotted") + 
  ylim(0, 10000) +
  scale_color_manual(values = color_tissue) +
  theme_classic() +
  theme(legend.position = "none",  
        plot.title = element_text(size =8),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(p)

p <- data_plot %>%
  ggplot(aes(x = `Fatty acid oxidation`, y = `Complex I`, color = Group, label = ratio)) +
  geom_point() +
  scale_color_manual(values = color_groups) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size =8),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(p)

## Calculate corrected ratios and mtPPS --------------------------------------------

all_pathways <- unique(data_pathways$Pathway)
data_long <- data_pathways %>%
  select(Tissue, Group, Pathway, score) %>%
  dplyr::rename(value = score)
pw_pairs <- data.frame(pw1 = rep(all_pathways, each = 149)) %>%
  group_by(pw1) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(.x, pw2 = unique(all_pathways)))) %>%
  unnest() %>%
  filter(!pw1==pw2) 
colnames(pw_pairs) <- c("pw1", "pw2")
data_added1 <- pw_pairs %>%
  ungroup() %>%
  dplyr::select(pw1) %>%
  dplyr::rename(Pathway = pw1) %>%
  inner_join(data_long, by = "Pathway") %>%
  unique() %>%
  dplyr::rename(value1 = value, Pathway1 = Pathway) %>%
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
  dplyr::rename(value2 = value, Pathway2 = Pathway) %>%
  dplyr::slice(rep(row_number(), 148)) %>%
  arrange(desc(Pathway2)) %>%
  unique() %>%
  group_by(Tissue, Group) %>%
  nest() %>%
  dplyr::rename(pw2 = data) 

data_ratios =full_join(data_added1, data_added2, by = c("Tissue", "Group")) %>%
  mutate(combined =  map2(pw1, pw2, ~ cbind(.x,  .y))) %>%
  dplyr::select(-c(pw1, pw2)) %>%
  unnest(cols = combined) %>%
  unique() %>%
  mutate(ratio = value1 / value2)%>%
  group_by(Pathway1, Pathway2) %>%
  mutate(average_ratio = mean(ratio)) %>%
  mutate(corrected_ratio = ratio / average_ratio)


mtPPS <- data_ratios %>%
  filter(!is.na(corrected_ratio)) %>%
  group_by(Tissue, Pathway1) %>%
  mutate(mtPPS = mean(corrected_ratio)) %>%
  select(Tissue, Group, Pathway1, mtPPS) %>%
  unique() %>%
  pivot_wider(names_from = "Pathway1", values_from = "mtPPS")

CI_FAO_ratios <- data_ratios %>% filter(Pathway1 %in% "Complex I", Pathway2 %in% "Fatty acid oxidation")



# Suppl FigS5D Plot Complex I corrected ratios --------------------------------------------------


## Cortex ------------------------------------------------------------------


cortex_CI <- data_ratios %>%
  filter(Pathway1 %in% "Complex I") %>%
  filter(Tissue %in% "Cerebral cortex") %>%
  group_by(Pathway1) %>%
  mutate(mtPPS = mean(corrected_ratio))
p<- cortex_CI %>% 
  ggplot(aes(label = Pathway2)) +
  geom_segment(aes(x=dplyr::desc(reorder(Pathway2,corrected_ratio)),
                   xend=dplyr::desc(reorder(Pathway2, corrected_ratio)),
                   y=0, yend = corrected_ratio), 
               size=0.6, alpha =0.7, color = "#1DB100") +
  labs(title= "Complex I ratios cortex") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = unique(cortex_CI$mtPPS), linetype = "dotted") +
  theme_classic() +
  theme(legend.position = "none",  
        plot.title = element_text(size =8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(p)

## Liver -------------------------------------------------------------------

liver_CI <- data_ratios %>%
  filter(Pathway1 %in% "Complex I") %>%
  filter(Tissue %in% "Liver")%>%
  group_by(Pathway1) %>%
  mutate(mtPPS = mean(corrected_ratio))

p<- liver_CI %>% 
  ggplot(aes(label = Pathway2)) +
  geom_segment(aes(x=dplyr::desc(reorder(Pathway2, corrected_ratio)),
                   xend=dplyr::desc(reorder(Pathway2, corrected_ratio)),
                   y=0, yend=corrected_ratio), 
               size=0.6, alpha =0.7, color = "#EF5FA7") +
  labs(title= "Complex I ratios liver") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = unique(liver_CI$mtPPS), linetype = "dotted") +
  theme_classic() +
  theme(legend.position = "none",  
        plot.title = element_text(size =8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(p)

# Suppl FigS5E Complex I all tissues ---------------------------------------------------

p <-  data_ratios %>%
  filter(Pathway1 %in% "Complex I") %>%
  group_by(Tissue) %>%
  mutate(mtPPS = mean(corrected_ratio)) %>%
  select(Tissue, Group, mtPPS) %>%
  mutate( `Complex I` =   mtPPS) %>%
  ggplot(aes(color = Group)) +
  geom_segment(aes(x=dplyr::desc(reorder(Tissue,  mtPPS)),
                   xend=dplyr::desc(reorder(Tissue,   mtPPS)),
                   y=0, yend= mtPPS), 
               size=1.3, alpha =0.7) +
  labs(title= "Complex I mitoPPS") +
  scale_color_manual(values = color_groups) + 
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none",  
        plot.title = element_text(size =8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(p)
