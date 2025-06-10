## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  03 March 2025 

## Notes: 
#  

## Clear environment 
rm(list = ls()) 

## Set seed 
set.seed(1234)

## Set functions: 

## set working directory 
setwd("/Users/panayiotislouca/Documents")

## load up packages: 

### core 
library(tidyverse)

# for network plots 
library(igraph)
library(ggraph)

# -------------------------------------------------------------------------- # 

# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 

# -------------------------------------------------------------------------- #  

##   Dataset ---- 
path <-
  file.path(
    "/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Analysis/data/phages_microbes_metabs_full_DATASET_SCFA_v2.csv"
  )
df <- read.csv(path) %>% as.data.frame()

# import phage - SCFA/BA results 
path <- file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Analysis/Phage_first_approach/combine_bile_acids_SCFA/combined_phage_metab_results.csv")
phage_metab_res <- read.csv(path)

phage_metab_res %>% count(depend_var)
phage_metab_res %>% pull(independ_var) %>% n_distinct()

phage_metab_res_sig <- phage_metab_res %>%
  filter(fdr_combined < 0.001)

# phage assoc metab names 
phage_assoc_metab_names = phage_metab_res_sig %>%
  pull(depend_var) %>%
  unique()

phage_metab_res_sig$independ_var %>% n_distinct()
metab_assoc_phage_names = phage_metab_res_sig %>%
  pull(independ_var) %>%
  unique()

# ----------------------------------------------------------------------------------------------------------------------- #  

df_plot <- df %>%
  select(all_of(phage_assoc_metab_names),
         all_of(metab_assoc_phage_names))

# -------------------------------------------------------------------------- #  

# Correlation Network Visualization 
phage_data <- df[ ,metab_assoc_phage_names]
metab_data <- df[ ,c(phage_assoc_metab_names)]

# -------------------------------------------------------------------------- #  

# Compute Spearman correlations between phages and bacteria 
cor_matrix <- cor(phage_data, metab_data, method = "spearman", use = "pairwise.complete.obs")

# Convert to long format and filter for moderate correlations abs(> 0.3) 
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Phage = Var1, Metabolite = Var2, Correlation = Freq) %>%
  filter(abs(Correlation) > 0.3) %>% 
  mutate(EdgeColour = ifelse(Correlation > 0, "Positive", "Negative"))

length(unique(cor_df$Phage))
length(unique(cor_df$Metabolite))

#  Compute p-values for these pairs 
cor_df <- cor_df %>%
  mutate(p_value = pmap_dbl(list(Phage, Metabolite), function(phg, metab) {
    x <- phage_data[[phg]]
    y <- metab_data[[metab]]
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  }))

summary(cor_df$p_value)
summary(cor_df$Correlation)

# -------------------------------------------------------------------------- #  

cor_df %>% count(Metabolite)

metab_mapping <- c( # create a metab mapping to clean names for plot 
  "Isoursodeoxycholate (serum)" = "serum_metab_M57577",
  "Butyrate (stool)" = "stool_SCFA_Butyrate",
  "Propionate (stool)" = "stool_SCFA_Propionate",
  "Deoxycholate (stool)" = "stool_metab_M1114",
  "Dehydrolithocholate (stool)" = "stool_metab_M31891",
  "7-Ketodeoxycholate (stool)" = "stool_metab_M31904",
  "Hyocholate (stool)" = "stool_metab_M34093",
  "Ursocholate (stool)" = "stool_metab_M57573",
  "Isoursodeoxycholate (stool)" = "stool_metab_M57577"
)

# Reverse the mapping so names match cor_df$Metabolite 
metab_lookup <- setNames(names(metab_mapping), metab_mapping)

# Create a dataframe for the mapping 
metab_lookup_df <- tibble(
  Metabolite = names(metab_lookup),
  Metabolite_label = unname(metab_lookup)
)

# Apply the mapping  
cor_df <- cor_df %>%
  left_join(metab_lookup_df, by = "Metabolite") %>%
  mutate(Metabolite = Metabolite_label) %>%
  select(-Metabolite_label)

cor_df %>% count(Metabolite)

# -------------------------------------------------------------------------- #  

# Create network and set attributes  
net <- graph_from_data_frame(cor_df[, c("Phage", "Metabolite", "Correlation")], directed = FALSE)

# Set edge weights and colors
E(net)$weight <- abs(cor_df$Correlation)
E(net)$color <- cor_df$EdgeColour

table(E(net)$color, cor_df$EdgeColour)

# Set node attributes 
V(net)$type <- ifelse(grepl("phage", V(net)$name), "Phage", "Metabolite")

# Label all bacteria and top 2% of other nodes
degree_vals <- degree(net)

V(net)$label <- gsub("_", " ", V(net)$name) # Store all names 


# -------------------------------------------------------------------------- #  

n_connections <- 100  # Threshold for labeling nodes 
Metabolite_nodes <- V(net)$name[V(net)$type == "Metabolite"]

# -------------------------------------------------------------------------- #  

# set seed 
set.seed(1234) 

# save plot 
path <- file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metabs_obesity/Write_up/Figures/phage_metabolite_network_plots/Phage_SCFA_bile_acids_ggraph_network_v2.png")
png(path, width = 16, height = 14, units = 'in', res = 300)

# Plot network 
ggraph(net, layout = "fr") +
  geom_edge_link0(aes(edge_alpha = weight, edge_colour = color, edge_width = weight), 
                  show.legend = FALSE) +
  scale_edge_alpha_continuous(range = c(0.2, 1)) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  geom_node_point(aes(color = type, size = degree_vals, 
                      shape = ifelse(degree_vals > quantile(degree_vals, 0.95), "triangle", "circle")), 
                  alpha = 0.9) +
  geom_node_text(aes(label = ifelse(label == "Isoursodeoxycholate (serum)" | degree_vals >= n_connections, label, "")),
                 size = 4,
                 repel = TRUE,
                 force = 8,
                 box.padding = 0,
                 point.padding = 0,
                 max.overlaps = 4500,
                 color = "black", fontface = "bold", segment.color = "grey80") +
  scale_edge_colour_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8"),  
                           name = "Correlation") +
  scale_color_manual(values = c("Phage" = "#6A3D9A", "Metabolite" = "#33A02C"),
                     name = "Node Type") +
  scale_shape_manual(values = c("circle" = 16, "triangle" = 17), guide = "none") +
  scale_size_area(max_size = 6) +
  guides(
    size = "none",
    color = guide_legend(
      title = "Node Type",
      override.aes = list(size = 6, shape = 16),
      keywidth = 1.5,
      keyheight = 1.2
    )
  ) +
  theme_void() +
  theme(
    # Panel and plot appearance 
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # Legend styling 
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5)        ) +
  labs(title = "Phage-SCFA/Bile Acid Correlation Network",
       subtitle = "Rho > 0.3; Red = Positive, Blue = Negative" 
       
  )

