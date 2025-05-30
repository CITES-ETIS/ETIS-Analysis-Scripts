################################################################################
# Perform cluster analysis and prepare heatmap and dendrogram                  #
#                                                                              #
################################################################################


#### Run cluster analysis ------------------------------------------------------

# Load cluster variables
cluster_vars <- read_csv(paste0("Processed Data/cluster_vars_", analysis_name, ".csv"))

# Remove countries with no seizures in or out during the CoP_Period
excl_ctys <- cluster_vars %>%
  filter(if_all(c(wt_in_1, sz_out_1, wt_out_1, wt_in_2, sz_out_2, wt_out_2), ~ .x == 0)) %>%
  pull(country)

cluster_vars <- cluster_vars %>%
  filter(!(country %in% excl_ctys))

# Log transform all variables
cluster_vars <- cluster_vars %>%
  mutate(across(where(is.numeric), ~ log(.x + 1)))

# Get variables for clustering
X <- select(cluster_vars, -country)

# Run agglomerative hierarchical clustering
cl <- agnes(X, stand = TRUE, method = "ward")  



#### Prepare dendrogram plot ---------------------------------------------------

# Convert clustering to dendrogram
dendrogram <- as.dendrogram(cl)

# Sort using dendsort and reverse order to put higher values on the left
dend.sorted <- rev(dendsort(dendrogram))

# Extract dendrogram data for plotting
dend.data <- dendro_data(dend.sorted)

# Dendrogram plot
dendrogram_plot <- ggplot(data = dend.data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_y_continuous(limits = c(0, max(dendro_data$segments$yend)), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(1-0.5,length(labels(dend.sorted))+0.5), expand = c(0, 0)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold", size=10),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),  
    axis.line.y = element_line(colour = "grey"),
    axis.ticks.y = element_line(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.margin = margin(t = 5, b = 2)  # Reduce bottom margin so that dendrogram leaves reach closer to heatmap columns when aligned
  )



#### Prepare heatmap plot ------------------------------------------------------

# Reorder cluster_vars countries by the order in the dendrogram, and rename NAM to NA
cty_order <- labels(dend.sorted)
cluster_vars$country[which(cluster_vars$country == "NAM")] <- "NA"
cluster_vars$country <- factor(cluster_vars$country, levels = unique(cluster_vars$country)[cty_order])

# Standardize variables and reformat for ggplot
heatmap_vars <- cluster_vars %>%
  mutate(across(where(is.numeric), ~ scale(.))) %>%
  pivot_longer(cols = -country, names_to = "variable", values_to = "value")
# Note the standardization here is scale which divides by sd, whereas the standardization in the agnes clustering function divides by mean absolute deviation. Negligible difference to visualization of variables

# Reverse order variables for heatmap y-axis display
heatmap_vars$variable <- factor(heatmap_vars$variable, levels = rev(c("lambda_1", "lambda_2", "lambda_3", "lambda_4", "lambda_5", "wt_in_1", "wt_in_2", "sz_out_1", "sz_out_2", "wt_out_1", "wt_out_2")))

# Heatmap plot
heatmap_plot <- ggplot(data = heatmap_vars) +
  geom_tile(aes(x = country, y = variable, fill = value)) +
  scale_fill_gradientn(name = "standardized\nvalue", colors = c("#F5ECE0", "#F0D8BA", "#EDA674", "#E56A3C", "#C2372C", "#7A1E14")) + 
  ylab(" ") + xlab("") +
  scale_y_discrete(expand = c(0, 0), labels = c("lambda_1" = "TI raw < 10 kg", "lambda_2" = "TI raw 10 - 100 kg", "lambda_3" = "TI raw \u2265 100 kg", "lambda_4" = "TI worked < 1 kg", "lambda_5" = "TI worked \u2265 1 kg", 
                                                "wt_in_1" = "wt-in < 500 kg", "sz_out_1" = "sz-out < 500 kg", "wt_out_1" = "wt-out < 500 kg", "wt_in_2" = "wt-in \u2265 500 kg", "sz_out_2" = "sz-out \u2265 500 kg", "wt_out_2" = "wt-out \u2265 500 kg")) +
  scale_x_discrete(expand = c(0, 0)) +  # expand = c(0, 0) removes padding between heatmap tiles and axis
  theme_bw() +
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=10, hjust=0.5, vjust=0.5 ,angle=90),  
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        axis.line = element_line(colour = "grey"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 0))      # Reduce top margin so that dendrogram leaves reach closer to heatmap columns when aligned


