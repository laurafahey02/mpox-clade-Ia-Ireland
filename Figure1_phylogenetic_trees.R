# -------------------------------
# Visualisation of MPXV phylogenetic trees:
# Panel A: Full tree highlighting NVRL cluster
# Panel B: Subtree with mutation rectangles
# -------------------------------

# Load required libraries
library(ggplot2)
library(ggtree)
library(dplyr)
library(ggnewscale)
library(cowplot)
library(tidyr)

# -----------------------------
# PANEL A: Full tree (Clade 1) highlighting NVRL cluster
# -----------------------------

# Read tree and prune outgroups
tree <- read.tree("957_seq.treefile")
tips_to_drop <- c(
  "hMpxV/Ireland/un-NVRL-12199/2025|EPI_ISL_19742488|2025-02-05",
  "KJ642617.1",
  "KJ642615.1"
)
tree_pruned <- drop.tip(tree, tips_to_drop)
p <- ggtree(tree_pruned)

# Collapse Clade 1b
collapsed_tree <- collapse(p, node = 1309, mode = 'max', fill = '#666666') +
  geom_cladelab(
    node = 1309,
    label = "Clade 1b",
    fontsize = 3,
    fontface = "bold",
    vjust = +1.4,
    offset.text = -0.00005
  )

# Identify descendant nodes
node1061_descendants <- offspring(p, 1222)$node
node1076_descendants <- offspring(p, 1241)$node
node1061_only <- setdiff(node1061_descendants, node1076_descendants)

# Annotate clusters
tree_data <- collapsed_tree$data %>%
  mutate(cluster = case_when(
    node %in% node1076_descendants ~ "NVRL cluster",
    node %in% node1061_only ~ "NVRL parent",
    TRUE ~ "Other"
  )) %>%
  mutate(cluster = factor(cluster, levels = c("Other", "NVRL parent", "NVRL cluster")))

# Legend tree
legend_tree <- collapsed_tree +
  geom_tippoint(
    data = tree_data,
    aes(fill = cluster, shape = cluster),
    size = 4,
    color = "black",
    stroke = 0.5
  ) +
  scale_shape_manual(
    values = c("Other" = 21, "NVRL parent" = 21, "NVRL cluster" = 23),
    labels = c(
      "Other" = "clade Ia (other)",
      "NVRL cluster" = "2024 Kinshasa outbreak - subcluster with Irish sequence",
      "NVRL parent" = "clade Ia - 2024 Kinshasa outbreak"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Other" = "#0072B2",
      "NVRL cluster" = "#1B9E77",
      "NVRL parent" = "#CC79A7"
    ),
    labels = c(
      "Other" = "clade Ia (other)",
      "NVRL cluster" = "2024 Kinshasa outbreak - subcluster with Irish sequence",
      "NVRL parent" = "clade Ia - 2024 Kinshasa outbreak"
    ),
    name = NULL
  ) +
  theme(
    legend.position = c(0.59, 0.65),
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  )

# Add bootstrap values
nodes_to_label <- c(959, 961, 962, 983, 984, 1188, 1189, 1218, 1221, 1222, 1241)

node_labels <- tree_data %>%
  filter(node %in% nodes_to_label) %>%
  mutate(hjust = case_when(
    node == 1241 ~ 1.9,
    node == 1222 ~ 1.8,
    TRUE ~ 1.2
  ))

nodes_labelled <- legend_tree +
  geom_text2(
    data = node_labels,
    aes(x = x, y = y, label = label, hjust = hjust),
    vjust = -0.7,
    size = 3.5,
    fontface = "bold"
  )

# Adjust ylim and add scale
pt_s <- nodes_labelled + ylim(-30, 959)

pt_d <- pt_s +
  geom_treescale(x = 0, y = 0, fontsize = 4, offset = -30, linesize = 1) +
  annotate("text", x = 0.00009, y = -24, label = "subs/site", size = 4)

# Clean labels
pt_d$data$label <- pt_s$data$label %>%
  gsub("^'(.*)'$", "\\1", .) %>%
  gsub("^hMPXV/Ireland/NVRL-12199/2025$", "Ireland/NVRL-12199/2025-02-05", .)

# Final plot for panel A
fig1a <- pt_d +
  geom_tiplab(
    data = subset(pt_d$data, label == "Ireland/NVRL-12199/2025-02-05") %>%
      mutate(y = y + 0.3),
    aes(label = label),
    fontface = "bold",
    size = 4,
    offset = 0.00001
  )


# -----------------------------
# PANEL B: Subtree showing mutation counts as rectangles
# -----------------------------

# Read subtree and mutation data
tree <- read.tree("subcluster.tree")
tree2 <- tree[[2]]
p <- ggtree(tree2)

tree_data <- p$data %>%
  mutate(cluster = case_when(
    node == 1 & isTip ~ "Other",
    node == 2 ~ "NVRL parent",
    isTip ~ "NVRL cluster",
    TRUE ~ NA_character_
  ))

# Read mutation counts
mut <- read.table("mutations.txt", header = TRUE)

# Merge mutation data
tree_data_mut <- tree_data %>%
  left_join(mut, by = "node")

# Prepare rectangles
tree_plot <- ggtree(tree_data_mut)
tree_data_df <- tree_plot$data

mutation_rects <- tree_data_df %>%
  filter(!is.na(APOBEC) & !is.na(Other)) %>%
  pivot_longer(cols = c("APOBEC", "Other"), names_to = "mutation_type", values_to = "count") %>%
  filter(count > 0) %>%
  group_by(node) %>%
  uncount(count) %>%
  mutate(offset = row_number())

# Define colours
mutation_colors <- c("APOBEC" = "#D55E00", "Other" = "#F0E442")

# Fill colours for tips
tree_data <- tree_data_df %>%
  mutate(fill_color = case_when(
    cluster == "Other" ~ "#0072B2",
    cluster == "NVRL cluster" ~ "#1B9E77",
    cluster == "NVRL parent" ~ "#CC79A7",
    TRUE ~ "#cccccc"
  ))

# Build plot
col_p <- p +
  geom_point2(
    data = tree_data,
    aes(x = x, y = y, fill = I(fill_color), shape = cluster),
    size = 4,
    color = "black",
    stroke = 0.5,
    show.legend = FALSE
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = mutation_rects,
    aes(x = branch + 0.000002 + offset * -0.0000015, y = y, fill = mutation_type),
    width = 0.0000009,
    height = 0.3,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = mutation_colors,
    name = "Mutation Types",
    labels = c(
      "APOBEC" = "APOBEC3 (TC→TT / GA→AA)",
      "Other" = "Other SNVs"
    )
  ) +
  scale_shape_manual(
    values = c("Other" = 21, "NVRL cluster" = 23, "NVRL parent" = 21),
    name = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.60),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.7, "lines")
  ) +
  annotate("text", x = 9.590e-05, y = 12.5, label = "Each rectangle = one mutation", size = 4, hjust = 0)

# Clean labels and highlight NVRL
col_p$data$label <- col_p$data$label %>%
  gsub("^'(.*)'$", "\\1", .) %>%
  gsub("^PP_002XE2K.1$", "Ireland/NVRL-12199/2025-02-05", .)

col_p$data$size <- ifelse(col_p$data$label == "Ireland/NVRL-12199/2025-02-05", 4, 3)
col_p$data$fontface <- ifelse(col_p$data$label == "Ireland/NVRL-12199/2025-02-05", "bold", "plain")

# Finalise labels and tree scale
p_pab <- col_p +
  geom_tiplab(aes(fontface = fontface, size = size), offset = 0.000001) +
  scale_size_identity()

p_s <- p_pab +
  geom_treescale(x = 0, y = 0, fontsize = 4, width = 5e-06, linesize = 1, offset = -0.8) +
  annotate("text", x = 0.000012, y = -0.65, label = "subs/site", size = 4) +
  annotate("text", x = 0.000030, y = -0.6, label = " (~ 1 SNV per genome)", size = 4)

p_x <- p_s + xlim(0, 0.00013)


# -----------------------------
# Combine Panels A and B into final figure
# -----------------------------

grid_plot <- plot_grid(fig1a, p_x, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1))

# Save final figure
ggsave("Figure1_combined_tree_panels.png", grid_plot, width = 10, height = 8, dpi = 350, bg = "white")