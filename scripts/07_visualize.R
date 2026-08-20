# 07_visualize.R
# Final visualizations: mutation hotspot plot + polished phylogenetic tree

library(ggplot2)
library(ape)

dir.create("figures", showWarnings = FALSE)

# ---- 1. Mutation hotspot entropy plot (protein-level) ----

hotspot_data <- read.csv("results/protein_position_entropy.csv")

# Flag top 20 positions so we can highlight them distinctly
top20 <- read.csv("results/top_20_protein_hotspots.csv")
hotspot_data$is_hotspot <- hotspot_data$codon_position %in% top20$codon_position

p <- ggplot(hotspot_data, aes(x = codon_position, y = entropy)) +
  geom_line(color = "grey60", linewidth = 0.3) +
  geom_point(data = subset(hotspot_data, is_hotspot),
             aes(x = codon_position, y = entropy),
             color = "firebrick", size = 2) +
  labs(
    title = "HIV-1 pol Amino Acid Variability Across 40 Public Sequences",
    subtitle = "Red points = top 20 mutation hotspots (highest Shannon entropy)",
    x = "Amino acid position (alignment-relative)",
    y = "Shannon entropy (variability)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/mutation_hotspot_plot.png", plot = p, width = 10, height = 5, dpi = 150)

cat("Saved figures/mutation_hotspot_plot.png\n")

# ---- 2. Polished phylogenetic tree ----

tree <- read.tree("results/hiv_pol_tree.nwk")
tree$tip.label <- sub("\\..*$", "", sub(" .*", "", tree$tip.label))
tree$tip.label <- sub("^([A-Z]+[0-9]+\\.[0-9]+).*", "\\1", tree$tip.label)

png("figures/phylogenetic_tree_final.png", width = 1100, height = 1100, res = 130)
plot(
  tree,
  type = "fan",              # circular layout - reads better with many tips
  cex = 0.55,
  no.margin = TRUE,
  edge.color = "grey40",
  tip.color = "steelblue4"
)
title("HIV-1 pol Neighbor-Joining Tree (40 public sequences)", cex.main = 0.9)
# Shorten tip labels to just the accession number (e.g. "PX471211.1")
dev.off()

cat("Saved figures/phylogenetic_tree_final.png\n")