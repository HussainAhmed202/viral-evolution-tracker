# 04_analysis.R
# Phylogenetic tree + mutation hotspot detection from aligned HIV-1 pol sequences

library(Biostrings)
library(DECIPHER)
library(ape)

aligned <- readRDS("data/aligned_hiv_pol.rds")

# ---- 1. Phylogenetic tree ----

# Compute pairwise distance matrix from the alignment
dist_matrix <- DistanceMatrix(aligned, correction = "JC69") # Jukes-Cantor

# Build a neighbor-joining tree
tree <- nj(as.dist(dist_matrix))  # typecasting because nj only accepts dist datatype objects

# Save tree to file (Newick format - standard, usable in other tools too)
write.tree(tree, file = "results/hiv_pol_tree.nwk")

# Save a plot of the tree
dir.create("figures", showWarnings = FALSE)
png("figures/phylogenetic_tree.png", width = 1000, height = 1000, res = 120)
plot(tree, type = "unrooted", cex = 0.5, main = "HIV-1 pol Neighbor-Joining Tree")
dev.off()

cat("Tree saved to results/hiv_pol_tree.nwk and figures/phylogenetic_tree.png\n")

# ---- 2. Mutation hotspot detection ----

# Get per-position base composition across all sequences
cm <- consensusMatrix(aligned, as.prob = TRUE)

# Keep only the 4 nucleotide rows (drop gaps/ambiguity codes for this calc)
cm_bases <- cm[c("A", "C", "G", "T"), ]

# Shannon entropy per position = measure of how "mixed"/variable that
# position is across our 40 sequences. High entropy = high variablity in code which means it could be a mutation hotspot.
shannon_entropy <- function(p) {
  p <- p[p > 0]  # avoid log(0)
  -sum(p * log2(p))
}

position_entropy <- apply(cm_bases, 2, shannon_entropy)  # 2 = apply across columns

# Build a results table: position + entropy score
hotspot_table <- data.frame(
  position = seq_along(position_entropy),
  entropy = position_entropy
)

# Sort to find the top hotspots
top_hotspots <- hotspot_table[order(-hotspot_table$entropy), ][1:20, ]  # order by entropy desc

write.csv(hotspot_table, "results/position_entropy.csv", row.names = FALSE)
write.csv(top_hotspots, "results/top_20_hotspots.csv", row.names = FALSE)

cat("Top 5 mutation hotspot positions:\n")
print(head(top_hotspots, 5))