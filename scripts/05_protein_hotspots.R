# 05_protein_hotspots.R
# Translate alignment to protein and find amino-acid-level mutation hotspots

library(Biostrings)

aligned <- readRDS("data/aligned_hiv_pol.rds")

# Translate the nucleotide alignment to protein (reading frame 1)
# if.fuzzy.codon="solve" handles ambiguous bases (R, Y, etc.) reasonably
protein_aligned <- translate(aligned, if.fuzzy.codon = "solve")

# Save protein alignment
writeXStringSet(protein_aligned, filepath = "data/aligned_hiv_pol_protein.fasta")

# Per-position amino acid entropy (same idea as before, protein level)
cm_aa <- consensusMatrix(protein_aligned, as.prob = TRUE)

shannon_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log2(p))
}

aa_entropy <- apply(cm_aa, 2, shannon_entropy)

hotspot_table_aa <- data.frame(
  codon_position = seq_along(aa_entropy),
  entropy = aa_entropy
)

top_hotspots_aa <- hotspot_table_aa[order(-hotspot_table_aa$entropy), ][1:20, ]

write.csv(hotspot_table_aa, "results/protein_position_entropy.csv", row.names = FALSE)
write.csv(top_hotspots_aa, "results/top_20_protein_hotspots.csv", row.names = FALSE)

cat("Top 10 amino acid mutation hotspots (codon position : entropy):\n")
print(head(top_hotspots_aa, 10))