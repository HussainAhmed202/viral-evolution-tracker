# 05_protein_hotspots.R
# Translate raw (unaligned) sequences to protein first, then align proteins,
# then find amino-acid-level mutation hotspots

library(Biostrings)
library(DECIPHER)

# Use the RAW sequences (no gaps yet) so reading frame is intact
raw_seqs <- readDNAStringSet("data/raw/hiv_pol_sequences.fasta")

# Trim each sequence to a length that's a multiple of 3
# (partial/incomplete codons at the end would otherwise break translation)
trim_to_codon <- function(s) {
  len <- length(s)
  trimmed_len <- len - (len %% 3)
  subseq(s, start = 1, end = trimmed_len)
}

raw_seqs_trimmed <- DNAStringSet(lapply(raw_seqs, trim_to_codon))

# Translate each sequence in its own frame
protein_seqs <- translate(raw_seqs_trimmed, if.fuzzy.codon = "solve")

# Now align the PROTEIN sequences directly (gaps only inserted now, post-translation)
protein_aligned <- AlignSeqs(protein_seqs)

writeXStringSet(protein_aligned, filepath = "data/aligned_hiv_pol_protein.fasta")

cm_aa <- consensusMatrix(protein_aligned, as.prob = TRUE)

shannon_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log2(p))
}

aa_entropy <- apply(cm_aa, 2, shannon_entropy)

hotspot_table_aa <- data.frame(
  aa_position = seq_along(aa_entropy),  # this should be amino acid position
  entropy = aa_entropy
)

top_hotspots_aa <- hotspot_table_aa[order(-hotspot_table_aa$entropy), ][1:20, ]

write.csv(hotspot_table_aa, "results/protein_position_entropy.csv", row.names = FALSE)
write.csv(top_hotspots_aa, "results/top_20_protein_hotspots.csv", row.names = FALSE)

cat("Top 10 amino acid mutation hotspots (codon position : entropy):\n")
print(head(top_hotspots_aa, 10))