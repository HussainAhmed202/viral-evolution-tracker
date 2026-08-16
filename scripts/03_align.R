# 03_align.R
# Multiple sequence alignment of HIV-1 pol sequences using DECIPHER

library(Biostrings)
library(DECIPHER)

# Load the raw sequences
seqs <- readDNAStringSet("data/raw/hiv_pol_sequences.fasta")

# Align them
aligned <- AlignSeqs(seqs, anchor = NA)

# Save the alignment (as FASTA, and as an R object for later steps)
writeXStringSet(aligned, filepath = "data/aligned_hiv_pol.fasta")
saveRDS(aligned, "data/aligned_hiv_pol.rds")  # serialize the R object. Will refer to it later when creating the phylogenetic tree 

cat("Alignment complete.\n")
cat("Aligned length (with gaps):", width(aligned)[1], "bp\n")
cat("Number of sequences aligned:", length(aligned), "\n")