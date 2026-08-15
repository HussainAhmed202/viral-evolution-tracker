# 02_inspect_sequences.R
# Load and sanity-check the fetched HIV-1 pol sequences

library(Biostrings)

seqs <- readDNAStringSet("data/raw/hiv_pol_sequences.fasta")

# How many sequences, and how long is each one?
length(seqs)
seqs

# Quick check: are all sequence lengths in a sane range for our target region?
summary(width(seqs))