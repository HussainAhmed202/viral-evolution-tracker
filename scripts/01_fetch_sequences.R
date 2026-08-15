# 01_fetch_sequences.R
# Fetch real HIV-1 pol gene (protease + RT) sequences from NCBI GenBank

library(rentrez)

# Search NCBI nucleotide database for HIV-1 pol gene sequences
# Filtering: HIV-1, "pol" gene, length restricted to roughly the
# protease + partial RT region (avoids full genomes / full pol CDS)
search_result <- entrez_search(
  db = "nuccore",
  term = "HIV-1[Organism] AND pol[Gene] AND 1000:3000[SLEN]",
  retmax = 40
)

# How many did we find, and what are the IDs?
print(search_result$count)
print(search_result$ids)

# Fetch the actual sequences (FASTA format) for the IDs we found
sequences_fasta <- entrez_fetch(
  db = "nuccore",
  id = search_result$ids,
  rettype = "fasta"
)

# Save raw FASTA to file
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
writeLines(sequences_fasta, "data/raw/hiv_pol_sequences.fasta")

cat("Saved", length(search_result$ids), "sequences to data/raw/hiv_pol_sequences.fasta\n")