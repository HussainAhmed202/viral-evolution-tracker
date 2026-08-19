# 06_map_to_hxb2.R
# Map our alignment's amino acid positions onto standard HXB2 reference numbering

library(Biostrings)
library(rentrez)

# ---- 1. Fetch HXB2 pol CDS from NCBI ----
# K03455 is the HXB2 reference genome. fasta_cds_na returns each coding
# sequence (CDS) separately as nucleotide FASTA, which lets us pull out
# just the pol gene without hand-parsing genome coordinates.
hxb2_cds <- entrez_fetch(
  db = "nuccore",
  id = "K03455.1",
  rettype = "fasta_cds_na"
)

writeLines(hxb2_cds, "data/raw/hxb2_all_cds.fasta")

# Load all CDS entries and find the one whose header mentions "pol"
all_cds <- readDNAStringSet("data/raw/hxb2_all_cds.fasta")
pol_entry <- all_cds[grepl("pol", names(all_cds), ignore.case = TRUE)]

cat("Found", length(pol_entry), "CDS entr(y/ies) matching 'pol':\n")
print(names(pol_entry))

# Take the first pol match, trim to a multiple of 3, translate
hxb2_pol_nt <- pol_entry[[1]]
trimmed_len <- width(hxb2_pol_nt) - (width(hxb2_pol_nt) %% 3)
hxb2_pol_nt <- subseq(hxb2_pol_nt, start = 1, end = trimmed_len)
hxb2_pol_protein <- translate(DNAStringSet(hxb2_pol_nt), if.fuzzy.codon = "solve")[[1]]

cat("HXB2 pol protein length:", length(hxb2_pol_protein), "aa\n")

# ---- 2. Build a consensus sequence from our own protein alignment ----
protein_aligned <- readAAStringSet("data/aligned_hiv_pol_protein.fasta")
cm_aa <- consensusMatrix(protein_aligned, as.prob = FALSE)

# Most common amino acid at each position = our consensus
our_consensus_chars <- apply(cm_aa, 2, function(col) names(which.max(col)))
our_consensus <- AAString(paste(our_consensus_chars, collapse = ""))

# ---- 3. Pairwise align our consensus to HXB2 reference protein ----
alignment <- pairwiseAlignment(
  pattern = our_consensus,
  subject = hxb2_pol_protein,
  type = "global-local"  # our sequence may be a fragment of HXB2's full pol
)

print(alignment)

# ---- 4. Build position mapping: our alignment position -> HXB2 position ----
our_aligned_seq <- as.character(pattern(alignment))
hxb2_aligned_seq <- as.character(subject(alignment))

our_pos <- 0
hxb2_pos <- start(subject(alignment)) - 1  # offset to true HXB2 start
mapping <- data.frame(our_position = integer(0), hxb2_position = integer(0))

for (i in seq_len(nchar(our_aligned_seq))) {
  our_char <- substr(our_aligned_seq, i, i)
  hxb2_char <- substr(hxb2_aligned_seq, i, i)
  
  if (our_char != "-") our_pos <- our_pos + 1
  if (hxb2_char != "-") hxb2_pos <- hxb2_pos + 1
  
  if (our_char != "-" && hxb2_char != "-") {
    mapping <- rbind(mapping, data.frame(our_position = our_pos, hxb2_position = hxb2_pos))
  }
}

write.csv(mapping, "results/position_mapping_to_hxb2.csv", row.names = FALSE)

# ---- 5. Apply mapping to our top hotspots ----
top_hotspots_aa <- read.csv("results/top_20_protein_hotspots.csv")
top_hotspots_mapped <- merge(top_hotspots_aa, mapping,
                             by.x = "aa_position", by.y = "our_position")
top_hotspots_mapped <- top_hotspots_mapped[order(-top_hotspots_mapped$entropy), ]

write.csv(top_hotspots_mapped, "results/top_hotspots_hxb2_mapped.csv", row.names = FALSE)

cat("\nTop hotspots mapped to HXB2 numbering:\n")
print(top_hotspots_mapped)