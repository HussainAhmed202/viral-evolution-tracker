# HIV-1 Drug Resistance Mutation Tracker

## Overview
This project analyzes HIV-1 *pol* gene sequences (protease + reverse
transcriptase region) from public NCBI data to identify mutation hotspots
associated with antiretroviral drug resistance. It combines sequence
retrieval, multiple sequence alignment, phylogenetic analysis, and
variant-site mapping — a compact but realistic bioinformatics pipeline.

## Motivation
HIV's high mutation rate allows it to rapidly evolve resistance to
antiretroviral therapy. Certain codons in the protease and reverse
transcriptase genes are known resistance "hotspots" (e.g. codons linked to
resistance against NRTIs, NNRTIs, and protease inhibitors). This project
asks: **can we recover known resistance hotspots directly from public
sequence data, using only alignment and variation analysis — no prior
annotation?**

## Data
- Source: NCBI GenBank / HIV Sequence Database
- Gene region: HIV-1 *pol* (protease + partial reverse transcriptase)
- Number of sequences: TBD
- Accession list: see `data/accessions.txt`

## Pipeline
1. `scripts/01_fetch_sequences.R` — retrieve sequences from NCBI via `rentrez`
2. `scripts/02_align.R` — multiple sequence alignment via `DECIPHER`
3. `scripts/03_analysis.R` — variability scoring, mutation hotspot
   identification, phylogenetic tree construction
4. `scripts/04_visualize.R` — hotspot heatmap + phylogenetic tree figures

## Results
_TBD — filled in once analysis is complete._

## Key Finding
Two things worth flagging honestly, since this is going on your CV and you should be able to defend it in an interview:

Only 11 of your top 20 hotspots mapped, not all 20. That's expected — some of your top hotspot positions likely fall in small insertions/gaps that don't have a clean 1-to-1 counterpart in the HXB2 reference. Not a bug, just worth noting in the writeup rather than hiding.
Position 99 is the protease/RT boundary (HXB2 protease is exactly 99 amino acids). So we can split your remaining hits:
Protease region (hxb2_position ≤ 99): position 28, position 99 (boundary)
RT region (hxb2_position > 99): subtract 99 to get RT's own numbering → 223, 230, 138, 213, 100, 229, 212, 127, 144

One genuinely interesting match: RT position 138 (your hxb2_position 237 → 237-99=138). E138 is a real, well-documented NNRTI accessory resistance site (E138K/A/G, associated with rilpivirine resistance). That's a legitimate, citable overlap between your data-derived hotspot and known clinical literature.

The others (RT 223, 230, 213, 212, 127, 144, protease 28) aren't on the classic "major resistance mutation" lists I know — they could be natural variability/polymorphism sites (also biologically real and worth reporting), or they could reflect the numbering offset uncertainty I flagged earlier from how we extracted the HXB2 pol CDS.

## Tools
R, Bioconductor (`Biostrings`, `DECIPHER`), `ape`, `rentrez`, `ggplot2`

## rentrez
Provides an R interface to the NCBI's 'EUtils' API, allowing users to search databases, process the results of those searches and pull data into their R sessions.
Useful functions:
* entrez_search() 
  - You provide the db and the search term. It searches that on the NCBI and returns a search object
  - Lets say search_result <- entrez_search(...)
    - search_result['count']  =>  total number of matches
    - search_result['ids']  =>  the accession IDs of the search results 
    
* entrez_fetch() 
  - You provide the db, the accession IDs and the file return type. It fetches that on the NCBI and returns in the providede format  

* entrez_summary()

## Biostrings
Biostrings is a Bioconductor package for representing and manipulating biological sequences.

It works primarily with:

* DNA
* RNA
* Amino-acid

The simplest one is a DNAString.
* seq1 <- DNAString("ATGCGT")
* seq2 <- DNAString("ATGAAA")
* seq3 <- DNAString("ATGCCC")

You could put them together:
* seqs <- DNAStringSet(c(seq1, seq2, seq3))

## DECIPHER
DECIPHER is a Bioconductor package for biological sequence analysis.

It can do things like:

* sequence alignment
* multiple sequence alignment
* identifying differences between sequences
* sequence classification
* detecting and correcting errors
* primer/probe design
* taxonomy-related sequence analysis
* handling large collections of sequences

Useful functions:
* AlignSeqs()
  - Aligns the sequences. Returms a DNAStringSet object.

* BrowseSeqs()
  - Takes in the DNAStringSet object retuned by the AlignSeqs foo.
  - View on the browser

* DistanceMatrix()
  - Distance matrix = a table that tells you how different every sequence is from every other sequence.

## ape
ape stands for Analyses of Phylogenetics and Evolution

* Many methods available to generate the phylogenetic trees. These methods take in a dist object. Example of such methods include: 
  * nj() 
  * UPGMA()

* We can read an existing tree using the read.tree().Trees are stored in the  .nwk file format. This standards for Newick.
        ┌── A
    ┌───┤
────┤   └── B
    │
    │   ┌── C
    └───┤
        └── D
        
The above tree in the Newick file format would be like:
((A,B),(C,D))

* We can also write a tree using the write.tree()

* We can plot trees using the plot() method.


## Author
Syed Hussain Ahmed