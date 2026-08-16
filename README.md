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
_TBD — one-sentence takeaway once we have results._

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

## Author
TBD