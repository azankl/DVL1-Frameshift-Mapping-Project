# This script recreates my patient's 1bp deletion in DVL1
# the results confirm that the deletion causes a W->G change at amino acid 532
# and a frameshift that leads to a premature stop codon 142 amino acids downstream.

# if (!require("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }

#BiocManager::install("Biostrings", ask = FALSE)

#if (!require("rentrez", quietly = TRUE)) {
#    install.packages("rentrez")
# }

library(Biostrings)
library(rentrez)

# Fetch the sequence
cds_record <- entrez_fetch(
    db = "nuccore",
    id = "NM_001330311.2",
    rettype = "fasta_cds_na"
)
writeLines(cds_record, "NM_001330311_CDS.fasta")

# Parse the CDS into a DNAString object for downstream analysis.
cds <- readDNAStringSet("NM_001330311_CDS.fasta")
cds_seq <- cds[[1]]
cds_seq

# Translate the CDS
protein <- translate(cds_seq)
protein

# Introduce a 1 bp deletion at position 1593 (1-based index)
fs_seq <- replaceAt(cds_seq, at = IRanges(1593, 1593), value = DNAStringSet(""))
fs_protein <- translate(fs_seq)
subseq(fs_protein, start = 532, end = start(matchPattern("*", fs_protein)))
