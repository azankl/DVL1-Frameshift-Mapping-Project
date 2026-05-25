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

# Fetch the sequences
NM_001330311_CDS_record <- entrez_fetch(
    db = "nuccore",
    id = "NM_001330311.2",
    rettype = "fasta_cds_na"
)
writeLines(NM_001330311_CDS_record, "NM_001330311_CDS.fasta")
NM_004421_CDS_record <- entrez_fetch(
    db = "nuccore",
    id = "NM_004421.2",
    rettype = "fasta_cds_na"
)
writeLines(NM_004421_CDS_record, "NM_004421_CDS.fasta")

# Parse the NM_001330311_CDS into a DNAString object for downstream analysis.
NM_001330311_CDS <- readDNAStringSet("NM_001330311_CDS.fasta")
NM_001330311_CDS_seq <- NM_001330311_CDS[[1]]
NM_001330311_CDS_seq

# Parse the NM_004421_CDS into a DNAString object for downstream analysis.
NM_004421_CDS <- readDNAStringSet("NM_004421_CDS.fasta")
NM_004421_CDS_seq <- NM_004421_CDS[[1]]
NM_004421_CDS_seq

# Translate the NM_001330311_CDS
NM_001330311_protein_seq <- translate(NM_001330311_CDS_seq)
NM_001330311_protein_seq

# Translate the NM_004421_CDS
NM_004421_protein_seq <- translate(NM_004421_CDS_seq)
NM_004421_protein_seq

# Introduce a 1 bp deletion at position 1593 (1-based index) in NM_001330311
del1593_seq <- replaceAt(
    NM_001330311_CDS_seq,
    at = IRanges(1593, 1593),
    value = DNAStringSet("")
)
del1593_NM_001330311_protein_seq <- translate(del1593_seq)
del1593_new_protein_seq <- subseq(
    del1593_NM_001330311_protein_seq,
    start = 532,
    end = start(matchPattern("*", del1593_NM_001330311_protein_seq))[1] - 1 # the [1] is needed because matchPattern returns a list of matches, and we want the first one, the -1 is needed to get the position of the last amino acid before the stop codon
)

# Introduce a 1 bp deletion at position 1508 (1-based index) in NM_004421
del1508_seq <- replaceAt(
    NM_004421_CDS_seq,
    at = IRanges(1508, 1508),
    value = DNAStringSet("")
)
del1508_NM_004421_protein_seq <- translate(del1508_seq)
del1508_new_protein_seq <- subseq(
    del1508_NM_004421_protein_seq,
    start = 503,
    end = start(matchPattern("*", del1508_NM_004421_protein_seq))[1] - 1
)
# library(pwalign)

# create a pairwise alignment of the two protein sequences
# aln <- pairwiseAlignment(
#     del1593_NM_001330311_protein_seq,
#     del1508_NM_004421_protein_seq,
#     substitutionMatrix = "BLOSUM62",
#     type = "global"
# )
# aln

# create a multiple sequence alignment of the two protein sequences
library(msa)
aln <- msa(
    AAStringSet(c(
        as.character(del1593_new_protein_seq),
        as.character(del1508_new_protein_seq)
    )),
    method = "ClustalW"
)
aln
msaPrettyPrint(
    aln,
    output = "pdf",
    shadingMode = "identical",
    shadingColors = "blues",
    showNames = "none"
)
