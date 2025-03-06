# CRISPR-Cas12a
Search the best gRNA for your DNA sequence.
Overview
This Python script identifies the top 5 gRNAs for the FGFR3 gene (NM_000142) using CRISPR/Cas12a, based on criteria like GC content, PAM sequence (TTTN), and avoidance of poly-nucleotide repeats. It fetches the FGFR3 mRNA sequence from NCBI, analyzes exonic regions (proxied via mRNA), and outputs results in a table format.


The script automatically fetches the FGFR3 mRNA sequence (NM_000142) from NCBI.
Enter your email for NCBI access (required for API usage).

5 best gRNAs for FGFR3 (NM_000142):  
gRNA                  PAM   Position     GC_Content  
GCACTACGCCGATGTCAACG  TTTT  123-142      50.0  
CAGGTACGTGCGCGACGTGA  TTTC  456-475      50.0  
... (additional rows)  
Notes
PAM Sequence : Targets TTTN (Cas12a-compatible).
GC Content : Filters sequences between 40–70%.
Exon Support : Uses mRNA (NM_000142) as a proxy for exonic regions. Modify the script to include genomic coordinates from Ensembl for precise exon filtering.
Off-Target Analysis : Basic checks included; use tools like Cas-OFFinder for advanced validation.

By MrEgorin. 
https://mregorin.github.io/
