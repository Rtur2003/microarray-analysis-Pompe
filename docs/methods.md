# Methods

This analysis reproduces the IOPD microarray workflow described in the manuscript.
Data were retrieved from GEO accession [GSE38680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38680).
Affymetrix U133 Plus 2.0 CEL files are RMA normalised using the **affy** package.
Differential expression is assessed with **limma** and Benjamini-Hochberg correction.
Functional enrichment uses **clusterProfiler** with GO and KEGG databases (`org.Hs.eg.db` and KEGG human).
miRNA targets can be queried via `multiMiR` (optional and requires network access).
No personal health information is included in this repository.
