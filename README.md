# Asthma Gene Expression Analysis (GSE152004)

Differential gene expression, pathway enrichment, and protein-protein
interaction network analysis of a public asthma transcriptomics dataset.

## Dataset
- **Source:** NCBI GEO, accession [GSE152004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004)
- **Platform:** Illumina
- **Samples:** 393 total — 257 asthma, 136 healthy control

## Workflow
1. **Differential expression** — GEO2R (limma), filtered to adj. p-value < 0.05 and |log2FC| > 1 → **33 DEGs**
2. **Functional enrichment** — `clusterProfiler` (R/Bioconductor): KEGG (9 significant pathways) and GO Biological Process (165 terms, simplified to 53)
3. **Pathway visualization** — `pathview`, fold-change mapped onto the KEGG Chemokine Signaling Pathway (hsa04062)
4. **PPI network** — STRING-db → Cytoscape, node degree analysis to identify hub genes

## Key finding
**IFNG (interferon-gamma)** was the dominant network hub (degree = 12), consistent
with a chemokine/leukocyte-recruitment signature (CXCL10, CXCL11, CXCL13, CCL20)
driving the enriched KEGG/GO pathways — a pattern consistent with known asthma
immunopathology. Full discussion in the report.

## Repo contents

| File | Description |
|---|---|
| `Asthma_DEG_Report.pdf` | Full write-up: methods, results, figures, discussion, limitations |
| `Asthma_DEG_Report_project_1.docx` | Same report, editable Word format |
| `analysis_pipeline.R` | Full R pipeline: gene ID conversion → KEGG/GO enrichment → pathview → hub gene calc |
| `asthma_deg_list_full.csv` | Full DEG table — 33 genes with padj, log2FC, baseMean, etc. |
| `deg_gene_entrez_mapping.csv` | All 33 DEGs with mapped Entrez IDs (used for enrichment) |
| `figures/kegg_barplot.png` | Top enriched KEGG pathways |
| `figures/go_dotplot.png` | Top 15 non-redundant enriched GO terms |
| `figures/pathway_diagram_hsa04062.png` | KEGG chemokine pathway with fold-change overlay |
| `figures/ppi_network_cytoscape.png` | STRING/Cytoscape interaction network, styled by hub degree |
| `figures/string_interactions.tsv` | Raw STRING interaction edge list |

## Tools used
R (clusterProfiler, org.Hs.eg.db, enrichplot, pathview) · GEO2R · STRING-db · Cytoscape
