## ===================================================================
## Asthma GEO Dataset — DEG, Pathway Enrichment & PPI Network Analysis
## Dataset: GSE152004 (GEO, Illumina platform, 393 samples)
##          257 asthma / 136 control
## ===================================================================
##
## NOTE ON REPRODUCIBILITY:
## The differential expression step (GEO2R group assignment + DEG
## export) was performed manually through NCBI's web-based GEO2R tool,
## not scripted in R. If you want a fully scripted, one-click pipeline,
## replace Step 1 below with a GEOquery + limma script that pulls
## GSE152004 directly and reproduces the same group comparison.
## The steps from DEG list onward (Step 2 onward) are exactly what
## was run in this project.
## ===================================================================


## ---- Step 0: Install required packages (run once) -----------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
                        "enrichplot", "DOSE", "pathview"))


## ---- Step 1: DEG list (from GEO2R) ---------------------------------
## Exported from NCBI GEO2R after defining "Asthma" (n=257) vs
## "Control" (n=136) groups for GSE152004.
## DEG cutoff applied: adj.P.Val < 0.05  AND  |log2FoldChange| > 1
## -> 33 genes passed this threshold.

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

deg <- read.csv("asthma_deg_list_full.csv")
head(deg)
str(deg)

## Expected columns based on this project's GEO2R export:
## GeneID | padj | pvalue | lfcSE | stat | log2FoldChange | baseMean
## | Symbol | Description | Significance


## ---- Step 2: Convert gene symbols to Entrez IDs --------------------
gene_list <- deg$Symbol
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
length(gene_list)   # 33 in this project

gene_conversion <- bitr(gene_list,
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)

nrow(gene_conversion)   # 33/33 mapped successfully in this project
entrez_genes <- gene_conversion$ENTREZID


## ---- Step 3: KEGG pathway enrichment --------------------------------
kegg_result <- enrichKEGG(gene = entrez_genes,
                           organism = "hsa",
                           pvalueCutoff = 0.05)

nrow(as.data.frame(kegg_result))   # 9 significant pathways
head(as.data.frame(kegg_result))

png("figures/kegg_barplot.png", width = 1000, height = 800)
barplot(kegg_result, showCategory = 9)
dev.off()


## ---- Step 4: GO Biological Process enrichment -----------------------
go_result <- enrichGO(gene = entrez_genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

nrow(as.data.frame(go_result))   # 165 significant terms

## Reduce redundancy (near-duplicate GO terms)
go_simplified <- simplify(go_result, cutoff = 0.7,
                           by = "p.adjust", select_fun = min)
nrow(as.data.frame(go_simplified))   # 53 non-redundant terms

png("figures/go_dotplot.png", width = 1000, height = 800)
dotplot(go_simplified, showCategory = 15)
dev.off()

## Optional: relationship map between enriched terms
png("figures/go_emapplot.png", width = 1000, height = 800)
emapplot(pairwise_termsim(go_result), showCategory = 15)
dev.off()


## ---- Step 5: Pathway-level visualization (fold-change overlay) -----
library(pathview)

fc_values <- deg$log2FoldChange
names(fc_values) <- gene_conversion$ENTREZID[match(deg$Symbol,
                                                     gene_conversion$SYMBOL)]

## Chemokine signaling pathway - top KEGG hit in this project
pathview(gene.data = fc_values,
          pathway.id = "hsa04062",
          species = "hsa")
## Produces: hsa04062.pathview.png in the working directory


## ---- Step 6: Export gene list for STRING ---------------------------
## Take deg$Symbol (33 genes) to string-db.org:
##   STRING -> "Multiple proteins" -> paste gene list ->
##   organism: Homo sapiens -> Search
## Export interactions as TSV ("string_interactions.tsv"),
## included in figures/ for reference.


## ---- Step 7: Hub gene identification (degree centrality) ------------
## Network was imported into Cytoscape from string_interactions.tsv
## (File -> Import -> Network from File), analyzed via
## Tools -> NetworkAnalyzer -> Analyze Network (undirected),
## and styled by degree (node size + color) with a force-directed layout.
##
## Equivalent degree calculation done in R for reference/reproducibility:

interactions <- read.delim("figures/string_interactions.tsv",
                            stringsAsFactors = FALSE)

library(dplyr)
degree_table <- c(interactions$node1, interactions$node2) |>
  table() |>
  sort(decreasing = TRUE)

head(degree_table, 10)
## Top hub in this project: IFNG (degree = 12), followed by
## CXCL10 and CD163 (degree = 9 each)

## ===================================================================
## End of pipeline. See Asthma_DEG_Report.pdf for full write-up,
## discussion, and biological interpretation.
## ===================================================================
