# Analysis of hMSC-derived chondrocytes in response to stimulation with SARS-CoV-2 spike proteins

## Description of Project
Previous studies have detected the presence of SARS-CoV-2 viral proteins in bronchial cartilage chondrocytes. We hypothesised that the leakage of viral proteins to local joint tissue was due to virus-induced endothelial dysfunction. Studies have also shown upregulation of endothelin-1 (ET-1), the most potent vasoconstrictor, in COVID patients. We are investigating the direct effect of SARS-CoV-2 spike protein (SP) and the host response to chondrocytes. 

Human mesenchymal stem cells (hMSCs)-differentiated chondrocytes were treated with either full-length SARS-CoV-2 spike protein (SP) only or a combination of spike protein, neutralising antibody to S1 and endothelin-1 (SAE) to mimic viral insult and host response respectively. RNA sequencing was performed to compared the change in transcriptome in control, SP, and SAE. All samples were processed in the same batch. Default quality control parameters were used. 

## Groupings
1. Control - human mesenchymal stem cells-differentiated chondrocytes (hMSCs)
2. ⁠SP - hMSCs + 5 ug/mL full-length SARS-CoV-2 spike protein (SP), 24h
3. ⁠SAE - hMSCs + 5 ug/mL SP + 5 ug/mL neutralising antibody to S1 + 100nM ET-1, 24h
Three replicates for each group. Samples are collected for bulk RNA seq.

## Bioinformatics Analysis
Raw sequencing reads were first filtered for adapter and low-quality sequences, retaining only reads with a length of 40bp or more. Low-quality sequences were identified as reads with over 5% unknown bases (“N”) or reads with more than 50% of bases having a quality value of 11 or lower. Subsequently, the filtered reads were aligned to the GRCh38 reference genome using STAR 2.7.8 (1). RNA expression levels were quantified by RSEM 1.2.31 (2). The raw sequencing data was deposited in the Gene Expression Omnibus (GSE267009).

Two independent researchers performed differential gene expression analysis, henceforth referred to as 'Original Analysis' and 'New Analysis'. 

In Original Analysis, the default EBSeq 1.18.0 (3) was used to conduct differential expression analysis. DEGs with FDR less than 0.05 were included in KEGG enrichment analysis using GAGE 2.52.0 (4). Compared with Control, the top 20 upregulated KEGG pathways in the SP group, and all upregulated pathways with p-value less than 0.05 in the SAE group were plotted with ggplot2 3.5.0. Rich ratio is calculated by the number of DEGs in a given pathway divided by the total number of genes in this pathway. Details of KEGG enrichment analysis can be found in 'source codes_kegg.Rmd'. Results for each comparison and the generated figures are included in 'Original Analysis Results.zip'.

Genes involved in cellular senescence and cytokine-cytokine interaction were studied. Selected DEGs of interest in these pathways were plotted with gplots 3.1.3.1. Normalized gene expression using TPM was presented in heatmaps, with adjusted p-value. Colours were scaled with respect to each row to emphasize the inter-sample difference. Row Z-Score was produced by subtracting the mean value and then divided by the standard deviation of each row. Details of KEGG enrichment analysis can be found in 'source codes_heatmap.Rmd'. The generated figures are included in 'Original Analysis Results.zip'.

In New Analysis, differential expression analysis was performed using limma 3.58.1 in R. Comparisons were made between the Ctrl and SP groups, Ctrl and SAE groups, and SP and SAE groups. Differentially expressed genes and transcripts were identified based on a threshold of p-value < 0.05. KEGG pathway and GO functional enrichment analyses of the differentially expressed genes (p< 0.05, Fold Change >1) were conducted using ClusterProfiler 4.0.5 in R70. Enrichment significance was assessed using the hypergeometric test, and p-values were adjusted using the Benjamini-Hochberg (BH) method to obtain q-values. Pathways and GO terms with q<0.05 were considered significantly enriched. Heatmaps depicting the expression levels of genes of interest across samples were generated using ggplot2 3.3.5 in R, with unadjusted p-values. 

## RNAseq Data
| Sample Information | SRA Accession No. | GEO Accession No. |
| ----------- | ----------- |----------- |
| Control, replicate 1| SRR28789040 | GSM8258185 |
| Control, replicate 2| SRR28789039 | GSM8258186 |
| Control, replicate 3| SRR28789038 | GSM8258187 |
| Treated with 5 ug/mL full-length SARS-CoV-2 spike protein, replicate 1 | SRR28789037 | GSM8258188 |
| Treated with 5 ug/mL full-length SARS-CoV-2 spike protein, replicate 2 | SRR28789036 | GSM8258189 |
| Treated with 5 ug/mL full-length SARS-CoV-2 spike protein, replicate 3 | SRR28789035 | GSM8258190 |
| Treated with a combination of SARS-CoV-2 spike protein (5 ug/mL), the neutralising antibody to S1 region (5 ug/mL), and 100nM endothelin-1, replicate 1 | SRR28789034 | GSM8258191 |
| Treated with a combination of SARS-CoV-2 spike protein (5 ug/mL), the neutralising antibody to S1 region (5 ug/mL), and 100nM endothelin-1, replicate 2 | SRR28789033 | GSM8258192 |
| Treated with a combination of SARS-CoV-2 spike protein (5 ug/mL), the neutralising antibody to S1 region (5 ug/mL), and 100nM endothelin-1, replicate 3 | SRR28789032 | GSM8258193 |

## Directory


## References
1. A. Dobin et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21 (2013).
2. B. Li, C. N. Dewey, RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. Bmc Bioinformatics 12,  (2011).
3. N. Leng et al., EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments (vol 29, pg 1035, 2013). Bioinformatics 29, 2073-2073 (2013).
4. W. Luo et al., GAGE: generally applicable gene set enrichment for pathway analysis. Bmc Bioinformatics 10, 161 (2009).
