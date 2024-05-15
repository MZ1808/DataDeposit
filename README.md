# Analysis of hMSC-derived chondrocytes in response to stimulation with SARS-CoV-2 spike proteins

## Description of Project
Previous studies have detected the presence of SARS-CoV-2 viral proteins in bronchial cartilage chondrocytes. We hypothesised that the leakage of viral proteins to local joint tissue was due to virus-induced endothelial dysfunction. Studies have also shown upregulation of endothelin-1 (ET-1), the most potent vasoconstrictor, in COVID patients. We are investigating the direct effect of SARS-CoV-2 spike protein (SP) and the host response to chondrocytes. Human mesenchymal stem cells (hMSCs)-differentiated chondrocytes were treated with either full-length SARS-CoV-2 spike protein (SP) only or a combination of spike protein, neutralising antibody to S1 and endothelin-1 (SAE) to mimic viral insult and host response respectively. RNA sequencing was performed to compared the change in transcriptome in control, SP and SAE. All samples were processed in the same batch. Default quality control parameters were used. 

## Groupings
1. Control - human mesenchymal stem cells-differentiated chondrocytes (hMSCs)
2. ⁠SP - hMSCs + 10ug full-length SARS-CoV-2 spike protein (SP), 24h
3. ⁠SAE - hMSCs + 10ug SP + 10ug neutralising antibody to S1 + 100nM ET-1, 24h
Three replicates for each group. Samples are collected for bulk RNA seq.

## Bioinformatics Analysis
Raw sequencing reads were first filtered for adapter and low-quality sequences, retaining only reads with a length of 40bp or more. Low-quality sequences were identified as reads with over 5% unknown bases (“N”) or reads with more than 50% of bases having a quality value of 11 or lower.
Subsequently, the filtered reads were aligned to the GRCh38 reference genome using STAR 2.7.8 (1). RNA expression levels were quantified by RSEM 1.2.31 (2). The default EBSeq 1.18.0 (3) was used to conduct differential expression analysis. The raw sequencing data was deposited in the Gene Expression Omnibus (GSE267009).

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

## References
1.​A. Dobin et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21 (2013).
2.​B. Li, C. N. Dewey, RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. Bmc Bioinformatics 12,  (2011).
3.​N. Leng et al., EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments (vol 29, pg 1035, 2013). Bioinformatics 29, 2073-2073 (2013).
