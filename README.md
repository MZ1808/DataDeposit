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
