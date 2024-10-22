## Nono deficiency compromises TET1 chromatin association and impedes neuronal differentiation of mouse embryonic stem cells

Wenjing Li<sup>1,2,4</sup>, Violetta Karwacki-Neisius<sup>3,4,\*</sup>, Chun Ma<sup>1</sup>, Li Tan<sup>1</sup>, Yang Shi<sup>3</sup>, Feizhen Wu<sup>1,\*</sup>, Yujiang Geno Shi<sup>2,*</sup>

1 Laboratory of Epigenetics, Institutes of Biomedical Sciences, Fudan University, Shanghai, 200032, China
2 Endocrinology Division, Brigham and Women’s Hospital, Harvard Medical School, 221 Longwood Avenue, Boston, MA, 02115, USA
3 Division of Newborn Medicine and Program in Epigenetics, Children’s Hospital, Harvard Medical School, 300 Longwood Avenue, Boston, MA, 02115, USA

4 The authors wish it to be known that, in their opinion, the first 2 authors should be regarded as joint First Authors.

\* Correspondence should be addressed to Yujiang Geno Shi. Tel: 1-(617) 525-8097; Fax: (617) 582-6193; Email: yujiang_shi@hms.harvard.edu
\* Correspondence may also be addressed to Feizhen Wu (wufz@fudan.edu.cn) and Violetta Karwacki-Neisius (violetta.karwacki-neisius@childrens.harvard.edu).


### ABSTRACT
NONO is a DNA/RNA-binding protein, which plays a critical regulatory role during cell stage transitions of mouse embryonic stem cells (mESCs). However, its function in neuronal lineage commitment and the molecular mechanisms of its action in such processes are largely unknown. Here we report that NONO plays a key role during neuronal differentiation of mESCs. Nono deletion impedes neuronal lineage commitment largely due to a failure of up-regulation of specific genes critical for neuronal differentiation. Many of the NONO regulated genes are also DNA demethylase TET1 targeted genes. Importantly, re-introducing wild type NONO to the Nono KO cells, not only restores the normal expression of the majority of NONO/TET1 coregulated genes but also rescues the defective neuronal differentiation of Nono-deficient mESCs. Mechanistically, our data shows that NONO directly interacts with TET1 via its DNA binding domain and recruits TET1 to genomic loci to regulate 5-hydroxymethylcytosine levels. Nono deletion leads to a drastic dissociation of TET1 from chromatin and dysregulation of DNA hydroxymethylation of neuronal genes. Taken together, our findings reveal a key role and an epigenetic mechanism of action of NONO in regulation of TET1-targeted neuronal genes, offering new functional and mechanistic understanding of NONO in stem cell functions, lineage commitment and specification.

### Repository Description
The repository is to deposit the commands, scripts, and middle data that were used to analyze the mRNA-seq and ChIP-seq dataset.

### Commands
|  Item   | Data-type |Description  |
|  ----   | ----  | ----  |
|ChIP-seq_process_commands.sh|linux commands|trimming and mapping reads, call-peaks|  
|mRNA-seq_process_commands.sh|linux commands|mapping, analyze differential expression genes|

#### Scripts

|  Item   | Data-type |Description  |
|  ----   | ----  | ----  |
|Figure2A_script.R |R script|Draw Figure 2 A|
|Figure2B_script.R |R script|Draw Figure 2 B|
|Figure2C_script.R |R script|Draw Figure 2 C|
|Figure2A_E_F_G_H_I_J.R|R script|Draw Figure 2 A, E, F, G, H, I, and J|
|Figure3_F_WT_Nono_correlation.R |R script|Draw Figure 3 F|
|Figure5_A_B.R |R script|Draw Figure5 A and B|
|Figure5_C_D_E.R |R script|Draw Figure 5 C, D, and E|
|Supplementary_Figure1B.R|R script|Draw Supplementary Figure 1 B |
|Supplementary_Figure2A_B_C.R |R script|Draw Supplementary Figure2 A,B, and C|
|Supplementary_Figure3A_D.R |R script|Draw Supplementary Figure3 A and D|
|Supplementary_Figure3B_E.R |R script|Draw Supplementary Figure3 B and E|
|Supplementary_Figure3C_F.R |R script|Draw Supplementary Figure3 C and F|
|Supplementary_Figure6.R|R script|Draw Supplementary Figure 6 |
|FPKM2TPM.R|R script|Convert FPKM generated by cufflinks to TPM|
|signalplot|python script|Draw Fig3D,Fig4D-G,SupFig4BCEF, and SupFig5CD|

#### Gene Expression Dataset

|Item|Data-type|Description|
|  ----   |   ----  |     ----  |
|WT_D12_vs_D0.RData|R data|WT Day12/Day0|
|NonoKO_D12_vs_D0.RData|R data|NonoKO Day12/Day0|
|RE_D12_vs_D0.RData|R data|NonoKO+WT Day12/Day0|
|WT_D3_vs_D0.RData|R data|WT Day3/Day0|
|NonoKO_D3_vs_D0.RData|R data|NonoKO Day3/Day0|
|RE_D3_vs_D0.RData|R data|NonoKO+WT Day3/Day0|
|WT_D6_vs_D3.RData|R data|WT Day6/Day3|
|NonoKO_D6_vs_D3.RData|R data|NonoKO Day6/Day3|
|RE_D6_vs_D3.RData|R data|NonoKO+WT Day6/Day3|
|WT_D12_vs_D6.RData|R data|WT Day12/Day6|
|NonoKO_D12_vs_D6.RData|R data|NonoKO Day12/Day6|
|RE_D12_vs_D6.RData|R data|NonoKO+WT Day12/Day6|
|3NonoKO_D0_vs_3WT_D0.RData|R data|NonoKO/WT in Day0|

Note: Each R data is a list object, which contains items: exp, diffgenes, Enrich, UpGeneEnrich, DwGeneEnrich fpkmcutoff, foldchange, and pvalue. The exp is a matrix from gene_exp.diff generated by a program,cuffdiff, in cufflinks package. The 3rd column and 4th column in the exp are the equivalent of average of Fpkm of replicated WT samples and replicated KO samples, respectively. the diffgenes is a differential expression gene list at the cutoff of the foldchange, the fpkmcutoff, and the pvalue. The Enrich, UpGeneEnrich, and DwGeneEnrich are GO and KEGG enrichment analyses for all diffgenes, up-regulated, and down-regulated genes, respectively.


#### Other Dataset
|Item|Data-type|Description|
|  ----   |   ----  |     ----  |
|set1_Expression_TPM.xls|Excel file|set1 gene expression abundance|  
|set2_Expression_TPM.xls|Excel file|set2 gene expression abundance|
|WT_NonoKO_hMeDIP_at_Promoter2K.txt|Text| 5hmC density in promoters of WT and NonoKO cells|  
|WT_NonoKO_Tet1_at_Promoter2K.txt|Text| Tet1 density in promoters of WT and NonoKO cells| 

If you have any questions or comments on these scripts and data, please let us know. Thank You!
