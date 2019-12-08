#trim low-quality-reads and carrying-adapter-reads
trim_galore --illumina --paired Read1.fastq Read2.fastq

#mapping trimmed-reads onto genome (mm9)
tophat  -o output mm9_bowtie2Index trimmed_Read1.fastq trimmed_Read2.fastq
cd ./output
cufflinks accepted_hits.bam -G mm9_genes.gtf

#analyze differential expression genes
cuffdiff -o Cuffdiff_output  mm9_genes.gtf Control_samples Treatment_samples

