#trim low-quality-reads and carrying-adapter-reads
trim_galore --illumina --paired Read1.fastq Read2.fastq

#mapping trimmed-reads onto genome (mm9)
fq=sample_name
bowtie2  -k 1 -X 2000 -x mm9_Bowtie2Index -1 trimmed-reads.fastq  -2 trimmed-reads.fastq -S ${fq}tmp.sam
samtools view -Sbq 1 ${fq}tmp.sam > ${fq}tmp.bam
rm -rf ${fq}tmp.sam

#remove PCR-duplicated reads for mapped reads
samtools rmdup -sS ${fq}tmp.bam ${fq}tmp_rmdup.bam
rm -rf ${fq}tmp.bam 

#sort/index bam file
samtools sort -o ${fq}_rmdup_sorted.bam ${fq}tmp_rmdup.bam
samtools index ${fq}_rmdup_sorted.bam
rm -rf ${fq}tmp_rmdup.bam

#convert bam into bed file
bedtools bamToBed -i ${fq}_rmdup_sorted.bam > ${fq}.bed

#call narrow Peaks for TET1 ChIP-seq sample
macs2 callpeak -t ChIP.bed -c Input.bed -g mm -q 0.05 --name=sample-Name

#call broad Peaks for NONO ChIP-seq sample
macs2 callpeak -t ChIP.bed -c Input.bed -g mm -q 0.05 --name=sample-Name --broad
