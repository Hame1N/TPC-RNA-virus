

## Bowtie2 v2.2.5 and Samtools v1.15.1

### abundance caulation for vContigs, vOTUs, cox and rps3 
bowtie2-build final_rdrp_contig_emboss_fliter.fna  bt2/final_rdrp_contig_emboss_fliter
bowtie2 -p 6  -x bt2/final_rdrp_contig_emboss_fliter   -1 ${sample}unalign_fwd.fq.gz  -2 ${sample}unalign_rev.fq.gz  | samtools view -bS -@ 20 -q 30 | samtools sort -@ 20 >bt2bam/${sample}.sorted.bam 

### CoverM v0.6.1
coverm contig  --bam-files ${sample}.sorted.bam   -t 10  --methods trimmed_mean  -o process/${sample}_trimmed_mean.tsv  
coverm contig  --bam-files ${sample}.sorted.bam  -t 10  --methods count  -o process/${sample}_count.tsv   
coverm contig  --bam-files ${sample}.sorted.bam  -t 10  --methods tpm -o process/${sample}_tpm.tsv   
coverm contig  --bam-files ${sample}.sorted.bam  -t 10  --methods covered_fraction -o process/${sample}_coverage.tsv   

### abundance caulation for host MAGs
## 15,640 MAGs were collected from the glacier, lake sediment, wetland, and upland of the TP cryosphere
## dRep v3.4.5

dRep dereplicate MAG_99_final   -pa 0.99 -cm larger -p 10 -sa 0.95 -comp 50 -con 10 -g Sdrep95MAG/*.fa

BBmap v39.01   Samtools v1.15.1.

