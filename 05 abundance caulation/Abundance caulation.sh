
### abundance caulation for vContigs, vOTUs, cox and rps3 
## Bowtie2 v2.2.5 and Samtools v1.15.1

seq= {final_rdrp_contig_emboss_fliter.fna, final_rdrp_contig_emboss_fliter_vOTUs.fna, cox.fna, rps3.fna}


bowtie2-build ${seq}  bt2/${seq}
bowtie2 -p 6  -x bt2/${seq}  -1 ${sample}unalign_fwd.fq.gz  -2 ${sample}unalign_rev.fq.gz  | samtools view -bS -@ 20 -q 30 | samtools sort -@ 20 >bt2bam/${seq}_${sample}.sorted.bam 

### CoverM v0.6.1
coverm contig  --bam-files ${seq}_${sample}.sorted.bam   -t 10  --methods trimmed_mean  -o process/${seq}_${sample}_trimmed_mean.tsv  
coverm contig  --bam-files ${seq}_${sample}.sorted.bam  -t 10  --methods count  -o process/${seq}_${sample}_count.tsv   
coverm contig  --bam-files ${seq}_${sample}.sorted.bam  -t 10  --methods tpm -o process/${seq}_${sample}_tpm.tsv   
coverm contig  --bam-files ${seq}_${sample}.sorted.bam  -t 10  --methods covered_fraction -o process/${seq}_${sample}_coverage.tsv   

### abundance caulation for host MAGs
## 15,640 MAGs were collected from the glacier, lake sediment, wetland, and upland of the TP cryosphere
## dRep v3.4.5

dRep dereplicate MAG_99_final   -pa 0.99 -cm larger -p 10 -sa 0.95 -comp 50 -con 10 -g Sdrep95MAG/*.fa

cat MAG_99_final
## BBmap v39.01 and Samtools v1.15.1

bbmap.sh ref=All99MAG.fa in=${sample}.unalign_fwd.fq.gz in2=${sample}.unalign_rev.fq.gz xmtag=t ambiguous=random outm=bam_bbmap/${sample}.bam threads=80  -Xmx500g

samtools sort -@ 20 -q 30  ${sample}.bam >${sample}.sorted.bam
### CoverM v0.6.1
coverm contig  --bam-files ${sample}.sorted.bam  -t 10  --methods tpm -o process/${sample}_tpm.tsv   


