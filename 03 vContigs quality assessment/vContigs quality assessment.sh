
### checkv v1.01
## checkv for quality assessment 

checkv end_to_end  final_rdrp_contig.fna   -t 10 : -d /data02nfs/biodatabase/checkv_db/checkv-db-v1.5

### emboss v6.5.7.0
## assessment of RdRp domain completeness

# tanslate to amino acid (AA) sequences using different codon table and frames
for i in `seq 0 1 23` ;do  transeq final_rdrp_contig.fa  final_rdrp_contig.fa.all6_t${i}.faa -frame=6 -table ${i}  ;done 

# translated AA sequences were searched against the 65 RdRp HMM profiles from Wolf et al.

hmmsearch --cpu 20  -E 0.01 --domE 0.01 --incdomE 0.01 --incE 0.01 -o final_rdrp_contig.fa.all6_t${i}.txt --domtblout final_rdrp_contig.fa.all6_t${i}_dom.txt  --tblout final_rdrp_contig.fa.all6_t${i}_tb.txt /datanode02/zhangzh/database/rdrp/65RdRP.hmm final_rdrp_contig.fa.all6_t${i}.faa

# select optime codon table and frame
for i in *dom.txt ; do  grep -v "^#" ${i} |awk -v n="$i" '{print $1,$14+0,$3,$4,$6,$18,$19,n}' OFS="\t" |sort -t$'\t' -s -k1,1  -k2,2nr    | awk -F"\t" '!a[$1]++'  OFS="\t" ; done >> All_emboss.tsv
awk -F"\t" '$9=$7-$6+1{split($1,array,"_");$10=array[length(array)];sub("_[0-9]{1,2}$","",$1);print  }' All_emboss.tsv | sort   -s -k1,1  -k2,2nr -k9,9nr -k8,8n |awk  '!a[$1]++'  OFS="\t"   >All_emboss_result.tsv


### Bowtie2 v2.2.5 and Samtools v1.15.1
## mapping metagenomic reads to putative RNA vContigs and AMG for verificaiton

bowtie2-build final_rdrp_contig.fna    bt2/final_rdrp_contig
bowtie2 -p 6  -x bt2/final_rdrp_contig  -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz | samtools view -q 30 -F 0x08 -b -f 0x2 | samtools sort -@ 5 >bt2bam/${sample}.bam 

bowtie2-build final_rdrp_contig_amg.fna    bt2/final_rdrp_contig_amg
bowtie2 -p 6  -x bt2/final_rdrp_contig_amg  -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz | samtools view -q 30 -F 0x08 -b -f 0x2 | samtools sort -@ 5 >bt2bam/${sample}.bam 

## CoverM v0.6.1
# filtering contigs with coverage >=75% and depth <= 1x
coverm contig  --bam-files ${sample}.sorted.bam  -t 10  --methods covered_fraction -o process/${sample}_coverage.tsv 
coverm contig  --bam-files ${sample}.sorted.bam   -t 10  --methods trimmed_mean  -o process/${sample}_trimmed_mean.tsv    


