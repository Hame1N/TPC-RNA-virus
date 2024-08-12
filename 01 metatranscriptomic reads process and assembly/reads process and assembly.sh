
### rawdata link  PRJNA1105542

### Trimmomatic v0.39

## quality control of raw reads by Trimmonatic

trimmomatic PE -threads 4 "${sample}".R1.fq.gz "${sample}"R2.fq.gz "${sample}"R1_trimmed.fq.gz "${sample}"R1_trimmed_U.fq.gz "${sample}"R2_trimmed.fq.gz "${sample}"R2_trimmed_U.fq.gz ILLUMINACLIP:/data01nfs/user/liupf/common_files/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

### sortmeRNA v4.3.6

## removal of rRNA reads by sortmeRNA

sortmerna --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/rfam-5.8s-database-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/rfam-5s-database-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-arc-16s-id95.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-arc-23s-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-bac-16s-id90.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-bac-23s-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-euk-18s-id95.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-euk-28s-id98.fasta --fastx -a 20 -v --log --reads /datanode02/zhangzh/MetaT/02Cleandata/Trimgalore_out/${sample}_clean_R1.fq.gz  --reads /datanode02/zhangzh/MetaT/02Cleandata/Trimgalore_out/${sample}_clean_R2.fq.gz --aligned ${sample}.align --other ${sample}.unalign --paired_in --out2 --workdir /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${sample}  


### megahit  v1.29

## assmbly of unalign reads of sortmeRNA

megahit  -t 8 -o 03Assembly/megahit/${sample} -1 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${sample}/${sample}.unalign_fwd.fq.gz   -2 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${sample}/${sample}.unalign_fwd.fq.gz 

##pcik contigs over 1k length
source /apps/source/seqkit-2.2.0.sh
seqkit seq -m 1000 ${sample}_final.contigs.fa > ${sample}.contigs1k.fa



