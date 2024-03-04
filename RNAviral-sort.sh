
### md5检验
 md5sum -c MD5.txt
##软连接文件

\ls /data02nfs/Project/Metatranscriptomics/*/*gz |sed 's/^/ln -s /g' |sh

### 质控
#软件
source /datanode02/zhangzh/.apps/trimgalore.sh

time trim_galore --paired --quality 30 --length 100 --max_n 5 -j 8 --fastqc  -o 02Cleandata/Trimgalore_out 01Rawdata/${i}.R1.fq.gz 01Rawdata/${i}.R2.fq.gz 



### sortmeRNA v4.3.6
#软件
source /datanode02/zhangzh/.apps/sortmeRNA.sh

sortmerna --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/rfam-5.8s-database-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/rfam-5s-database-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-arc-16s-id95.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-arc-23s-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-bac-16s-id90.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-bac-23s-id98.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-euk-18s-id95.fasta --ref /data01nfs/user/liupf/software_lpf/sortmeRNA_db/silva-euk-28s-id98.fasta --fastx -a 20 -v --log --reads /datanode02/zhangzh/MetaT/02Cleandata/Trimgalore_out/${i}_clean_R1.fq.gz  --reads /datanode02/zhangzh/MetaT/02Cleandata/Trimgalore_out/${i}_clean_R2.fq.gz --aligned ${i}.align --other ${i}.unalign --paired_in --out2 --workdir /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${i}  




### megahit
#软件
source /apps/source/megahit-1.2.9.sh

megahit  -t 8 -o 03Assembly/megahit/${i} -1 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${i}/${i}.unalign_fwd.fq.gz   -2 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/${i}/${i}.unalign_fwd.fq.gz 

/datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/gz

#contig修改名称

sed 's/>/>2006QYC1_/g' /data01nfs/user/zhangxf/megahit_2006QYC1/final.contigs.fa -i

mv /data01nfs/user/zhangxf/megahit_2006QYC1/final.contigs.fa /data01nfs/user/zhangxf/megahit_2006QYC1/2006QYC1_final.contigs.fa

#挑选大于1k的contig
source /apps/source/seqkit-2.2.0.sh
seqkit seq -m 1000 2006QYC1_final.contigs.fa > 2006QYC1_final.contigs1k.fa

### prodigal
source /apps/source/prodigal-2.6.3.sh

prodigal -p meta -q -m -i 2006QYC1_final.contigs1k.fa -a 2006QYC1_final.contigs1k_prodigal.faa -d 2006QYC1_final.contigs1k_prodigal.fna -o 2006QYC1_final.contigs1k_prodigal.gff -f gff

###  rdrp鉴定

## hmmsearch v3.3.2
source /datanode02/zhangzh/.apps/orthofinder.sh

eval=0.01
THREADS=10
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_NeoRdRp.txt /data01nfs/user/qinfsh/NeoRdRp-HMM.v1.1.hmm ${output_dir}/${ini_name}_${out_su}.faa
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_all_virus.txt /data01nfs/user/qinfsh/all_virus_RdRP_profiles.hmm ${output_dir}/${ini_name}_${out_su}.faa
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_rdrpscan.txt /data01nfs/user/liupf/software_lpf/RdRp-scan/Profile_db_and_alignments/ALL_RdRp_scan_profiles.hmm ${output_dir}/${ini_name}_${out_su}.faa

## palmscan
/datanode03/zhangxf/palmscan-main/bin/palmscan --search_pp ${output_dir}/${ini_name}_${out_su}.faa   -rdrp -ppout ${output_dir}/${ini_name}_${out_su}_palmscan_pp.fa -report ${output_dir}/${ini_name}_${out_su}_palmscan_pp.txt -fevout ${output_dir}/${ini_name}_${out_su}_palmscan_pp.fev -hiconf -threads ${THREADS}

### emboss 确认完整度和codon table
#软件
source /datanode02/zhangzh/.apps/emboss.sh

for i in `seq 0 1 23` ;do echo " transeq ../final_rdrp_contig.fa  final_rdrp_contig.fa.all6_t${i}.faa -frame=6 -table ${i} " ;done 

source /datanode02/zhangzh/.apps/orthofinder.sh
hmmsearch --cpu 20  -E 0.01 --domE 0.01 --incdomE 0.01 --incE 0.01 -o hmmsearch/${i}.txt --domtblout hmmsearch/${i}_dom.txt  --tblout hmmsearch/${i}_tb.txt /datanode02/zhangzh/database/rdrp/65RdRP.hmm ${i} 


for i in *dom.txt ; do  grep -v "^#" ${i} |awk -v n="$i" '{print $1,$14+0,$3,$4,$6,$18,$19,n}' OFS="\t" |sort -t$'\t' -s -k1,1  -k2,2nr    | awk -F"\t" '!a[$1]++'  OFS="\t" >deal/${i}_filter; done
cat *filter >../All_emboss.tsv
awk -F"\t" '$9=$7-$6+1{split($1,array,"_");$10=array[length(array)];sub("_[0-9]{1,2}$","",$1);print  }' All_emboss.tsv | sort   -s -k1,1  -k2,2nr -k9,9nr -k8,8n |awk  '!a[$1]++'  OFS="\t"  > ../All_emboss_result.tsv


#不同codon table的核酸序列
source /apps/source/prodigal-2.6.3.sh
for i in `seq 0 1 33` ;do echo " prodigal  -q -m -i final2_rdrp_contig_emboss_fliter.fna  -a prodigal/final2_rdrp_contig_emboss_fliter_t${i}.faa -d prodigal/final2_rdrp_contig_emboss_fliter_t${i}.fna -o prodigal/final2_rdrp_contig_emboss_fliter_t${i}.gff -f gff -g ${i} " ;done 


### 物种鉴定


# rdrp level

source  /datanode02/zhangzh/.apps/orthofinder.sh
#RVMT rdrp database
diamond blastp --db /datanode02/zhangzh/database/RVMT_RdRP/RVMT_DB.dmnd  --query ../../final2_rdrp_contig_emboss_fliter_rdrp.faa --out final2_rdrp_contig_emboss_fliter_rdrp_RVMT.tsv -f 6 qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore  --sensitive -k 1 --max-target-seqs 1   -p 5
#TaraOcean rdrp database
diamond blastp --db /datanode02/zhangzh/database/TaraOcean/TaraOcean_DB.dmnd --query ../../final2_rdrp_contig_emboss_fliter_rdrp.faa --out final2_rdrp_contig_emboss_fliter_rdrp_TO.tsv -f 6 qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore  --sensitive -k 1 --max-target-seqs 1   -p 5

#paldb
source /apps/source/blast-2.5.sh
psiblast -query ${input_pal} -db /datanode02/zhangzh/database/palmdb-main/palmdb_serratus -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out ${output}_palmscan_psiblast_serratusDB.tsv  -num_threads ${THREADS} 
psiblast -query ${input_pal} -db /datanode02/zhangzh/database/palmdb-main/palmdb_name  -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out ${output}_palmscan_psiblast_nameDB.tsv  -num_threads ${THREADS} 


# Vcontig level 

sed 's/>/>RVMT_/g' /datanode02/zhangzh/database/RVMT_contig/RiboV1.6_Contigs.fasta >RiboV1.6_Contigs_rename.fasta
sed 's/>/>REFV_/g' /datanode02/zhangzh/database/refseq_viral230504/RNA_viral.fna >RNA_viral_rename.fna
sed 's/>/>TO_/g' /datanode02/zhangzh/database/TaraOcean/44779_RdRP_contigs.fna >44779_RdRP_contigs_rename.fna
sed 's/>/>TBPCV_/g' ../final2_rdrp_contig_emboss_fliter.fna  > final2_rdrp_contig_emboss_fliter_rename.fna

cat final2_rdrp_contig_emboss_fliter_rename.fna  RiboV1.6_Contigs_rename.fasta RNA_viral_rename.fna  44779_RdRP_contigs_rename.fna >TBPCV_TO_REFV_RVMT_conitg.fna

source /apps/source/blast-2.5.sh
makeblastdb -in TBPCV_TO_REFV_RVMT_conitg.fna  -dbtype nucl -out TBPCV_TO_REFV_RVMT_conitg
blastn -query TBPCV_TO_REFV_RVMT_conitg.fna   -db TBPCV_TO_REFV_RVMT_conitg -outfmt '6 std qlen slen'  -max_target_seqs 10000 -out blast/TBPCV_TO_REFV_RVMT_conitg_ani.tsv -num_threads 20

source /apps/source/checkv-1.0.1.sh
anicalc.py -i blast/TBPCV_TO_REFV_RVMT_conitg_ani.tsv  -o blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv 
aniclust.py --fna TBPCV_TO_REFV_RVMT_conitg.fna --ani blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv  --out blast/TBPCV_TO_REFV_RVMT_conitg_ani_clusters9090.tsv  --min_ani 90 --min_tcov 90 --min_qcov 0
aniclust.py --fna TBPCV_TO_REFV_RVMT_conitg.fna --ani blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv  --out blast/TBPCV_TO_REFV_RVMT_conitg_ani_clusters9075.tsv  --min_ani 90 --min_tcov 75 --min_qcov 0

# uclust & mcl 完整rdrp
cat final2_rdrp_contig_emboss_fliter_compl.faa  /datanode02/zhangzh/database/RVMT_RdRP/RVMT_final_complete_contig_emboss_win_rdrp.faa  /datanode02/zhangzh/database/TaraOcean/RdRp_footprints_Tara_Genbank_Wolf2020_centroids_50_percent_near_complete.faa >TBPCV_TO_RVMT_compl_rdrp.faa

/data01nfs/user/liupf/software_lpf/usearch_v11_dyc  --cluster_fast  TBPCV_TO_RVMT_compl_rdrp.faa  -id 0.50 -sort length  -centroids  TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa  -uc clusters.uc -clusters Final

source /apps/source/blast-2.5.sh
makeblastdb -in TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa -dbtype prot -out TBPCV_TO_RVMT_compl_rdrp_0.5uclust
blastp -query TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa -db TBPCV_TO_RVMT_compl_rdrp_0.5uclust  -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore'   -gapopen 9 -gapextend 1 -word_size 3  -out Pairwise_CVTP_TR_RVMT_rdrp_0.5uclust.tsv  -num_threads 5 

grep "#" -v Pairwise_CVTP_TR_RVMT_rdrp_0.5uclust.tsv |awk -F"\t" '{if($12==0) print $1,$2,"1"; else print $1,$2,$12}' OFS="\t" >MCL_input.tsv

source  /datanode02/zhangzh/.apps/mcl.sh

mcxload -abc MCL_input.tsv --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o MCL_out.mci -write-tab MCL_out.tab

for i in `seq 1.1 0.1 2` ;do echo " mcl MCL_out.mci  -I ${i} -use-tab MCL_out.tab " ;done 

## IMGV

blastn -query ../../final2_rdrp_contig_emboss_fliter.fna -db /datanode02/zhangzh/database/IMGV/IMGVR_all_nucleotides  -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final2_rdrp_contig_emboss_fliter_blastIMGV.tsv -num_threads 20


##Kraken
source /datanode02/zhangzh/.apps/kraken2.sh
for i in `cat nonGlacier.sample` ;do echo " kraken2 --db /datanode02/zhangzh/database/kraken2_database  --threads 10 --report 02Cleandata/reads_analysis/kraken/${i}.report  --report-minimizer-data --output 02Cleandata/reads_analysis/kraken/${i}.output  --paired 02Cleandata/sortmeRNA/gz/${i}.unalign_fwd.fq.gz -r 02Cleandata/sortmeRNA/gz/${i}.unalign_rev.fq.gz  " ;done 

source /datanode02/zhangzh/.apps/palmid.sh

for i in ./kraken/*report; do  n=$(basename $i .report) ;bracken -d /datanode02/zhangzh/database/kraken2_database  -i ${i} -o bracken/${n}.S.bracken -w bracken/${n}.bracken.report -r 150 -l S ;  done

for i in ./kraken/*report; do n=$(basename $i .report) ; python ~/minisoft/KrakenTools/kreport2mpa.py -r bracken/${n}.bracken.report   -o bracken/${n}.tax ;done

for i in bracken/*tax; do n=$(basename $i .tax); awk -F"|" '$NF~/k__/{print $0}' ${i} |sed "1i Kindom\t${n}" >kra_tax/${n}_kindom.tsv ;done &
for i in bracken/*tax; do n=$(basename $i .tax); awk -F"|" '$NF~/s__/{print $0}' ${i} |sed "1i Species\t${n}" >kra_tax/${n}_species.tsv ;done &




### 功能注释

/datanode02/zhangzh/database/CDD2022/All_CDD.hmm 
~/minisoft/interproscan/interproscan-5.62-94.0/data/pfam/35.0/pfam_a.hmm
~/minisoft/interproscan/interproscan-5.62-94.0/data/gene3d/4.3.0/gene3d_main.hmm 
/datanode02/zhangzh/database/ECOD2017/ecodf.hmm
/datanode02/zhangzh/database/LysDB/Merge/AllHMMs.hmm 

for i in  `cat /datanode02/zhangzh/database/db.list`  ;do n=$(basename $i .hmm)  ; echo "/datanode02/zhangzh/minisoft/hmm_to_db.sh 10 ${n}_hmm ../final2_rdrp_contig_emboss_fliter.faa  ${n} ${i}" ;done

/datanode02/zhangzh/minisoft/interproscan/interproscan-5.62-94.0/interproscan.sh -i ../final2_rdrp_contig_emboss_fliter.fna -t n  -appl MobiDBLite,PRINTS,Phobius,TMHMM  -goterms  -pa  -cpu 5

source /datanode2/zhangzh/.apps/metacerberus.sh
metacerberus.py --protein  ../final2_rdrp_contig_emboss_fliter.faa  --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out metacerberus_result_faa --cpus 10


source /datanode2/zhangzh/.apps/hhsuite.sh

hhblits -cpu 10 -i  ${i}  -d /datanode02/zhangzh/database/hhsuite_db/UniRef30/UniRef30_2023_02 -oa3m a3m/${i#*/}.a3m   -n 1 -e 0.001
hhsearch -cpu 10 -i a3m/${i#*/}.a3m  -d /datanode02/zhangzh/database/hhsuite_db/pfam/pfam  -o hh_r/${i#*/}.hhr   -n 1 -e 0.001 


source /datanode2/zhangzh/.apps/orthofinder.sh
diamond blastp --db ~/database/VFDB/VFDB_setB.dmnd  --query ../final2_rdrp_contig_emboss_fliter.faa   --more-sensitive --outfmt 6 sseqid qseqid evalue gapopen pident length qstart qend sstart send qlen slen  evalue bitscore --max-target-seqs 10 --threads 10  --out final2_rdrp_contig_emboss_fliter_vfdb_more_sensitive.tsv

source /apps/source/blast-2.13.sh
blastp  -query ../final2_rdrp_contig_emboss_fliter.faa -db ~/database/CARD/CARD_blastdb  -outfmt '7 std qlen slen'   -max_target_seqs 10  -out final2_rdrp_contig_emboss_fliter_CARD.tsv  -num_threads 10   


### 宿主

#cripser
#minced 从MAG检索crispr序列

source /datanode2/zhangzh/.apps/minced.sh
for i in GGG_4237MAG/* ;do echo "minced $PWD/${i} $PWD/Crispers/${i#*/}.cripsers  $PWD/Crispers/${i#*/}.gff" ;done

#vContig比对crisper库
#spacer序列提取
for i in `ls Crispers/*fa.cripsers  -lSrh |awk '$5+0!=0 {print $NF}'` ;do awk -F"\t" '{print $4}' ${i} |sed '/^$/d' |sed '/-/d'|sed '/Length/d' | sed "s:^:>${i#*/}\n:g" > Spacer/${i#*/} ;done
for i in `fgrep ">" -c Spacer/* |awk -F":" '$2>3{print $1}' ` ;do cat ${i} |seqkit replace  -p ".fa.cripsers" -r "_{nr}" --nr-width 5 >>All_Spacerover3.fa ;done 

#TBP库
makeblastdb -in Allspacer/All_soil_wetland_glacier.fa  -dbtype nucl -out All_soil_wetland_glacierdb

blastn -query  final_vConitg.fna  -db /datanode02/zhangzh/MetaT/011CryVirus/08Host/All_soil_wetland_glacierdb -dust no -word_size 7  -out All_RVMT_again_soil_wetland_glacierdb_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'  -num_threads 10

#iphop库
makeblastdb -in  /data01nfs/user/huangxy/database/iphop_db/Sept_2021_pub_rw/db/All_CRISPR_spacers_nr_clean.fna -dbtype nucl -out All_iphopdb

blastn -query  final_vConitg.fna   -db /datanode02/zhangzh/MetaT/011CryVirus/08Host/All_iphopdb -dust no -word_size 7  -out All_TaraOcean_again_iphop_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10

#

blastn -query  ~/database/RVMT_contig/RiboV1.6_Contigs.fasta   -db ../../08Host/All_soil_wetland_glacierdb -dust no -word_size 7  -out All_RVMT_again_soil_wetland_glacierdb_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'  -num_threads 10

blastn -query  ~/database/TaraOcean/44779_RdRP_contigs.fna   -db ../../08Host/All_soil_wetland_glacierdb -dust no -word_size 7  -out All_TaraOcean_again_soil_wetland_glacierdb_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'  -num_threads 10

blastn -query  ~/database/TaraOcean/44779_RdRP_contigs.fna  -db ../../08Host/All_iphopdb -dust no -word_size 7  -out All_TaraOcean_again_iphop_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10

blastn -query  ~/database/RVMT_contig/RiboV1.6_Contigs.fasta  -db ../../08Host/All_iphopdb -dust no -word_size 7  -out All_RVMT_again_iphop_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10

blastn -query  IMGVR_all_nucleotides.fna   -db ../../08Host/All_soil_wetland_glacierdb -dust no -word_size 7  -out All_IMGVR4_again_soil_wetland_glacierdb_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'  -num_threads 10

blastn -query  IMGVR_all_nucleotides.fna  -db ../../08Host/All_iphopdb -dust no -word_size 7  -out All_IMGVR_again_iphop_host_prediction.txt -max_target_seqs 20  -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10




#EVEs

/datanode02/zhangzh/database/NCBI_database/NT_202307

tblastn -query ../final2_rdrp_contig_emboss_fliter_rdrp_complete.faa -db /datanode02/zhangzh/database/NCBI_database/NT_202307/NT_202307    -out NT_tblastn_r.txt -max_target_seqs 20 -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10  -evalue  0.1

# host abundance

bowtie2-build Nr_All_81contigs_cox1_pro_final_clean_rmdup.fna  bt2/Nr_All_81contigs_cox1_pro_final_clean_rmdup.fna
bowtie2 -p 6  -x bt2/Nr_All_81contigs_cox1_pro_final_clean_rmdup.fna  -1 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/gz/1QY1.unalign_fwd.fq.gz  -2 /datanode02/zhangzh/MetaT/02Cleandata/sortmeRNA/gz/1QY1.unalign_rev.fq.gz | samtools view -bS | samtools sort -@ 5 >bt2sam/1QY1.unalign_fwd.fq.gz.sam.sorted

for i in bt2sam/* ;do echo "coverm contig  --bam-files ${i}  -t 5  --methods trimmed_mean  -o process/${i##*/}_trimmed_mean.tsv  " ;done
for i in bt2sam/* ;do echo "coverm contig  --bam-files ${i}  -t 5  --methods count  -o process/${i##*/}_count.tsv  " ;done 
for i in bt2sam/* ;do echo "coverm contig  --bam-files ${i}  -t 5  --methods tpm -o process/${i##*/}_tpm.tsv  " ;done 
for i in bt2sam/* ;do echo "coverm contig  --bam-files ${i}  -t 5  --methods covered_fraction -o process/${i##*/}_coverage_histogram.tsv  " ;done 

###Others

#RNA type

bowtie2-build ../final2_rdrp_contig_emboss_fliter_vOTU.fna  bt2/final2_rdrp_contig_emboss_fliter_vOTU
for i in  *unalign_fwd.fq.gz ;do echo "bowtie2 -p 6  -x bt2/final2_rdrp_contig_emboss_fliter_vOTU   -1 $PWD/${i}  -2 $PWD/${i%unalign_fwd.fq.gz}unalign_rev.fq.gz  >bt2sam/${i}.sam " ;done 

for i in bt2sam/*sam ;do echo " samtools view -b -f 99  ${i}  |  samtools sort -@ 5 > bt2sam_pos/${i##*/}_f99.sorted   " ;done
for i in bt2sam/*sam ;do echo " samtools view -b -f 147  ${i}  |  samtools sort -@ 5 > bt2sam_pos/${i##*/}_f147.sorted   " ;done  
for i in bt2sam/*sam ;do echo " samtools view -b -f 83  ${i}  |  samtools sort -@ 5 > bt2sam_neg/${i##*/}_f83.sorted   " ;done 
for i in bt2sam/*sam ;do echo " samtools view -b -f 163  ${i}  |  samtools sort -@ 5 > bt2sam_neg/${i##*/}_f163.sorted   " ;done 

for i in *sam ;do echo "  coverm contig  --bam-files  bt2sam_pos/${i}_f99.sorted  -t 3  --methods trimmed_mean  -o  ../coverm/${i}_coverm_f99.tsv  --min-read-percent-identity 90 --min-read-aligned-percent 75 -m count " ;done 
for i in *sam ;do echo "  coverm contig  --bam-files  bt2sam_pos/${i}_f147.sorted  -t 3  --methods trimmed_mean  -o  ../coverm/${i}_coverm_f147.tsv  --min-read-percent-identity 90 --min-read-aligned-percent 75 -m count " ;done >
for i in *sam ;do echo "  coverm contig  --bam-files  bt2sam_neg/${i}_f83.sorted  -t 3  --methods trimmed_mean  -o  ../coverm/${i}_coverm_f83.tsv  --min-read-percent-identity 90 --min-read-aligned-percent 75 -m count " ;done >
for i in *sam ;do echo "  coverm contig  --bam-files  bt2sam_neg/${i}_f163.sorted  -t 3  --methods trimmed_mean  -o  ../coverm/${i}_coverm_f163.tsv  --min-read-percent-identity 90 --min-read-aligned-percent 75 -m count " ;done >










