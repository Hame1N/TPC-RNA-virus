
### CD-HIT v4.8.1

# CD-HIT for vOTU clustering with a clustering threshold of 90% amino acid identity (AAI) and default parameters

cd-hit -i  final_rdrp_contig_emboss_fliter_rdrp.faa  -o  final_rdrp_contig_emboss_fliter_05.fna  -c  0.9  -d 0 

# CD-HIT for C90 VContigs clustering with a clustering threshold of 90% average nucleotide identity (ANI) and 80% alignment fraction (AF)

cd-hit -i  final_rdrp_contig_emboss_fliter.fna   -aS 0.8  -o  final_rdrp_contig_emboss_fliter_90.fna   -c  0.9  -d 0


### Markov Cluster Algorithm (MCL) clustering for taxonomy assignment 

### USEARCH v10.0.240
# complete rdrp domain sequences form TPC, RVMT and TO dataset were clustered by USEARCH
cat final_rdrp_contig_emboss_fliter_rdrp_compl.faa RVMT_rdrp_compl.faa  TO_rdrp_compl.faa >TBPCV_TO_RVMT_compl_rdrp.faa
/data01nfs/user/liupf/software_lpf/usearch_v11_dyc  --cluster_fast  TBPCV_TO_RVMT_compl_rdrp.faa  -id 0.50 -sort length  -centroids  TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa  -uc clusters.uc -clusters Final

### blastp v2.13
# all centroid sequences were cross-compared running all-against-all pairwise 
makeblastdb -in TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa -dbtype prot -out TBPCV_TO_RVMT_compl_rdrp_0.5uclust
blastp -query TBPCV_TO_RVMT_compl_rdrp_0.5uclust.faa -db TBPCV_TO_RVMT_compl_rdrp_0.5uclust  -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore'   -gapopen 9 -gapextend 1 -word_size 3  -out Pairwise_CVTP_TR_RVMT_rdrp_0.5uclust.tsv  -num_threads 5 

### MCL v14.137
grep "#" -v Pairwise_CVTP_TR_RVMT_rdrp_0.5uclust.tsv |awk -F"\t" '{if($12==0) print $1,$2,"1"; else print $1,$2,$12}' OFS="\t" >MCL_input.tsv
mcxload -abc MCL_input.tsv --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o MCL_out.mci -write-tab MCL_out.tab
mcl MCL_out.mci  -I 1.1 -use-tab MCL_out.tab -te 20


### rdrp level 
### blastp v2.13
# blast to reference databses for taxonomy assignment 
#RVMT rdrp database
diamond blastp --db /datanode02/zhangzh/database/RVMT_RdRP/RVMT_DB.dmnd  --query final_rdrp_contig_emboss_fliter_rdrp.faa --out final_rdrp_contig_emboss_fliter_rdrp_RVMT.tsv -f 6 qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore  --sensitive -k 1 --max-target-seqs 1   -p 5
#TaraOcean rdrp database
diamond blastp --db /datanode02/zhangzh/database/TaraOcean/TaraOcean_DB.dmnd --query final_rdrp_contig_emboss_fliter_rdrp.faa --out final_rdrp_contig_emboss_fliter_rdrp_TO.tsv -f 6 qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send evalue bitscore  --sensitive -k 1 --max-target-seqs 1   -p 5
#paldb
psiblast -query final_rdrp_contig_emboss_fliter_rdrp_palmscan_pp.faa -db /datanode02/zhangzh/database/palmdb-main/palmdb_serratus -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_rdrp_palmscan_psiblast_serratusDB.tsv  -num_threads 20 
psiblast -query final_rdrp_contig_emboss_fliter_rdrp_palmscan_pp.faa -db /datanode02/zhangzh/database/palmdb-main/palmdb_name  -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_rdrp_palmscan_psiblast_nameDB.tsv  -num_threads 20


# vContig level 
## blast to reference databaase of vContigs
cat final_rdrp_contig_emboss_fliter.fna  RiboV1.6_Contigs_rename.fasta RNA_viral_rename.fna  44779_RdRP_contigs_rename.fna >TBPCV_TO_REFV_RVMT_conitg.fna

makeblastdb -in TBPCV_TO_REFV_RVMT_conitg.fna  -dbtype nucl -out TBPCV_TO_REFV_RVMT_conitg
blastn -query TBPCV_TO_REFV_RVMT_conitg.fna   -db TBPCV_TO_REFV_RVMT_conitg -outfmt '6 std qlen slen'  -max_target_seqs 10000 -out blast/TBPCV_TO_REFV_RVMT_conitg_ani.tsv -num_threads 20
## checkv-1.0.1
anicalc.py -i blast/TBPCV_TO_REFV_RVMT_conitg_ani.tsv  -o blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv 
aniclust.py --fna TBPCV_TO_REFV_RVMT_conitg.fna --ani blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv  --out blast/TBPCV_TO_REFV_RVMT_conitg_ani_clusters9090.tsv  --min_ani 90 --min_tcov 90 --min_qcov 0
aniclust.py --fna TBPCV_TO_REFV_RVMT_conitg.fna --ani blast/TBPCV_TO_REFV_RVMT_conitg_aniclust.tsv  --out blast/TBPCV_TO_REFV_RVMT_conitg_ani_clusters9075.tsv  --min_ani 90 --min_tcov 75 --min_qcov 0

## megablast 
blastn -query final_rdrp_contig_emboss_fliter.fna -db ~/database/RVMT_contig/RiboV1.6_Contigs -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_blastRVMT.tsv -num_threads 20
blastn -query final_rdrp_contig_emboss_fliter.fna -db ~/database/TaraOcean/Tara_RdRP_contigs  -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_blastTO.tsv -num_threads 20
blastn -query final_rdrp_contig_emboss_fliter.fna -db /datanode02/zhangzh/database/IMGV/IMGVR_all_nucleotides  -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_blastIMGV.tsv -num_threads 20
blastn -query final_rdrp_contig_emboss_fliter.fna  -db ~/database/refseq_viral230504/refseq_RNAviral -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_blastrefviralRNA.tsv -num_threads 20
blastn -query final_rdrp_contig_emboss_fliter.fna -db ~/database/RVMT_contig/RiboV1.6_Contigs -use_index true -task megablast -outfmt '7 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1  -out final_rdrp_contig_emboss_fliter_blastRVMT.tsv -num_threads 20