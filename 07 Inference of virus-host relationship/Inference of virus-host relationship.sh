

#### CRISPR-spacer screening
### retrieve CRISPR-spacer sequence from 15,640 medium-high quality
## MinCED v0.4.2 and prodigal v2.6.3

prodigal -p meta -q -m -i ${MAG}.fna -a ${MAG}_prodigal.faa -d ${MAG}_prodigal.fna -o ${MAG}_prodigal.gff -f gff
minced ${MAG}.fna ${MAG}.cripsers  ${MAG}_prodigal.gff

### All spacers were queried for matches against viral contigs
## BLAST v2.13

makeblastdb -in All_crisper.fa  -dbtype nucl -out All_crisperdb
blastn -query  final_rdrp_contig_fliter.fa  -db All_crisperdb -perc_identity 95 -dust no -word_size 7 -out All_host_prediction.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend sstart send evalue bitscore'  -num_threads 10  -evalue 1e-10 


#### EVE based
## EVEs are genetic sequences derived from DNA or RNA viruses that have integrated into host genomes
# tblastn
tblastn -query final_rdrp_contig_emboss_fliter_rdrp_complete.faa -db NCBI_database/NT_202307/NT_202307    -out NT_tblastn_r.txt -max_target_seqs 20 -outfmt '7 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend qlen sstart send  slen evalue bitscore'  -num_threads 10  -evalue 1e-20 


#### zoonic_rnak detection
## zoonotic_rank with default parameters was applied to identify and prioritize potential human-infecting viruses

Rscript ~/minisoft/zoonotic_rank-main/Scripts/PredictNovel.R fasta  final_rdrp_contig_emboss_fliter.fa metadata.tsv  TPC_risk_re

#### Identification of Shine-Dalgarno (SD) sequences of RNA viruses
### The Shine-Dalgarno (SD) sequence is an RBS in bacterial and archaeal mRNA. The prevalence of SD in an RNA virus indicates it is more likely a prokaryotic RNA virus
# OSTIR v1.1.0 

ostir -i final_rdrp_contig_emboss_fliter.fa -o OSTIR -j 10
