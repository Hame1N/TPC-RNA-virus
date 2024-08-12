
#### Function annoation for vContigs
### ORFs based  default code 
## InterProScan v5.62-94.0
# database: MobiDBLite (v2.0),PRINTS (v.42.0), Phobius (v.1.0.1) and TMHMM (v.2.0c) 
interproscan.sh -i final_rdrp_contig_emboss_fliter.fna -t n  -appl MobiDBLite,PRINTS,Phobius,TMHMM  -goterms  -pa  -cpu 50

### ORFs based on optimized code
## HHsuite v3.3.0
# database: UniRef30_2020_03 database and Pfam 35
# seq= single ORFs protein sequence wiht optimized genetic code

hhblits -cpu 10 -i  ${seq}  -d UniRef30_2023_02 -oa3m a3m/${seq}.a3m   -n 1 -e 0.001
hhsearch -cpu 10 -i a3m/${seq}.a3m  -d  pfam  -o hh_r/${seq}.hhr   -n 1 -e 0.001 

## Domian-based annotation
## HMMER V3.3.2
# database:PFam 35, CDD v.3.19, Gene3D v4.3, ECOD 2017 release93, LysDB
# All_CDD.hmm 
# pfam_a.hmm
# gene3d_main.hmm 
# ecodf.hmm
# AllHMMs.hmm 

#!/bin/bash
## Aouhter: Zhangzh
## Date: 29,July,2023  

##################################################################################################
if [[ $# -eq 0 ]]; then
	echo ' 
   run hmmsearch against different database;     
   Arguments:
   #	Desc (suggestion)	
   1	Threads
   2	output directory
   3	Input (expected *.faa)
   4	suffix for output 
   5    hmm database (full path)

'
	exit
fi
##################################################################################################
THREADS=$1
output_dir=$2
input=$3
out_su=$4
hmm_db=$5

eval=0.01

ini_name=$(basename $input ".faa")
echo $ini_name "start hmmsearch"

if test -f "$output_dir"; then
   echo "$FILE exist"
else
   mkdir $output_dir 
fi

hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  -o ${output_dir}/${ini_name}_${out_su}_hmm_result.tsv --domtblout ${output_dir}/${ini_name}_${out_su}_domtblout_hmm_result.tsv --tblout ${output_dir}/${ini_name}_${out_su}_tabular_hmm_result.tsv ${hmm_db}  ${input}

echo "$hmm_db hmmserch finished "

echo "All done"

#### Dram annotation for host MAGs
## DRAM v1.4.6

DRAM.py annotate -i ${MAG}.fna --use_uniref --use_vogdb  -o Dram/${MAG}.annotation --threads 10
DRAM.py distill -i Dram/${MAG}.annotation/annotations.tsv  -o  Dram/${MAG}.annotation/dram_distll


