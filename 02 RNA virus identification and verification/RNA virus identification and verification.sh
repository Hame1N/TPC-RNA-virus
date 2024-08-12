

### prodigal v2.6.3
## Open reading Frame predication

prodigal -p meta -q -m -i  ${sample}.contigs1k.fa -a ${sample}.contigs1k_prodigal.faa -d ${sample}.contigs1k_prodigal.fna -o ${sample}.contigs1k_prodigal.gff -f gff

## palmscan

### RdRP domians identification

/datanode03/zhangxf/palmscan-main/bin/palmscan --search_pp ${sample}.contigs1k_prodigal.faa  -rdrp -ppout ${output_dir}/${ini_name}_${out_su}_palmscan_pp.faa  -report ${output_dir}/${ini_name}_${out_su}_palmscan_pp.txt -fevout ${output_dir}/${ini_name}_${out_su}_palmscan_pp.fev -hiconf -threads ${THREADS}

## hmmsearch v3.3.2
### RdRP domians identification accroding the three hmm file;  HMMs from the Tara Oceans RNA virus dataset NeoRdRp pipeline and the RVMT global RNA virome dataset were used for RdRp domain search.

eval=0.01
THREADS=20
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_NeoRdRp.txt /data01nfs/user/qinfsh/NeoRdRp-HMM.v1.1.hmm ${output_dir}/${ini_name}_${out_su}_palmscan_pp.faa
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_all_virus.txt /data01nfs/user/qinfsh/all_virus_RdRP_profiles.hmm ${output_dir}/${ini_name}_${out_su}_palmscan_pp.faa
hmmsearch --cpu ${THREADS} -E ${eval} --domE ${eval} --incdomE ${eval} --incE ${eval}  --tblout ${output_dir}/${ini_name}_${out_su}_hmmout_rdrpscan.txt /data01nfs/user/liupf/software_lpf/RdRp-scan/Profile_db_and_alignments/ALL_RdRp_scan_profiles.hmm ${output_dir}/${ini_name}_${out_su}_palmscan_pp.faa

## merge hmmsearch results and get fianl rdrp sequecens


