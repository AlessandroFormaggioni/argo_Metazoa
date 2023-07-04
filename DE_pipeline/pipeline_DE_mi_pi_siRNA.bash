#!/bin/bash

sp=DaRe
sp2=dar


#Remove all the reads that have been already assigned to miRNAs. 
bedtools intersect -v -a ../DaRe_sorted.bam -b ../mirdeep/mirna_results_11_11_2022_t_17_54_51/known_+_novel_pre.bed > ../DaRe_intersect.bam



for a in 20
do
b=$(( a-2 ))

#Select reads of $a length
sed -n '1~4s/^@/>/p;2~4p' ../${sp}_merged.fastq  | bioawk -v a=$a -c fastx '(length($seq) == a) {print ">"$name;print $seq}' - > ../split_fasta/${sp}_merged_${a}.fasta

#Collpase identical reads
collapse_reads_md.pl ../split_fasta/${sp}_merged_${a}.fasta $sp2 > ../split_fasta/${sp}_merged_collapsed_${a}.fasta


#overlapping_reads.py script; it looks for piRNA (10 bp overlap) and siRNA ($a-2 overlap) signature
snakemake -s Snakefile --cores 7 ${sp}_${a}_overlap_{10,$b}

#overlap_tab_fasta.py script; the script produces for each ${sp}_${a}_overlap_{10,$a-2} file a ${sp}_${a}_overlap_{10,$a-2}_table_fast file. The latter file is a tabular file that reports a overlapping pair for each line, reporting the coverage of each pair.   
snakemake -s Snakefile2 --cores 7 ${sp}_${a}_overlap_{10,${b}}_table_fast

if [ -d built ]; then
rm -r built
fi

mkdir built

#Selecting miRNAs of $a length
bioawk -v a=$a -c fastx '(length($seq) == a) {print ">"$name;print $seq}' known_+_novel_matured.fa > built/mirna_${a}.fa
mv ${sp}_${a}_overlap_*_table_fast.fasta built

#Indexing fasta files containing miRNAs, piRNAs and siRNAs of $a length
bowtie-build built/${sp}_${a}_overlap_10_table_fast.fasta built/${sp}_${a}_overlap_10_table_fast.fasta
bowtie-build built/${sp}_${a}_overlap_${b}_table_fast.fasta built/${sp}_${a}_overlap_${b}_table_fast.fasta
bowtie-build built/mirna_${a}.fa built/mirna_${a}.fa

#Mapping $a long reads against $a long miRNAs, piRNAs and siRNAs
bowtie  -v 0 built/mirna_${a}.fa  -f ../split_fasta/${sp}_merged_${a}.fasta  -S mirna_${a}_mapped
bowtie  -v 0 built/${sp}_${a}_overlap_10_table_fast.fasta  -f ../split_fasta/${sp}_merged_${a}.fasta  -S pirna_${a}_mapped
bowtie  -v 0 built/${sp}_${a}_overlap_${b}_table_fast.fasta  -f ../split_fasta/${sp}_merged_${a}.fasta  -S sirna_${a}_mapped

done

#Producing the smallRNA-count table for each sample (library).
#for miRNA

for a in `awk '(substr($1,1,1) == ">") {print$1}' known_+_novel_matured.fa  | sed 's/>//g' | sort | uniq`; do 

	gp=$(grep -w $a  mirna_ALL_mapped.sam)

	sp1=$(echo "$gp" | grep -c "SRR11713497")
	sp2=$(echo "$gp" | grep -c "SRR11713498")
	sp3=$(echo "$gp" | grep -c "SRR11713499")
	sp4=$(echo "$gp"| grep -c "SRR11713500")
	sp5=$(echo "$gp"| grep -c "SRR11713501")
	sp6=$(echo "$gp"| grep -c "SRR11713502")


	echo -e "${a}\t${sp1}\t${sp2}\t${sp3}\t${sp4}\tmiRNA" >> miRNA_read_count.tab

done



#for piRNA

for b in `seq 20 30`; do 

	samtools view pirna_${b}_mapped.bam > pirna_${b}_mapped.sam

	for a in `awk '(substr($1,1,1) != ">") {print$1}' AnGa_${b}_overlap_10_table_fast.fasta | sort | uniq`; do 
		gp=$(grep -w $a  pirna_${b}_mapped.sam)

		sp1=$(echo "$gp" | grep -c "SRR11713497")
		sp2=$(echo "$gp" | grep -c "SRR11713498")
		sp3=$(echo "$gp" | grep -c "SRR11713499")
		sp4=$(echo "$gp"| grep -c "SRR11713500")
		sp5=$(echo "$gp"| grep -c "SRR11713501")
		sp6=$(echo "$gp"| grep -c "SRR11713502")

		echo -e "${a}\t${sp1}\t${sp2}\t${sp3}\t${sp4}\t${sp5}\t${sp6}\tpiRNA" >> piRNA_read_count.tab

	done
done


#for siRNA

for b in `seq 20 30`; do 

	c=$(( b-2 ))
	samtools view sirna_${b}_mapped.bam > sirna_${b}_mapped.sam

	for a in `awk '(substr($1,1,1) != ">") {print$1}' AnGa_${b}_overlap_${c}_table_fast.fasta | sort | uniq`; do 
		gp=$(grep -w $a  sirna_${b}_mapped.sam)

		sp1=$(echo "$gp" | grep -c "SRR11713497")
		sp2=$(echo "$gp" | grep -c "SRR11713498")
		sp3=$(echo "$gp" | grep -c "SRR11713499")
		sp4=$(echo "$gp"| grep -c "SRR11713500")
		sp5=$(echo "$gp"| grep -c "SRR11713501")
		sp6=$(echo "$gp"| grep -c "SRR11713502")

		echo -e "${a}\t${sp1}\t${sp2}\t${sp3}\t${sp4}\t${sp5}\t${sp6}\tsiRNA" >> siRNA_read_count.tab

	done
done

#For piRNA and siRNA we need a script that sums up reads that differ of a few bases but shares the same location on the genome 

python piRNA_readcount_x_sample.py
python siRNA_readcount_x_sample.py


cat miRNA_read_count.tab piRNA_readcount_x_sample.tab siRNA_readcount_x_sample.tab > ALL_readcount.tab

