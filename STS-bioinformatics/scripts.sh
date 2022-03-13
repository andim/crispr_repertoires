# Input:
# g : GenBank accession number
# f : genomic sequence .fna filename

###################################
### generate CRISPR info and spacer .fna files ###
###################################

# construct table of CRISPR arrays and fasta file of spacers for each genbank
# since strings in .csv files have "" around them, use both ',' and '"' as field separators
while read g; do
awk -F',' -v g="\"$g\"" ' $1==g ' CRISPRCasdb/crisprdata_4.csv > crisprs/${g}_CRISPR.csv
awk -F',' -v g="\"$g\"" ' $1==g ' CRISPRCasdb/spacerdata_4.csv | awk -F'[,"]' '{ printf ">%s|%s|%s|%s|%s\n%s\n" ,$2,$5,$7,$8,length($10),$10 }' > crisprs/${g}_spacers.fna
done < gbs_cc.txt

###################################
### BLAST ###
###################################

# BLAST spacers against own genome
makeblastdb -in $f -parse_seqids -blastdb_version 5 -dbtype nucl
blastn -task blastn-short -evalue 1 -num_threads 16 -perc_identity 60 -word_size 4 -reward 5 -penalty -4 -ungapped -outfmt 10 -db $f -query crisprs/${g}_spacers.fna -out blast/${g}_blast.csv

###################################
### find prophages ###
###################################

# run VirSorter2 on NCBI genomes
virsorter run -w $g.out -i $f --min-length 1500 -j 16 all

# check for presence of prophage
for f in $(find . -name final-viral-score.tsv | sort); do
if [ $(cat $f | wc -l) -gt 1 ]; then
	echo 1 >> prophage.txt
else
	echo 0 >> prophage.txt
fi
done

###################################
### match up CRISPRCasdb and NCBI ###
###################################

# compile NCBI headers and CRISPRCasdb seqdata together for each sequence of genomes
while read g; do
grep $g CRISPRCasdb/seqdata.csv | tr -d '"' | awk -F',' '{OFS = "\t"}{print $4,$5,$6,$7,$8}' | sort -nr -k3 > A.txt
f=$(find genomes/$g/ncbi_dataset/data/$g -name *.fna)
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' $f | tr -d ">" | paste -sd "\t\n" | sort -nr -t$'\t' -k2 | awk -F' ' '{OFS = "\t"}{print $1,$0}' | paste - A.txt > seqinfo/${g}_seqinfo.tsv
done < dl.txt
