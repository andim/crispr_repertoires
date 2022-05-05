Phage sequences were downloaded as follows:

Go to https://www.ncbi.nlm.nih.gov/labs/virus/vssi/

Select top three bacterial hosts and their viruses (E coli, Salmonella enterica 28901,  P aeruginosa 287)
download complete nucleotide sequences as fasta

For reproducibility the accession numbers for all the phage genomes are provided in the file phage_accession_numbers.txt

cluster using 
mmseqs easy-cluster phages_taxid28901.fasta filtered/phages_taxid28901 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0
