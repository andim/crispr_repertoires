# 1) CRISPR loci (evidence level 4)
SELECT str.genbank, seq.id AS seqid, cl.start AS locus_start, cl.length AS locus_length, cl.orientation, r.sequence,
	   cl.drconservation, cl.spacerconservation, cl.potentialorientation, cl.evidencelevelreeval, cl.blastscore
FROM crisprlocus cl
LEFT JOIN region r ON cl.drconsensus = r.id
LEFT JOIN sequence seq ON cl.sequence = seq.id
LEFT JOIN strain str ON seq.strain = str.id
WHERE cl.evidencelevel = 4
ORDER BY str.taxon, str.genbank, seq.length DESC, cl.start;

# 2) spacers (category = 3), and 3) direct repeats (category = 1)
SELECT str.genbank, seq.id AS seqid, cl.start AS locus_start, cr.start, r.sequence
FROM crisprlocus_region cr
LEFT JOIN region r ON cr.region = r.id
LEFT JOIN crisprlocus cl ON cr.crisprlocus = cl.id
LEFT JOIN sequence seq ON cl.sequence = seq.id
LEFT JOIN strain str ON seq.strain = str.id
WHERE cl.evidencelevel = 4 AND r.category = 3
ORDER BY str.taxon, str.genbank, seq.length DESC, cl.start, cr.start;

# 4) cas loci
SELECT str.genbank, seq.id AS seqid, cc.start AS locus_start, cc.length AS locus_length, cc.class AS subtype
FROM clustercas cc
LEFT JOIN sequence seq ON cc.sequence = seq.id
LEFT JOIN strain str ON seq.strain = str.id
ORDER BY str.taxon, str.genbank, seq.length DESC, cc.start;

# 5) cas genes
SELECT str.genbank, seq.id AS seqid, cc.start AS locus_start, cg.gene, cg.start AS gene_start, cg.length AS gene_length, cg.orientation AS gene_orientation
FROM clustercas_gene cg
LEFT JOIN clustercas cc ON cg.clustercas = cc.id
LEFT JOIN sequence seq ON cc.sequence = seq.id
LEFT JOIN strain str ON seq.strain = str.id
ORDER BY str.taxon, str.genbank, seq.length DESC, cc.start, cg.start;

# 6) genome sequences
SELECT str.taxon AS taxid, str.genbank, str.refseq, seq.id AS seqid, seq.category, seq.length AS seq_length, seq.ncount, seq.description
FROM strain str
LEFT JOIN sequence seq ON str.id = seq.strain
ORDER BY str.taxon, str.genbank, seq.length DESC;

# 7) taxa
SELECT *
FROM taxon
ORDER BY id;

####

# get column names
SELECT *
FROM INFORMATION_SCHEMA.COLUMNS
WHERE TABLE_NAME = N'crisprlocus';

# look at spacers (r.category = 3) or drs (r.category = 1)
SELECT *
FROM crisprlocus_region cr
LEFT JOIN region r ON cr.region = r.id
WHERE r.category = 1;

# look at cas clusters ordered by length
SELECT str.genbank, cc.class AS clustercas_class, cc.start AS clustercas_start, cc.length AS clustercas_length
FROM clustercas cc
LEFT JOIN sequence seq ON cc.sequence = seq.id
LEFT JOIN strain str ON seq.strain = str.id
ORDER BY cc.length DESC;

# look at orientation of cas genes and CRISPR loci
SELECT cg.orientation, cl.orientation, cl.potentialorientation
FROM crisprlocus cl
INNER JOIN clustercas cc ON cl.sequence = cc.sequence
INNER JOIN clustercas_gene cg ON cg.clustercas = cc.id
WHERE cg.orientation = 1;

####

# bash command for loading .sql file
psql -U postgres -d CRISPRCasdb_2021 < 20210121_ccpp.sql