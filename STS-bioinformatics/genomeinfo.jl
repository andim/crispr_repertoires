#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

# @show ARGS  # put any Julia code here

##################################################
# 2021_11_03,11_11
# Compile genome-level information including overall level of self-targeting

##### on RONIN, pwd: crispr-genome-analysis
# in Dropbox, pwd: genome-analysis

# usage: julia genomeinfo.jl
##################################################

using DelimitedFiles

seqdata = readdlm("../CRISPRCasdb/2021/seqdata.csv",',')[2:end,:] # sequence information from CRISPRCasdb
casdata = readdlm("../CRISPRCasdb/2021/casdata.csv",',')[2:end,:] # cas locus information from CRISPRCasdb
acrdata = readdlm("AcrDB/operon_statistics.csv",',')[2:end,:] # AcrDB
refseq0_acrdata = [split(r,'.')[1] for r in acrdata[:,2]] # ignore RefSeq version number
prophagedata = readdlm("prophage.tsv",'\t') # VirSorter2 results for 100 genomes

##################################################
### Run through genomes and read array, sequence, and blast output info ###
##################################################

# # read genome
# g = "GCA_000006175.2"
# # g = ARGS[1]

# read genomes
genomes = readdlm("dl.txt")

for g in genomes
    # 2) RefSeq
    refseq = seqdata[findfirst(x->x==g, seqdata[:,2]),3]

    # 3) cas subtype(s)
    subtypes = sort(unique(casdata[findall(x->x==g, casdata[:,1]),5]))
    ststring = join(subtypes,'|')

    # 4) has acr
    if refseq != "" && (ind_acr = findfirst(x->x==split(refseq,'.')[1], refseq0_acrdata)) != nothing # ignore version number
        has_acr = acrdata[ind_acr,3] == "yes" ? 1 : 0
    else
        has_acr = -1
    end

    # 5) has prophage
    if (ind_ph = findfirst(x->x==g, prophagedata[:,1])) != nothing
        has_ph = prophagedata[ind_ph,2]
    else
        has_ph = -1
    end

    # 6-9) no. of spacers at STS levels 0-3
    spacerinfo = readdlm("spacerinfo/$(g)_spacerinfo.csv",',') # spacer-level information of genome
    numspcats_g = zeros(Int,4)
    for i in 1:4
        numspcats_g[i] = count(x->x==i-1,spacerinfo[:,16])
    end

    # 10) STS level of genome
    if numspcats_g[4] > 0
        ST_g = 3
    elseif numspcats_g[3] > 0
        ST_g = 2
    # elseif numspcats_g[2]/numspcats_g[1] > .2
    elseif numspcats_g[2] >= 20
        ST_g = 1
    else
        ST_g = 0
    end

    ##################################################
    ### Append to file ###
    ##################################################
    open("genomeinfo.csv", "a") do io
        write(io, "$g,$refseq,$ststring,$has_acr,$has_ph,$(numspcats_g[1]),$(numspcats_g[2]),$(numspcats_g[3]),$(numspcats_g[4]),$ST_g\n")
    end
end