#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

# @show ARGS  # put any Julia code here

##################################################
# 2021_08_17, 10_16, 10_29-30
# Compile spacer information for each genome, including level of self-targeting

# on RONIN, pwd: crispr-genome-analysis
# in Dropbox, pwd: genome-analysis

# usage: julia spacerinfo.jl
##################################################

##################################################
### Read spacer, dr, and Vink data ###
##################################################
using DelimitedFiles

spacerdata = readdlm("../CRISPRCasdb/2021/spacerdata_4.csv",',')
drdata = readdlm("../CRISPRCasdb/2021/drdata_4.csv",',')
Vinkdata = readdlm("Vinkdata.csv",',')[2:end,1:13]

#=
Parse Vinkdata
=#
seqIDs = String[]
spinds_seqs = Array{Int}[]

for i in 1:size(Vinkdata,1)
    IDs = strip.(split(Vinkdata[i,4],'|'))
    for id in IDs
        if in(id,seqIDs) # accession no. already present
            ind = findfirst(x->x==id,seqIDs)
            push!(spinds_seqs[ind],i)
        else # new accession no.
            push!(seqIDs,id)
            push!(spinds_seqs,[i])
        end
    end
end


##################################################
### Functions ###
##################################################
#=
Find the BLAST output row containing the closest match outside a CRISPR array to save
=#
function find_closestmatch(blastout, inds_s, arrayinfo, seqinfo)
    ind_tosave = 0
    len_max = 0
    for ind in inds_s
        m = blastout[ind,:]

        # 1) get no. of matches in alignment
        len_aln = m[4] - m[5]

        # 2) check if it's longer than the currently saved longest alignment
        if len_aln > len_max

            # 3) check if it's outside a CRISPR array
            if inarray(string(m[2]),sort([m[9],m[10]]),arrayinfo,seqinfo) == false
                len_max = len_aln
                ind_tosave = ind
            end
        end
    end

    return len_max, ind_tosave
end

#=
Returns true if alignment overlaps with *any* CRISPR array, false otherwise
Input:
    seq : accession no. of matching sequence
    pos : [start,end] of alignment on matching sequence (start < end)
    arrayinfo : information about CRISPR arrays in genome (rows of crisprdata_4.csv)
    seqinfo : information about sequences in genome
=#
function inarray(seq::String, pos::Array{Int}, arrayinfo, seqinfo)
    seqid = seqinfo[findfirst(x->x==seq,seqinfo[:,1]),4] # for identifying accession no. with sequence in CRISPRCasdb
    for i in 1:size(arrayinfo,1)
        if seqid == arrayinfo[i,2] && # array is on matching sequence
            length(intersect(pos[1]:pos[2],arrayinfo[i,3]:arrayinfo[i,3]+arrayinfo[i,4]-1)) > 0 # alignment overlaps with array
            ### println("spacer match to CRISPR array")
            return true
        end
    end
    return false
end

#=
Compute spacer category (0, 1, 2, or 3)
=#
function get_spcat(l_s::Int,len_max::Int)::Int
    if len_max == 0 || l_s - len_max > 10 # no self-target outside of CRISPR arrays (or with >10 mismatches)
        return 0
    elseif l_s - len_max > 3 # self-target with 4-10 mismatches
        return 1
    elseif l_s - len_max > 0 # self-target with 1-3 mismatches
        return 2
    else # exact self-target
        return 3
    end
end


##################################################
### Run through genomes and read array, sequence, and blast output info ###
##################################################

# # read genome
# g = "GCA_000006175.2"
# g = ARGS[1]

# read genomes
genomes = readdlm("dl.txt")

for g in genomes
    arrayinfo = readdlm("crisprs/$(g)_CRISPR.csv",',') # information about CRISPR arrays in genome (rows of crisprdata_4.csv)
    seqinfo = readdlm("seqinfo/$(g)_seqinfo.tsv",'\t') # information about sequences in genome
    # blastout = readdlm("test-STS/test.csv",',')
    blastout = readdlm("blast/$(g)_blast.csv",',')

    #=
    Parse spacer info (BLAST output col. 1); sid = [genbank, seqid, locus_start, spacer_start, spacer_length]
    =#
    sids_blastout = [split(m,'|') for m in blastout[:,1]]


    ##################################################
    ### Run through CRISPR arrays and compile spacer info ###
    ##################################################
    for i in 1:size(arrayinfo,1)
        # array information
        seqid_a = strip(arrayinfo[i,2],'"') # seqid
        start_a = arrayinfo[i,3] # array start
        drcons_a = strip(arrayinfo[i,6],'"') # dr consensus
        accnum_a = seqinfo[findfirst(x->x==seqid_a,seqinfo[:,4]),1] # sequence accession no.
        Vind_a = findfirst(x->x==accnum_a, seqIDs) # index in Vink seqIDs

        # all BLAST output rows corresponding to array
        inds_ainblastout = findall(x->x[2]==seqid_a && parse(Int,x[3])==start_a, sids_blastout)

        # all spacer IDs in array
        sids_a = unique(sids_blastout[inds_ainblastout])

        # find corresponding rows in spacerdata & drdata
        inds_s0 = findall(i->spacerdata[i,1]==g && spacerdata[i,2]==seqid_a && spacerdata[i,3]==start_a, 1:size(spacerdata,1))
        inds_d0 = findall(i->drdata[i,1]==g && drdata[i,2]==seqid_a && drdata[i,3]==start_a, 1:size(drdata,1))

        # run through each spacer in array
        for j in 1:length(sids_a)
            sid = sids_a[j]

            # spacer information
            start_s = parse(Int,sid[4])
            l_s = parse(Int,sid[5])
            seq_s = spacerdata[inds_s0[j],5]
            dr1_s = drdata[inds_d0[j],5]
            dr2_s = drdata[inds_d0[j+1],5]

            # obtain spacer information from Vinkdata
            # Vind_s = spinds_seqs[Vind_a][findfirst(x->x==seq_s, Vinkdata[spinds_seqs[Vind_a],2])]
            # subtype_s = Vinkdata[Vind_s,5]
            # orient_s = Vinkdata[Vind_s,12]
            # PAM_s = Vinkdata[Vind_s,13]
            
            if Vind_a != nothing && (Vind = findfirst(x->x==seq_s, Vinkdata[spinds_seqs[Vind_a],2])) != nothing
                Vind_s = spinds_seqs[Vind_a][Vind]
                subtype_s = Vinkdata[Vind_s,5]
                orient_s = Vinkdata[Vind_s,12]
                PAM_s = Vinkdata[Vind_s,13]
            else
                orient_s = ""
                subtype_s = ""
                PAM_s = ""
            end

            # all BLAST output rows corresponding to spacer
            inds_s = findall(x->x[2]==sid[2] && x[3]==sid[3] && x[4]==sid[4], sids_blastout)
            
            len_max, ind_tosave = find_closestmatch(blastout, inds_s, arrayinfo, seqinfo)

            c = get_spcat(l_s,len_max)


            ##################################################
            ### Output files ###
            ##################################################
            # append to spacer info file (1 line per spacer)
            open("spacerinfo/$(g)_spacerinfo.csv", "a") do io
                write(io, "$g,$i,$accnum_a,$start_a,$drcons_a,$subtype_s,$orient_s,$PAM_s,$j,$start_s,$l_s,$seq_s,$dr1_s,$dr2_s,$len_max,$c\n")
            end

            # append to new BLAST output file (1 line per spacer)
            open("blast1/$(g)_blast1.csv", "a") do io
                if c > 0
                    arr = string.(blastout[ind_tosave,:])
                    str = string([arr[i]*"," for i in 1:11]...,arr[end])
                    write(io, "$str\n")
                else
                    write(io, "$(blastout[inds_s[1],1]),,,,,,,,,,,\n")
                end
            end
        end
    end
end