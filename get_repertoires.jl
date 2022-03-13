#=
Return the array IDs and list of spacers for all arrays of each gb.

Output:
    arrayids : for each strain in gbs, a list of CRISPR arrays in the format [gb, seqid, locus_start]
    spacerarrays : for each strain in gbs, a list of spacers in each array
=#
function get_spacerarrays(gbs::Array)
    arrayids = Array[]
    spacerarrays = Array[]

    for g in gbs
        arrayids_g = Array[]
        spacerarrays_g = Array{String}[]

        for c in genbank2crisprloci[g]
            push!(arrayids_g,c)
            push!(spacerarrays_g,crisprlocus2spacers[c])
        end

        push!(arrayids,arrayids_g)
        push!(spacerarrays,spacerarrays_g)
    end

    return arrayids, spacerarrays
end


#=
Return the list of all spacers for each gb.
=#
function get_spacers(gbs::Array)
    arrayids, spacerarrays = get_spacerarrays(gbs)

    return [vcat(s...) for s in spacerarrays]
end


#=
Return the number of gbs, repertoire size, and all spacer lengths for each gb.
Output:
    num_gbs : length of gbs
    repsizes : number of spacers in repertoire
    spacerlens : lengths of spacers in repertoire
=#
function get_repertoires(gbs::Array; verbose=1)
    spacers = get_spacers(gbs)
    
    num_gbs = length(gbs)
    if verbose==1; println("$num_gbs repertoires returned."); end

    repsizes = Int[length(s) for s in spacers]
    spacerlens = Array{Int}[length.(s) for s in spacers]

    return num_gbs, repsizes, spacerlens
end
