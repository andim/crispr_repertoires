#=
Get the frequency of occurrence of each contiguous block of unique elements in an array.

Note: This returns 3 (not 2) values for [a,b,a].
=#
function get_distrib_sortedarr(arr::Array)
    vals = []
    freqs = Int[]

    val_curr = arr[1]
    num_curr = 1
    for i in 2:length(arr)
        if arr[i] != val_curr
            push!(vals,val_curr)
            push!(freqs,num_curr)
            val_curr = arr[i]
            num_curr = 1
        else
            num_curr += 1
        end
    end
    push!(vals,val_curr)
    push!(freqs,num_curr)

    return vals, freqs
end

#=
Get the cumulative frequencies.

Note: length(cum) == length(freqs)+1
=#
function get_cum(freqs::Array)
    cum = zeros(Int,length(freqs)+1)
    for i in 2:length(freqs)+1
        cum[i] = cum[i-1] + freqs[i-1]
    end
    return cum
end

function get_cum_sortedarr(arr::Array)
    vals, freqs = get_distrib_sortedarr(arr)
    return vals, get_cum(freqs)
end