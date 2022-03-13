#=
Equal number of elements in each bin.
=#
function bin_equalheight(X::Array;binsize,mergelast::Bool=false)
	inds = sortperm(X)

	inds_bybin = Array{Int}[]

	a = 1
	b = binsize
	while b < length(inds)
		push!(inds_bybin,inds[a:b])
		a += binsize
		b += binsize
	end

	if mergelast == false
		push!(inds_bybin,inds[a:end])
	else
		inds_bybin[end] = vcat(inds_bybin[end],inds[a:end])
	end

	return inds_bybin
end