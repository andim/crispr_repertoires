#=
Inputs:
	bs = 50 : no. of strains in each bin
	ml = false : whether to merge last bin with previous one
=#

include("linreg.jl")
include("bin.jl")

lin(x;m=0,c=0) = m*x + c # linear function

function generate_scaling_plot(repsizes, mspacerlens; bs=50, ml=false, alldata=true,
	x="N",
	ytext="mean spacer length",
	cas=true,
	pvalue=0,
	fs=12,
	saveplot=false)
	
	# plot all data or not
	if alldata == true
		figsize = (6.4,7); ylim1 = 25; ylim2 = 56.5; ytickvals=25:5:55
	else
	    figsize = (6.4,4.8); ylim1 = 31; ylim2 = 37.5; ytickvals = 31:37
    end

    # x-axis is N or L
    if x == "N"
    	xrange = collect(0:0.1:3); xlim1 = 0; xlim2 = 3; textx = 2.7
    	xtick_vals = log10.(vcat(collect(1:10),collect(20:10:100),collect(200:100:1000)))
		xtick_labels = string.(vcat(["1"],fill("",8),["10"],fill("",8),["100"],fill("",8),["1000"]))
		xtext="spacer repertoire size"
    elseif x == "L"
    	xrange = collect(log10(7e5):0.01:log10(1.7e7)); xlim1 = log10(7e5); xlim2 = log10(1.7e7); textx=7.09
    	xtick_vals = log10.(vcat(collect(7e5:1e5:1e6),collect(2e6:1e6:1e7)))
    	xtick_labels = vcat(fill("",3),["\$10^6\$"],fill("",8),["\$10^7\$"])
    	xtext = "genome length"
    end

	### linear regression ###

	m, c, dm, dc = linreg(x=log10.(repsizes), y=mspacerlens)

	m1 = round(m, digits=2)
	dm1 = round(dm, digits=2)

	propconst = log(10)/m
	dpropconst = log(10)*dm/m^2

	propconst1 = round(propconst, digits=1)
	dpropconst1 = round(dpropconst, digits=1);

	if cas == false
		slopetext = "fit, slope = $m1 \$\\pm\$ $dm1"
	else
		slopetext = "fit, \$\\ln\\,$x = ($propconst1 \\pm $dpropconst1)\\ l_s\$ + const."
	end

	### bin strains by repsize, and find mean and standard error in bins ###
	# Sort repsizes and partition them into bins of size bs.
	# If `mergelast=true`, the last bin (which may contain < bs strains) is merged with the penultimate one.
	# Strains with equal repsize are randomly shuffled and may fall in either bin if they lie on the boundary between 2 bins.

	inds_inbins = bin_equalheight(repsizes, binsize=bs, mergelast=ml)

	mean_logN_inbins = [mean(log10.(repsizes[inds])) for inds in inds_inbins]
	stderr_logN_inbins = [std(log10.(repsizes[inds]))/sqrt(length(inds)) for inds in inds_inbins]

	mean_l_inbins = [mean(mspacerlens[inds]) for inds in inds_inbins]
	stderr_l_inbins = [std(mspacerlens[inds])/sqrt(length(inds)) for inds in inds_inbins];

	### correlation coefficient ###

	rho = round(cor(log10.(repsizes),mspacerlens), digits=2)

	### plot binned data and linear fit ###

	fig, ax = subplots(figsize=figsize, dpi=200)

	# plot mean and standard error of binned data
	errorbar(mean_logN_inbins, mean_l_inbins, linestyle="none", marker="o", markersize=1,
	    xerr=stderr_logN_inbins, yerr=stderr_l_inbins, elinewidth=1, capsize=0)

	# plot linear fit
	plot(xrange, lin.(xrange, m=m, c=c))


	# plot unbinned data
	if alldata == true
		scatter(log10.(repsizes), mspacerlens, color="green", s=0.4, alpha=0.5)

		if pvalue == 0
			pvaluetext = "<10^{-6}"
		else
			pvaluetext = string("=",round(pvalue, digits=2))
		end
		text(textx,49.8,"Pearson's \$r = $rho,\\ p$pvaluetext\$", fontsize=fs)

		legendtext = [slopetext, "individual datapoints", "data (binned)"]
    	legendloc = (0.9,0.83)
    	filename = "fit_alldata.svg"
    else
    	legendtext = [slopetext, "data (binned)"]
        legendloc = "upper left"
        filename = "fit.svg"
	end

	xlabel(xtext, fontsize=fs)
	xlim(xlim1,xlim2)
	xticks(xtick_vals, xtick_labels, fontsize=fs)
	
	ylabel(ytext, fontsize=fs)
	ylim(ylim1,ylim2)
	yticks(ytickvals, fontsize=fs)

	legend(legendtext, fontsize=fs, loc=legendloc)

	ax.spines["top"].set_visible(false)
	ax.spines["right"].set_visible(false)

	if saveplot==true
		savefig(filename,format="svg")
	end

	return rho, m, c, dm, dc, propconst1, dpropconst1, inds_inbins, mean_logN_inbins, mean_l_inbins, stderr_logN_inbins, stderr_l_inbins
end
