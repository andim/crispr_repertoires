using Statistics

#=
Linear regression

Note: var() and cov() use corrected=true (1/(N-1) rather than 1/N) by default.

Reference: `Data Analysis: A Bayesian Tutorial` by D. Sivia, Section 3.5.1
=#
function linreg(;x,y,dy=0)
    if dy==0 # dy is NOT provided
        m = cov(x,y)/var(x)
        c = mean(y) - m*mean(x)
        N = length(x)
        S2 = (N-1)/(N-2)*var(m*x .+ c .- y)
        dm = sqrt(S2/N/var(x,corrected=false))
        dc = sqrt(S2*mean(x.^2)/N/var(x,corrected=false))

    else # dy is provided
        α = sum(x.^2 ./ dy.^2)
        β = sum(1    ./ dy.^2)
        γ = sum(x    ./ dy.^2)
        p = sum(x.*y ./ dy.^2)
        q = sum(y    ./ dy.^2)
        
        denom = α*β - γ^2
        m = (β*p - γ*q)/denom
        c = (α*q - γ*p)/denom
        dm = sqrt(β/denom)
        dc = sqrt(α/denom)
    end
    
    return m, c, dm, dc
end