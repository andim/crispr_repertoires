# p-value for one-sided permutation test on linear regression slope
function pvalue_linreg(m0, x, y; num::Int = Int(1e6))
    p = 0
    for i in 1:num
        m, c, dm, dc = linreg(x=shuffle(x), y=y)
        if m > m0
            p += 1
        end
    end
    
    return p/num
end

# p-value for one-sided permutation test on Pearson's r
function pvalue_Pearson(rho0, x, y; num::Int = Int(1e6))
    p = 0
    for i in 1:num
        rho = cor(shuffle(x),y)
        if rho > rho0
            p += 1
        end
    end
    
    return p/num
end
