#=
We compute the probability, $p_{\text{self}}$, of a match between a spacer of length $l_t = l_s + l_p$ and a host genome of length $L$. Both sequences are assumed to be random and uncorrelated, and have equal nucleotide usage.

* $l_s$ = spacer length
* $l_p$ = PAM length
* $L$ = host genome length

We consider the following model of cross-reactivity: within a length-$l_t$ alignment between spacer and host genome, up to $k_{\text{fix}}$ mismatches are tolerated at fixed positions, and up to $k_{\text{var}}$ mismatches are tolerated at any other position.

* $k_{\text{fix}}$ = no. of fixed-position mismatches
* $k_{\text{var}}$ = no. of variable-position mismatches

For a given alignment, the matching probability, $p_m$, is
$$\begin{aligned}
p_m(l_t,k_{\text{fix}},k_{\text{var}}) &= 4^{-l_t + k_{\text{fix}}} \sum_{i=0}^{k_{\text{var}}}\binom{l_t - k_{\text{fix}}}{i} 3^{i}\\
&\approx 4^{-l_t + k_{\text{fix}}} \binom{l_t - k_{\text{fix}}}{k_{\text{var}}} 3^{k_{\text{var}}},
\end{aligned}$$
since the largest term tends to dominate the sum.

Then, the probability that there exists an alignment somewhere in the host genome is
$$\begin{aligned}
p_{\text{self}}(L,l_t,k_{\text{fix}},k_{\text{var}}) &= 1 - [1 - p_m(l_t,k_{\text{fix}},k_{\text{var}})]^L\\
&\approx L p_m(l_t,k_{\text{fix}},k_{\text{var}}),
\end{aligned}$$
assuming $p_m \ll 1$.

Below, we define the following functions:

|                           | all terms of sum | largest term of sum |
|---------------------------|:----------------:|:-------------------:|
| not taking $p_m \ll 1$    | `F_self()`       | `p_self1()`         |
| taking $p_m \ll 1$        |                  | `p_self()`          |
=#

using SpecialFunctions

logbinom(l,k) = loggamma(l+1) - loggamma(l-k+1) - loggamma(k+1)

k_tot(k_fix,k_var,l_t) = k_fix + (logbinom(l_t-k_fix,k_var) + k_var*log(3))/log(4)

p_m(l_t,k_fix,k_var) = 4.0^(- l_t + k_tot(k_fix,k_var,l_t))

p_self(L,l_t,k_fix,k_var) = L * p_m(l_t,k_fix,k_var)

p_self1(L,l_t,k_fix,k_var) = 1 - exp(L * log(1 - p_m(l_t,k_fix,k_var)))

#=
$F_{\text{self}}(k)$ is also the CDF of the number of alignments with up to $k$ variable-position mismatches. Its PDF is
$$p_{nn}(k) = F_{\text{self}}(k) - F_{\text{self}}(k-1).$$

Note: `F_self()` and `p_nn()` require integer `k_var`, whereas `p_self()` and `p_self1()` allow non-integer `k_var` as input.
=#

F_self(L,l_t,k_fix,k_var::Int) = 1 - exp(L * log(1 - sum([p_m(l_t,k_fix,k) for k in 0:k_var])))

p_nn(L,l_t,k_fix,k_var::Int) = F_self(L,l_t,k_fix,k_var) - F_self(L,l_t,k_fix,k_var-1)