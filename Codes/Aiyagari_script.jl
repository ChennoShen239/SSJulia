using Plots
include("HA_fcns.jl")
include("HA_setup.jl")


compute_ss!(ss,par,grid);



# make the classic Aiyagari (1994) figure
Kgrid = KCompleteMarkets .+ collect(0:0.25:1.5) * (ss[:K] - KCompleteMarkets);
Ksupply = [K+checkK!(K,par,grid) for K in Kgrid];
Rgrid = [getWR(par,K)[2] for K in Kgrid];
plot(Ksupply,Rgrid,label="Capital supply",ylabel="R",xlabel="K",width=3);
p = plot!([Kgrid;48],[Rgrid;getWR(par,48)[2]],xlim=[40;60], label="Capital demand",legend=:bottomright,width=3)
display(p)
savefig("../Figures/AiyagariSD.png")

# plot the distribution of savings
normfactor = [1;diff(grid.a)];
density = reshape(ss[:D],length(grid.a),:) ./normfactor;
p = plot(grid.a,density,xlabel="Assets (before interest)",ylabel="Density", label=["Low endow" "High endow"],width=3,xlim=[-3;250])
display(p)

# Compute the cumulative distribution functions (CDFs) for the two groups

# First, recover each group’s probability mass by multiplying the density by the integration weights
prob_low  = density[:, 1] .* normfactor   # probability mass for low endowment group
prob_high = density[:, 2] .* normfactor   # probability mass for high endowment group

# Compute the CDFs as cumulative sums of the probability masses
cdf_low  = cumsum(prob_low)
cdf_high = cumsum(prob_high)

# Plot the CDFs for both groups
p_cdf = plot(grid.a, cdf_low, label="Low endow CDF", xlabel="Assets", ylabel="CDF", lw=3)
plot!(p_cdf, grid.a, cdf_high, label="High endow CDF", lw=3)
display(p_cdf)
savefig("../Figures/Aiyagari_CDF.png")

# Define a helper function to compute the Gini coefficient from asset levels and their probability weights
function gini_coefficient(a, weights)
    μ = sum(a .* weights)    # Mean asset level
    # Compute the pairwise absolute differences matrix between asset levels
    diff_matrix = abs.(a .- a')
    # Form the outer product of weights (each element gives the joint probability weight)
    weight_matrix = weights * weights'
    return sum(diff_matrix .* weight_matrix) / (2 * μ)
end

# Compute the Gini coefficient for the entire population
total_prob = prob_low .+ prob_high  # Combined probability mass for all individuals
gini_total = gini_coefficient(grid.a, total_prob)

println("Gini coefficient (Total population): ", gini_total)
