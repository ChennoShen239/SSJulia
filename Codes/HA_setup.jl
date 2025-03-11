
using Parameters
using Roots

# parameters
@with_kw struct Par 
    α = 0.4
    β = 0.98
    γ = 2.0
    δ = 0.02
    ρ = 0.95
end
par = Par();

# see the website for the whole process 
struct HAGrids
    a::Array{Float64,1}
    na::Int64
    e::Array{Float64,1}
    ne::Int64
    Πe::Array{Float64,2} 
end   

function HAGrids()

    amin = 0.0  # borrowing constraint
    #exogenous transition matrix
    #has form cols = future state, rows = current  state
    Ne = 2
    lamw = 0.6 # job finding prob
    sigma = 0.2 # job separation prob
    exogTrans = [1-lamw  lamw; sigma  1.0-sigma] 
    # a 2-state income process

    #labor endowments
    endow = [1.0;2.5]
    Lbar = sum(endow .* stationarydistribution(exogTrans'))
    endow = endow / Lbar
    Lbar = 1.0

    #grid for savings
    agridmin, agridmax, agridsize = amin, 250, 251
    agrid = range(agridmin^(0.25), agridmax^(0.25), length= agridsize).^4
    # this makes the gird more dense when a is low
    return HAGrids(agrid,agridsize,endow,Ne,exogTrans)
end
grid = HAGrids()

# functions for marginal utility and its inverse
uPrime = (par,c) -> c.^(-par.γ);
uPrimeInv = (par,up) -> up.^(-1.0/par.γ);

"""
    getWR(par, K; Z=1)

Compute the wage and return on capital based on the input parameters and capital level.

# Arguments
- par: A parameter structure that includes model parameters, at least α (capital share) and δ (depreciation rate).
- K  : The level of aggregate capital.
- Z  : Total factor productivity (default is 1).

# Returns
A two-element array [W, R] where:
- W is the wage computed as (1-α)*Z*(K/Lbar)^α,
- R is the return on capital computed as α*Z*(K/Lbar)^(α-1) + 1 - δ.

Note: Lbar is set to 1 within the function.
"""
function getWR(par, K; Z=1)
    @unpack α, δ = par
    Lbar = 1.0
    W = (1-α)*Z.*(K / Lbar).^α
    R = α*Z.*(K / Lbar).^(α-1) .+ 1 .- δ
    return [W, R]
end
#initial guess of consumption function and initialize steady state dictionary
#KCompleteMarkets =  find_zero(K -> par.β * getWR(par,K)[2] - 1, (40,60)) ;
# the equilibirium captial stock in comeplete market that statisfies β(1+r)=1
KCompleteMarkets = fzero(K -> par.β *getWR(par,K)[2]-1,45.63711005516378)

W,R = getWR(par,48);
c0 = (R-1)/R * grid.a .+ W*grid.e';
ss = Dict{Symbol,Any}(:c =>c0)


function aggregator(g,c,D)
    return D' * reshape(g,:);
end


"""
    EGMStepBack

    One step of the endogenous grid method.  Given the marginal value of assets
    Va_p in t+1 and prices at date t Xt, we produce the marginal value of assets at t 
    and the saving and consumption policy rules at t.
"""
function EGMStepBack(Va_p,Xt,par,grid)
    @unpack β, γ = par
    W,R = Xt

    uc = β * Va_p * grid.Πe';  # marginal utility today given Va' on grid of b'  
    # actually we take the expectation using Πe'
    cimplied = uPrimeInv(par, uc);   # convert to consumption
    
    currentassets = (grid.a .+ cimplied .- grid.e'*W)/R; # current assets from budget constraint
    
    g = zeros(size(currentassets)); 
    
    for j = 1:grid.ne
        g[:,j],_ = interp(currentassets[:,j],grid.a,grid.a,right_extrap=true); # interpolate savings rule onto grid
    end
    
    g = max.(g,grid.a[1])  # impose the borrowing constraint
    c = (R*grid.a .+ W*grid.e' .-  g); # consumption from the budget constraint
    Va = R*uPrime(par,c); # envelope condition
    return Va, g, c
end



"""
    checkK!(K, par, grid)

Evaluate the steady state capital accumulation condition by computing the difference between 
the implied capital stock from the endogenous grid method and the current guess K.

This function overwrites the global consumption guess (c0) with the result of SolveEGM. It then 
computes the implied capital stock (K_implied) by aggregating the savings policy and the stationary 
distribution derived from the grid's transition matrix. The returned value (K_implied - K) is used 
to iteratively adjust the guess for the steady state capital stock.

# Arguments
- K: The current capital stock guess.
- par: Model parameters.
- grid: The grid structure used in the endogenous grid method.

# Returns
- The difference between the implied capital stock and the guess K.
"""
function checkK!(K, par, grid)
    _, g, c0[:] = SolveEGM(c0, getWR(par, K), par, grid)
    K_implied = aggregator(g, c0, stationarydistribution(maketransmat(g, grid)))
    return K_implied - K
end


function compute_ss!(ss,par,grid)
    print("Computing steady state...")
    KLarge = 46.5;
    
    ss[:K] = fzero(kk -> checkK!(kk,par,grid),46.5)
   # ss[:K] = find_zero(kk -> checkK!(kk,par,grid), (KCompleteMarkets+0.03, KLarge))
    ss[:X] = getWR(par,ss[:K]);  # store the steady state prices in :X
    ss[:Va], ss[:g], ss[:c] = SolveEGM(ss[:c],ss[:X],par,grid);  #compute the ss decision rules 
    ss[:D] = stationarydistribution(maketransmat(ss[:g],grid)); #compute the ss distribution
    ss[:agg] = aggregator(ss[:g],ss[:c],ss[:D]); # store the ss aggregates
    ss[:Nagg] = length(ss[:agg])  # it is convenient to have the size of the aggregates
    println("done")
end