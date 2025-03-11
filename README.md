# HetAgentsv2

This repository contains a collection of heterogeneous agent macroeconomic models implemented in Julia, designed to accompany the [notes on heterogeneous agent macro](https://alisdairmckay.com/Notes/HetAgentsV2).

## Project Overview

HetAgentsv2 provides Julia implementations of several important macroeconomic models, including:

- **Aiyagari Model**: Studies savings behavior and wealth distribution in incomplete markets
- **HANK Model**: Heterogeneous Agent New Keynesian model, used to analyze monetary policy effects with household heterogeneity
- **RBC Model**: Real Business Cycle model, examines economic fluctuations
- **MIT Shock Analysis**: Studies dynamic effects of policy shocks

## Project Structure

The project is divided into two main parts:

### Notes/
Contains theoretical notes and mathematical derivations for each model:
- `Aiyagari.md` - Theoretical foundations of the Aiyagari model
- `HANK.md` - Heterogeneous Agent New Keynesian model
- `RBC.md` - Real Business Cycle model
- `EGM.md` - Endogenous Grid Method
- `MIT.md` - MIT shock analysis
- `SSJac.md` - Steady state and Jacobian calculation methods

### Codes/
Contains Julia implementations of all models:
- **Core Files**:
  - `HA_setup.jl` - Basic setup for heterogeneous agent models
  - `HA_fcns.jl` - Core functions for heterogeneous agent models
  - `HANK_setup.jl` - Basic setup for HANK models

- **Model Implementations**:
  - `Aiyagari_script.jl` - Aiyagari model
  - `HANK_script.jl` - HANK model
  - `RBC.jl` - Real Business Cycle model
  - `MIT.jl` - MIT shock related code

- **Method Implementations**:
  - `EGM_script.jl` - Endogenous Grid Method
  - `SSJac.jl` - Steady state and Jacobian calculations
  - `Transition_SSJ.jl` - Code for handling transition dynamics

- **Auxiliary Tools**:
  - `_aux/ModelUtils.jl` - Model utility functions

## Usage

### Environment Setup
Ensure you have Julia installed (version 1.6 or higher recommended) and the necessary dependency packages:

```julia
using Pkg
Pkg.add(["Parameters", "Plots", "ForwardDiff", "LinearAlgebra", "Roots", "SparseArrays", "Symbolics"])
```

### Running Models
In the Julia REPL, you can run each model as follows:

```julia
# Run Aiyagari model
include("Codes/Aiyagari_script.jl")

# Run RBC model
include("Codes/RBC.jl")

# Run HANK model
include("Codes/HANK_script.jl")
```

### Customizing Parameters
The parameter structures in each model file can be modified as needed, for example:

```julia
@with_kw struct Par 
    α = 0.4    # Capital share
    β = 0.98   # Discount factor
    γ = 2.0    # Risk aversion coefficient
    δ = 0.02   # Capital depreciation rate
    ρ = 0.95   # Exogenous shock persistence
end
```

## Core Features

- **Steady State Solutions**: Calculate steady-state equilibria for various models
- **Policy Function Calculation**: Solve household optimal decisions using the Endogenous Grid Method (EGM)
- **Transition Matrix Construction**: Construct asset distribution transition matrices using the `maketransmat` function
- **Impulse Response Analysis**: Calculate dynamic responses of the economy to exogenous shocks
- **Visualization**: Plot policy functions, asset distributions, impulse responses, etc.

## Contributions and Feedback

Feedback and contributions via Issues or Pull Requests are welcome. For any questions, please refer to the related theoretical notes or contact the project maintainers.

## License

This project is licensed under the MIT License. 
