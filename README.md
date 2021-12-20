# extrinsic_CLE
Extrinsic chemical Langevin Equation for biochemical systems in dynamic environments

# Description
The extrinsic chemical Langevin equation (E-CLE) is a computationally efficient simulation framework for modelling inherently stochastic biochemical reaction networks embedded in fluctuationg environments, or then biochemical networks subject to both intrinsic and extrinsic noise.  The code is written and implemented in The Julia Programming Language.

The E-CLE is a physically motivated extension of the chemical Langevin equation.  Specifically, it is a continuous approximation to the Chemical Master Equation (CME) with time-varying propensities. Noise is incorporated at the level of the CME, and can account for the full dynamics of the exogenous noise process, irrespective of timescales and their mismatches.  

We use the Euler-Maruyama (EM) method for simulating the E-CLE.  Other implementations are  possible, however, we focus on the EM method as it is the standard method for simulating the CLE.  The implementation takes as an input a pre-simulated time series of the exogenous process (I), typically either (i) an approximate simulation of a stochastic differential equation or (ii) an exact simulation of a chemical reaction network.  

# Examples
The script [e-cle_functions.jl](https://github.com/theosysbio/extrinsic_CLE/blob/main/e-cle_functions.jl) contains the core functions to simulate the E-CLE and the script [e-cle_working_example.jl](https://github.com/theosysbio/extrinsic_CLE/blob/main/e-cle_working_example.jl) implements a minimal working example of a non-linear gene regulatory model where the transcription rate is varied according to an exogenous Cox-Ingersoll-Ross process.

# Literature
For more details on the method and additional examples please refer to the associated paper: [The chemical Langevin equation for biochemical systems in dynamic environments]()
