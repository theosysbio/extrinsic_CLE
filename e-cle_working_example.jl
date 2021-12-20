################################################################################
# THIS SCRIPT SOLVES THE NON-LINEAR MODEL: NEGATIVE FEEDBACK LOOP WITH TIME-
# VARYING TRANSCRIPTION RATE K(t), WHERE K(t) VARIES ACCORDING TO A CIR PROCESS
# BY WAY OF THE E-CLE AND EXTRANDE
################################################################################
# Import packages
using DifferentialEquations, Distributions, Plots, Colors
# Include functions
in_path = "/Users/mcoomer/phd/papers/paper3/code/"
include(in_path * "e-cle_functions.jl")
################################################################################
# Declare exogenous process (I) model & simulation parameters
κ,a,Θ = (1.25,1.0,10)          # Parameters used in Fig.5(d.) where a varies
u0 = complex(1.0)               # Initial # of RNA polymerase molecules
dt = 0.01                       # Time discritization
tspan_I = 7*10^5.0              # Timespan of problem
# Declare primary process model parameters
k2,k3,k4 = (0.0001,1,1)
u1 = complex(0.0)               # Initial # of unbound molecules
u2 = complex(10.0)              # Initial # of proteins
# Declare E-CLE simulation parameters
dt = 0.01                       # Time discritization
tspan = 3*10^5                  # Timespan of problem
burn = 2000                     # Burin-in period
Scle = 1000                     # Sample every Scle time points
# Declare Extrande simulation parameters
T = 6*10^5.0                    # Timespan of problem
Sext = 50                       # Sample every Sext time points
################################################################################
# 1. Generate the exogenouse noise process (I)
I = abs.(last.(solve_CIR(u0,dt,tspan_I,κ,a,Θ)))
# 2. Solve the E-CLE where the first vector of the solution is u1 (# unbound
# molecules over time) and the second vector u2 (# of protein molecues over time).
# Values are in complex space since we employ the complex CLE.
sol_ecle = solve_e_cle(u1,u2,dt,tspan,I,k2,k3,k4)
# Run burn-in and collect every Scle protein samples
protein_ecle = map(x->x[2],sol_ecle[2][burn:Scle:end])
# Compute the first two moments
mean_protein_ecle = abs.(sum(protein_ecle)/length(protein_ecle))
var_protein_ecle = abs.(sum(((protein_ecle).-(mean(protein_ecle))).^2)/length(protein_ecle))
# # Plot the E-CLE distribution
# histogram(protein_ecle,normed=true,label="E-CLE",line=(1),color=:orangered2,
#     xlims=(0,100), xlabel="Protein", ylabel="Probability p(X2)")

# 3. Solve the same system using the Extrande method
sol_ext = extrande(burn,Sext,k2,k3,k4,T,I)
protein_ext = Int.(sol_ext[1])
# Compute the first two moments
mean_protein_ext = mean(protein_ext)
var_protein_ext = var(protein_ext)
# # Normalise data for plotting
# samplesize=length(protein_ext)
# max_copy_number = 100                       # sets max copy number to be plotted
# copylist1 = [0 for i=1:max_copy_number]      # for use in scaling
# for n in 1:samplesize
#     cn = protein_ext[n] + 1
#     if cn <=max_copy_number
#         copylist1[cn] += 1.0
#     end
# end
# copylist1 = copylist1 / samplesize
# sum(copylist1)
# n=0:1:100
# scatter!(copylist1,color=:darkturquoise)
################################################################################
