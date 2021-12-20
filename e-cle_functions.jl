################################################################################
# THIS SCRIPT CONTAINS THE NECESSARY FUNCTIONS TO RUN THE FILE:
# "e-cle_working_example.jl"
################################################################################

# 1. The E-CLE describing the "Non-linear model: negative feedback loop with
# time-varying transcription rate".  See eq. 21 of the article.

# u1 = # molecules in unbound state; u2 = # protein molecules; k1 = time-varying
# transcription rate K(t); k2 = sb; k3 = su; k4 = δ
function neg_feedback1(u1,u2,k1,k2,k3,k4)
    return -k2*u1*u2 + k3*(1-u1)
end
function neg_feedback2(u1,u2,k1,k2,k3,k4)
    return k1*u1 - k2*u1*u2 + k3*(1-u1) - k4*u2
end

function lang_noise1(u1,u2,k1,k2,k3,k4)
    return sqrt(k1*u1)
end
function lang_noise2(u1,u2,k1,k2,k3,k4)
    return -sqrt(k2*u1*u2)
end
function lang_noise3(u1,u2,k1,k2,k3,k4)
    return sqrt(k3*(1-u1))
end
function lang_noise4(u1,u2,k1,k2,k3,k4)
    return -sqrt(k4*u2)
end

# 2. Custom solver for the E-CLE (Euler-Maruyama method)

# u1 = # molecules in unbound state; u2 = # protein molecules;
# dt = time discretization; tspan = timespan for the problem, I = exogenouse
# input process; k2 = sb; k3 = su; k4 = δ
function solve_e_cle(u1,u2,dt,tspan,I,k2,k3,k4)
    dW(dt) = rand(Normal(0,sqrt(dt)))
    gene = Vector{Tuple{Float64,Complex{Float64}}}(undef,length(collect(0:dt:tspan)))
    protein = Vector{Tuple{Float64,Complex{Float64}}}(undef,length(collect(0:dt:tspan)))
    iter = 0
    @inbounds for (i,t) in enumerate(0:dt:tspan)
        iter+=1
        K_new=I[iter]
        noise1=lang_noise1(u1,u2,K_new,k2,k3,k4)*dW(dt)
        noise2=lang_noise2(u1,u2,K_new,k2,k3,k4)*dW(dt)
        noise3=lang_noise3(u1,u2,K_new,k2,k3,k4)*dW(dt)
        noise4=lang_noise4(u1,u2,K_new,k2,k3,k4)*dW(dt)
        u1=u1+neg_feedback1(u1,u2,K_new,k2,k3,k4)*dt+noise2+noise3
        u2=u2+neg_feedback2(u1,u2,K_new,k2,k3,k4)*dt+noise1+noise2+noise3+noise4
        gene[i] = (t,u1)
        protein[i] = (t,u2)
    end
    return gene, protein
end

# 3. Extrande algorithm [Voliotis,2016] adpated to The Julia Programming Language

# burn = burn-in period; S = collect every S samples to avoid sample correlation;
# T = timespan for the problem, I = exogenouse input process; k2 = sb; k3 = su;
# k4 = δ
function extrande(burn,S,k2,k3,k4,T,I)
    samples=[]
    # 1: INITIALIZE TIME AND NETWORK STATE
    t = 0
    s = [0.0 ; 1.0]
    stoichiometry = [0 -1 1 0 ; 1 -1 1 -1]
    u = rand()
    max = maximum(I)
    n_old = s[2]
    count = 0
    while t < T
        while t >= S*burn                # only collect every S samples to avoid sample correlation
            burn += 1
            push!(samples,n_old)
        end
        n_old = s[2]
        count += 1
    # 3: DETERMINE PROPENSITY BOUND = UPPER BOUND ON TOTAL PROPENSITY
        L = T-t
        B = sum([max*s[1],s[1]*s[2]*k2,k3*(1-s[1]),s[2]*k4])
    # 4: GENERATE REACTION TIME
        τ = rand(Exponential(1/B))
    # 5: REJECT - STATE OF THE NETWORK REMAINS UNCHANGED
        if τ > L
    # 6: UPDATE TIME - TIME ADVANCES BY L; NO UPDATE TO STATE OF NETWORK
            t = t + L
    # 7: ELSE
        else
    # 8: UPDATE TIME - TIME ADVANCES BY τ
            t = t + τ
    # 9: UPDATE PROPENSITIES THAT DEPEND ON I(t) (K1)
            index=Int(round((t*100)+1))
            a = (I[index]*s[1],(s[1]*s[2]*k2),k3*(1-s[1]),(s[2]*k4))
            propensitySum = sum(a)
    # 10: GENEARTE UNIFORMALLY DISTRIBUTED RANDOM VARIABLE AND MULTIPLY BY B
            u = rand()
            r = u*B
    # 11: ACCEPT IF CONDITION IS MET
            if propensitySum >= r
    # 12: PICK REACTION TO FIRE
                reaction = choose_reaction(a,u,B)
    # 13: UPDATE STATE: ON/OFF AND # PROTEINS
                s = s + stoichiometry[:,reaction]
            end
        end
    end
    return samples, count
end

function choose_reaction(a, u, B)
    if sum(a) <= B*u
        return "error"
    else
        j=1
        jthsum = a[1]
        while jthsum < B*u
            j += 1
            jthsum += a[j]
        end
        return j
    end
end

# 4. SDE describing the exogenous input process I.  Here we use the CIR model of
# transcription that evolves according to eq. 22 of the article
function drift(u,κ,a,Θ)
    return a*Θ - κ*u
end

function diffusion(u,κ,Θ)
    dW(dt)=rand(Normal(0,sqrt(dt)))
    return sqrt(2*κ*Θ*u)*dW(dt)
end

# 5. Custom solver for the CIR process (Euler-Maruyama method)

# u = initial # RNA polymerase molecules; dt = time discretization; tspan =
# timespan for the problem; κ = degradation rate, a = production rate; θ = scale
# parameter.  See article for detailed description of model parameters.
function solve_CIR(u,dt,tspan,κ,a,Θ)
    RNA = Vector{Tuple{Float64,Complex{Float64}}}(undef,length(collect(0:dt:tspan)))
    @inbounds for (i,t) in enumerate(0:dt:tspan)
        u=u+drift(u,κ,a,Θ)*dt+diffusion(u,κ,Θ)
        RNA[i] = (t,u)
    end
    return RNA
end
################################################################################
