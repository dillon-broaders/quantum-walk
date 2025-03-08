using Pkg
Pkg.activate(".")
using LinearAlgebra, Plots, StatsBase  

# total steps in quantum walk
nsteps = 200 

# main loop: simulating quantum walk over nsteps
# finding the shannon entropy after each step
steps, shannon = [], [] 
for N in 0:nsteps-1

    # walks spreads left and right at each step giving P possible positions
    P = 2 * N + 1

    # define coin states
    coin0 = [1, 0] 
    coin1 = [0, 1] 

    # coin flip operator C_hat modelled via Hadamard
    C_hat = (1 / sqrt(2)) * [1 1; 1 -1]

    # construct the shift operators
    # forw shifts forward 
    # backw shifts backward
    forw = circshift(Matrix(I, P, P), (1, 0))
    backw = circshift(Matrix(I, P, P), (-1, 0))

    # define conditional shift operator S_hat
    # it shifts forward if coin state is C11
    # It shifts backward if coin state is C00
    C00 = coin0 * coin0'
    C11 = coin1 * coin1'
    S_hat = kron(C11, forw) + kron(C00, backw)
    
    # define full unitary operator W 
    W = S_hat * kron(C_hat, I(P))
    
    # walker starts at position N
    posn0 = zeros(P)
    posn0[N + 1] = 1

    # initial state is superposition of coin states at position N
    psi0 = kron((coin0 + 1im * coin1) / sqrt(2), posn0)

    # state after N steps
    psiN = W^N * psi0
    
    # loop calculates the position probability distribution
    prob = Float64[]
    for k in 1:P
        posn = zeros(P)
        posn[k] = 1
        M_hat_k = kron(I(2), posn * posn')
        proj = M_hat_k * psiN
        push!(prob, real(proj' * proj))
    end
    
    # saving entropy as function of steps walked
    push!(steps, N)
    push!(shannon, entropy(prob, 2))
end

# plotting results
plot(steps, shannon, label="Entropy", xlabel="Steps", ylabel="Randomness", title="Randomness in system as a function of qw steps", color=:blue, frame=:box)
savefig("quantum_walk_entropy.png")

