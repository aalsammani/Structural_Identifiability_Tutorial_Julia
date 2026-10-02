# Section 5.5: the discrete ambiguity of the SIWR model observed through incidence is not a
# plain exchange of two parameters. Writing z = beta_W*W and kappa = beta_W*xi, the
# identifiable functions are beta_I, gamma + mu, gamma*mu and beta_I*mu + kappa. Hence
#     (gamma, mu, kappa)  ->  (mu, gamma, kappa + beta_I*(mu - gamma))
# preserves all of them. This script constructs the matching initial state and checks by
# simulation that both parameter sets produce the same incidence y = S*(beta_I*I + z).
using OrdinaryDiffEq
include("common.jl")

function siwr_scaled!(du, u, p, t)
    S, I, z = u
    beta_I, gamma, mu, kappa = p
    F = beta_I * I + z                 # force of infection
    du[1] = -S * F
    du[2] =  S * F - gamma * I
    du[3] =  kappa * I - mu * z
end

incidence(sol, t) = (u = sol(t); u[1] * (sol.prob.p[1] * u[2] + u[3]))

"Initial state of the exchanged parameter set that reproduces the output of (p, u0)."
function exchanged_initial_state(p, u0)
    beta_I, gamma, mu, kappa = p
    S0, I0, z0 = u0
    F0 = beta_I * I0 + z0
    dF0 = beta_I * (S0 * F0 - gamma * I0) + (kappa * I0 - mu * z0)
    # With the exchanged parameters (gamma2, mu2, kappa2) = (mu, gamma, kappa + beta_I*(mu-gamma)),
    # I = (F' - beta_I*S*F + mu2*F) / (kappa2 + beta_I*(mu2 - gamma2)) = (F' - beta_I*S*F + gamma*F) / kappa.
    I0b = (dF0 - beta_I * S0 * F0 + gamma * F0) / kappa
    return [S0, I0b, F0 - beta_I * I0b]
end

# beta_I, gamma, mu, kappa; chosen so that all quantities of the exchanged set stay positive
pA = [0.3, 0.5, 0.25, 0.2]
pB = [pA[1], pA[3], pA[2], pA[4] + pA[1] * (pA[3] - pA[2])]
u0A = [0.99, 0.01, 0.02]
u0B = exchanged_initial_state(pA, u0A)

tspan = (0.0, 60.0)
solA = solve(ODEProblem(siwr_scaled!, u0A, tspan, pA), Tsit5(); reltol = 1e-12, abstol = 1e-12)
solB = solve(ODEProblem(siwr_scaled!, u0B, tspan, pB), Tsit5(); reltol = 1e-12, abstol = 1e-12)

ts = range(tspan...; length = 601)
yA = [incidence(solA, t) for t in ts]
yB = [incidence(solB, t) for t in ts]

header("SIWR exchange symmetry checked by simulation")
println("Parameter set A (beta_I, gamma, mu, beta_W*xi): ", pA)
println("Parameter set B (beta_I, gamma, mu, beta_W*xi): ", pB)
println("Initial state A (S, I, beta_W*W): ", u0A)
println("Initial state B (S, I, beta_W*W): ", round.(u0B; sigdigits = 6))
println("max |y_A - y_B| on [0, 60]:  ", maximum(abs.(yA .- yB)))
println("max y_A on [0, 60]:          ", maximum(yA))
println("max |I_A - I_B| on [0, 60]:  ", maximum(abs.([solA(t)[2] - solB(t)[2] for t in ts])))
