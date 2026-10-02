# Wall-clock times of the local and global tests for the case-study models.
# Each call is run once to trigger compilation and then timed on a second run.
# Times depend on the machine and on Julia's compilation state; they are reported only to
# illustrate relative costs, not as a benchmark.
include("common.jl")

models = [
    ("Exponential decay", @ODEmodel(x'(t) = -k * x(t), y(t) = x(t))),
    ("SIR, y = beta*S*I", @ODEmodel(
        S'(t) = -beta * S(t) * I(t),
        I'(t) =  beta * S(t) * I(t) - gamma * I(t),
        R'(t) =  gamma * I(t),
        y(t)  =  beta * S(t) * I(t))),
    ("PK, y = x2", @ODEmodel(
        x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
        x2'(t) =  a21 * x1(t) - a12 * x2(t),
        y(t)   =  x2(t))),
    ("Viral dynamics, y = V", @ODEmodel(
        T'(t) = s - d_T * T(t) - beta * T(t) * V(t),
        I'(t) = beta * T(t) * V(t) - delta * I(t),
        V'(t) = p * I(t) - c * V(t),
        y(t)  = V(t))),
    ("SIWR, y = incidence", @ODEmodel(
        S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
        I'(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
        W'(t) =  xi * I(t) - mu * W(t),
        R'(t) =  gamma * I(t),
        y(t)  =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t))),
    ("SEIR-H, y1 = eta*I, y2 = beta*S*I", @ODEmodel(
        S'(t) = -beta * S(t) * I(t),
        E'(t) =  beta * S(t) * I(t) - sigma * E(t),
        I'(t) =  sigma * E(t) - (gamma + eta) * I(t),
        H'(t) =  eta * I(t) - rho * H(t),
        R'(t) =  gamma * I(t) + rho * H(t),
        y1(t) =  eta * I(t),
        y2(t) =  beta * S(t) * I(t))),
]

header("Timing of assess_local_identifiability and assess_identifiability (seconds)")
println(rpad("Model", 38), rpad("local", 12), "global")
for (name, ode) in models
    assess_local_identifiability(ode; loglevel = SI_LOGLEVEL)
    assess_identifiability(ode; loglevel = SI_LOGLEVEL)
    _, t_loc = timed(() -> assess_local_identifiability(ode; loglevel = SI_LOGLEVEL))
    _, t_glob = timed(() -> assess_identifiability(ode; loglevel = SI_LOGLEVEL))
    println(rpad(name, 38), rpad(string(round(t_loc; sigdigits = 2)), 12),
            round(t_glob; sigdigits = 2))
end
println()
println("CPU threads: ", Sys.CPU_THREADS, " | ", Sys.cpu_info()[1].model)
