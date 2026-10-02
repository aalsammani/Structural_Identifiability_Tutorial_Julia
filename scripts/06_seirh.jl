# Section 5.6: SEIR model with hospitalization.
include("common.jl")
Random.seed!(1)

seirh = @ODEmodel(
    S'(t) = -beta * S(t) * I(t),
    E'(t) =  beta * S(t) * I(t) - sigma * E(t),
    I'(t) =  sigma * E(t) - (gamma + eta) * I(t),
    H'(t) =  eta * I(t) - rho * H(t),
    R'(t) =  gamma * I(t) + rho * H(t),
    y1(t) =  eta * I(t)
)

header("SEIR-H, y1 = eta*I: global identifiability")
report(assess_identifiability(seirh; loglevel = SI_LOGLEVEL))

header("SEIR-H, y1 = eta*I: identifiable functions")
report_functions(find_identifiable_functions(seirh; loglevel = SI_LOGLEVEL))

seirh_two = @ODEmodel(
    S'(t) = -beta * S(t) * I(t),
    E'(t) =  beta * S(t) * I(t) - sigma * E(t),
    I'(t) =  sigma * E(t) - (gamma + eta) * I(t),
    H'(t) =  eta * I(t) - rho * H(t),
    R'(t) =  gamma * I(t) + rho * H(t),
    y1(t) =  eta * I(t),
    y2(t) =  beta * S(t) * I(t)
)

header("SEIR-H, y1 = eta*I, y2 = beta*S*I: global identifiability")
report(assess_identifiability(seirh_two; loglevel = SI_LOGLEVEL))

header("SEIR-H, y1 = eta*I, y2 = beta*S*I: identifiable functions")
report_functions(find_identifiable_functions(seirh_two; loglevel = SI_LOGLEVEL))

# Fixing the incubation rate from clinical data (here sigma = 1/5 per day).
seirh_two_sigma = set_parameter_values(seirh_two, Dict(var(seirh_two, "sigma") => 1 // 5))

header("SEIR-H, y1 = eta*I, y2 = beta*S*I, sigma fixed at 1/5: global identifiability")
report(assess_identifiability(seirh_two_sigma; loglevel = SI_LOGLEVEL))

# Hospital occupancy added as a third output.
seirh_three = @ODEmodel(
    S'(t) = -beta * S(t) * I(t),
    E'(t) =  beta * S(t) * I(t) - sigma * E(t),
    I'(t) =  sigma * E(t) - (gamma + eta) * I(t),
    H'(t) =  eta * I(t) - rho * H(t),
    R'(t) =  gamma * I(t) + rho * H(t),
    y1(t) =  eta * I(t),
    y2(t) =  beta * S(t) * I(t),
    y3(t) =  H(t)
)

header("SEIR-H, y1 = eta*I, y2 = beta*S*I, y3 = H: global identifiability")
report(assess_identifiability(seirh_three; loglevel = SI_LOGLEVEL))
