# Section 5.2: SIR model with incidence observed.
include("common.jl")
Random.seed!(1)

sir = @ODEmodel(
    S'(t) = -beta * S(t) * I(t),
    I'(t) =  beta * S(t) * I(t) - gamma * I(t),
    R'(t) =  gamma * I(t),
    y(t)  =  beta * S(t) * I(t)
)

header("SIR, y = beta*S*I: global identifiability")
report(assess_identifiability(sir; loglevel = SI_LOGLEVEL))

header("SIR, y = beta*S*I: identifiable functions")
report_functions(find_identifiable_functions(sir; loglevel = SI_LOGLEVEL))
