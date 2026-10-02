# Section 7.7 (limitations): rationalizing a non-rational model by adding a state.
# Generalized growth model C' = r*C^alpha with incidence y = C' observed
# (Liyanage, Chowell, Pogudin, Tuncer, Viruses 17:496, 2025, Section 3).
include("common.jl")
Random.seed!(1)

# Lifting 1: x = r*C^alpha. Then C' = x and x' = alpha*x^2/C, with y = x.
ggm_lift1 = @ODEmodel(
    C'(t) = x(t),
    x'(t) = alpha * x(t)^2 / C(t),
    y(t)  = x(t)
)

header("GGM, lifting x = r*C^alpha: global identifiability")
report(assess_identifiability(ggm_lift1; loglevel = SI_LOGLEVEL))

# r = x/C^alpha is not rational in (x, C, alpha); its identifiability follows from that of
# alpha, C and x, which are reported above.

# Lifting 2: w = C^alpha. Then C' = r*w and w' = alpha*r*w^2/C, with y = r*w.
ggm_lift2 = @ODEmodel(
    C'(t) = r * w(t),
    w'(t) = alpha * r * w(t)^2 / C(t),
    y(t)  = r * w(t)
)

header("GGM, lifting w = C^alpha: global identifiability")
report(assess_identifiability(ggm_lift2; loglevel = SI_LOGLEVEL))

header("GGM, lifting w = C^alpha: identifiable functions (with states)")
report_functions(find_identifiable_functions(ggm_lift2; with_states = true,
                                             loglevel = SI_LOGLEVEL))
