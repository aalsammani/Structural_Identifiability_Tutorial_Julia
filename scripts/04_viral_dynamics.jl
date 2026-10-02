# Section 5.4: within-host viral dynamics.
include("common.jl")
Random.seed!(1)

viral = @ODEmodel(
    T'(t) = s - d_T * T(t) - beta * T(t) * V(t),
    I'(t) = beta * T(t) * V(t) - delta * I(t),
    V'(t) = p * I(t) - c * V(t),
    y(t)  = V(t)
)

header("Viral dynamics, y = V: global identifiability")
report(assess_identifiability(viral; loglevel = SI_LOGLEVEL))

header("Viral dynamics, y = V: identifiable functions")
report_functions(find_identifiable_functions(viral; loglevel = SI_LOGLEVEL))

# The known initial condition must be a state of the same model object, hence the lookup
# by name (a variable T created by a later @ODEmodel call would belong to a different
# polynomial ring and would raise "parents do not match").
T_viral = var(viral, "T")

header("Viral dynamics, y = V, T(0) known (experimental keyword known_ic)")
report(assess_identifiability(viral; known_ic = [T_viral], loglevel = SI_LOGLEVEL))

header("Viral dynamics, y = V, T(0) known: identifiable functions")
report_functions(find_identifiable_functions(viral; known_ic = [T_viral],
                                             loglevel = SI_LOGLEVEL))

# Two observables: viral load and target cells.
viral_two_outputs = @ODEmodel(
    T'(t) = s - d_T * T(t) - beta * T(t) * V(t),
    I'(t) = beta * T(t) * V(t) - delta * I(t),
    V'(t) = p * I(t) - c * V(t),
    y1(t) = V(t),
    y2(t) = T(t)
)

header("Viral dynamics, y1 = V, y2 = T: global identifiability")
report(assess_identifiability(viral_two_outputs; loglevel = SI_LOGLEVEL))

header("Viral dynamics, y1 = V, y2 = T: identifiable functions")
report_functions(find_identifiable_functions(viral_two_outputs; loglevel = SI_LOGLEVEL))
