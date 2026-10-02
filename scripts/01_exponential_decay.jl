# Section 5.1: exponential decay, and the scaled-observation variant used in the text.
include("common.jl")
Random.seed!(1)

decay = @ODEmodel(
    x'(t) = -k * x(t),
    y(t)  =  x(t)
)

header("Exponential decay, y = x: local identifiability")
report(assess_local_identifiability(decay; loglevel = SI_LOGLEVEL))

header("Exponential decay, y = x: global identifiability")
report(assess_identifiability(decay; loglevel = SI_LOGLEVEL))

header("Exponential decay, y = x: identifiable functions")
report_functions(find_identifiable_functions(decay; loglevel = SI_LOGLEVEL))
