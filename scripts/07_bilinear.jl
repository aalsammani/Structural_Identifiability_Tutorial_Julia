# Section 6.1: the bilinear toy model dx/dt = -p*q*x.
include("common.jl")
Random.seed!(1)

bilinear = @ODEmodel(
    x'(t) = -p * q * x(t),
    y(t)  =  x(t)
)

header("Bilinear model, y = x: global identifiability")
report(assess_identifiability(bilinear; loglevel = SI_LOGLEVEL))

header("Bilinear model, y = x: identifiable functions")
report_functions(find_identifiable_functions(bilinear; loglevel = SI_LOGLEVEL))
