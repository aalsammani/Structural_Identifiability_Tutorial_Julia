# Section 5.3: two-compartment pharmacokinetic model.
include("common.jl")
Random.seed!(1)

# Peripheral compartment observed, known input u(t) into the central compartment.
pk = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    y(t)   =  x2(t)
)

header("PK, y = x2, input u: global identifiability")
report(assess_identifiability(pk; loglevel = SI_LOGLEVEL))

header("PK, y = x2: symmetric combinations")
report(assess_identifiability(pk; funcs_to_check = [a01 + a12, a01 * a12],
                              loglevel = SI_LOGLEVEL))

header("PK, y = x2: identifiable functions")
report_functions(find_identifiable_functions(pk; loglevel = SI_LOGLEVEL))

# Central compartment observed, input with unknown gain b.
pk_central = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + b * u(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    y(t)   =  x1(t)
)

header("PK, y = x1, input b*u: global identifiability")
report(assess_identifiability(pk_central; loglevel = SI_LOGLEVEL))

header("PK, y = x1, input b*u: identifiable functions")
report_functions(find_identifiable_functions(pk_central; loglevel = SI_LOGLEVEL))

# No input: peripheral and central compartment observed in turn.
pk_noinput_x2 = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    y(t)   =  x2(t)
)
pk_noinput_x1 = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    y(t)   =  x1(t)
)

header("PK, no input, y = x2: global identifiability and identifiable functions")
report(assess_identifiability(pk_noinput_x2; loglevel = SI_LOGLEVEL))
report_functions(find_identifiable_functions(pk_noinput_x2; loglevel = SI_LOGLEVEL))

header("PK, no input, y = x1: global identifiability and identifiable functions")
report(assess_identifiability(pk_noinput_x1; loglevel = SI_LOGLEVEL))
report_functions(find_identifiable_functions(pk_noinput_x1; loglevel = SI_LOGLEVEL))
