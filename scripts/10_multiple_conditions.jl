# Section 7.7 (limitations): pooling several experimental conditions by model augmentation.
# Each condition receives its own copy of the states (with its own unknown initial
# conditions); parameters that are shared across conditions appear once, and
# condition-specific parameters receive their own names.
include("common.jl")
Random.seed!(1)

# Two repeated experiments with the same design (all parameters shared).
pk_repeat = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u1(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    z1'(t) = -(a01 + a21) * z1(t) + a12 * z2(t) + u2(t),
    z2'(t) =  a21 * z1(t) - a12 * z2(t),
    yA(t)  =  x2(t),
    yB(t)  =  z2(t)
)

header("PK, two conditions, all parameters shared: global identifiability")
report(assess_identifiability(pk_repeat; funcs_to_check = [a01, a12, a21],
                              loglevel = SI_LOGLEVEL))

# Second condition with a condition-specific elimination rate b01
# (for example, the same subject studied with and without an elimination inhibitor).
pk_two_conditions = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u1(t),
    x2'(t) =  a21 * x1(t) - a12 * x2(t),
    z1'(t) = -(b01 + a21) * z1(t) + a12 * z2(t) + u2(t),
    z2'(t) =  a21 * z1(t) - a12 * z2(t),
    yA(t)  =  x2(t),
    yB(t)  =  z2(t)
)

header("PK, two conditions, condition-specific elimination b01: global identifiability")
report(assess_identifiability(pk_two_conditions; funcs_to_check = [a01, b01, a12, a21],
                              loglevel = SI_LOGLEVEL))

header("PK, two conditions, condition-specific elimination b01: identifiable functions")
report_functions(find_identifiable_functions(pk_two_conditions; loglevel = SI_LOGLEVEL))
