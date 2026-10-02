# Section 2.6: coefficients of input-output equations need not be identifiable.
# Example 2.14 of Hong, Ovchinnikov, Pogudin and Yap (CPAM 73:1831-1879, 2020):
#     x' = 0,  y1 = x,  y2 = a*x + a^2.
# The input-output equations are y1' = 0 and y2 - a*y1 - a^2 = 0, whose coefficients
# include a. From a single experiment, however, a solves a^2 + x(0)*a - y2 = 0, which has
# two roots, so a is locally but not globally identifiable.
include("common.jl")
Random.seed!(1)

ex214 = @ODEmodel(
    x'(t)  = 0,
    y1(t)  = x(t),
    y2(t)  = a * x(t) + a^2
)

header("Hong et al. Example 2.14: assess_identifiability (warnings shown)")
report(assess_identifiability(ex214; loglevel = SI_LOGLEVEL))

header("Hong et al. Example 2.14: find_identifiable_functions (warnings shown)")
report_functions(find_identifiable_functions(ex214; loglevel = SI_LOGLEVEL))

header("Hong et al. Example 2.14: local test, single-experiment (:SE) and multi-experiment (:ME)")
report(assess_local_identifiability(ex214; funcs_to_check = [a], type = :SE,
                                    loglevel = SI_LOGLEVEL))
# With type = :ME the function returns a tuple; its second entry is the replica counter num_exp.
me_result, n_experiments = assess_local_identifiability(ex214; funcs_to_check = [a], type = :ME,
                                                        loglevel = SI_LOGLEVEL)
report(me_result)
println("  second return value (num_exp): ", n_experiments)
