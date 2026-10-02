# Section 3.4: SBML file -> SBMLImporter.jl -> Catalyst ReactionSystem -> StructuralIdentifiability.jl.
# The SBML file models/enzyme_kinetics.xml encodes the mechanism E + S <-> C -> E + P.
using SBMLImporter
using Catalyst
include("common.jl")
Random.seed!(1)
println("SBMLImporter.jl ", pkgversion(SBMLImporter), " | Catalyst.jl ", pkgversion(Catalyst))

sbml_file = joinpath(@__DIR__, "..", "models", "enzyme_kinetics.xml")
# SBMLImporter imports stoichiometries symbolically, which Catalyst's conservation-law
# elimination does not support; the analysis is therefore run with remove_conserved = false.
# The compartment volume `cell` is imported as a parameter and is declared known.
rn, callbacks = load_SBML(sbml_file)

println()
println("Imported ReactionSystem:")
println(rn)
println("Reaction-rate equations:")
for eq in equations(ode_model(rn))
    println("  ", eq)
end

header("SBML enzyme model, product P measured: global identifiability")
report(assess_identifiability(rn; measured_quantities = [:P], known_p = [:cell], remove_conserved = false, loglevel = SI_LOGLEVEL))

header("SBML enzyme model, product P measured: identifiable functions")
report_functions(find_identifiable_functions(rn; measured_quantities = [:P], known_p = [:cell], remove_conserved = false,
                                             loglevel = SI_LOGLEVEL))

header("SBML enzyme model, substrate S and product P measured: global identifiability")
report(assess_identifiability(rn; measured_quantities = [:S, :P], known_p = [:cell], remove_conserved = false, loglevel = SI_LOGLEVEL))

header("SBML enzyme model, complex C measured: global identifiability")
report(assess_identifiability(rn; measured_quantities = [:C], known_p = [:cell], remove_conserved = false, loglevel = SI_LOGLEVEL))

header("SBML enzyme model, complex C measured: identifiable functions")
report_functions(find_identifiable_functions(rn; measured_quantities = [:C], known_p = [:cell], remove_conserved = false,
                                             loglevel = SI_LOGLEVEL))
