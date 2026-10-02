# Shared helpers for the companion scripts of
# "A Tutorial on Symbolic Structural Identifiability Analysis of ODE Models in Julia".
#
# Every script activates the environment in the parent directory (Project.toml + Manifest.toml),
# so it can be run from any working directory with
#     julia --project=.. scripts/<name>.jl      (from inside scripts/)
# or  julia --project=. scripts/<name>.jl        (from the package root).

using StructuralIdentifiability
using Logging
using Random

# Warnings (for example the Wronskian warning about multi-experiment identifiability) are
# still printed; the routine progress messages of StructuralIdentifiability.jl are suppressed.
const SI_LOGLEVEL = Logging.Warn

"Print a section header."
function header(title)
    println()
    println("="^72)
    println(title)
    println("="^72)
end

"Print an identifiability or observability dictionary, one entry per line."
function report(result)
    for (k, v) in result
        println("  ", rpad(string(k), 28), " => ", v)
    end
end

"Print a vector of identifiable functions, one per line."
function report_functions(funcs)
    for f in funcs
        println("  ", f)
    end
end

"Return the parameter or state of `ode` whose name is `name` (a String); states may be given as `\"T\"` or `\"T(t)\"`."
function var(ode, name::AbstractString)
    matches = filter(v -> string(v) in (name, name * "(t)"), vcat(ode.parameters, ode.x_vars))
    isempty(matches) && error("no parameter or state named $name")
    return only(matches)
end

"Return `(value, seconds)` for a timed call, printing nothing."
timed(f) = (t0 = time(); v = f(); (v, time() - t0))

println("Julia ", VERSION, " | StructuralIdentifiability.jl ",
        pkgversion(StructuralIdentifiability))
