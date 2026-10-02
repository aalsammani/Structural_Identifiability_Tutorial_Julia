# Run every companion script in a fresh Julia process and store its output in
# expected_output/<script>.txt. Usage, from the package root:
#     julia --project=. run_all.jl
using Pkg
Pkg.instantiate()

const ROOT = @__DIR__
const SCRIPTS = filter(f -> occursin(r"^\d\d_.*\.jl$", f), readdir(joinpath(ROOT, "scripts")))
const OUTDIR = joinpath(ROOT, "expected_output")
mkpath(OUTDIR)

for script in SCRIPTS
    out = joinpath(OUTDIR, replace(script, ".jl" => ".txt"))
    cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no $(joinpath(ROOT, "scripts", script))`
    print(rpad(script, 36))
    t = @elapsed ok = open(out, "w") do io
        success(pipeline(cmd; stdout = io, stderr = io))
    end
    println(ok ? "ok" : "FAILED", "  (", round(t; digits = 1), " s)")
end
