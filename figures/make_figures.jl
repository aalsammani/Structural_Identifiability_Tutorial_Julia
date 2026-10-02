# Regenerates the data figures of the manuscript as PDF files:
#   fig_verdict_matrix.pdf    (verdict matrix; computed here with StructuralIdentifiability.jl)
#   fig_nonidentifiability.pdf (numerical demonstrations; computed here with OrdinaryDiffEq.jl)
# Usage, from the package root:  julia --project=. figures/make_figures.jl
using StructuralIdentifiability
using OrdinaryDiffEq
using Plots
using LaTeXStrings
using Logging
using Random

gr()
default(; fontfamily = "Computer Modern")
const OUT = @__DIR__
Random.seed!(1)

# ---------------------------------------------------------------------------------------
# Verdict matrix
# ---------------------------------------------------------------------------------------
rows = [
    (L"\mathrm{Exp.\ decay},\ y=x", @ODEmodel(x'(t) = -k * x(t), y(t) = x(t))),
    (L"\mathrm{SIR},\ y=\beta SI", @ODEmodel(
        S'(t) = -beta * S(t) * I(t), I'(t) = beta * S(t) * I(t) - gamma * I(t),
        R'(t) = gamma * I(t), y(t) = beta * S(t) * I(t))),
    (L"\mathrm{PK},\ y=x_2", @ODEmodel(
        x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
        x2'(t) = a21 * x1(t) - a12 * x2(t), y(t) = x2(t))),
    (L"\mathrm{PK},\ y=x_1,\ \mathrm{input}\ bu", @ODEmodel(
        x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + b * u(t),
        x2'(t) = a21 * x1(t) - a12 * x2(t), y(t) = x1(t))),
    (L"\mathrm{Viral},\ y=V", @ODEmodel(
        T'(t) = s - d_T * T(t) - beta * T(t) * V(t), I'(t) = beta * T(t) * V(t) - delta * I(t),
        V'(t) = p * I(t) - c * V(t), y(t) = V(t))),
    (L"\mathrm{Viral},\ y=(V,T)", @ODEmodel(
        T'(t) = s - d_T * T(t) - beta * T(t) * V(t), I'(t) = beta * T(t) * V(t) - delta * I(t),
        V'(t) = p * I(t) - c * V(t), y1(t) = V(t), y2(t) = T(t))),
    (L"\mathrm{SIWR},\ y=\mathrm{incidence}", @ODEmodel(
        S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
        I'(t) = beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
        W'(t) = xi * I(t) - mu * W(t), R'(t) = gamma * I(t),
        y(t) = beta_I * S(t) * I(t) + beta_W * S(t) * W(t))),
    (L"\mathrm{SIWR},\ y=(\mathrm{incidence},W)", @ODEmodel(
        S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
        I'(t) = beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
        W'(t) = xi * I(t) - mu * W(t), R'(t) = gamma * I(t),
        y1(t) = beta_I * S(t) * I(t) + beta_W * S(t) * W(t), y2(t) = W(t))),
    (L"\mathrm{SIWR},\ y=I", @ODEmodel(
        S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
        I'(t) = beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
        W'(t) = xi * I(t) - mu * W(t), R'(t) = gamma * I(t), y(t) = I(t))),
    (L"\mathrm{SEIR\!-\!H},\ y=\eta I", @ODEmodel(
        S'(t) = -beta * S(t) * I(t), E'(t) = beta * S(t) * I(t) - sigma * E(t),
        I'(t) = sigma * E(t) - (gamma + eta) * I(t), H'(t) = eta * I(t) - rho * H(t),
        R'(t) = gamma * I(t) + rho * H(t), y1(t) = eta * I(t))),
    (L"\mathrm{SEIR\!-\!H},\ y=(\eta I,\beta SI)", @ODEmodel(
        S'(t) = -beta * S(t) * I(t), E'(t) = beta * S(t) * I(t) - sigma * E(t),
        I'(t) = sigma * E(t) - (gamma + eta) * I(t), H'(t) = eta * I(t) - rho * H(t),
        R'(t) = gamma * I(t) + rho * H(t), y1(t) = eta * I(t), y2(t) = beta * S(t) * I(t))),
    (L"\mathrm{Bilinear},\ y=x", @ODEmodel(x'(t) = -p * q * x(t), y(t) = x(t))),
]

columns = ["k", "beta", "gamma", "a01", "a12", "a21", "b", "p", "delta", "c", "s", "d_T",
           "beta_I", "beta_W", "xi", "mu", "sigma", "eta", "rho", "q"]
column_labels = [L"k", L"\beta", L"\gamma", L"a_{01}", L"a_{12}", L"a_{21}", L"b", L"p",
                 L"\delta", L"c", L"s", L"d_T", L"\beta_I", L"\beta_W", L"\xi", L"\mu",
                 L"\sigma", L"\eta", L"\rho", L"q"]

code = Dict(:globally => 3, :locally => 2, :nonidentifiable => 1)   # 0 = absent
M = zeros(Int, length(rows), length(columns))
for (i, (_, ode)) in enumerate(rows)
    res = assess_identifiability(ode; funcs_to_check = ode.parameters, loglevel = Logging.Warn)
    for (par, verdict) in res
        j = findfirst(==(string(par)), columns)
        j === nothing && error("parameter $par missing from the column list")
        M[i, j] = code[verdict]
    end
end
println("Verdict matrix (0 absent, 1 NI, 2 L, 3 G):")
display(M)

colors = [RGB(0.87, 0.87, 0.87), RGB(0.75, 0.22, 0.17), RGB(0.85, 0.65, 0.05), RGB(0.18, 0.55, 0.34)]
letters = Dict(0 => "–", 1 => "NI", 2 => "L", 3 => "G")
nr, nc = size(M)
pm = plot(; xlim = (0.5, nc + 0.5), ylim = (0.5, nr + 0.5), yflip = true, legend = false,
          framestyle = :none, size = (1000, 500), left_margin = 48Plots.mm,
          top_margin = 4Plots.mm, bottom_margin = 1Plots.mm)
for i in 1:nr, j in 1:nc
    c = colors[M[i, j] + 1]
    plot!(pm, Shape([j - 0.46, j + 0.46, j + 0.46, j - 0.46], [i - 0.44, i - 0.44, i + 0.44, i + 0.44]);
          fillcolor = c, linecolor = c)
    annotate!(pm, j, i, text(letters[M[i, j]], M[i, j] == 0 ? 9 : 10,
                             M[i, j] == 0 ? :gray40 : :white, :center))
end
for (j, lab) in enumerate(column_labels)
    annotate!(pm, j, 0.1, text(lab, 12, :black, :center))
end
for (i, (lab, _)) in enumerate(rows)
    annotate!(pm, 0.3, i, text(lab, 11, :black, :right))
end
legend_items = [(3, "G  globally identifiable", 1.0), (2, "L  locally only", 7.6),
                (1, "NI  non-identifiable", 11.8), (0, "–  absent", 17.2)]
for (v, lab, x0) in legend_items
    plot!(pm, Shape([x0 - 0.35, x0 + 0.35, x0 + 0.35, x0 - 0.35],
                    [nr + 1.05, nr + 1.05, nr + 1.65, nr + 1.65]); fillcolor = colors[v + 1],
          linecolor = :black)
    annotate!(pm, x0 + 0.6, nr + 1.35, text(lab, 11, :black, :left))
end
plot!(pm; ylim = (-0.3, nr + 1.8), xlim = (0.5, nc + 0.6))
savefig(pm, joinpath(OUT, "fig_verdict_matrix.pdf"))
println("wrote fig_verdict_matrix.pdf")

# ---------------------------------------------------------------------------------------
# Numerical demonstrations
# ---------------------------------------------------------------------------------------
default(; fontfamily = "Computer Modern", titlefontsize = 11, guidefontsize = 11,
        tickfontsize = 9, legendfontsize = 8, linewidth = 2.2, framestyle = :box)

# (A) bilinear model: four (p, q) pairs with pq = 6
tt = range(0, 2; length = 300)
pA = plot(; xlabel = L"t", ylabel = L"x(t)", title = L"(A)\ \dot x=-pqx,\ pq=6")
for ((p, q), ls) in zip([(2, 3), (6, 1), (0.5, 12), (1.5, 4)], [:solid, :dash, :dot, :dashdot])
    plot!(pA, tt, exp.(-p * q .* tt); linestyle = ls, label = latexstring("(p,q)=($p,$q)"))
end

# (B) two-compartment model, a01 <-> a12 exchange, bolus x1(0) = 10, x2(0) = 0
function pk!(du, u, p, t)
    a01, a12, a21 = p
    du[1] = -(a01 + a21) * u[1] + a12 * u[2]
    du[2] = a21 * u[1] - a12 * u[2]
end
solB1 = solve(ODEProblem(pk!, [10.0, 0.0], (0.0, 10.0), (1.0, 2.0, 0.5)), Tsit5(); reltol = 1e-12, abstol = 1e-12)
solB2 = solve(ODEProblem(pk!, [10.0, 0.0], (0.0, 10.0), (2.0, 1.0, 0.5)), Tsit5(); reltol = 1e-12, abstol = 1e-12)
tB = range(0, 10; length = 400)
y1 = [solB1(t)[2] for t in tB]; y2 = [solB2(t)[2] for t in tB]
dB = maximum(abs.(y1 .- y2))
println("Panel B: max |x2_A - x2_B| = ", dB)
pB = plot(tB, y1; xlabel = L"t", ylabel = L"y(t)=x_2(t)", label = L"(a_{01},a_{12},a_{21})=(1,2,0.5)",
          title = L"(B)\ \mathrm{PK},\ a_{01}\leftrightarrow a_{12}")
plot!(pB, tB, y2; linestyle = :dash, label = L"(a_{01},a_{12},a_{21})=(2,1,0.5)")

# (C) after reparameterization kappa = pq: distinct kappa, distinct outputs
tC = range(0, 5; length = 300)
pC = plot(; xlabel = L"t", ylabel = L"x(t)", title = L"(C)\ \dot x=-\kappa x")
for (κ, ls) in zip([0.5, 1, 2, 4], [:solid, :dash, :dot, :dashdot])
    plot!(pC, tC, exp.(-κ .* tC); linestyle = ls, label = latexstring("\\kappa=$κ"))
end

# (D) SIWR with incidence observed: (gamma, mu, beta_W*xi) -> (mu, gamma, beta_W*xi + beta_I*(mu-gamma))
function siwr!(du, u, p, t)
    S, I, z = u
    bI, g, m, κ = p
    F = bI * I + z
    du[1] = -S * F; du[2] = S * F - g * I; du[3] = κ * I - m * z
end
p1 = [0.3, 0.5, 0.25, 0.2]
p2 = [p1[1], p1[3], p1[2], p1[4] + p1[1] * (p1[3] - p1[2])]
u1 = [0.99, 0.01, 0.02]
u2 = [u1[1], u1[2] + (p1[2] - p1[3]) * u1[3] / p1[4], u1[3] * p2[4] / p1[4]]
sD1 = solve(ODEProblem(siwr!, u1, (0.0, 60.0), p1), Tsit5(); reltol = 1e-12, abstol = 1e-12)
sD2 = solve(ODEProblem(siwr!, u2, (0.0, 60.0), p2), Tsit5(); reltol = 1e-12, abstol = 1e-12)
tD = range(0, 60; length = 600)
inc(s, t) = (u = s(t); u[1] * (s.prob.p[1] * u[2] + u[3]))
yD1 = [inc(sD1, t) for t in tD]; yD2 = [inc(sD2, t) for t in tD]
dD = maximum(abs.(yD1 .- yD2))
println("Panel D: max |y_A - y_B| = ", dD, "   (max y = ", maximum(yD1), ")")
println("Panel D: parameter sets ", p1, " and ", p2, "; initial states ", u1, " and ", u2)
pD = plot(tD, yD1; xlabel = L"t", ylabel = L"y(t)=S(\beta_I I+\beta_W W)",
          label = L"(\gamma,\mu,\beta_W\xi)=(0.5,0.25,0.2)",
          title = L"(D)\ \mathrm{SIWR},\ \gamma\leftrightarrow\mu")
plot!(pD, tD, yD2; linestyle = :dash, label = L"(\gamma,\mu,\beta_W\xi)=(0.25,0.5,0.125)")

pall = plot(pA, pB, pC, pD; layout = (2, 2), size = (900, 640), left_margin = 4Plots.mm,
            bottom_margin = 3Plots.mm)
savefig(pall, joinpath(OUT, "fig_nonidentifiability.pdf"))
println("wrote fig_nonidentifiability.pdf")
