# Section 5.5: SIWR model of environmentally transmitted disease.
include("common.jl")
Random.seed!(1)

siwr_incidence = @ODEmodel(
    S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
    I'(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
    W'(t) =  xi * I(t) - mu * W(t),
    R'(t) =  gamma * I(t),
    y(t)  =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t)
)

siwr_incidence_W = @ODEmodel(
    S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
    I'(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
    W'(t) =  xi * I(t) - mu * W(t),
    R'(t) =  gamma * I(t),
    y1(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t),
    y2(t) =  W(t)
)

siwr_prevalence = @ODEmodel(
    S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
    I'(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
    W'(t) =  xi * I(t) - mu * W(t),
    R'(t) =  gamma * I(t),
    y(t)  =  I(t)
)

siwr_prevalence_W = @ODEmodel(
    S'(t) = -beta_I * S(t) * I(t) - beta_W * S(t) * W(t),
    I'(t) =  beta_I * S(t) * I(t) + beta_W * S(t) * W(t) - gamma * I(t),
    W'(t) =  xi * I(t) - mu * W(t),
    R'(t) =  gamma * I(t),
    y1(t) =  I(t),
    y2(t) =  W(t)
)

for (label, ode) in [("y = incidence", siwr_incidence),
                     ("y1 = incidence, y2 = W", siwr_incidence_W),
                     ("y = I (prevalence)", siwr_prevalence),
                     ("y1 = I, y2 = W", siwr_prevalence_W)]
    header("SIWR, $label: global identifiability")
    report(assess_identifiability(ode; loglevel = SI_LOGLEVEL))
    header("SIWR, $label: identifiable functions")
    report_functions(find_identifiable_functions(ode; loglevel = SI_LOGLEVEL))
end
