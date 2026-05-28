using MAGEMin_C
using Plots

# ---------------------------------------------------------------------------
# Reproduce Fig. 6a of Sun & Yao (2026): H2O–CO2 saturation isobars and
# isopleths for a low-K rhyolite at 900°C.
# ---------------------------------------------------------------------------

# 1. MAGEMin equilibrium to obtain the anhydrous melt composition
data    = Initialize_MAGEMin("ig", verbose=-1, solver=0)
P_ref, T_ref = 2.5, 900.0
Xoxides = ["SiO2","TiO2","Al2O3","Cr2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O"]
X       = [77.75, 0.09, 12.56, 0.0, 0.71, 0.01, 0.05, 0.45, 3.21, 5.17, 5.0]
out     = single_point_minimization(P_ref, T_ref, data, X=X, Xoxides=Xoxides, sys_in="wt")
Finalize_MAGEMin(data)

# 2. Generate isobars (fixed P_total, sweep X_H2O in fluid)
P_isobars_kbar = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
n_pts   = 300
xH2O_v  = range(0.0, 1.0, length=n_pts)   # X_H2O_fluid: 0 = pure CO2, 1 = pure H2O

isobar_curves = [begin
    H2O_wt = Float64[]; CO2_wt = Float64[]
    P_bar  = P_kbar * 1000.0
    for xH2O in xH2O_v
        P_H2O = xH2O * P_bar
        P_CO2 = (1.0 - xH2O) * P_bar
        s_h2o, s_co2 = MAGEMin_C.volatile_saturation_SY26(out;
                            P_H2O=P_H2O, P_CO2=P_CO2, P_total=P_bar)
        push!(H2O_wt, isnan(s_h2o) ? 0.0 : s_h2o)
        push!(CO2_wt, isnan(s_co2) ? 0.0 : s_co2 / 1e4)   # ppm → wt%
    end
    (H2O_wt, CO2_wt, P_kbar)
end for P_kbar in P_isobars_kbar]

# 3. Generate isopleths (fixed X_H2O_fluid, sweep P)
xH2O_isopleths = [0.2, 0.5, 0.8]
P_sweep_bar    = range(5.0, 3200.0, length=300)

isopleth_curves = [begin
    H2O_wt = Float64[]; CO2_wt = Float64[]
    for P_bar in P_sweep_bar
        P_H2O = xH2O * P_bar
        P_CO2 = (1.0 - xH2O) * P_bar
        s_h2o, s_co2 = MAGEMin_C.volatile_saturation_SY26(out;
                            P_H2O=P_H2O, P_CO2=P_CO2, P_total=P_bar)
        push!(H2O_wt, isnan(s_h2o) ? 0.0 : s_h2o)
        push!(CO2_wt, isnan(s_co2) ? 0.0 : s_co2 / 1e4)
    end
    (H2O_wt, CO2_wt, xH2O)
end for xH2O in xH2O_isopleths]

# 4. Plot
plt = plot(
    xlabel   = "Melt H₂O (wt%)",
    ylabel   = "Melt CO₂ (wt%)",
    title    = "(a) Rhyolite (900 °C)",
    xlims    = (0, 8.5),
    ylims    = (0, 0.1675),
    legend   = :topright,
    framestyle = :box,
    size     = (480, 420),
)

# isobars
for (i, (H2O, CO2, P_kbar)) in enumerate(isobar_curves)
    lbl = i == 1 ? "Isobar (Pₛₐₜ)" : ""
    plot!(plt, H2O, CO2, color=:steelblue, linewidth=1.5, label=lbl)
    # label at the bottom (H2O-rich, x-axis end) of each isobar
    annotate!(plt, H2O[end-2] + 0.2, 0.003,
              text("$(P_kbar)", :center, 7, :steelblue))
end

# isopleths
for (i, (H2O, CO2, xH2O)) in enumerate(isopleth_curves)
    lbl = i == 1 ? "Isopleth (X_H₂O)" : ""
    plot!(plt, H2O, CO2, color=:gray50, linestyle=:dash, linewidth=1.2, label=lbl)
    # label near the high-P end
    annotate!(plt, H2O[end] + 0.05, CO2[end],
              text("$(xH2O) mol", :left, 7, :gray50))
end

display(plt)
savefig(plt, joinpath(@__DIR__, "SY26_rhyolite_saturation.png"))
println("Saved → examples/SY26_rhyolite_saturation.png")
