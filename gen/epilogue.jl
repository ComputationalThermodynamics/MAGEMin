#
# START OF EPILOGUE
#

# stable phases
struct SS_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    Comp_wt::Vector{Cdouble}
    compVariables::Vector{Cdouble}
    compVariablesNames::Vector{String}
    emNames::Vector{String}
    emFrac::Vector{Cdouble}
    emFrac_wt::Vector{Cdouble}
    emChemPot::Vector{Cdouble}
    emComp::Vector{Vector{Float64}}
    emComp_wt::Vector{Vector{Float64}}
end

function Base.convert(::Type{SS_data}, a::stb_SS_phases) 
    return SS_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.entropy, a.enthalpy, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                                    unsafe_wrap( Vector{Cdouble},        a.Comp,             a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.Comp_wt,          a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.compVariables,    a.n_xeos),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.compVariablesNames,a.n_xeos)),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.emNames,          a.n_em)),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac,           a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac_wt,        a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emChemPot,        a.n_em),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp, a.n_em),  a.nOx),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp_wt, a.n_em),  a.nOx)   )
end

# metastable phases
struct mSS_data
    ph_name::String
    info::String
    ph_id::Cint
    n_xeos::Cint
    n_em::Cint
    G_Ppc::Cdouble
    DF_Ppc::Cdouble
    comp_Ppc::Vector{Cdouble}
    p_Ppc::Vector{Cdouble}
    mu_Ppc::Vector{Cdouble}
    xeos_Ppc::Vector{Cdouble}
end

function Base.convert(::Type{mSS_data}, a::mstb_SS_phases) 
    return  mSS_data(   unsafe_string(a.ph_name),
                        unsafe_string(a.info),
                        a.ph_id, a.n_xeos, a.n_em, a.G_Ppc, a.DF_Ppc,
                        unsafe_wrap( Vector{Cdouble},        a.comp_Ppc,           a.nOx),
                        unsafe_wrap( Vector{Cdouble},        a.p_Ppc,              a.n_em),
                        unsafe_wrap( Vector{Cdouble},        a.mu_Ppc,             a.n_em),
                        unsafe_wrap( Vector{Cdouble},        a.xeos_Ppc,           a.n_xeos)    )
end

# pure phases
struct PP_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    entropy::Cdouble
    enthalpy::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    Comp_wt::Vector{Cdouble}
end

function Base.convert(::Type{PP_data}, a::stb_PP_phases) 
    return PP_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.entropy, a.enthalpy, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                    unsafe_wrap(Vector{Cdouble},a.Comp, a.nOx),
                    unsafe_wrap(Vector{Cdouble},a.Comp_wt, a.nOx))
end





#
# END OF EPILOGUE
#