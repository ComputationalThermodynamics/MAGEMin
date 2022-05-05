#
# START OF EPILOGUE
#




struct SS_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
    compVariables::Vector{Cdouble}
    emNames::Vector{String}
    emFrac::Vector{Cdouble}
    emChemPot::Vector{Cdouble}
    emComp::Vector{Vector{Float64}}
end



function Base.convert(::Type{SS_data}, a::stb_SS_phases) 
    return SS_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                                    unsafe_wrap( Vector{Cdouble},        a.Comp,             a.nOx),
                                    unsafe_wrap( Vector{Cdouble},        a.compVariables,    a.n_xeos),
                    unsafe_string.( unsafe_wrap( Vector{Ptr{Int8}},      a.emNames,          a.n_em)),
                                    unsafe_wrap( Vector{Cdouble},        a.emFrac,           a.n_em),
                                    unsafe_wrap( Vector{Cdouble},        a.emChemPot,        a.n_em),
      unsafe_wrap.(Vector{Cdouble}, unsafe_wrap( Vector{Ptr{Cdouble}},   a.emComp, a.n_em),  a.nOx)   )
end

struct PP_data
    f::Cdouble
    G::Cdouble
    deltaG::Cdouble
    V::Cdouble
    alpha::Cdouble
    cp::Cdouble
    rho::Cdouble
    bulkMod::Cdouble
    shearMod::Cdouble
    Vp::Cdouble
    Vs::Cdouble
    Comp::Vector{Cdouble}
end

function Base.convert(::Type{PP_data}, a::stb_PP_phases) 
    return PP_data(a.f, a.G, a.deltaG, a.V, a.alpha, a.cp, a.rho, a.bulkMod, a.shearMod, a.Vp, a.Vs,
                    unsafe_wrap(Vector{Cdouble},a.Comp, a.nOx))
end





#
# END OF EPILOGUE
#