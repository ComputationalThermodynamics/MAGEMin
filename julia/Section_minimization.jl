"""
    SECTION_XY

    Diagram axes `(x, y)` of each section type used by [`section_point`](@ref),
    [`multi_point_section`](@ref) and [`section_ranges`](@ref):
    `:ps` → (S, P), `:pv` → (V, P), `:ts` → (T, S), `:tv` → (T, V), `:sv` → (S, V), `:ph` → (H, P), `:sh` → (S, H), `:vh` → (V, H),
    `:rh` → (ρ, H).

    Units: P [kbar], T [°C], S [kJ/K per mol of bulk oxides] (`out.entropy[1]`),
    V [cm³ per mol of bulk oxides] (`out.V`), H [kJ per mol of bulk oxides] (`out.enthalpy[1]`) and
    ρ [kg/m³] (`out.rho`, system density). For a fixed bulk ρ = 1000 M / V with M the molar mass of the bulk
    [g per mol of oxides], so a ρ–H target is solved as the V–H target V = 1000 M / ρ.
"""
const SECTION_XY = Dict(:ps => (:S, :P), :pv => (:V, :P), :ts => (:T, :S), :tv => (:T, :V), :sv => (:S, :V), :ph => (:H, :P), :sh => (:S, :H), :vh => (:V, :H), :rh => (:rho, :H))

"""
    section_info

    Per-point diagnostics returned by [`section_point`](@ref) and [`multi_point_section`](@ref).

    - `status`: 0 converged; 1 maximum iterations; 2 target outside the image of the P–T box
      (including a step that leaves the box with `max_clamp = 0`, and a solution on the box bound
      with `bound_tol`); 3 first-order jump not resolved; 4 repeated clamping on the box (`:sv`,
      `:sh`, `:vh`; more than `max_clamp` clamped steps); 6 point minimization failed.
    - `n_solve`: point minimizations; `n_newton`: section Newton calls.
    - `latent`: the result is a coexistence of assemblages (inside a first-order transition)
      that no single point minimization returns; `out.mSS_vec` then holds the minimization's
      solution phases for use as initial guess.
    - `dfree`: secant derivative of the target property with respect to the free variable near
      the solution (1-D sections; NaN when not available). Its sign identifies the branch of
      non-monotonic sections (`:pv`, `:ts`).
    - `residual`: largest relative mismatch of the returned state to the target(s).
"""
struct section_info
    status      :: Int64
    n_solve     :: Int64
    n_newton    :: Int64
    latent      :: Bool
    dfree       :: Float64
    residual    :: Float64
end

section_value(out, c::Symbol) = c === :P ? out.P_kbar : c === :T ? out.T_C : c === :S ? out.entropy[1] : c === :H ? out.enthalpy[1] : c === :rho ? out.rho : out.V

function section_targets(section::Symbol, x::Float64, y::Float64)
    haskey(SECTION_XY, section) || error("Unknown section $(section); use one of $(collect(keys(SECTION_XY)))")
    ax      = SECTION_XY[section]
    get_c(c) = ax[1] === c ? x : ax[2] === c ? y : NaN
    return (P = get_c(:P), T = get_c(:T), S = get_c(:S), V = get_c(:V), H = get_c(:H), rho = get_c(:rho))
end

section_assemblage(out) = sort(out.ph)

function section_state_union(a, wa::Float64, b, wb::Float64; dx::Float64 = 5e-2, nfree::Int64 = 0)
    sa, sb  = equilibrium_state(a), equilibrium_state(b)
    ss_ids  = copy(sa.ss_ids)
    ss_n    = wa .* sa.ss_n
    xs      = Vector{Vector{Float64}}()
    off     = 0
    for q in eachindex(a.SS_vec)
        k   = length(a.SS_vec[q].compVariables)
        push!(xs, sa.ss_x[off+1:off+k]); off += k
    end
    off     = 0
    for q in eachindex(b.SS_vec)
        k   = length(b.SS_vec[q].compVariables)
        xb  = sb.ss_x[off+1:off+k]; off += k
        f   = findfirst(i -> ss_ids[i] == sb.ss_ids[q] && length(xs[i]) == k && maximum(abs.(xs[i] .- xb)) < dx, eachindex(ss_ids))
        if isnothing(f)
            push!(ss_ids, sb.ss_ids[q]); push!(ss_n, wb*sb.ss_n[q]); push!(xs, xb)
        else
            n0      = ss_n[f]
            n1      = wb*sb.ss_n[q]
            xs[f]   = (n0 .* xs[f] .+ n1 .* xb) ./ max(n0 + n1, eps())
            ss_n[f] = n0 + n1
        end
    end
    pp_ids  = copy(sa.pp_ids)
    pp_n    = wa .* sa.pp_n
    for q in eachindex(sb.pp_ids)
        f = findfirst(==(sb.pp_ids[q]), pp_ids)
        if isnothing(f)
            push!(pp_ids, sb.pp_ids[q]); push!(pp_n, wb*sb.pp_n[q])
        else
            pp_n[f] += wb*sb.pp_n[q]
        end
    end
    nmax    = length(a.oxides) + nfree
    while length(ss_ids) + length(pp_ids) > nmax
        ks  = isempty(ss_n) ? 0 : argmin(ss_n)
        kp  = isempty(pp_n) ? 0 : argmin(pp_n)
        if kp == 0 || (ks > 0 && ss_n[ks] <= pp_n[kp])
            deleteat!(ss_ids, ks); deleteat!(ss_n, ks); deleteat!(xs, ks)
        else
            deleteat!(pp_ids, kp); deleteat!(pp_n, kp)
        end
    end
    return (ss_ids = ss_ids, ss_x = reduce(vcat, xs; init = Float64[]), ss_n = ss_n, pp_ids = pp_ids, pp_n = pp_n,
            gam = (wa .* sa.gam .+ wb .* sb.gam) ./ (wa + wb))
end

"""
    section_point(section, x, y, gv, z_b, DB, splx_data; Pbox, Tbox, seeds=nothing, P0=nothing, T0=nothing, scales=nothing, buffer_n=0.0, rm_list=nothing, W=nothing, scp=0, name_solvus=true, maxit=16, max_clamp=nothing, bound_tol=nothing, log_P=true, tol_G=1e-7, ...)

    Equilibrium at fixed diagram coordinates `(x, y)` of a thermodynamic section (see
    [`SECTION_XY`](@ref)), solved on one MAGEMin handle (`gv, z_b, DB, splx_data`, e.g. one thread
    of a `MAGEMin_Data` whose bulk rock is already defined). The free pressure and/or temperature
    are searched inside `Pbox = (Pmin, Pmax)` [kbar] and `Tbox = (Tmin, Tmax)` [°C].

    Each iteration runs the section Newton step ([`pwm_run_section_from_state`](@ref)) from the
    current stable state to the target, then one point minimization at the solved P, T. The Newton
    state is accepted when its Gibbs energy is not above the minimization's (relative tolerance
    `tol_G`). Otherwise the Newton step is repeated on the union of both assemblages, which
    resolves targets inside a first-order transition. When the Newton step cannot reach the
    target, the free variables are moved by point minimizations with secant/bisection (1-D
    sections) or Broyden (`:sv`, `:sh`, `:vh`) steps, and the Newton step is tried again from the new state.

    Start: the closest of `seeds` (vector of `gmin_struct`, e.g. converged neighbours) in the
    diagram coordinates scaled by `scales = (dx, dy)`, whose `mSS_vec` are also used as initial
    guess; otherwise a point minimization at `(P0, T0)` (default: box centre), with the fixed
    coordinate set.

    With `log_P` (default), 1-D sections with free pressure (`:ts`, `:tv`) step in log P, where
    the volume of fluid-rich assemblages (V ∝ 1/P at low pressure) is close to linear. `log_P_2d`
    extends this to the 2-D sections (`:sv`, `:sh`, `:vh`, `:rh`): Broyden steps and finite-difference
    Jacobians are then taken in (log P, T), which recovers low-pressure fluid-rich targets but spends
    more iterations on unreachable ones (use it to re-solve isolated failures).

    Box bounds: by default a step leaving the box is clamped to the bound (1-D: the bound is
    evaluated once, a second step past it gives status 2; `:sv`/`:sh`/`:vh`: more than 4 clamped steps give
    status 4). `max_clamp = k` stops after `k` clamped steps instead (`0`: as soon as a step leaves the
    box). `bound_tol` rejects a converged state whose free P or T lies within `bound_tol` × box width of
    a bound (status 2).

    TC, SB and BR databases are supported (section Newton step).

    `:rh` (ρ–H) targets are converted to V–H targets with the bulk molar mass `M = ρ V / 1000` of the closest
    seed (or of a point minimization at `(P0, T0)`), and solved as `:vh`.

    Returns `(out, info)`: `out::gmin_struct` at the solved `P_kbar`, `T_C`, and
    `info::section_info`.
"""
function section_point( section     :: Symbol,
                        x           :: Float64,
                        y           :: Float64,
                        gv, z_b, DB, splx_data;
                        Pbox        :: NTuple{2,Float64},
                        Tbox        :: NTuple{2,Float64},
                        seeds       = nothing,
                        P0          = nothing,
                        T0          = nothing,
                        scales      = nothing,
                        buffer_n    :: Float64 = 0.0,
                        rm_list     = nothing,
                        W           = nothing,
                        scp         :: Int64 = 0,
                        name_solvus :: Bool = true,
                        maxit       :: Int64 = 16,
                        max_clamp   :: Union{Nothing, Int64} = nothing,
                        bound_tol   :: Union{Nothing, Float64} = nothing,
                        log_P       :: Bool = true,
                        log_P_2d    :: Bool = false,
                        tol_G       :: Float64 = 1e-7,
                        warm        :: Bool = false,
                        verbose     :: Bool = false,
                        tol_point   :: Float64 = 1e-8,
                        seismic_cor = false, aspect_ratio = 0.3, seismic_water = 0, shallow_correction = false,
                        fluid_as_melt = false, anelastic_cor = false)

    if section === :rh
        probe = if !isnothing(seeds) && !isempty(seeds)
            sc0 = isnothing(scales) ? (1.0, 1.0) : scales
            seeds[argmin([((s.rho - x)/sc0[1])^2 + ((s.enthalpy[1] - y)/sc0[2])^2 for s in seeds])]
        else
            point_wise_minimization(isnothing(P0) ? 0.5*(Pbox[1] + Pbox[2]) : Float64(P0), isnothing(T0) ? 0.5*(Tbox[1] + Tbox[2]) : Float64(T0),
                                    gv, z_b, DB, splx_data; buffer_n = buffer_n, rm_list = rm_list, W = W, scp = scp, name_solvus = name_solvus)
        end
        M   = probe.rho * probe.V / 1000.0
        Vt  = 1000.0 * M / x
        sc  = isnothing(scales) ? nothing : (scales[1] * Vt / x, scales[2])
        return section_point(:vh, Vt, y, gv, z_b, DB, splx_data;
                             Pbox, Tbox, seeds, P0, T0, scales = sc, buffer_n, rm_list, W, scp, name_solvus, maxit, max_clamp,
                             bound_tol, log_P, log_P_2d, tol_G, warm, verbose, tol_point, seismic_cor, aspect_ratio, seismic_water,
                             shallow_correction, fluid_as_melt, anelastic_cor)
    end

    tg      = section_targets(section, x, y)
    fixP    = !isnan(tg.P)
    fixT    = !isnan(tg.T)
    rows    = [r for r in (:S, :V, :H) if !isnan(getfield(tg, r))]
    free    = [c for c in (:P, :T) if isnan(getfield(tg, c))]
    lgP     = log_P && Pbox[1] > 0.0 && (free == [:P] || (log_P_2d && :P in free))
    tP(v)   = lgP ? log(max(v, 1e-300)) : v
    uP(u)   = lgP ? exp(u) : u
    lo      = Dict(:P => tP(Pbox[1]), :T => Tbox[1])
    hi      = Dict(:P => tP(Pbox[2]), :T => Tbox[2])
    wid     = Dict(:P => tP(Pbox[2]) - tP(Pbox[1]), :T => Tbox[2] - Tbox[1])
    eps_th  = Dict(:P => 1e-5, :T => 1e-4)
    target  = Dict(:S => tg.S, :V => tg.V, :H => tg.H)
    rscale  = Dict(:S => max(abs(tg.S), 1e-3), :V => max(abs(tg.V), 1e-3), :H => max(abs(tg.H), 1e-3))
    n_solve = 0
    n_newt  = 0
    Gs      = LibMAGEMin.mSS_data[]
    gof(o)  = isempty(o.mSS_vec) ? Gs : o.mSS_vec

    solve(P, T, G) = begin
        n_solve += 1
        ig = warm && !(isnothing(G) || isempty(G))
        point_wise_minimization(max(P, 0.001), T, gv, z_b, DB, splx_data;
                                Gi = ig ? G : nothing, iguess = ig,
                                buffer_n = buffer_n, rm_list = rm_list, W = W, scp = scp, name_solvus = name_solvus,
                                seismic_cor = seismic_cor, aspect_ratio = aspect_ratio, seismic_water = seismic_water,
                                shallow_correction = shallow_correction, fluid_as_melt = fluid_as_melt, anelastic_cor = anelastic_cor)
    end
    newton(P, T, st) = begin
        n_newt += 1
        point_wise_section_newton_from_state(max(P, 0.001), T, gv, z_b, DB, splx_data, st;
                                             section = section, S = tg.S, V = tg.V, H = tg.H, W = W, name_solvus = name_solvus)
    end
    resid(o)    = [ (section_value(o, r) - target[r]) / rscale[r] for r in rows ]
    rnorm(o)    = maximum(abs.(resid(o)))
    fval(o, c)  = c === :P ? tP(section_value(o, c)) : section_value(o, c)
    theta(o)    = [ fval(o, c) for c in free ]
    inbox(o)    = all(c -> lo[c] - 1e-9 <= fval(o, c) <= hi[c] + 1e-9, free)
    pt_of(o)    = (fixP ? tg.P : o.P_kbar, fixT ? tg.T : o.T_C)
    on_fixed(o) = (!fixP || abs(o.P_kbar - tg.P) <= 1e-9*max(1.0, abs(tg.P))) && (!fixT || abs(o.T_C - tg.T) <= 1e-9*max(1.0, abs(tg.T)))
    stable_ok(c)  = c.status < 3
    no_higher(o, c) = o.G_system <= c.G_system + tol_G*max(1.0, abs(c.G_system))
    lower(o, c)     = o.G_system <  c.G_system - tol_G*max(1.0, abs(c.G_system))
    function contains(o, c)
        po = section_assemblage(o)
        for p in section_assemblage(c)
            k = findfirst(==(p), po)
            isnothing(k) && return false
            deleteat!(po, k)
        end
        return true
    end

    hist    = NamedTuple{(:th, :r, :out), Tuple{Vector{Float64}, Vector{Float64}, Any}}[]
    record!(o) = push!(hist, (th = theta(o), r = resid(o), out = o))

    function dfree_est()
        length(free) == 1 || return NaN
        k = length(hist)
        k >= 2 || return NaN
        dth = section_value(hist[k].out, free[1]) - section_value(hist[k-1].out, free[1])
        return dth != 0 ? (hist[k].r[1] - hist[k-1].r[1])*rscale[rows[1]]/dth : NaN
    end
    on_bound(o) = !isnothing(bound_tol) && any(c -> min(fval(o, c) - lo[c], hi[c] - fval(o, c)) <= bound_tol*wid[c], free)
    finish(o, st, lat) = ((!isnothing(seeds) && any(sd -> sd === o, seeds)) ? deepcopy(o) : o,
                          section_info((st == 0 && on_bound(o)) ? 2 : st, n_solve, n_newt, lat, dfree_est(), rnorm(o)))

    function newton_accept(P, T, st, G)
        o = newton(P, T, st)
        verbose && !(o.status == 0 && inbox(o) && rnorm(o) < 1e-6) && println("  section Newton from P $(round(P, digits = 6)) T $(round(T, digits = 4)) failed: status $(o.status) P $(round(o.P_kbar, digits = 6)) T $(round(o.T_C, digits = 4)) res $(rnorm(o))")
        (o.status == 0 && inbox(o) && rnorm(o) < 1e-6) || return nothing, nothing
        c = solve(o.P_kbar, o.T_C, G)
        verbose && println("  section Newton: P $(round(o.P_kbar, digits = 6)) T $(round(o.T_C, digits = 4)) G $(round(o.G_system, digits = 8)) $(join(section_assemblage(o), " ")) | check status $(c.status) G $(round(c.G_system, digits = 8)) $(join(section_assemblage(c), " ")) | res $(rnorm(c))")
        stable_ok(c) || return nothing, nothing
        if (no_higher(o, c) && contains(o, c)) || lower(o, c)
            same = section_assemblage(c) == section_assemblage(o)
            same && rnorm(c) < 1e-6 && return (c, false), c
            append!(empty!(o.mSS_vec), c.mSS_vec)
            return (o, !same), c
        end
        return nothing, (o, c)
    end

    if !isnothing(seeds) && !isempty(seeds)
        ax  = SECTION_XY[section]
        sc  = isnothing(scales) ? (1.0, 1.0) : scales
        d   = [ ((section_value(s, ax[1]) - x)/sc[1])^2 + ((section_value(s, ax[2]) - y)/sc[2])^2 for s in seeds ]
        cur = seeds[argmin(d)]
        Gs  = reduce(vcat, [s.mSS_vec for s in seeds]; init = LibMAGEMin.mSS_data[])
    else
        Ps  = fixP ? tg.P : (isnothing(P0) ? 0.5*(Pbox[1] + Pbox[2]) : Float64(P0))
        Ts  = fixT ? tg.T : (isnothing(T0) ? 0.5*(Tbox[1] + Tbox[2]) : Float64(T0))
        cur = solve(Ps, Ts, nothing)
        stable_ok(cur) || return finish(cur, 6, false)
    end
    own(o) = isnothing(seeds) || !any(sd -> sd === o, seeds)
    if on_fixed(cur) && own(cur)
        record!(cur)
    end

    n_clamp = 0
    Jb      = nothing

    for it in 1:maxit
        Ps, Ts      = pt_of(cur)
        if own(cur) && on_fixed(cur) && stable_ok(cur) && !isempty(cur.mSS_vec) && rnorm(cur) < tol_point
            return finish(cur, 0, false)
        end
        verbose && println("section iteration $(it): current P $(round(cur.P_kbar, digits = 6)) T $(round(cur.T_C, digits = 4)) $(join(section_assemblage(cur), " ")) res $(rnorm(cur))")
        res, other  = newton_accept(Ps, Ts, equilibrium_state(cur), gof(cur))
        isnothing(res) || return finish(res[1], 0, res[2])

        if !isnothing(other)
            o, c        = other
            res2, c2    = newton_accept(o.P_kbar, o.T_C, section_state_union(o, 0.5, c, 0.5; nfree = length(free)), gof(c))
            isnothing(res2) || return finish(res2[1], 0, res2[2])
            record!(c)
            cur = c
            if !isnothing(c2) && c2 isa Tuple
                record!(c2[2])
                cur = c2[2]
            end
            if own(cur) && on_fixed(cur) && rnorm(cur) < tol_point
                return finish(cur, 0, false)
            end
        end

        if !on_fixed(cur)
            cur = solve(Ps, Ts, gof(cur))
            stable_ok(cur) || return finish(cur, 6, false)
            record!(cur)
            continue
        end

        r0  = resid(cur)
        th0 = theta(cur)

        if length(free) == 1
            f       = free[1]
            neg     = [k for k in eachindex(hist) if hist[k].r[1] < 0]
            pos     = [k for k in eachindex(hist) if hist[k].r[1] > 0]
            if !isempty(neg) && !isempty(pos)
                ka, kb  = argmin(p -> abs(hist[p[1]].th[1] - hist[p[2]].th[1]), [(a, b) for a in neg, b in pos])
                a, b    = hist[ka].th[1], hist[kb].th[1]
                if abs(a - b) < 1e-2*wid[f] && section_assemblage(hist[ka].out) != section_assemblage(hist[kb].out)
                    ra, rb  = hist[ka].r[1], hist[kb].r[1]
                    lam     = ra/(ra - rb)
                    thm     = a + lam*(b - a)
                    Pm, Tm  = f === :P ? (uP(thm), Ts) : (Ps, thm)
                    st      = section_state_union(hist[ka].out, 1.0 - lam, hist[kb].out, lam; nfree = length(free))
                    res3, _ = newton_accept(Pm, Tm, st, vcat(gof(hist[ka].out), gof(hist[kb].out)))
                    isnothing(res3) || return finish(res3[1], 0, res3[2])
                end
                abs(a - b) < eps_th[f] && return finish(cur, 3, false)
                ra, rb  = hist[ka].r[1], hist[kb].r[1]
                lb, ub  = min(a, b), max(a, b)
                inside(t) = isfinite(t) && lb + 0.02*(ub - lb) < t < ub - 0.02*(ub - lb)
                k1, k2  = length(hist), length(hist) - 1
                ts      = (k2 >= 1 && hist[k1].th[1] != hist[k2].th[1]) ?
                          hist[k1].th[1] - hist[k1].r[1]*(hist[k1].th[1] - hist[k2].th[1])/(hist[k1].r[1] - hist[k2].r[1]) : NaN
                tf      = a - ra*(b - a)/(rb - ra)
                thn     = inside(ts) ? ts : inside(tf) ? tf : 0.5*(a + b)
            else
                if length(hist) >= 2 && hist[end].th[1] != hist[end-1].th[1]
                    der = (hist[end].r[1] - hist[end-1].r[1])/(hist[end].th[1] - hist[end-1].th[1])
                else
                    h   = 1e-2*wid[f]*(th0[1] + 1e-2*wid[f] > hi[f] ? -1 : 1)
                    Pn, Tn = f === :P ? (uP(th0[1] + h), Ts) : (Ps, th0[1] + h)
                    cf  = solve(Pn, Tn, gof(cur))
                    stable_ok(cf) || return finish(cur, 6, false)
                    record!(cf)
                    der = (resid(cf)[1] - r0[1])/h
                    cur = cf
                    sign(resid(cf)[1]) != sign(r0[1]) && continue
                    r0, th0 = resid(cf), theta(cf)
                end
                (isfinite(der) && der != 0) || (der = sign(r0[1])/wid[f])
                thn = th0[1] + clamp(-r0[1]/der, -0.25*wid[f], 0.25*wid[f])
            end
            if !(lo[f] <= thn <= hi[f])
                bnd      = thn < lo[f] ? lo[f] : hi[f]
                k_b      = findlast(k -> abs(hist[k].th[1] - bnd) <= 1e-12*max(1.0, abs(bnd)), eachindex(hist))
                isnothing(k_b) || return finish(hist[k_b].out, 2, false)
                (!isnothing(max_clamp) && n_clamp >= max_clamp) && return finish(cur, 2, false)
                thn      = bnd
                n_clamp += 1
            end
            Pn, Tn = f === :P ? (uP(thn), Ts) : (Ps, thn)
        else
            if isnothing(Jb)
                Jb = zeros(2, 2)
                for (j, f) in enumerate(free)
                    h   = 1e-2*wid[f]*(th0[j] + 1e-2*wid[f] > hi[f] ? -1 : 1)
                    Pn  = f === :P ? uP(th0[j] + h) : cur.P_kbar
                    Tn  = f === :T ? th0[j] + h : cur.T_C
                    cf  = solve(Pn, Tn, gof(cur))
                    stable_ok(cf) || return finish(cur, 6, false)
                    Jb[:, j] = (resid(cf) .- r0) ./ h
                end
            elseif length(hist) >= 2
                s   = hist[end].th .- hist[end-1].th
                yv  = hist[end].r .- hist[end-1].r
                sum(abs2, s) > 0 && (Jb .+= ((yv .- Jb*s) * s') ./ sum(abs2, s))
            end
            step = try
                -(Jb \ r0)
            catch
                return finish(cur, 1, false)
            end
            all(isfinite, step) || return finish(cur, 1, false)
            clamped = false
            for (j, f) in enumerate(free)
                step[j] = clamp(step[j], -0.25*wid[f], 0.25*wid[f])
                v       = th0[j] + step[j]
                if !(lo[f] <= v <= hi[f])
                    step[j] = clamp(v, lo[f], hi[f]) - th0[j]
                    clamped = true
                end
            end
            n_clamp += clamped
            n_clamp > something(max_clamp, 4) && return finish(cur, 4, false)
            Pn = uP(th0[findfirst(==(:P), free)] + step[findfirst(==(:P), free)])
            Tn = th0[findfirst(==(:T), free)] + step[findfirst(==(:T), free)]
        end

        c = solve(Pn, Tn, gof(cur))
        stable_ok(c) || return finish(cur, 6, false)
        if length(free) == 2 && section_assemblage(c) != section_assemblage(cur)
            res4, _ = newton_accept(c.P_kbar, c.T_C, section_state_union(cur, 0.5, c, 0.5; nfree = length(free)), gof(c))
            isnothing(res4) || return finish(res4[1], 0, res4[2])
        end
        record!(c)
        cur = c
    end
    return finish(cur, 1, false)
end


function section_seed_grid(MAGEMin_db, Pbox, Tbox, n; kwargs...)
    Pv  = collect(range(Pbox[1], Pbox[2], length = n))
    Tv  = collect(range(Tbox[1], Tbox[2], length = n))
    P   = vec([p for p in Pv, t in Tv])
    T   = vec([t for p in Pv, t in Tv])
    return multi_point_minimization(P, T, MAGEMin_db; kwargs...)
end

"""
    section_ranges(section, MAGEMin_db; Pbox, Tbox, n=7, X=nothing, B=nothing, Xoxides=String[], sys_in="mol", test=0, rm_list=nothing, W=nothing, progressbar=false)

    Diagram ranges of a thermodynamic section ([`SECTION_XY`](@ref)): the bounding box of the image
    of the P–T box, from point minimizations on an `n × n` grid of the box. P and T axes keep the box
    limits.

    Returns `(xrange, yrange, grid)`, with `grid` the vector of grid results (reusable as `grid` in
    [`multi_point_section`](@ref)).
"""
function section_ranges(section     :: Symbol,
                        MAGEMin_db  :: MAGEMin_Data;
                        Pbox        :: NTuple{2,Float64},
                        Tbox        :: NTuple{2,Float64},
                        n           :: Int64 = 7,
                        X = nothing, B = nothing, Xoxides = String[], sys_in = "mol", test = 0,
                        rm_list = nothing, W = nothing, scp = 0, progressbar = false)
    Bv      = isnothing(B) ? nothing : fill(Float64(B isa Number ? B : B[1]), n*n)
    grid    = section_seed_grid(MAGEMin_db, Pbox, Tbox, n; X = X, B = Bv, Xoxides = Xoxides, sys_in = sys_in, test = test,
                                rm_list = rm_list, W = W, scp = scp, name_solvus = true, progressbar = progressbar)
    ok      = filter(o -> o.status < 3, grid)
    isempty(ok) && error("section_ranges: no point minimization of the $(n)×$(n) grid converged in Pbox = $(Pbox), Tbox = $(Tbox)")
    rng(c)  = c === :P ? Pbox : c === :T ? Tbox : (minimum(section_value.(ok, c)), maximum(section_value.(ok, c)))
    ax      = SECTION_XY[section]
    return rng(ax[1]), rng(ax[2]), grid
end

function section_retry_seeds(section::Symbol, x::Float64, y::Float64, pool, tried, Pbox, Tbox; n_max::Int64 = 6)
    tg      = section_targets(section, x, y)
    fixc    = isnan(tg.P) ? :T : :P
    prop    = !isnan(tg.S) ? :S : !isnan(tg.H) ? :H : :V
    span    = fixc === :P ? Pbox[2] - Pbox[1] : Tbox[2] - Tbox[1]
    target  = getfield(tg, prop)
    fixv    = getfield(tg, fixc)
    cand    = [o for o in pool if !any(t -> t === o, tried)]
    isempty(cand) && return cand
    dfix    = [abs(section_value(o, fixc) - fixv)/span for o in cand]
    near    = cand[dfix .<= minimum(dfix) + 1e-9 + 0.2]
    above   = sort([o for o in near if section_value(o, prop) >= target], by = o -> section_value(o, prop) - target)
    below   = sort([o for o in near if section_value(o, prop) <  target], by = o -> target - section_value(o, prop))
    out     = Any[]
    for k in 1:max(length(above), length(below))
        k <= length(below) && push!(out, below[k])
        k <= length(above) && push!(out, above[k])
        length(out) >= n_max && break
    end
    return out
end

function section_foreach(f::Function, n::Int)
    next = Threads.Atomic{Int}(0)
    @threads :static for w in 1:Threads.nthreads()
        id = Threads.threadid()
        while true
            k = Threads.atomic_add!(next, 1) + 1
            k > n && break
            f(k, id)
        end
    end
    return nothing
end

"""
    multi_point_section(section, x, y, MAGEMin_db; Pbox, Tbox, seeds=nothing, grid=nothing, seed_grid=7, xrange=nothing, yrange=nothing, X=nothing, B=nothing, Xoxides=String[], sys_in="mol", test=0, rm_list=nothing, W=nothing, scp=0, name_solvus=true, max_clamp=:auto, bound_tol=1e-6, log_P=true, log_P_2d=false, retry_passes=2, retry_maxit=16, retry_max_clamp=:auto, progressbar=true, callback_fn=nothing, callback_int=1, ...)

    Threaded [`section_point`](@ref) over the targets `(x[i], y[i])` of a thermodynamic section
    ([`SECTION_XY`](@ref)), with one MAGEMin handle per thread as in
    [`multi_point_minimization`](@ref). The bulk rock follows `multi_point_minimization`: one
    composition `X` (with `Xoxides`, `sys_in`) for all points, one composition per target
    (`X::Vector{Vector{Float64}}`), or the built-in `test`. `B` is the buffer value, one for all points
    or one per target. With one composition per target, the internal P–T seed grid is not built:
    pass `seeds` (e.g. points computed at each target's bulk) or `grid`; targets without seeds start
    from a point minimization at the box centre. Failed targets are retried from their nearest
    converged targets, in `(x, y)` or in the coordinates `retry_xy = (u, v)` (one value per target, e.g.
    the composition axis of an S–X diagram). The retries (`retry_passes`, each attempt from one neighbour,
    all attempts in parallel) use at most `retry_maxit` iterations and stop after `retry_max_clamp` steps
    clamped on the search box. With `max_clamp = :auto` (default, also for `retry_max_clamp`),
    2-D sections (`:sv`, `:sh`, `:vh`, `:rh`) leave a target blank as soon as a step leaves the box
    (`max_clamp = 0`, see [`section_point`](@ref)), while 1-D sections keep the clamp to the bound,
    which is evaluated once (`nothing`): 1-D properties such as V(P) of a fluid-rich bulk are very
    non-linear, and a secant step from a distant seed often overshoots the box although the target
    is inside. Solutions on the box bound are rejected (`bound_tol = 1e-6`).

    Seeds per target: `seeds[i]` (vector of `gmin_struct`, e.g. the converged neighbours of a refined
    mesh point). Without `seeds`, or for an empty `seeds[i]`, the closest node of a
    `seed_grid × seed_grid` P–T grid of the box is used (`grid` reuses the grid of [`section_ranges`](@ref); `seed_Pbox`/`seed_Tbox` set the grid's box when
    it should differ from the search box, default `Pbox`/`Tbox`). `xrange`/`yrange` scale the distances
    (default: the image of the grid).

    Returns `(out, info)`: `Vector{gmin_struct}` and `Vector{section_info}`.

    Examples
    --------
    ```julia
    data    = Initialize_MAGEMin("mp", verbose = false)
    xr, yr, grid = section_ranges(:pv, data; Pbox = (2.0, 10.0), Tbox = (500.0, 800.0), test = 0)
    V       = collect(range(xr..., length = 5))
    P       = fill(6.0, 5)
    out, info = multi_point_section(:pv, V, P, data; Pbox = (2.0, 10.0), Tbox = (500.0, 800.0), test = 0, grid = grid)
    ```
"""
function multi_point_section(section     :: Symbol,
                             x           :: AbstractVector{Float64},
                             y           :: AbstractVector{Float64},
                             MAGEMin_db  :: MAGEMin_Data;
                             Pbox        :: NTuple{2,Float64},
                             Tbox        :: NTuple{2,Float64},
                             seeds       = nothing,
                             grid        = nothing,
                             seed_grid   :: Int64 = 7,
                             seed_Pbox   = nothing,
                             seed_Tbox   = nothing,
                             retry_xy    = nothing,
                             retry_passes:: Int64 = 2,
                             retry_maxit :: Int64 = 16,
                             retry_max_clamp = :auto,
                             max_clamp   = :auto,
                             bound_tol   :: Union{Nothing, Float64} = 1e-6,
                             log_P       :: Bool = true,
                             log_P_2d    :: Bool = false,
                             xrange      = nothing,
                             yrange      = nothing,
                             X           = nothing,
                             B           = nothing,
                             Xoxides     = String[],
                             sys_in      :: String = "mol",
                             test        :: Int64 = 0,
                             rm_list     = nothing,
                             W           = nothing,
                             scp         :: Int64 = 0,
                             name_solvus :: Bool = true,
                             maxit       :: Int64 = 16,
                             warm        :: Bool = false,
                             progressbar :: Bool = true,
                             callback_fn :: Union{Nothing, Function} = nothing,
                             callback_int:: Int64 = 1,
                             seismic_cor = false, aspect_ratio = 0.3, seismic_water = 0, shallow_correction = false,
                             fluid_as_melt = false, anelastic_cor = false)

    @assert length(x) == length(y) "x and y must have the same length"
    haskey(SECTION_XY, section) || error("Unknown section $(section); use one of $(collect(keys(SECTION_XY)))")
    isnothing(seeds) || @assert length(seeds) == length(x) "seeds must have one entry per target"
    per_pt  = X isa AbstractVector && !isempty(X) && first(X) isa AbstractVector
    per_pt && @assert length(X) == length(x) "X must have one bulk composition per target"
    bufs    = isnothing(B) ? fill(0.0, length(x)) : B isa Number ? fill(Float64(B), length(x)) :
              length(B) == length(x) ? Float64.(B) : fill(Float64(B[1]), length(x))
    buffer  = bufs[1]

    need_pool = isnothing(seeds) || any(isempty, seeds)
    if need_pool && isnothing(grid) && !per_pt
        grid = section_seed_grid(MAGEMin_db, something(seed_Pbox, Pbox), something(seed_Tbox, Tbox), seed_grid; X = X, B = fill(buffer, seed_grid^2), Xoxides = Xoxides,
                                 sys_in = sys_in, test = test, rm_list = rm_list, W = W, scp = scp, name_solvus = name_solvus,
                                 progressbar = false)
    end

    if !per_pt
        for i in 1:Threads.maxthreadid()
            MAGEMin_db.gv[i] = isnothing(X) ? use_predefined_bulk_rock(MAGEMin_db.gv[i], test, MAGEMin_db.db) :
                                              define_bulk_rock(MAGEMin_db.gv[i], X, Xoxides, sys_in, MAGEMin_db.db)
        end
    end

    ax      = SECTION_XY[section]
    one_d   = any(c -> c === :P || c === :T, ax)
    max_clamp       === :auto && (max_clamp = one_d ? nothing : 0)
    retry_max_clamp === :auto && (retry_max_clamp = one_d ? nothing : 0)
    pool    = isnothing(grid) ? nothing : filter(o -> o.status < 3, grid)
    span(v) = (isempty(v) || maximum(v) == minimum(v)) ? 1.0 : maximum(v) - minimum(v)
    axspan(c, r) = !isnothing(r) ? r[2] - r[1] : c === :P ? Pbox[2] - Pbox[1] : c === :T ? Tbox[2] - Tbox[1] :
                   !isnothing(pool) ? span(section_value.(pool, c)) :
                   (!isnothing(seeds) && any(!isempty, seeds)) ? span(section_value.(reduce(vcat, seeds), c)) : 1.0
    scales  = (axspan(ax[1], xrange), axspan(ax[2], yrange))

    n       = length(x)
    outs    = Vector{gmin_struct{Float64, Int64}}(undef, n)
    infos   = Vector{section_info}(undef, n)
    progr   = progressbar ? Progress(n, desc = "Computing $(n) section points...") : nothing
    count   = Threads.Atomic{Int}(0)

    section_foreach(n) do i, id
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        per_pt && (gv = define_bulk_rock(gv, X[i], Xoxides, sys_in, MAGEMin_db.db))

        sd          = isnothing(seeds) ? pool : seeds[i]
        if isnothing(seeds) || isempty(sd)
            if isnothing(pool)
                sd  = nothing
            else
                d   = [ ((section_value(s, ax[1]) - x[i])/scales[1])^2 + ((section_value(s, ax[2]) - y[i])/scales[2])^2 for s in pool ]
                sd  = pool[partialsortperm(d, 1:min(2, length(d)))]
            end
        end

        solve_i(sds) = section_point(section, x[i], y[i], gv, z_b, DB, splx_data;
                                    Pbox = Pbox, Tbox = Tbox, seeds = sds, scales = scales, buffer_n = bufs[i],
                                    rm_list = rm_list, W = W, scp = scp, name_solvus = name_solvus, maxit = maxit, warm = warm,
                                    max_clamp = max_clamp, bound_tol = bound_tol, log_P = log_P, log_P_2d = log_P_2d,
                                    seismic_cor = seismic_cor, aspect_ratio = aspect_ratio, seismic_water = seismic_water,
                                    shallow_correction = shallow_correction, fluid_as_melt = fluid_as_melt, anelastic_cor = anelastic_cor)
        out, info   = solve_i(sd)
        if info.status != 0 && !isnothing(pool) && one_d
            n0, k0  = info.n_solve, info.n_newton
            for alt in section_retry_seeds(section, x[i], y[i], pool, something(sd, []), Pbox, Tbox)
                out2, info2 = solve_i([alt])
                n0 += info2.n_solve; k0 += info2.n_newton
                if info2.status == 0
                    out, info = out2, info2
                    break
                end
            end
            info = section_info(info.status, n0, k0, info.latent, info.dfree, info.residual)
        end
        outs[i]     = out
        infos[i]    = info

        k = Threads.atomic_add!(count, 1) + 1
        isnothing(progr) || next!(progr)
        if !isnothing(callback_fn) && mod(k, callback_int) == 0
            callback_fn(k, n, time())
        end
    end
    isnothing(progr) || finish!(progr)

    ru, rv  = isnothing(retry_xy) ? (x, y) : (Float64.(retry_xy[1]), Float64.(retry_xy[2]))
    rsc     = isnothing(retry_xy) ? scales : (span(ru), span(rv))
    tried   = [Set{Int}() for _ in 1:n]
    for pass in 1:retry_passes
        bad  = [i for i in 1:n if infos[i].status != 0]
        good = [i for i in 1:n if infos[i].status == 0]
        (isempty(bad) || isempty(good)) && break
        jobs = Tuple{Int,Int,Int}[]
        for b in bad
            d = [((ru[g] - ru[b])/rsc[1])^2 + ((rv[g] - rv[b])/rsc[2])^2 for g in good]
            for (r, g) in enumerate(good[partialsortperm(d, 1:min(3, length(d)))])
                g in tried[b] && continue
                push!(jobs, (b, g, r))
                push!(tried[b], g)
            end
        end
        isempty(jobs) && break
        res = Vector{Tuple{gmin_struct{Float64, Int64}, section_info}}(undef, length(jobs))
        rprog = progressbar ? Progress(length(jobs), desc = "Retrying $(length(bad)) unresolved point(s), pass $(pass)...") : nothing
        section_foreach(length(jobs)) do k, id
            b, g, _     = jobs[k]
            gv          = MAGEMin_db.gv[id]
            z_b         = MAGEMin_db.z_b[id]
            DB          = MAGEMin_db.DB[id]
            splx_data   = MAGEMin_db.splx_data[id]
            per_pt && (gv = define_bulk_rock(gv, X[b], Xoxides, sys_in, MAGEMin_db.db))
            res[k]      = section_point(section, x[b], y[b], gv, z_b, DB, splx_data;
                                        Pbox = Pbox, Tbox = Tbox, seeds = [outs[g]], scales = scales, buffer_n = bufs[b],
                                        rm_list = rm_list, W = W, scp = scp, name_solvus = name_solvus, maxit = retry_maxit,
                                        max_clamp = retry_max_clamp, bound_tol = bound_tol, log_P = log_P, log_P_2d = log_P_2d, warm = warm,
                                        seismic_cor = seismic_cor, aspect_ratio = aspect_ratio, seismic_water = seismic_water,
                                        shallow_correction = shallow_correction, fluid_as_melt = fluid_as_melt, anelastic_cor = anelastic_cor)
            isnothing(rprog) || next!(rprog)
        end
        isnothing(rprog) || finish!(rprog)
        by_pt = Dict{Int,Vector{Int}}()
        for (k, (b, _, _)) in enumerate(jobs)
            push!(get!(by_pt, b, Int[]), k)
        end
        for (b, ks) in by_pt
            sort!(ks, by = k -> jobs[k][3])
            kb     = something(findfirst(k -> res[k][2].status == 0, ks), 0)
            used   = kb == 0 ? ks : ks[1:kb]
            kb     = kb == 0 ? 0 : ks[kb]
            n0     = infos[b].n_solve + sum(k -> res[k][2].n_solve, used)
            k0     = infos[b].n_newton + sum(k -> res[k][2].n_newton, used)
            if kb > 0
                o, inf     = res[kb]
                outs[b]    = o
                infos[b]   = section_info(inf.status, n0, k0, inf.latent, inf.dfree, inf.residual)
            else
                infos[b]   = section_info(infos[b].status, n0, k0, infos[b].latent, infos[b].dfree, infos[b].residual)
            end
        end
    end
    isnothing(callback_fn) || callback_fn(n, n, time())

    return outs, infos
end
