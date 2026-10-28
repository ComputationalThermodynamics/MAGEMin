#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
using Test
using MAGEMin_C
using MAGEMin_C: MAGEMIN_DATABASES, db_infos, ss_infos, harvest_db_infos, SOLVUS_FAMILY,
                 get_mineral_name,
                 harvest_oxide_list, harvest_default_dataset, get_oxide_list,
                 harvest_phase_oxide_support

"""
    _solvus_sweep(db, ss, levels, vars)

    Every name `get_mineral_name` returns for `(db, ss)` over the full grid of `levels`
    applied to the compositional variables `vars`. A function, not inline test-body code:
    a `@testset` block executes in global scope, where the loop variables are untyped
    globals and the same sweep runs ~25x slower.
"""
function _solvus_sweep(db::String, ss::String, levels::NTuple{N,Float64}, vars::NTuple{M,Int}) where {N,M}
    names = Set{String}()
    x     = zeros(12)
    sv    = (compVariables = x,)
    for k in 0:(N^M - 1)
        m = k
        for v in vars
            x[v] = levels[m % N + 1]
            m    = m ÷ N
        end
        push!(names, get_mineral_name(db, ss, sv))
    end
    return names
end

const _REGEN_HINT = "run `make USE_MPI=0 lib && julia --project=. gen/generate_db_infos.jl` " *
                    "to regenerate julia/db_infos_generated.jl"

@testset "database registry vs compiled library" begin

    @testset "library version" begin
        v = get_MAGEMin_version()
        @test v isa String && !isempty(v)
        @test get_MAGEMin_version() === v

        data = Initialize_MAGEMin("mp", verbose = false)
        data = use_predefined_bulk_rock(data, 0)
        out  = point_wise_minimization(4.0, 400.0, data)
        Finalize_MAGEMin(data)
        @test out.MAGEMin_ver == v
    end

    @testset "Warr (2021) mapping covers the Stixrude phases" begin
        # every sb solution phase carries a citation tag (sb_gss_function.c) and a CSV row
        for reg in MAGEMIN_DATABASES
            reg.research_group == "sb" || continue
            inf = retrieve_solution_phase_information(reg.db_name)
            tag = uppercase(reg.db_name)
            for s in inf.data_ss
                @test s.ss_fName == s.ss_name * "_" * tag
                @test get_Warr_name(s.ss_fName) != s.ss_fName
                # `mw` deliberately differs: bare it is the THERMOCALC buffer, tagged it is
                # the Stixrude ferropericlase solution phase - see the explicit checks below
                s.ss_name == "mw" && continue
                @test get_Warr_name(s.ss_fName) == get_Warr_name(s.ss_name)
            end
        end
        # bare `mw` stays the THERMOCALC magnetite-wuestite buffer, the tagged one is ferropericlase
        @test get_Warr_name("mw")      == "mw"
        @test get_Warr_name("mw_SB11") == "Fper"
        # names containing an underscore that is not a citation tag are untouched
        @test get_Warr_name("H_SUCCINa") == "H_SUCCINa"
    end

    @testset "SOLVUS_FAMILY matches get_mineral_name" begin
        solvus_family = filter(kv -> is_db(kv.first[1]), SOLVUS_FAMILY)

        for ((db, ss), fam) in solvus_family
            @test is_db(db)
            inf = retrieve_solution_phase_information(db)
            @test ss in inf.ss_name || ss in [s.ss_fName for s in inf.data_ss]
            @test !isempty(fam) && allunique(fam)
        end

        levels = (-0.2, 0.005, 0.2, 0.55, 1.2)
        vars   = (1, 2, 3, 4, 6, 7, 8)

        seen = Dict((db, ss) => _solvus_sweep(db, ss, levels, vars) for ((db, ss), _) in solvus_family)

        for ((db, ss), fam) in solvus_family
            extra = setdiff(seen[(db, ss)], fam)
            isempty(extra) || error("SOLVUS_FAMILY[(\"$db\", \"$ss\")] is missing $(sort(collect(extra))) " *
                                    "- get_mineral_name can return it")
            @test isempty(extra)
        end
    end

    @testset "registry is well formed" begin
        names = [d.db_name for d in MAGEMIN_DATABASES]
        @test length(unique(names)) == length(names)
        @test all(d -> d.research_group in ("tc", "sb", "gh", "br"), MAGEMIN_DATABASES)
        for rg in ("tc", "sb", "gh", "br")
            codes = [d.EM_database for d in MAGEMIN_DATABASES if d.research_group == rg]
            @test length(unique(codes)) == length(codes)
        end
        for d in MAGEMIN_DATABASES
            @test is_db(d.db_name)
            @test get_db_label(d.db_name)      == d.db_info
            @test get_research_group(d.db_name) == d.research_group
            @test has_dataset_choice(d.db_name) == (d.research_group == "tc")
        end
        @test get_db_list() == names
        @test get_db_list(research_group = "sb") == ["sb11", "sb21", "sb24"]
        @test_throws ErrorException get_db("not_a_database")
    end

    for reg in MAGEMIN_DATABASES
        dtb = reg.db_name

        @testset "$dtb" begin
            live = harvest_db_infos(dtb)
            ref  = retrieve_solution_phase_information(dtb)

            @test ref.db_name    == live.db_name
            @test ref.db_info    == live.db_info
            @test ref.db_dataset == live.db_dataset
            @test ref.dataset_opt == live.dataset_opt

            @test ref.ss_name == live.ss_name  || error("$dtb: solution-phase list drifted, $_REGEN_HINT")
            @test ref.data_pp == live.data_pp  || error("$dtb: pure-phase list drifted, $_REGEN_HINT")
            @test length(ref.data_ss) == length(live.data_ss)

            if length(ref.data_ss) == length(live.data_ss)
                for (i, ss) in enumerate(live.data_ss)
                    for f in fieldnames(ss_infos)
                        got, want = getfield(ref.data_ss[i], f), getfield(ss, f)
                        got == want || error("$dtb: data_ss[$i] ($(ss.ss_name)).$f drifted, $_REGEN_HINT\n" *
                                             "  committed: $got\n  library:   $want")
                        @test got == want
                    end
                end
            end

            ox = harvest_oxide_list(dtb)
            @test get_oxide_list(dtb) == ox || error("$dtb: get_oxide_list drifted from gv.ox\n" *
                                                     "  get_oxide_list: $(get_oxide_list(dtb))\n  gv.ox:          $ox")

            sup_live = harvest_phase_oxide_support(dtb)
            sup_ref  = get_phase_oxide_support(dtb)
            @test Set(keys(sup_ref)) == Set(keys(sup_live)) ||
                  error("$dtb: get_phase_oxide_support phase list drifted, $_REGEN_HINT")
            for ph in sort(collect(keys(sup_live)))
                sup_ref[ph] == sup_live[ph] ||
                    error("$dtb: get_phase_oxide_support[\"$ph\"] drifted, $_REGEN_HINT\n" *
                          "  committed: $(sort(collect(sup_ref[ph])))\n" *
                          "  library:   $(sort(collect(sup_live[ph])))")
                @test sup_ref[ph] == sup_live[ph]
            end

            @test Set(keys(sup_live)) == Set(vcat(live.ss_name, live.data_pp))

            if reg.research_group == "tc"
                @test harvest_default_dataset(dtb) == reg.db_dataset
                @test reg.db_dataset in available_TC_ds
                @test all(d -> d in available_TC_ds, reg.dataset_opt)
            end
        end
    end
end
