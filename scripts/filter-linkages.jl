using Pkg; Pkg.activate(joinpath(@__DIR__, ".."), io=devnull);

using StaticArrays
using StatsBase
using LinearAlgebra
using ArgParse
# using Distributions
using TOML
using Arrow
using DataFrames
using Logging
using CSV
using ProgressMeter

import heliolinc_core as hlc
import HelioLinC as hl
import Kepler

# """
# Cometary formated orbit elements. Assumes ecliptic reference frame
# """
# struct ComOrb{T}
#     q::T
#     e::T
#     i::T
#     Om::T
#     w::T
#     t::T
#     tp::T
# end

# """ get comentary elemetns from pos, vel
# """
# function ComOrb(pos, vel, t0, gm; R = I)

# end

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "path to config file"
        arg_type = String
        # default = joinpath(@__DIR__, "../defaults/config.toml")

        "--hypothesis"
        help = "path to hypothesis grid file. Overrides config value"
        arg_type = String

        "--out-path"
        help = "output directory path"

        "--out-stem"
        help = "output file stem"

        "--var-cut-scale"
        # default = Inf
        arg_type = Float64
        help     = "limit on variance proportional to the clustering radius to reject clusters"

        "--var-cut-const"
        # default = Inf
        arg_type = Float64
        help     = "constant limit on variance to reject clusters, units of distance"

        # "--deduplicate"
        # action = :store_true
        # help   = "assign at most one cluster for any detection"
        "--deduplicate-times"
        action = :store_true
        help   = "remove from clusters tracklets with overlapping times, keeping the one closest to the clusters mean state"

        "--deduplicate-linkages"
        action = :store_true
        help   = "remove clusters with the same set of detections, keeping only the one with the lowest cluster variance"
        
        "--max-rms"
        # default = Inf
        arg_type = Float64
        help     = "maximum rms of the residuals generated using a clusters mean state, before performing differential correction."

        "--iod-samples"
        arg_type = Int
        help     = "maximum number of samples to try from a cluster for IOD"

        # "--fit"
        # action = :store_truec

        "--reject-time-duplicates"
        action = :store_true

        "--od-steps"
        # default = 0
        arg_type = Int
        help     = "maximum number of differential correction steps to perform"

        "--min-nights"
        arg_type = Int
        help     = "minimum number of nights for a valid linkage"

        "--print-log"
        action = :store_true

        "--cluster-radius"
        arg_type = Float64

        "--timescale"
        arg_type = Float64

        "--n-threads"
        arg_type = Int
    end
    return parse_args(ARGS, s)
end

function (@main)(_args)
    args = parse_cli()
    config = !isnothing(args["config"]) ? TOML.parsefile(args["config"]) : Dict()

    # out_path = !isnothing(args["out-path"]) ? args["out-path"] : haskey(config, "output") && haskey(config["output"], "path") ? config["output"]["path"] : "."
    out_path = resolve_program_parameters(args["out-path"], config, ["output", "path"], "./out")
    if !args["print-log"]
        log_io = open(joinpath(out_path, "$(args["out-stem"])_linkages.log"), "w")
        logger = SimpleLogger(log_io)
        global_logger(logger)
    end

    # resolve command line vs config settings
    cluster_radius  = resolve_program_parameters(args["cluster-radius"],  config, ["clustering", "radius"],          1e-3)
    var_cut_scale   = resolve_program_parameters(args["var-cut-scale"],   config, ["clustering", "var_cut_scale"],   Inf)
    var_cut_const   = resolve_program_parameters(args["var-cut-const"],   config, ["clustering", "var_cut_const"],   Inf)
    # dedupe_times   = resolve_program_parameters(args["dedupe-times"],    config, ["clustering", "dedupe_times"],    false)
    # dedupe_links   = resolve_program_parameters(args["dedupe-linkages"], config, ["clustering", "dedupe_linkages"], false)
    min_nights      = resolve_program_parameters(args["min-nights"],      config, ["clustering", "min_nights"],      3)
    # od_steps       = resolve_program_parameters(args["od-steps"],        config, ["clustering", "od_steps"],        0)
    # iod_samples     = resolve_program_parameters(args["iod-samples"], config, ["clustering", "iod_samples"], 100)
    # max_rms        =
    hyp_path = resolve_program_parameters(args["hypothesis"], config, ["hypothesis", "file"], nothing)

    n_threads = !isnothing(args["n-threads"]) ? args["n-threads"] : Threads.nthreads()

    # var_cut_scale  = !isnothing(args["var-cut-scale"])   ? args["var-cut-scale"]   : haskey(config, "clustering") && haskey(config["clustering"], "var_cut_scale")    ? config["clustering"]["var_cut_scale"]    : Inf
    # var_cut_const  = !isnothing(args["var-cut-const"])   ? args["var-cut-const"]   : haskey(config, "clustering") && haskey(config["clustering"], "var_cut_const")    ? config["clustering"]["var_cut_const"]    : Inf
    # dedupe_times   = !isnothing(args["dedupe-times"])    ? args["dedupe-times"]    : haskey(config, "clustering") && haskey(config["clustering"], "dedupe_times")     ? config["clustering"]["dedupe_times"]     : true
    # dedupe_links   = !isnothing(args["dedupe-linkages"]) ? args["dedupe-linkages"] : haskey(config, "clustering") && haskey(config["clustering"], "dedupe_linkages")  ? config["clustering"]["dedupe_linkages"]  : true
    # max_rms        = !isnothing(args["max-rms"])         ? args["max-rms"]         : haskey(config, "clustering") && haskey(config["clustering"], "max_rms")          ? config["clustering"]["max_rms"]          : Inf
    # od_steps       = !isnothing(args["od-steps"])        ? args["od-steps"]        : haskey(config, "clustering") && haskey(config["clustering"], "od_steps")         ? config["clustering"]["od_steps"]         : 0
    # # iod_tol       = 
    # min_nights     = !isnothing(args["min-nights"])     ? args["min-nights"]     : haskey(config, "clustering") && haskey(config["clustering"], "min-nights")    ? config["clustering"]["min-nights"]    : 3
    max_var = min(var_cut_const, cluster_radius*var_cut_scale)^2

    @info "variance cut off         : $max_var"
    # @info "time scale               : $timescale"
    @info "minimum number of nights : $min_nights"
    # @info "maximum linkage initial rms : $max_rms"
    # @info "number of OD steps          : $od_steps"

    det_cols = config["detections"]["columns"]
    obs_cols = config["observer"]["columns"] 

    # load the detections
    det_files = config["detections"]["files"]
    det_tbl   = Arrow.Table(det_files) |> DataFrame

    ra     = det_tbl[:, det_cols["ra"]]
    dec    = det_tbl[:, det_cols["dec"]]
    times  = det_tbl[:, det_cols["time"]]
    angles = [SVector{3}(hlc.sphere2cart(a, b; deg = true)) for (a, b) in zip(ra, dec)]
    # nights = det_tbl[:, det_cols["night"]]

    weights = [SMatrix{2, 2}((1/deg2rad(sigma), 0.0, 0.0, 1/deg2rad(sigma))) for sigma in det_tbl[:, "astrometricSigma_deg"]]

    t0 = minimum(times)
    tf = maximum(times)
    t_ref = 0.5(t0 + tf)

    # load the observer
    observer = hl.load_observer_interp(config["observer"]["file"], obs_cols)
    obs      = [observer(t) for t in times]

    # load the hypotheses
    # hyp_table  = CSV.read(hyp_path, DataFrame; delim = ' ', ignorerepeated = true)
    # have to denormalize the 3rd component
    # hypotheses = collect(Iterators.map(x -> hlc.Hypothesis(
    #     x["r(au)"], 
    #     x["rdot(AU/day)"], 
    #     -hl.GM_sun/x["r(au)"]^2*x["mean_accel"],
    #     t_ref,
    #     hl.GM_sun, 
    #     hl.C_au_day,
    #     ), eachrow(hyp_table)))


    # clusters_files = if isdir(joinpath(out_path, "clusters"))
    #     readdir(joinpath(out_path, "clusters"); join = true)
    # else
    #     [joinpath(out_path, "clusters.arr")]
    # end

    track_files = sort(filter(s -> split(s, '_')[end] == "tracks.csv", readdir(joinpath(out_path, "clusters"))))
    orbit_files = sort(filter(s -> split(s, '_')[end] == "orbits.csv", readdir(joinpath(out_path, "clusters"))))

    @info "  reducing clusters to linkages"
    prog = Progress(length(track_files); enabled = args["print-log"])
    # linkages = Dict{Vector{Int}, Tuple{Float64, SVector{3, Float64}, SVector{3, Float64}}}()
    _linkages = Vector{Any}(nothing, length(track_files))
    # Threads.@threads for i in 1:length(clusters_files)
    #     f = clusters_files[i]
        # clusters_tbl = Arrow.Table(f) |> DataFrame

    Threads.@threads for (i, (lfile, ofile)) in collect(enumerate(zip(track_files, orbit_files)))
        ltbl = CSV.read(joinpath(out_path, "clusters", lfile), DataFrame; ntasks = 1)
        otbl = CSV.read(joinpath(out_path, "clusters", ofile), DataFrame; ntasks = 1)

        # if nrow(otbl) == 0
        #     continue
        # end

        _linkages[i] = Dict{Vector{Int}, Tuple{Float64, SVector{3, Float64}, SVector{3, Float64}}}()

        orbits = groupby(otbl, "linkage")
        for ((l,), _ltbl) in pairs(groupby(ltbl, "linkage"))
            # link = _ltbl[:, "detection"]
            d1 = _ltbl[:, "detection_1"]
            d2 = _ltbl[:, "detection_2"]
            link = unique(vcat(d1, d2))
            orb  = orbits[(l,)][1, :]

            # if length(unique()) 

            if !haskey(_linkages[i], link) || _linkages[i][link][1] > orb["scatter"]
                _linkages[i][link] = (
                    orb["scatter"],
                    SVector{3}(orb["x"], orb["y"], orb["z"]),
                    SVector{3}(orb["vx"], orb["vy"], orb["vz"])
                )
            end
        end

        # clusters_tbl[!, "cluster_index"] .= 0
        # _linkages[i] = 
        
        # process_clusters(
        #     clusters_tbl,
        #     angles,
        #     times,
        #     obs,
        #     hypotheses;
        #     n_threads = n_threads,
        #     min_nights = min_nights,
        #     reject_duplicate_times = true,
        # )

        # linkages = merge_linkages((linkages, _linkages))
        next!(prog)
    end
    finish!(prog)

    filter!(x -> !isnothing(x), _linkages)
    linkages = merge_linkages(_linkages)
    @info "  $(length(linkages)) linkages to fit orbits"

    # check for sufficient coverage
    filter!(x -> hlc.is_discoverable(sort(det_tbl[x[1], "fieldJD_TDB"]); dt_min = 5/60/24, dt_max = 2/24, n_nights = min_nights), linkages)

    # check for time duplicates
    if args["reject-time-duplicates"]
        filter!(x -> length(unique(det_tbl[x[1], "fieldJD_TDB"])) == length(x[1]), linkages)
    end

    # fit the orbits
    linkages = kepler_fits(
        linkages, angles, weights, times, obs, t_ref; 
        gm = hl.GM_sun, c = hl.C_au_day, n_threads = n_threads, progress = args["print-log"],
        max_iter = 100, mu = 1000.0
    )

    # assign each detection to it's best linkage
    linkages = assign_unique_linkages(linkages)

    # recheck for sufficient coverage
    filter!(x -> hlc.is_discoverable(sort(det_tbl[x[1], "fieldJD_TDB"]); dt_min = 5/60/24, dt_max = 2/24, n_nights = min_nights), linkages)

    # refit remaining
    linkages = kepler_fits(
        linkages, angles, weights, times, obs, t_ref; 
        gm = hl.GM_sun, c = hl.C_au_day, n_threads = n_threads, progress = args["print-log"],
        max_iter = 100, mu = 1000.0
    )

    # check for duplicates

    # write out the remaining linkages
    stem = isnothing(args["out-stem"]) ? "linkages" : "$(args["out-stem"])_linkages"
    write_linkages(linkages, stem, out_path, t_ref) # retain all linkages for analysis, reject based on fit in post

    if !args["print-log"]
        close(log_io)
    end
end

function resolve_program_parameters(arg, config, config_key, default)
    if !isnothing(arg)
        arg
    elseif haskey(config, config_key)
        config[config_key]
    else
        default
    end
end

function resolve_program_parameters(arg, config, config_keys::Union{AbstractVector, Tuple}, default)
    if !isnothing(arg)
        arg
    else
        d = config
        for k in config_keys
            if haskey(d, k)
                d = d[k]
            else
                return default
            end
        end

        return d
    end
end

function get_state(ang1, ang2, dt1, dt2, o1, o2, ref_r, ref_dr, ref_ddr, branch, gm, c)
    r1 = ref_r + ref_dr*dt1 + 0.5ref_ddr*dt1^2
    r2 = ref_r + ref_dr*dt2 + 0.5ref_ddr*dt2^2

    rho1 = hlc.sphere_intersect_branch(o1, ang1, r1; neg_branch = branch)
    rho2 = hlc.sphere_intersect_branch(o2, ang2, r2; neg_branch = branch)

    lt1 = dt1 - rho1/c
    lt2 = dt2 - rho2/c

    pos1 = o1 + rho1*ang1
    pos2 = o2 + rho2*ang2
    pos, vel, dt = hlc.lambert_central_diff(pos1, pos2, lt1, lt2)
    pos, vel = Kepler.propagate(pos, vel, -dt, gm)
    # pos, vel = try 
    #     Kepler.propagate(pos, vel, -dt, hypi.gm)
    # catch e
    #     println("r  : $ref_r")
    #     println("dr : $ref_dr")
    #     println("ddr: $ref_ddr")
    #     println("t0 : $(hypi.t0)")
    #     println("pos: $pos")
    #     println("vel: $vel")
    #     println("dt : $dt")
    #     println("gm : $(hypi.gm)")
    #     println()
    #     throw(e)
    # end

    # assume no ddr0 correction for time being
    # assume orbit bounds from the clustering
    return pos, vel
end

function get_cluster_stats!(vpos, vvel, det1, det2, branches, angles, times, obs, hypi)
    for (i, (j1, j2, branch)) in enumerate(zip(det1, det2, branches))
        pos, vel = get_state(
            angles[j1],
            angles[j2],
            times[j1] - hypi.t0,
            times[j2] - hypi.t0,
            obs[j1],
            obs[j2], 
            hypi.r0,
            hypi.dr0,
            hypi.ddr0,
            branch, 
            hl.GM_sun, 
            hl.C_au_day
            )

        # assume no ddr0 correction for time being
        # assume orbit bounds from the clustering
        vpos[i] = pos
        vvel[i] = vel
    end

    n = length(det1)
    spos = std(view(vpos, 1:n))
    svel = std(view(vvel, 1:n))
    return spos, svel
end

# function process_cluster!(cluster_pos, cluster_vel, linkages, det1, det2, branches, angles, times, nights, obs, hyp; max_var = Inf, min_nights = 3, reject_duplicate_times = true)
function process_cluster!(cluster_pos, cluster_vel, linkages, det1, det2, branches, angles, times, obs, hyp; max_var = Inf, min_nights = 3, reject_duplicate_times = true)
    gm = hl.GM_sun
    c  = hl.C_au_day

    n_cluster = length(det1)
    spos, _ = get_cluster_stats!(cluster_pos, cluster_vel, det1, det2, branches, angles, times, obs, hyp)

    s2 = dot(spos, spos)
    if s2 > max_var
        return nothing
    end
    # passes mean state filter, pass linkage on to iod
    link = sort(unique((vcat(det1, det2))))

    l_times  = times[link]
    # l_nights = nights[link]

    ul_times  = unique(l_times)
    if (
           length(ul_times) < 2min_nights
        # || (reject_duplicate_times && length(ul_times) < length(l_times))
        || !hlc.is_discoverable(sort(ul_times); n_nights = min_nights)
    )

        return nothing
    end

    # just use the cluster mean as the iod guess
    guess_pos = mean(view(cluster_pos, 1:n_cluster))
    guess_vel = mean(view(cluster_vel, 1:n_cluster))

    # deduplicate. record the lowest iod X2 and the initial guess for a linkage
    scatter = norm(spos)
    if (!haskey(linkages, link) || linkages[link][1] < scatter) # && (1 - cdf(Chisq(n_meas - 6), X2) > 0.25)
        linkages[link] = (scatter, guess_pos, guess_vel)
    end

    return nothing
end

function process_clusters(clusters_tbl, angles, times, obs, hypotheses; n_threads = 1, progress = false, max_var = Inf, min_nights = 3, reject_duplicate_times = true)
    grps = groupby(clusters_tbl, [:cluster, :hypothesis])
    grp_nrows = combine(grps, nrow)
    n_max = maximum(grp_nrows[:, :nrow])

    if progress
        @info "read in $(length(grps)) clusters, $(nrow(clusters_tbl)) total tracklets ($(Base.summarysize(clusters_tbl) * 1e-6) Mb)"
    end

    # prog = Progress(length(grps); enabled = progress)

    proc_structs = Vector{Any}(undef, n_threads)
    Threads.@threads for iproc in 1:n_threads
        proc_structs[iproc] = (
            cluster_pos = Vector{SVector{3, Float64}}(undef, n_max),
            cluster_vel = Vector{SVector{3, Float64}}(undef, n_max),
            linkages    = Dict{Vector{Int}, Tuple{Float64, SVector{3, Float64}, SVector{3, Float64}}}()
        )

        for ci in iproc:n_threads:length(grps)
            cluster     = grps[ci]
            det1        = cluster[:, "det1"]
            det2        = cluster[:, "det2"]
            branches    = cluster[:, "branch"]
            hyp         = hypotheses[cluster[1, "hypothesis"]]

            process_cluster!(
                proc_structs[iproc].cluster_pos, 
                proc_structs[iproc].cluster_vel, 
                proc_structs[iproc].linkages, 
                det1, det2, branches, 
                angles, times, obs, 
                hyp; 
                max_var = max_var,
                min_nights = min_nights,
                reject_duplicate_times = reject_duplicate_times,
                )
            # next!(prog)
        end
    end
    # finish!(prog)

    # merge linkage sets
    linkages = merge_linkages(Iterators.map(x -> x.linkages, proc_structs))

    if progress
        @info "$(length(linkages)) linkages to fit orbits"
    end

    return linkages
end

function merge_linkages(chunks)
    linkages = Dict{Vector{Int}, Tuple{Float64, SVector{3, Float64}, SVector{3, Float64}}}()
    for chunk in chunks
        for (l, v) in pairs(chunk)
            if !haskey(linkages, l) || linkages[l][1] > v[1]
                linkages[l] = v
            end
        end
    end
    return linkages
end

function write_linkages(linkages, stem, out_path, t_ref; prog = false)
    linkage_file = open(joinpath(out_path, "$stem.csv"), "w")
    states_file  = open(joinpath(out_path, "$(stem)_orbits.csv"), "w")

    write(states_file,  "linkage,t,x,y,z,vx,vy,vz,rX2\n")
    write(linkage_file, "linkage,detection\n")

    # flatten the linkages, filtering on pval

    # for (i, (c, state)) in enumerate(zip(clusters, orbit_estimates))
    progress_bar = Progress(length(linkages); enabled = prog, desc = "writing out linkages")
    for (i, (linkage, (pval, pos, vel))) in enumerate(pairs(linkages))
        # write the initial orbit guess
        # heliolinc cluster vectors are in ICRF, 
        # elems = ComOrb(pos, vel, t_ref, hl.GM_sun, hl.ICRF)
        write(states_file, "$i,$t_ref,"*join(pos, ',')*','*join(vel, ',')*",$pval"*"\n")

        # write the linkage detections
        for j in linkage
            write(linkage_file, "$i,$j\n")
        end
        next!(progress_bar)
    end
    finish!(progress_bar)

    close(linkage_file)
    close(states_file)
end

function kepler_fits(linkages, angles, weights, times, obs, t0; gm = hl.GM_sun, c = hl.C_au_day, n_threads = 1, progress = true, e_max = 1.0, tol = 1e-8, max_iter = 2000, mu = 1000.0, dt_tol = 1/86400)
    # Threads.@threads 
    prog = Progress(length(linkages), desc = "performing preliminary orbit fits..."; enabled = progress)

    new_linkages = [empty(linkages) for _ in 1:Threads.nthreads()]

    l_keys = collect(keys(linkages))
    # for (i, (link, (_, pos, vel))) in enumerate(pairs(linkages))
    Threads.@threads for ith in 1:n_threads
        for i in ith:n_threads:length(linkages)
            link = l_keys[i]

            _, pos, vel = linkages[link]

            # try to fit a keplerian orbit
            # if the fit exceeds the specified eccentricity, set the score to 0
            # otherwise, set the score to the probability of drawing a larger X^2 value for the measurements given the fit
            # l_angles  = angles[link]
            # l_weights = weights[link]
            l_times   = times[link]
            # l_obs     = obs[link]

            obs2fit = [(
                time        = times[d],
                angles      = angles[d],
                position    = obs[d],
                weight      = weights[d]  
            ) for d in link]

            # obs2fit = [Kepler.angles_observation(angles[l])]

            # pv, X2, _ = hlc.fit_kepler_analytic(l_angles, l_weights, l_times, l_obs, pos, vel, t0, gm, c; tol = tol, max_iter = max_iter, mu = mu)
            pv, resids = try
                xv0 = Kepler.Cartesian(pos, vel, t0, gm)
                _pv = Kepler.fit_kepler(obs2fit, xv0; c = c, rel_tol = 1e-8, max_iter = 10)[1]
                _resids = Iterators.map(o -> Kepler.compute_weighted_residuals(o, _pv, c), obs2fit)
                _pv, _resids
                # hlc.fit_kepler_analytic(l_angles, l_weights, l_times, l_obs, pos, vel, t0, gm, c; tol = tol, max_iter = max_iter, mu = mu)
            catch err
                # throw(err)
                @info "kepler fit failure caught at linkage $i ($err)"
                continue
            end

            X2 = reduce(+, Iterators.map(r -> dot(r, r), resids))
            
            # check for duplicates TODO: replace with outlier rejection scheme?
            # mask = trues(l_times)
            # for (i, j) in Iterators.product(eachindex(l_times), eachindex(l_times))
            #     if i < j && mask[i] && mask[j] && abs(l_times[i] - l_times[j]) < dt_tol 
            #         mask[norm(l_weights[i] * resids[i]) < norm(l_weights[j] * resids[j]) ? i : j] = false
            #     end
            # end

            # if !all(mask)
            #     link      = link[mask]
            #     l_times   = times[link]
            #     l_obs     = obs[link]

            #     pv, X2, _ = try 
            #         # hlc.fit_kepler_analytic(l_angles, l_weights, l_times, l_obs, pos, vel, t0, gm, c; tol = tol, max_iter = max_iter, mu = mu)
            #         Kepler.fit_kepler()
            #     catch err
            #         @info "kepler fit failure caught at linkage $i"
            #         continue
            #     end
            # end

            # posf = SVector{3}(pv[1], pv[2], pv[3])
            # velf = SVector{3}(pv[4], pv[5], pv[6])

            # q, e, inc, Om, w, tp = Kepler.cometary(hl.ICRF2ECL*posf, hl.ICRF2ECL*velf, t0, gm)
            nu   = 2*length(link) - 6
            X2  /= nu 
            # pval = e < e_max ? 1 - cdf(Chisq(nu), X2) : 0.0
            new_linkages[ith][link] = (X2, pv.position, pv.velocity)

            next!(prog)
        end
    end
    finish!(prog)

    final_linkages = empty(linkages)
    for chunk in new_linkages
        for (k, v) in pairs(chunk)
            final_linkages[k] = v
        end
    end

    return final_linkages
end

function assign_unique_linkages(linkage_dict)
    # dets2links = Dict{Int, Tuple{Int, Float64}}()
    linkages = collect(keys(linkage_dict))
    vals = collect(values(linkage_dict))

    # rmap = [0 for _ in unique(Iterators.flatten(linkages))]
    rmap = Dict(i => 0 for i in unique(Iterators.flatten(linkages)))
    # @showprogress desc = "mapping detections to linakges" 
    for (l, link) in enumerate(linkages)
        # map each detection to the linkage which claims it with the best fit
        for i in link
            if rmap[i] == 0 || vals[l][1] < vals[rmap[i]][1]
                rmap[i] = l
            end
        end
    end
    # reconstruct new linkages
    # new_linkages = empty(linkages)
    new_linkages = [Int[] for _ in linkages]
    # @showprogress desc = "uniquely assigning detections to linkages" 
    for (d, l) in pairs(rmap)
        push!(new_linkages[l], d)
    end

    out = Dict(link => val for (link, val) in zip(new_linkages, vals))
    return out
end

# should test based on residual fits?
# function deduplicate_times(linkages, times; tol = 1/86400)
#     for (link, val) in pairs(linkages)
#         ltimes = times[link]
#         perm   = sortperm(ltimes)
#         slink  = link[perm]
#         stimes = ltimes[perm]

#         has_duplicates = false
#         for i in eachindex(stimes)[2:end]
#             dt = stimes[i] - stimes[i-1]
#             if dt < tol
#                 has_duplicates = true
#             end
#         end

#         if !has_duplicates
#             return
#         end

#         i = 1
#         j = 2
#         copies = []
#         while j < length(slinks)
#             while j < length(slink) && times[j] - times[i] < tol
#                 j += 1
#             end

#             for k in (i):(j-1)
#                 # for clink in copies
#             end
            
#         end

#         for clink in copies

#         end
#     end
# end
