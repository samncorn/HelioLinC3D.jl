using Pkg; Pkg.activate(joinpath(@__DIR__, ".."), io=devnull);

using ArgParse
using TOML

using CSV
using Arrow
using DataFrames

using LoggingExtras
using Printf

using StaticArrays
using ProgressMeter

using Dates
using InlineStrings

import heliolinc_core as hlc
import HelioLinC as hl

const default_out_path = "hl_out"
const default_cluster_radius = 0.001
const default_max_rms        = Inf
const default_velocity_scaling = 1.0
const default_dt_min = 60.0 / 86400.0
const default_dt_max = 2.0 / 24.0
const default_e_min  = 0.0
const default_e_max  = Inf
const default_i_min  = 0.0
const default_i_max  = 90.0
const default_v_max  = 2.0
const default_r_min  = 0.0
const default_n_threads = 1
const default_n_proc   = 1

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "path to config file"
        arg_type = String
        default = joinpath(@__DIR__, "../defaults/config.toml")

        "--hypothesis"
        help = "path to hypothesis grid file. Overrides config value"
        arg_type = String

        "--r"
        help = "hypothesis radial distance. overrides config value"
        arg_type = Float64

        "--dr"
        help = "hypothesis radial velocity. overrides config calue"
        arg_type = Float64

        "--out"
        help = "path to create the ouput directory at. overrides config value"
        arg_type = String
        default = default_out_path

        "--v-max"
        help = "maximum angular apparent motion for valid tracklets (deg/day)"
        arg_type = Float64
        default = default_v_max

        "--e-min"
        help = "minimum eccentricity. overrides config."
        arg_type = Float64
        default = default_e_min

        "--e-max"
        help = "maximum eccentricity. overrides config."
        arg_type = Float64
        default = default_e_max

        "--i-min"
        help = "minimum inclination (degrees). overrides config."
        arg_type = Float64
        default = default_i_min

        "--i-max"
        help = "maximum inclination (degrees). overrides config."
        arg_type = Float64
        default = default_i_max

        "--r-min"
        help = "minimum perihelion distance. overrides config."
        arg_type = Float64
        default = default_r_min

        "--dt-min"
        help = "minimum time seperation for a tracklet (days). overrides config. useful if detections on the same image may not agree by a tolerable amount"
        arg_type = Float64
        default = default_dt_min

        "--dt-max"
        help = "maximum time seperation for a tracklet (days). overrides config."
        arg_type = Float64
        default = default_dt_max

        "--cluster-radius"
        help = "dbscan clutsering parameter (AU). overrides config."
        arg_type = Float64
        default = default_cluster_radius

        "--velocity-scaling"
        help = "value to scale velocity vector by"
        arg_type = Float64
        default  = default_velocity_scaling

        "--progress-bar"
        help = "whether to display a progress bar over hypotheses."
        nargs = 0

        "--n-threads"
        help = "number of threads to give each hypothesis (NOTE: Parallel hypotheses will EACH get this many threads). overrides config."
        arg_type = Int
        default = default_n_threads

        "--n-proc"
        help = "number of hypothesis to run in parallel. overrides config."
        arg_type = Int
        default = default_n_proc

        "--clustering"
        help = "vectors to cluster states in. either "
        arg_type = String
        default = "rv"

        "--model"
        help = "heliocentric distance model"
        arg_type = String
        default = "quadratic"

        "--max-pts"
        help = "maximum number of points an arrow can have in its neighborhood before terminating the hypothesis"
        arg_type = Int
        default  = typemax(Int)

        "--min-pts"
        help = "minimum number of points to generate a cluster"
        arg_type = Int
        default  = 3

        "--max-rms"
        help = "max rms error for a cluster to pass inspection. same units as clustering radius"
        arg_type = Float64
        default = default_max_rms

        "--debug"
        nargs = 0

    end
    return parse_args(ARGS, s)
end

function (@main)(_args)
    args   = parse_cli()
    config = TOML.parsefile(args["config"])

    out_path = haskey(config["output"], "path") && args["out"] == default_out_path ? config["output"]["path"] : args["out"]
    mkpath(out_path)

    # set up the logger(s)
    open(joinpath(out_path, "heliolinc.log"), "w") do log_io
    
    log_lvl = args["debug"] ? Logging.Debug : Logging.Info
    logger = MinLevelLogger(FileLogger(log_io), log_lvl)

    det_cols = config["detections"]["columns"]
    obs_cols = config["observer"]["columns"]

    with_logger(logger) do
        cluster_radius = haskey(config["clustering"], "radius") && args["cluster-radius"] == default_cluster_radius ? config["clustering"]["radius"] : args["cluster-radius"]
        vel_scale      = haskey(config["clustering"], "velocity_scaling") && args["velocity-scaling"] == default_velocity_scaling ? config["clustering"]["velocity_scaling"] : args["velocity-scaling"]

        n_proc    = haskey(config["parallel"], "hypotheses") && args["n-proc"] == default_n_proc ? config["parallel"]["hypotheses"] : args["n-proc"]
        n_threads  = haskey(config["parallel"], "dbscan") && args["n-threads"] == default_n_threads ? config["parallel"]["dbscan"] : args["n-threads"]
        
        dt_min = haskey(config["tracklet"], "dt_min") && args["dt-min"] == default_dt_min ? config["tracklet"]["dt_min"] : args["dt-min"]
        dt_max = haskey(config["tracklet"], "dt_max") && args["dt-max"] == default_dt_max ? config["tracklet"]["dt_max"] : args["dt-max"] 
        e_min  = haskey(config["tracklet"], "e_min")  && args["e-min"]  == default_e_min  ? config["tracklet"]["e_min"]  : args["e-min"]
        e_max  = haskey(config["tracklet"], "e_max")  && args["e-max"]  == default_e_max  ? config["tracklet"]["e_max"]  : args["e-max"] 
        i_min  = haskey(config["tracklet"], "i_min")  && args["i-min"]  == default_i_min  ? config["tracklet"]["i_min"]  : args["i-min"]
        i_max  = haskey(config["tracklet"], "i_max")  && args["i-max"]  == default_i_max  ? config["tracklet"]["i_max"]  : args["i-max"] 
        r_min  = haskey(config["tracklet"], "r_min")  && args["r-min"]  == default_i_min  ? config["tracklet"]["r_min"]  : args["r-min"]
        v_max  = haskey(config["tracklet"], "v_max")  && args["v-max"]  == default_v_max  ? config["tracklet"]["v_max"]  : args["v-max"]

        cluster_params = hl.Cluster_Params(cluster_radius, args["max-rms"],vel_scale, args["min-pts"], args["max-pts"], n_threads)
        # cluster_params = hl.Cluster_Params(cluster_radius, Inf, vel_scale, args["min-pts"], args["max-pts"], n_threads)
        arrow_limits = hl.Arrow_Limits(r_min, e_min, e_max, i_min, i_max, dt_min, dt_max)

        state_clustering = if args["clustering"] == "rv"
            hlc.KeplerianPropagation
        elseif args["clustering"] == "he"
            hlc.MomentumEccentricity
        else
            @error "improperly specifed clustering method : $(args["clustering"])"
            exit(1)
        end

        distance_model = if args["model"] == "linear"
            hlc.LinearDistance
        elseif args["model"] == "quadratic"
            hlc.QuadraticDistance
        elseif args["model"] == "analytic"
            hlc.AnalyticDistance
        else   
            @error "improperly specifed distance model : $(args["model"])"
            exit(1)
        end

        @info "clustering radius            : $(cluster_radius)"
        @info "minimum tracklet time        : $(dt_min)"
        @info "maximum tracklet time        : $(dt_max)"
        @info "maximum tracklet velocity    : $(v_max)"
        @info "minimum arrow eccentricity   : $(e_min)"
        @info "maximum arrow eccentricity   : $(e_max)"
        @info "minimum arrow inclination    : $(i_min)"
        @info "maximum arrow inclination    : $(i_max)"
        @info "minimum arrow perihelion     : $(r_min)"
        @info "tracklet maximum on sky rate of motion: $(v_max) (deg/time)"
        @info "tracklet maximum time separation: $(dt_max)"
        @info "propagating distance with $(distance_model) model"
        @info "clustering with $(state_clustering) method"

        det_files = config["detections"]["files"]
        det_df = Arrow.Table(det_files) |> DataFrame
        @info "loading $(nrow(det_df)) detections"

        # sort!(det_df, [det_cols["id"]])
        t0 = minimum(det_df[:, det_cols["time"]])
        tf = maximum(det_df[:, det_cols["time"]])
        t_ref = 0.5(t0 + tf)
        @info @sprintf "reference time: %.6f" t_ref

        # instantiate hypotheses
        if !isnothing(args["r"]) && !isnothing(args["dr"])
            hypotheses = [hl.Hypothesis(hl.GM_sun, args["r"], args["dr"], t_ref)]
            @info "hypothesis: r: $(args["r"]) dr: $(args["dr"])"
        elseif !isnothing(args["r"]) || !isnothing(args["dr"])
            throw("single hypothesis value provided for r or dr, but not both")
        else
            # hypo_file = resolve_setting(config["hypothesis"]["file"], args["hypothesis"], "hypothesis file")
            hypo_file = haskey(args, "hypothesis") ? args["hypothesis"] : config["hypothesis"]["file"]
            hypothesis_df = CSV.read(hypo_file, DataFrame)
            hypotheses    = map(row -> hl.Hypothesis(hl.GM_sun, row[:r], row[:dr], t_ref), eachrow(hypothesis_df))
            @info "hypothesis file: $(hypo_file)"
            @info "number of hypothesis points  : $(length(hypotheses))"
        end

        # load the observer
        observer_file = config["observer"]["file"]
        obs_df = hl.read_horizon(observer_file)
        observer = hl.observer_interpolant(
             obs_df[:, obs_cols["t"]],
            -obs_df[:, obs_cols["x"]],
            -obs_df[:, obs_cols["y"]],
            -obs_df[:, obs_cols["z"]],
            -obs_df[:, obs_cols["dx"]],
            -obs_df[:, obs_cols["dy"]],
            -obs_df[:, obs_cols["dz"]],
            )

        # convert dataframe columns to detections structs
        detections = collect(Iterators.map(row -> det_row_map(row, observer, det_cols["time"], det_cols["ra"], det_cols["dec"]), eachrow(det_df)))

        prog = Progress(length(hypotheses), enabled = args["progress-bar"])
        chunksize = (length(hypotheses) ÷ n_proc) + 1 
        @info "distributing hypotheses across $(n_proc) processes of $(chunksize) points"
        
        start_time = now()
        # construct tracklets
        @info "generating tracklets ..."
        tracklets = hlc.find_tracklets(
            detections, 
            deg2rad(v_max),
            dt_min, dt_max;
            n_threads = Threads.nthreads(),
            progress = args["progress-bar"]
        )
        @info "found $(length(tracklets)) potential tracklets"
        Threads.@threads for i_proc in 1:n_proc
            for i in i_proc:n_proc:length(hypotheses)
                hyp = hypotheses[i]
                next!(prog)
                # try 
                    run_hypotheses_csv(
                        out_path, 
                        tracklets, 
                        detections, 
                        det_df, 
                        hyp, 
                        arrow_limits, 
                        cluster_params, 
                        "$(i).csv",
                        state_clustering = state_clustering,
                        distance_model   = distance_model,
                        )
                # catch e
                #     if isa(e, DomainError)
                #         @info "r:$(hyp.r), dr:$(hyp.dr), exceeded max points ($(args["max-pts"])) in cluster limit"
                #     else
                #         throw(e)
                #     end
                #     continue
                # end
            end
        end

        end_time = now()
        dt = end_time - start_time
        @info "completed in $(canonicalize(dt))"
    end
    end # io
end

function run_hypotheses_csv(path, tracklets, detections, det_df, hyp, arrow_limits, cluster_params, f_name; state_clustering = hlc.KeplerianPropagation, distance_model = hlc.QuadraticDistance)
    t_start = now()
    clusters = hlc.cluster_detections(
        detections, tracklets,
        hyp.r, hyp.dr, hyp.t, hyp.gm,
        cluster_params.cluster_radius;
        scaling = cluster_params.vel_scale,
        state_clustering = state_clustering,
        state_estimate   = hlc.CentralDifference1, 
        distance_model   = distance_model, 
        c = hl.C_au_day,
        r_min = arrow_limits.r_min,
        e_min = arrow_limits.e_min,
        e_max = arrow_limits.e_max,
        n_threads = cluster_params.n_threads,
        ddr0_tol = 1e-8,
        min_pts = cluster_params.min_pts,
        max_pts = cluster_params.max_pts,
        max_rms = cluster_params.max_rms,
    )
    t_stop = now()
    if !isnothing(clusters) 
        dets, arrows, times, labels = clusters
        @info "r:$(hyp.r), dr:$(hyp.dr), found $(length(labels)) in $(canonicalize(t_stop - t_start))"

        out_tbl = DataFrame(
            :cluster    => labels,
            :source_id1 => [t[1] for t in dets],
            :source_id2 => [t[2] for t in dets],
            :t   => times,
            :x   => [x[1] for x in arrows],
            :y   => [x[2] for x in arrows],
            :z   => [x[3] for x in arrows],
            :vx  => [x[4] for x in arrows],
            :vy  => [x[5] for x in arrows],
            :vz  => [x[6] for x in arrows],
            :t0  => hyp.t,
            :r0  => hyp.r,
            :dr0 => hyp.dr,
        )
        CSV.write(joinpath(path, f_name), out_tbl)
    end
end

# function run_hypotheses_arrow(out_stem, detections, det_df, hypotheses, arrow_limits, cluster_params, prog)
#     open(Arrow.Writer, out_stem*".arr") do io
#         offset = 0
#         for hyp in hypotheses
#             clusters = hl.cluster_hypothesis(detections, hyp, arrow_limits, cluster_params)
#             if !isnothing(clusters)
#                 write_clusters_arrow(io, det_df, clusters, hyp; offset=offset)
#                 offset += length(clusters)
#             end
            
#             next!(prog)
#         end
#     end
# end

# function write_clusters_arrow(io, det_df, clusters, hyp; offset=0)
#     cluster_ids  = Int[]
#     cluster_dets = Int[]
#     for (i, dets) in enumerate(clusters)
#         append!(cluster_ids, [i + offset for _ in dets])
#         append!(cluster_dets, dets)
#     end

#     rows = det_df[cluster_dets, ["ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]
#     rows[:, :cluster] = cluster_ids
#     rows[:, :r]  .= hyp.r
#     rows[:, :dr] .= hyp.dr
#     rows[:, :t]  .= hyp.t

#     Arrow.write(io, rows)
# end

""" 
Just write out tracklet i, j and arrow ranges
"""
function write_clusters_csv(io, clusters, hyp; offset = 0)
    for (c, cluster) in enumerate(clusters)
        trks, arrows, times = cluster
        for ((i, j), arrow, t) in zip(trks, arrows, times)
            println(io, @sprintf("%i,%i,%i,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,", c, i, j, t, arrow...)*"$(hyp.t),$(hyp.r),$(hyp.dr)")
        end
    end
end

# """ the detections dataframe should have unnecessary columns removed externaly
# """
# function write_clusters_csv(io, det_df, clusters, hyp; offset=0)
#     cluster_ids  = Int[]
#     cluster_dets = Int[]
#     for (i, dets) in enumerate(clusters)
#         append!(cluster_ids, [i + offset for _ in dets])
#         append!(cluster_dets, dets)
#     end

#     rows = det_df[cluster_dets, ["sourceID", "ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]
#     rows[:, :cluster] = cluster_ids
#     rows[:, :r]  .= hyp.r
#     rows[:, :dr] .= hyp.dr
#     rows[:, :t]  .= hyp.t
#     # open(out_file, "a") do io 
#     # CSV.write(io, rows)
#     for row in eachrow(rows)
#         println(io, "$(row["cluster"]),$(row["sourceID"]),$(row["ObjID"]),$(row["FieldID"]),$(row["RA_deg"]),$(row["Dec_deg"]),$(row["fieldJD_TDB"]),$(row["r"]),$(row["dr"]),$(row["t"])")
#         # @printf 
#     end
#     # end
# end

function det_row_map(row, observer_interp, t_col, ra_col, dec_col)
    t = row[t_col]
    ra  = deg2rad(row[ra_col])
    dec = deg2rad(row[dec_col])
    # ang  = hlc.sphere2cart(ra, dec)
    ang = SVector{2}(ra, dec)
    pos = observer_interp(t)
    return hlc.AnglesMeas(t, ang, pos)
end

# function cluster_iod(det_df, cluster, hypothesis)

# end