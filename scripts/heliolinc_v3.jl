using Pkg; Pkg.activate(joinpath(@__DIR__, ".."));

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

import HelioLinC as hl

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "path to config file"
        arg_type = String
        default = joinpath(@__DIR__, "../defaults/config.toml")
        # required = true

        # moved to list of files names (night selection logic offloaded to the user/external config script)
        # "--detections"
        # help = "path to detections directory. requires specific "
        # arg_type = String
        # required = true

        # now specified in config file
        # "--detections-list"
        # help = "text file specifying files to process as detections, in addition to any "

        # moved to config
        # "--observer"
        # help = "path to observer ephem file"
        # arg_type = String

        "--hypothesis"
        help = "path to hypothesis grid file. Overrides config value"
        arg_type = String
        # required = true

        "--r"
        help = "hypothesis radial distance. overrides config value"
        arg_type = Float64
        # required = true

        "--dr"
        help = "hypothesis radial velocity. overrides config calue"
        arg_type = Float64
        # required = true

        # specified by input files now
        # "--night"
        # arg_type = Int
        # required = true

        "--out"
        help = "path to create the ouput directory at. overrides config value"
        arg_type = String
        default = "out"

        # "--n-start"
        # arg_type = Int
        # default = 1

        # "--n-stop"
        # arg_type = Int
        # default = 1

        # managed by input files now
        # "--window"
        # arg_type = Int
        # default = 14

        # "--v-max"
        #     help = "maximum cartesian speed for valid tracklets (AU/day)"
        #     arg_type = Float64
        #     default = 0.018

        "--e-min"
        help = "minimum eccentricity. overrides config."
        arg_type = Float64
        default = 0.0

        "--e-max"
        help = "maximum eccentricity. overrides config."
        arg_type = Float64
        default = 1.0

        "--q-min"
        help = "minimum perihelion distance. overrides config."
        arg_type = Float64
        default = 0.4

        "--Q-max"
        help = "maximum aphelion distance. overrides config."
        arg_type = Float64
        default = Inf

        # "--a-max"
        # help = "maximum semimajor axis to consider"
        # arg_type = Float64
        # default = 5.0

        "--dt-max"
        help = "maximum time seperation for a tracklet (days). overrides config."
        arg_type = Float64
        default = 8.0 / 24.0

        "--cluster-radius"
        help = "dbscan clutsering parameter (AU). overrides config."
        arg_type = Float64
        default = 0.0014

        # "--loglvl"
        # help = "logging level for output"
        # arg_type = Int
        # default = 0

        "--progress-bar"
        help = "whether to display a progress bar over hypotheses."
        nargs = 0
        # default = false

        "--n-chunks"
        help = "number of chunks for clustering parallelization. overrides config."
        arg_type = Int
        default = 1

        "--n-proc"
        help = "number of hypothesis to run in parallel. overrides config."
        arg_type = Int
        default = 1

    end
    return parse_args(ARGS, s)
end

function (@main)(_args)
    args = parse_cli()
    config = TOML.parsefile(args["config"])

    out_path = mkpath(resolve_setting(config["output"]["path"], args["out"], "output path"))

    # set up the logger(s)
    # log_io = stdout
    log_io = open(joinpath(out_path, "log.txt"), "w+")
    logger = ConsoleLogger(log_io)

    det_cols = config["detections"]["columns"]
    obs_cols = config["observer"]["columns"]

    with_logger(logger) do
        # r_ref  = args["r"]
        # dr_ref = args["dr"]
        # hypothesis = args["hypothesis"]

        # window = args["window"]
        # @info "clustering window            : $(window)"

        cluster_radius = resolve_setting(config["clustering"]["radius"], args["cluster-radius"], "clustering radius")
        @info "clustering radius            : $(cluster_radius)"

        # dt_max = args["dt-max"]
        dt_max = resolve_setting(config["tracklet"]["dt_max"], args["dt-max"], "maximimum tracklet timespan")
        @info "maximum tracklet time        : $(dt_max)"

        # e_min = args["e-min"]
        e_min = resolve_setting(config["tracklet"]["e_min"], args["e-min"], "minimum eccentricity"; default = 0.0)
        @info "minimum arrow eccentricity   : $(e_min)"

        e_max = resolve_setting(config["tracklet"]["e_max"], args["e-max"], "maximum eccentricity"; default = 1.0)
        @info "maximum arrow eccentricity   : $(e_max)"

        # q_min = args["q-min"]
        q_min = resolve_setting(config["tracklet"]["q_min"], args["q-min"], "minimum periapse"; default = 0.0)
        @info "minimum arrow perihelion     : $(q_min)"

        # Q_max = args["Q-max"]
        if !haskey(config["tracklet"], "Q_max") 
            config["tracklet"]["Q_max"] = nothing
        end
        Q_max = resolve_setting(config["tracklet"]["Q_max"], args["Q-max"], "maximum apoapse"; default = Inf)
        @info "maximum arrow aphelion       : $(Q_max)"

        n_chunks = resolve_setting(config["parallel"]["dbscan"], args["n-chunks"], "n clustering chunks"; default = 1)

        # wrap orbit limits intop one struct to pass around
        arrow_limits = hl.Arrow_Limits(Q_max, q_min, e_min, e_max, dt_max)

        # wrap clustering settings into a struct
        vel_scale = config["clustering"]["velocity_scaling"]
        cluster_params = hl.Cluster_Params(cluster_radius, vel_scale, 3, n_chunks)

        # det_dir   = config["detections"]["directory"]
        det_files = config["detections"]["files"]
        # det_files = map(f -> joinpath(det_dir, f), det_files)
        det_df = Arrow.Table(det_files) |> DataFrame
        #
        # sort!(det_df, [det_cols["id"]])
        t0 = minimum(det_df[:, det_cols["time"]])
        tf = maximum(det_df[:, det_cols["time"]])
        t_ref = 0.5(t0 + tf)
        @info "reference time: $(t_ref)"

        # instantiate hypotheses
        if !isnothing(args["r"]) && !isnothing(args["dr"])
            hypotheses = [hl.Hypothesis(hl.GM_sun, args["r"], args["dr"], t_ref)]
            @info "hypothesis: r: $(args["r"]) dr: $(args["dr"])"
        elseif !isnothing(args["r"]) || !isnothing(args["dr"])
            throw("single hypothesis value provided for r or dr, but not both")
        else
            hypo_file = resolve_setting(config["hypothesis"]["file"], args["hypothesis"], "hypothesis file")
            hypothesis_df = CSV.read(hypo_file, DataFrame)
            hypotheses    = map(row -> hl.Hypothesis(hl.GM_sun, row[:r], row[:dr], t_ref), eachrow(hypothesis_df))
            @info "hypothesis file: $(hypo_file)"
            @info "number of hypothesis points  : $(length(hypotheses))"
        end

        # load the observer
        # observer_file = resolve_setting(config["observer"]["file"], args["observer"], "observer file")
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
        det_row_map = irow -> hl.Detection(
            irow[1],
            irow[2][det_cols["time"]],
            hl.sphere2cart(
                deg2rad(irow[2][det_cols["ra"]]),
                deg2rad(irow[2][det_cols["dec"]]),
                ),
            observer(irow[2][det_cols["time"]])
            )
        detections = collect(Iterators.map(det_row_map, enumerate(eachrow(det_df))))

        prog = Progress(length(hypotheses), enabled = args["progress-bar"])

        n_proc = resolve_setting(config["parallel"]["hypotheses"], args["n-proc"], "hypotheses processes"; default = 1)
        chunksize = (length(hypotheses) ÷ n_proc) + 1
        chunks = Iterators.partition(hypotheses, chunksize)

        @info "distributing hypotheses across $(n_proc) chunks of $(chunksize) points" 
        
        start_time = now()

        # for task, point to an output file and spinup
        tasks = []
        for (i, chunk) in enumerate(chunks)
            # out_file = joinpath(out_path, "$(i).csv")
            out_file = joinpath(out_path, "$(i)")
            if isfile(out_file)
                rm(out_file)
            end
            task = Threads.@spawn run_hypotheses_csv(out_file*".csv", detections, det_df, chunk, arrow_limits, cluster_params, prog)
            # task = Threads.@spawn run_hypotheses_arrow(out_file, detections, det_df, chunk, arrow_limits, cluster_params, prog)
            push!(tasks, task)
        end

        for task in tasks
            wait(task)
        end

        end_time = now()

        dt = (end_time - start_time).value / 1000

        # @info "completed in $(dt) seconds"

        sec = dt % 60.0
        minutes = floor(Int, dt / 60) % 60
        hours = floor(Int, dt / 3600) % 60

        @info @sprintf "completed in %02dH %02dM %02dS" hours minutes sec
    end
end

function run_hypotheses_csv(out_file, detections, det_df, hypotheses, arrow_limits, cluster_params, prog)
    open(out_file, "w") do io
        # @printf(io, "%-9s %-9s %-7s %-8s %-8s %-14s %-7s %-7s %-7s\n", 
        #     "sourceid", "objid", "fieldid", "ra", "dec", "jd", "r", "dr", "tref",
        #     )
        # "$(row["cluster"]),
        # $(row["sourceID"]),
        # $(row["ObjID"]),
        # $(row["FieldID"]),
        # $(row["RA_deg"]),
        # $(row["Dec_deg"]),
        # $(row["fieldJD_TDB"]),
        # $(row["r"]),
        # $(row["dr"]),
        # # $(row["t"])")
        println(io, "cluster,source_id,obj_id,field_id,ra_deg,dec_deg,field_jd_tdb,r,dr,t_ref")
        offset = 0
        for hyp in hypotheses
            clusters = hl.cluster_hypothesis(detections, hyp, arrow_limits, cluster_params)
            if !isnothing(clusters)
                write_clusters_csv(io, det_df, clusters, hyp; offset=offset)
                offset += length(clusters)
            end

            next!(prog)
        end
    end
end

function run_hypotheses_arrow(out_stem, detections, det_df, hypotheses, arrow_limits, cluster_params, prog)
    open(Arrow.Writer, out_stem*".arr") do io
        offset = 0
        for hyp in hypotheses
            clusters = hl.cluster_hypothesis(detections, hyp, arrow_limits, cluster_params)
            if !isnothing(clusters)
                write_clusters_arrow(io, det_df, clusters, hyp; offset=offset)
                offset += length(clusters)
            end
            
            next!(prog)
        end
    end
end

function write_clusters_arrow(io, det_df, clusters, hyp; offset=0)
    cluster_ids  = Int[]
    cluster_dets = Int[]
    for (i, dets) in enumerate(clusters)
        append!(cluster_ids, [i + offset for _ in dets])
        append!(cluster_dets, dets)
    end

    rows = det_df[cluster_dets, ["sourceID", "ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]
    rows[:, :cluster] = cluster_ids
    rows[:, :r]  .= hyp.r
    rows[:, :dr] .= hyp.dr
    rows[:, :t]  .= hyp.t

    Arrow.write(io, rows)

    # for (i, dets) in enumerate(clusters)
    #     # should really source the column names from a config
    #     rows = det_df[dets, ["sourceID", "ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]

    #     # record cluster and hypothesis info 
    #     rows[:, :cluster] .= i + offset
    #     rows[:, :r]  .= hyp.r
    #     rows[:, :dr] .= hyp.dr
    #     rows[:, :t]  .= hyp.t

    #     Arrow.write(io, rows)
    # end
end

""" the detections dataframe should have unnecessary columns removed externaly
"""
function write_clusters_csv(io, det_df, clusters, hyp; offset=0)
    cluster_ids  = Int[]
    cluster_dets = Int[]
    for (i, dets) in enumerate(clusters)
        append!(cluster_ids, [i + offset for _ in dets])
        append!(cluster_dets, dets)
    end

    rows = det_df[cluster_dets, ["sourceID", "ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]
    rows[:, :cluster] = cluster_ids
    rows[:, :r]  .= hyp.r
    rows[:, :dr] .= hyp.dr
    rows[:, :t]  .= hyp.t
    # open(out_file, "a") do io 
    # CSV.write(io, rows)
    for row in eachrow(rows)
        println(io, "$(row["cluster"]),$(row["sourceID"]),$(row["ObjID"]),$(row["FieldID"]),$(row["RA_deg"]),$(row["Dec_deg"]),$(row["fieldJD_TDB"]),$(row["r"]),$(row["dr"]),$(row["t"])")
        # @printf 
    end
    # end
end

"""
prioiritize cli arg value over config file setting. If neither, and no deault is provided, exit with error.
"""
function resolve_setting(config_val, arg_val, name; default=nothing)
    !isnothing(arg_val)     && return arg_val
    !isnothing(config_val)  && return config_val
    !isnothing(default)     && return default
    throw("$(name) not provided")
end