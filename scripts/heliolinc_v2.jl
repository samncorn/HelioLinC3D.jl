using Pkg; Pkg.activate(joinpath(@__DIR__, ".."));

using ArgParse
using TOML

using CSV
using Arrow
using DataFrames

using LoggingExtras
using Printf

using DataInterpolations
using StaticArrays
using ProgressMeter

using Dates

import HelioLinC as hl

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "path to config file"
        arg_type = String
        default = joinpath(@__DIR__, "../defaults/config.toml")
        # required = true

        "--detections"
        help = "path to detections file. Can be a directory"
        arg_type = String
        required = true

        "--observer"
        help = "path to observer ephem file"
        arg_type = String
        required = true

        "--hypothesis"
        help = "path to hypothesis grid file"
        arg_type = String
        # required = true

        "--r"
        help = "hypothesis radial distance"
        arg_type = Float64
        # required = true

        "--dr"
        help = "hypothesis radial velocity"
        arg_type = Float64
        # required = true

        "--night"
        arg_type = Int
        required = true

        "--out"
        help = "path to create the ouput directory at. The input file name is injected as the file stem."
        arg_type = String
        default = "out"

        # "--n-start"
        # arg_type = Int
        # default = 1

        # "--n-stop"
        # arg_type = Int
        # default = 1

        "--window"
        arg_type = Int
        default = 14

        # "--v-max"
        #     help = "maximum cartesian speed for valid tracklets (AU/day)"
        #     arg_type = Float64
        #     default = 0.018

        "--e-min"
        help = "minimum eccentricity"
        arg_type = Float64
        default = 0.0

        "--e-max"
        help = "maximum eccentricity"
        arg_type = Float64
        default = 1.0

        "--q-min"
        help = "minimum perihelion distance"
        arg_type = Float64
        default = 0.4

        "--Q-max"
        help = "maximum aphelion distance"
        arg_type = Float64
        default = Inf

        # "--a-max"
        # help = "maximum semimajor axis to consider"
        # arg_type = Float64
        # default = 5.0

        "--dt-max"
        help = "maximum time seperation for a tracklet (days)"
        arg_type = Float64
        default = 8.0 / 24.0

        "--cluster-radius"
        help = "dbscan clutsering parameter (AU)"
        arg_type = Float64
        default = 0.0014

        # "--loglvl"
        # help = "logging level for output"
        # arg_type = Int
        # default = 0

        "--progress-bar"
        help = "whether to display a progress bar over hypotheses"
        nargs = 0
        # default = false

        "--n-chunks"
        help = "number of chunks for clustering parallelization"
        arg_type = Int
        default = 1

        "--p"
        help = "number of hypothesis to run in parallel"
        arg_type = Int
        default = 1

    end
    return parse_args(ARGS, s)
end

function (@main)(_args)
    args = parse_cli()
    config = TOML.parsefile(args["config"])

    out_path = mkpath(args["out"])

    # set up the logger(s)
    # log_io = stdout
    log_io = open(joinpath(out_path, string(args["night"])*".log"), "w+")
    logger = ConsoleLogger(log_io)

    # hard code cols atm. should eventually move to a config file
    det_cols = config["detections"]["columns"]
    obs_cols = config["observer"]["columns"]

    with_logger(logger) do

        # r_ref  = args["r"]
        # dr_ref = args["dr"]
        # hypothesis = args["hypothesis"]

        window = args["window"]
        @info "clustering window            : $(window)"

        cluster_radius = args["cluster-radius"]
        @info "clustering radius            : $(cluster_radius)"

        dt_max = args["dt-max"]
        @info "maximum tracklet time        : $(dt_max)"

        e_min = args["e-min"]
        @info "minimum arrow eccentricity   : $(e_min)"

        e_max = args["e-max"]
        @info "maximum arrow eccentricity   : $(e_max)"

        q_min = args["q-min"]
        @info "minimum arrow perihelion     : $(q_min)"

        Q_max = args["Q-max"]
        @info "maximum arrow aphelion       : $(Q_max)"

        # wrap orbit limits intop one struct to pass around
        arrow_limits = hl.Arrow_Limits(Q_max, q_min, e_min, e_max, dt_max)

        # wrap clustering settings into a struct
        vel_scale = config["clustering"]["velocity_scaling"]
        cluster_params = hl.Cluster_Params(cluster_radius, vel_scale, 3, args["n-chunks"])

        nf = args["night"]
        n0 = nf - window

        night_files = Iterators.filter(f -> n0 <= parse(Int, split(f, '.')[1]) <= nf, readdir(args["detections"]))
        night_files = Iterators.map(f -> joinpath(args["detections"], f), night_files)

        det_df = Arrow.Table(collect(night_files)) |> DataFrame
        sort!(det_df, [det_cols["id"]])
        t0 = minimum(det_df[:, det_cols["time"]])
        tf = maximum(det_df[:, det_cols["time"]])
        t_ref = 0.5(t0 + tf)

        # id0 = det_df[1, det_cols[:id]] # ? doesn't get used?

        # instantiate hypotheses
        if !isnothing(args["hypothesis"]) && (!isnothing(args["r"]) || !isnothing(args["dr"]))
            throw("ERROR: both --hypothesis and --r or --dr specified, exiting...")
        end

        if isnothing(args["hypothesis"])
            hypotheses = [hl.Hypothesis(hl.GM_sun, args["r"], args["dr"], t_ref)]
            # hypothesis_pairs = [(1, (args["r"], args["dr"]))]
            # @info "number of hypothesis points  : $(length(hypothesis_pairs))"
            @info "hypothesis: r: $(args["r"]) dr: $(args["dr"])"
        else
            hypothesis_df    = CSV.read(args["hypothesis"], DataFrame)
            hypotheses       = map(row -> hl.Hypothesis(hl.GM_sun, row[:r], row[:dr], t_ref), eachrow(hypothesis_df))
            # hypothesis_pairs = collect(enumerate(hypothesis))
            @info "number of hypothesis points  : $(length(hypotheses))"
        end

        # load the observer
        obs_df = hl.read_horizon(args["observer"])
        observer = hl.observer_interpolant(
            obs_df[:, obs_cols["t"]],
            obs_df[:, obs_cols["x"]],
            obs_df[:, obs_cols["y"]],
            obs_df[:, obs_cols["z"]],
            obs_df[:, obs_cols["dx"]],
            obs_df[:, obs_cols["dy"]],
            obs_df[:, obs_cols["dz"]],
            )

        # convert dataframe columns to detections structs
        detections = Iterators.map(
            irow -> hl.Detection(
                # row[det_cols[:id]],
                irow[1],
                irow[2][det_cols["time"]],
                hl.sphere2cart(
                    deg2rad(irow[2][det_cols["ra"]]),
                    deg2rad(irow[2][det_cols["dec"]]),
                    ),
                observer(irow[2][det_cols["time"]])
                ),
            enumerate(eachrow(det_df))
            )
        detections = collect(detections)

        prog = Progress(length(hypotheses), enabled = args["progress-bar"])

        chunksize = (length(hypotheses) ÷ args["p"]) + 1
        chunks = Iterators.partition(hypotheses, chunksize)

        @info "distributing hypotheses across $(args["p"]) chunks of $(chunksize) points" 
        
        start_time = now()

        # for task, point to an output file and spinup
        tasks = []
        for (i, chunk) in enumerate(chunks)
            out_file = joinpath(out_path, "$(i).csv")
            if isfile(out_file)
                rm(out_file)
            end
            task = Threads.@spawn run_hypotheses_csv(out_file, detections, det_df, chunk, arrow_limits, cluster_params, prog)
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
    for hyp in hypotheses
        clusters = hl.cluster_hypothesis(detections, hyp, arrow_limits, cluster_params)
        if !isnothing(clusters)
            write_clusters_csv(out_file, det_df, clusters, hyp)
        end

        next!(prog)
    end
end

function run_hypotheses_arrow(out_stem, detections, det_df, hypotheses, arrow_limits, cluster_params, prog)
    open(Arrow.Writer, out_stem*".arr") do io
        offset = 0
        for hyp in hypotheses
            clusters = hl.cluster_hypothesis(detections, hyp, arrow_limits, cluster_params)
            write_clusters_arrow(io, det_df, clusters, hyp; offset=offset)

            offset += length(clusters)
            next!(prog)
        end
    end
end

function write_clusters_arrow(io, det_df, clusters, hyp; offset=0)
    for (i, dets) in enumerate(clusters)
        # should really source the column names from a config
        rows = det_df[dets, ["sourceID", "ObjID", "FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB"]]

        # record cluster and hypothesis info 
        rows[:, :cluster] .= i + offset
        rows[:, :r]  .= hyp.r
        rows[:, :dr] .= hyp.dr
        rows[:, :t]  .= hyp.t

        Arrow.write(io, rows)
    end
end

""" the detections dataframe should have unnecessary columns removed externaly
"""
function write_clusters_csv(out_file, det_df, clusters, hyp)
    open(out_file, "a") do io 
        for dets in clusters
            rows = det_df[dets, :]
            # write(io, "r  : $(hyp.r)\n")
            # write(io, "dr : $(hyp.dr)\n")
            # write(io, "t  : $(hyp.t)\n")

            @printf(io, "r  : %f\n", hyp.r)
            @printf(io, "dr : %f\n", hyp.dr)
            @printf(io, "t  : %f\n", hyp.t)
            # write(io, names(rows))
            # CSV.write(io, rows; append=true)

            @printf(io, "%-9s %-9s %-7s %-8s %-8s %-14s\n", "sourceid", "objid", "fieldid", "ra", "dec", "jd")
            for row in eachrow(rows)
                # for val in ("sourceID", "ObjID","FieldID", "RA_deg", "Dec_deg", "fieldJD_TDB")
                # # for val in names(det_df)
                #     write(io, string(row[val])*',')
                # end
                @printf(io, "%-9i %-9s %-7i %8.3f %8.3f %14.6f\n", row["sourceID"], row["ObjID"], row["FieldID"], row["RA_deg"], row["Dec_deg"], row["fieldJD_TDB"])
                # write(io, "\n")
            end

            write(io, "\n")
        end
    end
end