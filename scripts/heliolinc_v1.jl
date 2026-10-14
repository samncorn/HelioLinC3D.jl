using Pkg;
Pkg.activate(@__DIR__);

using ArgParse
using LoggingExtras
using TOML
using Printf
using CSV
using Arrow
using DataFrames
using DataInterpolations
using StaticArrays
using ProgressMeter

import HelioLinC as hl

function (@main)(_args)
    args = parse_cli()
    config = TOML.parsefile(args["config"])

    logger = ConsoleLogger(Logging.Info)

    out_dir = mkpath(args["out"])

    n0     = args["n-start"]
    nf     = args["n-stop"]
    window = args["window"]

    # detections_df = read_detections(args["detections"])
    observer      = load_observer(args["observer"], config["observer"]["columns"])
    # nights        = generate_nights(detections_df, config["detections"]["columns"]) # array of (t0, tf) (or dict(night => (t0, tf))?)
    hypotheses    = [(row["r"], row["dr"]) for row in eachrow(CSV.read(args["hypothesis"], DataFrame))]

    nights = Iterators.map(f -> parse(Int, split(f, '.')[1]), readdir(args["detections"]))
    nights = Iterators.filter(n -> n0 <= n <= nf, nights)

    prog = args["progress-bar"]

    with_logger(logger) do
        @info "clustering radius: $(args["cluster-radius"])"
        @info "$(length(hypotheses)) hypothesis pairs"
        # @info "processing night(s): "
        # for (n, tf) in sort(pairs(nights), by=x -> x[1]) # distributed here
        for n in nights
            out_path = mkpath(joinpath(out_dir, lpad(n, 4, '0')))

            in_files = Iterators.filter(f -> n - window <= parse(Int, split(f, '.')[1]) <= n, readdir(args["detections"]))
            in_paths = Iterators.map(f -> joinpath(args["detections"], f), in_files)
            detections_df = Arrow.Table(collect(in_paths)) |> DataFrame

            t0 = maximum(detections_df[:, config["detections"]["columns"]["t_jd"]]) - 2400000.5
            tf = maximum(detections_df[:, config["detections"]["columns"]["t_jd"]]) - 2400000.5
            t_ref = 0.5(t0 + tf)

            detections = collect(Iterators.map(
                row -> hl.Detection{UInt,Float64,Float64,Float64}(
                    row[config["detections"]["columns"]["id"]],
                    row[config["detections"]["columns"]["t_jd"]] - 2400000.5,
                    hl.sphere2cart(
                        deg2rad(row[config["detections"]["columns"]["ra"]]), 
                        deg2rad(row[config["detections"]["columns"]["dec"]])
                    ),
                    observer(row[config["detections"]["columns"]["t_jd"]] - 2400000.5),
                ), 
                eachrow(detections_df)
            ))

            # detections, (t1, tf) = load_detections(detections_df, observer, tf - window, tf, config["detections"]["columns"])
            # t0 = 0.5(t1 + tf)

            # n_clusters = 0
            p = Progress(length(hypotheses); desc=@sprintf("processing night % 4d ref mjd: %.4f", n, t0))
            Threads.@threads for (r_ref, dr_ref) in hypotheses 
                clusters = hl.cluster_detections(
                    detections,
                    r_ref,
                    dr_ref,
                    t_ref,
                    args["dt-max"],
                    args["e-min"],
                    args["e-max"],
                    args["q-min"],
                    args["Q-max"],
                    args["cluster-radius"],
                    # config["gm"]
                    hl.GM_sun
                ) # multithreaded here (inside clustering algorithm)
        
                if !isnothing(clusters)
                    cluster_ids = [x[1] for x in clusters] #.+ n_clusters
                    cluster_dets_1 = [x[2][1] for x in clusters]
                    cluster_dets_2 = [x[2][2] for x in clusters]

                    # just csvs for the moment
                    out_file = joinpath(out_path, "$(r_ref)_$(dr_ref).csv")
                    open(out_file, "w+") do io
                        header = "r_ref, dr_ref, cluster, det1, det2"
                        write(io, header)
                        CSV.write(io,
                            DataFrame(
                                :r_ref => r_ref,
                                :dr_ref => dr_ref,
                                :cluster => cluster_ids,
                                :source_id1 => cluster_dets_1,
                                :source_id2 => cluster_dets_2
                            ),
                            append=true
                        )
                    end
                    # n_clusters += length(unique(clusters))
                end
                prog && next!(p) #, showvalues=[("found clusters", n_clusters)])
            end # for r, dr
        end # for nights
    end # withlogger
end # main

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
        help = "path to config file"
        arg_type = String
        required = true

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
        required = true

        "--out"
        help = "path to create the ouput directory at. The input file name is injected as the file stem."
        arg_type = String
        default = "out"

        "--n-start"
        arg_type = Int
        default = 1

        "--n-stop"
        arg_type = Int
        default = 1

        "--window"
        arg_type = Float64
        default = 14.0

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

        "--loglvl"
        help = "logging level for output"
        arg_type = Int
        default = 0

        "--progress-bar"
        help = "whether to display a progress bar over hypotheses"
        nargs = 0
        # default = false

    end
    return parse_args(ARGS, s)
end

function read_detections(path; kwargs...)
    if isfile(path)
        return Arrow.Table(path; kwargs...) |> DataFrame
    elseif isdir(path)
        files = [joinpath(path, f) for f in readdir(path)]
        return Arrow.Table(files; kwargs...) |> DataFrame
    end
end

function load_observer(path, cols)
    df = hl.read_horizon(path)
    t = [row[cols["t"]] .- 2400000.5 for row in eachrow(df)]

    x = CubicHermiteSpline(df[:, cols["dx"]], df[:, cols["x"]], t)
    y = CubicHermiteSpline(df[:, cols["dy"]], df[:, cols["y"]], t)
    z = CubicHermiteSpline(df[:, cols["dz"]], df[:, cols["z"]], t)

    return s -> SVector{3}((x(s), y(s), z(s)))
end

function generate_nights(df, cols)
    T1 = eltype(df[:, cols["night"]])
    T2 = eltype(df[:, cols["t_jd"]])

    nights = Dict{T1,T2}()
    for (k, chunk) in pairs(groupby(df, cols["night"]))
        n = k[1]
        tf = maximum(chunk[:, cols["t_jd"]]) - 2400000.5
        nights[n] = tf
    end
    return nights
end

function generate_hypotheses(r_bounds, r_step, dr_bounds, dr_step)
    r = r_bounds[1]:r_step:r_bounds[2]
    dr = dr_bounds[1]:dr_step:dr_bounds[2]

    return [(ri, drj) for ri in r for drj in dr]
end