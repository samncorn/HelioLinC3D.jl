using Pkg; Pkg.activate(joinpath(@__DIR__, ".."), io=devnull);

using ArgParse
using TOML

using CSV
using Arrow
using HDF5
using DataFrames
using Statistics
using LinearAlgebra
using LoggingExtras
using Printf
using StaticArrays
using ProgressMeter
using Dates
# using InlineStrings

import HelioLinC as hl
import heliolinc_core as hlc

# const default_out_path = "./heliolinc_output"
# const default_cluster_radius   = 0.001
# const default_max_rms          = Inf
# const default_time_scale = 1.0
# const default_dt_min    = 5 * 60.0 / 86400.0 # five minute minimum tests well
# const default_dt_max    = 2.0 / 24.0
# const default_e_min     = 0.0
# const default_e_max     = Inf
# const default_i_min     = 0.0
# const default_i_max     = 90.0
# const default_v_max     = 2.0
# const default_r_min     = 0.0
# # const default_rho_min   = 0.01
# const default_n_threads = 1
# const default_n_proc    = 1
# const default_max_pts   = typemax(Int)
# const default_min_pts   = 3
# const default_angle_iters = 0

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

        "--exclude"
        help = "path to exclusion file, which specifies detections to skip"
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
        # default = default_out_path

        "--v-max"
        help = "maximum angular apparent motion for valid tracklets (deg/day)"
        arg_type = Float64
        # default = default_v_max

        "--v-min"
        help = "minimum angular apparent motion for valid tracklets (deg/day)"
        arg_type = Float64

        "--e-min"
        help = "minimum eccentricity. overrides config."
        arg_type = Float64
        # default = default_e_min

        "--e-max"
        help = "maximum eccentricity. overrides config."
        arg_type = Float64
        # default = default_e_max

        "--i-min"
        help = "minimum inclination (degrees). overrides config."
        arg_type = Float64
        # default = default_i_min

        "--i-max"
        help = "maximum inclination (degrees). overrides config."
        arg_type = Float64
        # default = default_i_max

        "--r-min"
        help = "minimum perihelion distance. overrides config."
        arg_type = Float64
        # default = default_r_min

        "--range-min"
        help = "minimum range to accept a tracklet projection"
        arg_type = Float64
        # default = 0.01

        "--dt-min"
        help = "minimum time seperation for a tracklet (days). overrides config. useful if detections on the same image may not agree by a tolerable amount"
        arg_type = Float64
        # default = default_dt_min

        "--dt-max"
        help = "maximum time seperation for a tracklet (days). overrides config."
        arg_type = Float64
        # default = default_dt_max

        "--cluster-radius"
        help = "dbscan clutsering parameter (AU). overrides config."
        arg_type = Float64
        # default = default_cluster_radius

        "--timescale"
        help = "value to scale velocity vector by"
        arg_type = Float64
        # default  = default_velocity_scaling

        "--progress-bar"
        help = "whether to display a progress bar over hypotheses."
        nargs = 0

        "--n-threads"
        help = "number of threads to give each hypothesis (NOTE: Parallel hypotheses will EACH get this many threads). overrides config."
        arg_type = Int
        # default = default_n_threads

        "--n-proc"
        help = "number of hypothesis to run in parallel. overrides config."
        arg_type = Int
        # default = default_n_proc

        # "--clustering"
        # help = "vectors to cluster states in. either "
        # arg_type = String
        # default = "rv"

        # "--model"
        # help = "heliocentric distance model"
        # arg_type = String
        # default = "quadratic"

        "--max-pts"
        help = "maximum number of points an arrow can have in its neighborhood before terminating the hypothesis"
        arg_type = Int
        # default  = typemax(Int)

        "--min-pts"
        help = "minimum number of points to generate a cluster"
        arg_type = Int
        # default  = 3

        "--max-rms"
        help = "max rms error for a cluster to pass inspection. same units as clustering radius"
        arg_type = Float64
        # default = default_max_rms

        "--debug"
        nargs = 0

        # "--dedupe"
        # help   = "run cluster deduplication"
        # action = :store_true

        "--angle-iters"
        help = "number of tracklet angluar iterations to use for estimating the ddr term"
        # default = default_angle_iters
    end
    return parse_args(ARGS, s)
end

function (@main)(_args)
    args   = parse_cli()
    config = TOML.parsefile(args["config"])

    out_path = resolve_program_parameters(args["out"], config, ["output", "path"], "./out")
    mkpath(out_path)
    mkpath(joinpath(out_path, "clusters"))

    for f in readdir(joinpath(out_path, "clusters"); join = true)
        rm(f)
    end

    # file_stem = "csv"
    file_stem = "arr"

    # set up the logger(s)
    open(joinpath(out_path, "heliolinc.log"), "w") do log_io
        log_lvl = args["debug"] ? Logging.Debug : Logging.Info
        logger = MinLevelLogger(FileLogger(log_io), log_lvl)

        det_cols = config["detections"]["columns"]
        obs_cols = config["observer"]["columns"]

        with_logger(logger) do
            cluster_radius  = resolve_program_parameters(args["cluster-radius"], config, ["clustering", "radius"],    1e-3)
            max_pts   = resolve_program_parameters(args["max-pts"], config, ["clustering", "max_pts"], typemax(Int))
            min_pts   = resolve_program_parameters(args["min-pts"], config, ["clustering", "min_pts"], 3)
            n_proc    = resolve_program_parameters(args["n-proc"],    config, ["parallel", "hypotheses"], 1)
            n_threads = resolve_program_parameters(args["n-threads"], config, ["parallel", "dbscan"],     1)
            dt_min    = resolve_program_parameters(args["dt-min"],    config, ["tracklet", "dt_min"],    5 * 60.0 / 86400.0)
            dt_max    = resolve_program_parameters(args["dt-max"],    config, ["tracklet", "dt_max"],    2.0/24.0)
            e_min     = resolve_program_parameters(args["e-min"],     config, ["tracklet", "e_min"],     0.0)
            e_max     = resolve_program_parameters(args["e-max"],     config, ["tracklet", "e_max"],     1.0)
            i_min     = resolve_program_parameters(args["i-min"],     config, ["tracklet", "i_min"],     0.0)
            i_max     = resolve_program_parameters(args["i-max"],     config, ["tracklet", "i_max"],     90.0)
            r_min     = resolve_program_parameters(args["r-min"],     config, ["tracklet", "r_min"],     0.0)
            v_min     = resolve_program_parameters(args["v-min"],     config, ["tracklet", "v_min"],     0.05)
            v_max     = resolve_program_parameters(args["v-max"],     config, ["tracklet", "v_max"],     3.0)
            rho_min   = resolve_program_parameters(args["range-min"], config, ["tracklet", "range_min"], 0.01)

            # hyp_path = !isnothing(args["hypothesis"]) ? args["hypothesis"] : config["hypothesis"]["file"]
            hyp_path = resolve_program_parameters(args["hypothesis"], config, ["hypothesis", "file"], nothing)

            # rms_max = haskey(config["tracklet"], "max_rms") && args["max-rms"] == default_v_max  ? config["tracklet"]["max_rms"] : args["max-rms"]
            # angle_iters = haskey(config["tracklet"], "angle_iters") && args["angle-iters"] == default_angle_iters  ? config["tracklet"][ "angle_iters"] : args["angle-iters"]
            rms_max     = resolve_program_parameters(args["max-rms"], config, ["tracklet", "max_rms"], Inf)
            angle_iters = resolve_program_parameters(args["angle-iters"], config, ["tracklet", "angle_iters"], 0)


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
            @info "minimum number of cluster points: $(min_pts)"
            @info "maximum number of cluster points: $(max_pts)"

            det_files = config["detections"]["files"]
            det_df = Arrow.Table(det_files) |> DataFrame

            ra     = det_df[:, det_cols["ra"]]
            dec    = det_df[:, det_cols["dec"]]
            times  = det_df[:, det_cols["time"]]
            angles = [SVector{3}(hlc.sphere2cart(a, b; deg = true)) for (a, b) in zip(ra, dec)]
            # nights = 

            @info "loading $(nrow(det_df)) detections"
            # sort!(det_df, [det_cols["id"]])
            t0 = minimum(times)
            tf = maximum(times)
            t_ref = 0.5(t0 + tf)
            @info @sprintf "reference time: %.6f" t_ref
            timescale = resolve_program_parameters(args["timescale"], config, ["clustering", "timescale"], (tf - t0)/4) # default value from timespan
            # instantiate hypotheses
            hyp_table  = CSV.read(hyp_path, DataFrame; delim = ' ', ignorerepeated = true)
            # println(hyp_table)
            # have to denormalize the 3rd copmp
            hypotheses = collect(Iterators.map(x -> hlc.Hypothesis(
                x["r(au)"], 
                x["rdot(AU/day)"], 
                -hl.GM_sun/x["r(au)"]^2*x["mean_accel"],
                t_ref,
                hl.GM_sun, 
                hl.C_au_day,
                ), eachrow(hyp_table)))

            cluster_params = hlc.ClusteringParams(cluster_radius, timescale, n_threads, rms_max, min_pts, max_pts)
            orbit_bounds   = hlc.OrbitBounds(r_min, e_min, e_max, rho_min, i_min, i_max)

            # hyp_start, hyp_end = resolve_program_parameters(nothing, config, ["hypothesis", "range"], (1, nrow(hypotheses)))
            # hypotheses = hypotheses[hyp_start:hyp_end]

            # load the observer
            observer = hl.load_observer_interp(config["observer"]["file"], obs_cols)
            obs = try
                [observer(t) for t in times]
            catch _
                throw("loading observer failed. time range = [$t0, $tf], config = $(args["config"])")
            end
            
            # convert dataframe columns to detections structs

            # exclude path is not configurable via config file, only command line
            # exclude_path = args["exclude"]
            exclude_path = resolve_program_parameters(args["exclude"], config, ["detections", "exclude_path"], nothing)
            @info "loading preattributed detections to exclude from $(exclude_path)..."

            exclude = if isnothing(exclude_path)
                nothing
            else
                # exclude_tbl = CSV.read(exclude_path, 
                exc = BitVector([false for _ in eachrow(det_df)])
                open(exclude_path, "r") do exc_io
                    for l in eachline(exc_io)
                        ll = parse(Int, l)
                        exc[ll] = true
                    end
                end
                @info "$(sum(exc)) preattributed detections excluded"
                exc
            end

            

            prog = Progress(length(hypotheses), enabled = args["progress-bar"])
            chunksize = (length(hypotheses) ÷ n_proc) + 1 
            @info "distributing hypotheses across $(n_proc) processes of $(chunksize) points"
            
            start_time = now()
            # construct tracklets
            @info "generating tracklets ..."
            # tracklets = hlc.find_tracklets(
            #     detections, 
            #     deg2rad(v_max),
            #     dt_min, dt_max;
            #     n_threads = Threads.nthreads(),
            #     progress = args["progress-bar"]
            # )
        
            tracklets = if isnothing(exclude)
                hlc.find_tracklets(angles, times, deg2rad(v_min), deg2rad(v_max), dt_min, dt_max; n_threads = Threads.nthreads(), progress = args["progress-bar"])
            else
                hlc.find_tracklets(angles, times, deg2rad(v_min), deg2rad(v_max), dt_min, dt_max, exclude; n_threads = Threads.nthreads(), progress = args["progress-bar"])
            end

            @info "found $(length(tracklets)) potential tracklets"

            ch = Channel{Int}(length(hypotheses))
            for i in 1:length(hypotheses)
                put!(ch, i)
            end
            # Threads.@threads for i_proc in 1:n_proc
            @sync for i_proc in 1:n_proc
                Threads.@spawn begin

                    thread_log_io = open(joinpath(out_path, "$(i_proc).log"), "w") 
                    thread_logger = MinLevelLogger(FileLogger(thread_log_io), log_lvl)

                    # writer = open(Arrow.Writer, joinpath(out_path, "clusters-$i_proc.$(file_stem)"))

                    # i_tbls = []
                    with_logger(thread_logger) do
                        while isready(ch)
                            # i_labels = []
                            # i_arrows = []
                            # i_trkids = []
                            # label_offset = 0
                            # for i in i_proc:n_proc:length(hypotheses)
                            i = try 
                                take!(ch)
                            catch take_err
                                continue
                            end
                            hyp = hypotheses[i]
                            @info "hypothesis r: $(hyp.r0), dr: $(hyp.dr0), ddr $(hyp.ddr0)"

                            # ====== #
                            labels, arrows, trkids, _ = hlc.heliolinc(
                            # linkages, states, scatters = hlc.heliolinc(
                                tracklets, 
                                angles, 
                                times, 
                                obs, 
                                hyp, 
                                cluster_params, 
                                orbit_bounds; 
                                n_threads = n_threads, 
                                R2ECL = hl.ICRF2ECL, 
                                angle_iters = angle_iters,
                            )

                            # replace at some point with custom table struct
                            det1 = [tracklets[ij][1] for ij in trkids]
                            det2 = [tracklets[ij][2] for ij in trkids]
                            # ======= #

                            df = DataFrame(hypothesis = i, cluster = labels, det1 = det1, det2 = det2, states = arrows)
                            # #---End---#

                            lfile = open(joinpath(out_path, "clusters", "$(i)_tracks.csv"), "w")
                            write(lfile, "linkage,detection_1,detection_2\n")

                            ofile = open(joinpath(out_path, "clusters", "$(i)_orbits.csv"), "w")
                            write(ofile, "linkage,t,x,y,z,vx,vy,vz,scatter\n")

                            # n_clusters = 0
                            clusters   = groupby(df, ["hypothesis", "cluster"])
                            @info "$(length(clusters)) clusters before validation"
                            for (l, ctbl) in enumerate(clusters)
                                link = vcat(ctbl[:, :det1], ctbl[:, :det2])
                                unique!(link)

                                ltimes = times[link]
                                perm   = sortperm(ltimes)
                                link   = link[perm]
                                ltimes = ltimes[perm]
                                # lra    = ra[link]
                                # ldec   = dec[link]
                                
                                # if !hlc.is_discoverable(ltimes, lra, ldec; n_nights = min_pts, dt_min = dt_min, dt_max = dt_max, v_min = v_min, v_max = v_max)
                                if !hlc.is_discoverable(ltimes; n_nights = min_pts, dt_min = dt_min, dt_max = dt_max)
                                    continue
                                end

                                # n_clusters += 1
                                for row in eachrow(ctbl)
                                    write(lfile, "$l,$(row["det1"]),$(row["det2"])\n")
                                end

                                orbit = mean(ctbl[:, "states"])
                                sigma = norm(std(ctbl[:, "states"]))

                                write(ofile, "$l,$t_ref,$(orbit[1]),$(orbit[2]),$(orbit[3]),$(orbit[4]),$(orbit[5]),$(orbit[6]),$sigma\n")
                            end
                            # 

                            n_clusters = countlines(joinpath(out_path, "clusters", "$(i)_orbits.csv")) - 1
                            @info "$(n_clusters) linkages after enforcing timespan coverage"

                            close(lfile)
                            close(ofile)

                            next!(prog)
                            if !isready(ch) # we have exhausted the queue, and can close, releasing any threads waiting
                                close(ch)
                            end
                        end
                    end
                    close(thread_log_io)
                end
                
            end
            finish!(prog)

            # why is this?
            # final_file = joinpath(out_path, "clusters.$(file_stem)")
            # if isfile(final_file) 
            #     # have to check for a pre-existing file with this name (ie overwriting a previour run)
            #     # otherwise it gets caught by the temprary file search and deleted
            #     rm(final_file)
            # end

            end_time = now()
            dt = end_time - start_time
            @info "completed in $(canonicalize(dt))"
        end
    end # io
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
