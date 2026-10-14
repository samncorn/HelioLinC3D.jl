# function read(path; kwargs...)
#     # TODO: if directory or vector, check formats. 
#     # if all the same, pass to vector read for that format
#     # otherwise read individually, and vcat
#     suffix = split(path, ".")[end]

#     if suffix ∈ ("arr", "arrow")
#         return read_arrow(path; kwargs...)
#     elseif suffix == "csv"
#         return read_CSV(path; kwargs...)
#     elseif suffix == "db"
#         return read_SQL(path; kwargs...)
#     end

#     throw("ERROR: IO.read: unsupported file type $(suffix) passed")
# end

# function get_hypotheses(args, config, t_ref, gm, c; ddr0 = 0.0)
#     if haskey(args, "r") && haskey(args, "dr") && !isnothing(args["r"]) && !isnothing(args["dr"])
#         hypotheses = [hlc.Hypothesis(args["r"], args["dr"], ddr0, t_ref, gm, c)]
#         @info "hypothesis: r: $(args["r"]) dr: $(args["dr"])"
#     elseif (haskey(args, "r") && !isnothing(args["r"])) || (haskey(args, "dr") && !isnothing(args["dr"]))
#         throw("single hypothesis value provided for r or dr, but not both")
#     else
#         hypo_file = haskey(args, "hypothesis") && !isnothing(args["hypothesis"]) ? args["hypothesis"] : config["hypothesis"]["file"]
#         @info "hypothesis file: $(hypo_file)"
#         hypothesis_df = CSV.read(hypo_file, DataFrame)
#         hypotheses    = map(row -> hlc.Hypothesis(row[:r], row[:dr], ddr0, t_ref, gm, c), eachrow(hypothesis_df))
#         @info "number of hypothesis points  : $(length(hypotheses))"
#     end
#     return hypotheses
# end

""" 
find the lines in a horizons vector file corresponding to the data
"""
function get_lines(f)
    start = 0
    stop  = 0
    for (i, l) in enumerate(eachline(f))
        if length(l) >= 5 && l[1:5] == "\$\$SOE" 
            start += i + 1
        end

        if length(l) >= 5 && l[1:5] == "\$\$EOE" 
            stop += i - 1
        end
    end
    return start, stop
end

function read_horizon(f)
    soe, eoe = get_lines(f)
    return CSV.read(f, DataFrame; skipto = soe, limit = eoe-soe, header = soe-3, stripwhitespace=true, select=1:11, ntasks=1)
end

function write_hdf5(path, tbl)
    h5open(path, "w") do fid
        for n in names(tbl)
            write(fid, n, tbl[:, n])
        end
    end
end

function load_observer_interp(file, obs_cols)
    obs_df = read_horizon(file)
    return observer_interpolant(
         obs_df[:, obs_cols["t"]],
        -obs_df[:, obs_cols["x"]],
        -obs_df[:, obs_cols["y"]],
        -obs_df[:, obs_cols["z"]],
        -obs_df[:, obs_cols["dx"]],
        -obs_df[:, obs_cols["dy"]],
        -obs_df[:, obs_cols["dz"]],
    )
end

# """ returns a table like structure (dictionary) to access the contents
# """
# function read_hdf(path)

# end

# function read_CSV(path; kwargs...)
#     if isfile(path)
#         return CSV.read(path, DataFrame; kwargs...)
#     elseif isdir(path)
#         paths = [joinpath(path, f) for f in readdir(path)]
#         return CSV.read(paths, DataFrame; kwargs...)
#     end
#     throw
# end

# """ 
# Should materilaize as a memmory mapped table
# """
# function read_arrow(path::String)
#     if isfile(path)
#         return ar.Table(path) |> DataFrame
#     elseif isdir(path)
#         paths = [joinpath(path, f) for f in readdir(path)]
#         return ar.Table(paths) |> DataFrame
#     end
#     throw("ERROR: read_arrow: Bad path")
# end

# function read_arrow(paths::Vector{String}; kwargs...)
#     return ar.Table(paths; kwargs...) |> DataFrame
# end

# """ Probably shouldn't use. Parquet2 documentation indicates parquet to be a suboptimal format
# """
# function read_parquet(path)
#     ds = Parquet2.Dataset(path)
#     if isdir(path)
#         for f in Parquet2.filelist(ds)
#             append!(ds, f)
#         end
#     end
#     df = DataFrame(ds; copycols=false,)
#     return df
# end

# function read_SQL(path; select=["*"], table_name="observations")
#     con = SQLite.DB(path)
#     cols = join(select, ", ")
#     return DBInterface.execute(con, "SELECT $(cols) FROM $(table_name)") |> DataFrame
# end

# function get_chunk_times(path; time_col = :fieldMJD_TDB)
#     if isdir(path)
#         return sort([get_nights(joinpath(path, f)) for f in readdir(path)], by = x -> x[1])
#     else
#         df = read(path)
#         return (minimum(df[:, time_col]), maximum(df[:, time_col]))
#         # if split(path, ".")[end] == "arr"
#         #     # can get intervals straight from the metadata
#         #     mdata = Arrow.getmetadata(Arrow.Table(path))
#         #     return (mdata)
#         # else
#         #     # have to read the actual data
#         #     df = read(path)
#         #     return (minimum(df[:, time_col]), maximum(df[:, time_col]))
#         # end
#     end
# end