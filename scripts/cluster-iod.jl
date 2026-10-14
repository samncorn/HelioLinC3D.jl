using Pkg; Pkg.activate(joinpath(@__DIR__, ".."), io=devnull);

using ArgParse
using CSV
using DataFrames
using LoggingExtras
using Printf
using StaticArrays
using ProgressMeter

function (@main)(_args)
    args = parse_cli()
    clusters = CSV.read(args["clusters"], DataFrame)
end

function parse_cli()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--clusters"
        help = "path to cluster file"
        arg_type = String
        required = true
    end
    return parse_args(ARGS, s)
end