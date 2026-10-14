module HelioLinC

# IO
using CSV
import Arrow as ar
using SQLite
using DataFrames

# math
using LinearAlgebra
# using DifferentialEquations
using DataInterpolations
using Rotations
using StatsBase
# using Statistics
using Random
using dbscan
using Kepler
# using heliolinc_core
import heliolinc_core as hlc
# data structures
import NearestNeighbors as nn

# misc
using Printf
using LoggingExtras
using StaticArrays
using ProgressMeter
using Reexport

struct MaximumIterationError 
    n :: Int
end

# @reexport using heliolinc_core

include("HelioLinC/constants.jl")
# include("HelioLinC/config.jl")
include("HelioLinC/IO.jl")
# include("HelioLinC/least_squares.jl")
# include("HelioLinC/herget.jl")
include("HelioLinC/interp.jl")
# include("HelioLinC/linkages.jl")

end # module HelioLinC
