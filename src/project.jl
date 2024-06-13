
using DataFrames, DataFramesMeta, CSV
using Glob
using ProgressMeter
using GenomicFeatures
using BigWig
using Plots, StatsPlots, Measures
using Statistics
#using Clustering,
using ClusterOrderTools

include("chromatlas.jl")
include("regions.jl")
include("plots.jl")
theme(:wong2)

### function to get base dir of project whether working in src or notebooks
function getprojectdir()
    d = pwd()
    if (basename(d) == "notebooks") || (basename(d) == "src")
        return dirname(d)
    else 
        return d ## assume working in project dir
    end
end