module MCES

using JuMP
using Gurobi
using MathOptInterface
using Graphs

const MOI  = MathOptInterface
abstract type AbstractMCESModel end

abstract type Config end

include("graph.jl")
include("readFile.jl")
include("modelisation.jl")
include("marencoModel.jl")
#include("symmetricMarencoModel.jl")
include("bahienseModel.jl")
include("degreeModel.jl")
include("reducedModel1.jl")
include("reducedModel2.jl")
include("modelAD.jl")
include("reducedAD.jl")

export Config, ConfigBMPS, ConfigML, ConfigR1, ConfigR2,ConfigD,ConfigAD, ConfigRAD
export run_config_in_thread, run_config, load_graph_ML, load_graph_ARG, clean_graph!, free_loops!
end
