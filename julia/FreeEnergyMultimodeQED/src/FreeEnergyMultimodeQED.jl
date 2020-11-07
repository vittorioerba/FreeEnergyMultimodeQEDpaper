module FreeEnergyMultimodeQED

ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux"

using Random                    # utils
using LinearAlgebra             # model
using CPLEX, Dates, Ipopt, JuMP # optimiziers

export *
# export fractionGenerator
# export gradMagnetizationEntropy, gradModelFreeEnergy, magnetizationEntropy, modelDimension, modelEntropy, modelFreeEnergy, modelInteraction
# export modelFindMinimum

include("utils.jl")      # utility functions
include("model.jl")      # model-related functions
include("optimizers.jl") # optimizers to find minima of the free energy of the model



end
