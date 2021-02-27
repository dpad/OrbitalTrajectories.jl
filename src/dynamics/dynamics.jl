module Dynamics
    using OrbitalTrajectories.SpiceUtils

    using SPICE
    using Plots
    using Unitful
    using DifferentialEquations
    using DiffEqBase
    using DiffEqSensitivity
    using SciMLBase
    using LinearAlgebra
    using ModelingToolkit
    using NLsolve
    using LineSearches
    using ForwardDiff
    using FiniteDiff
    using Memoize
    using SimpleTraits
    using StaticArrays
    using ProgressMeter
    using RecipesBase
    using PhysicalConstants.CODATA2014: NewtonianConstantOfGravitation

    # Required to ensure that we can precompile ODEFunctions (which needs to be
    # done in our own cache).
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)

    #----------------#
    # ABSTRACT TYPES #
    #----------------#

    # Dynamical models
    abstract type Abstract_DynamicalModel end
    abstract type Abstract_ModelODEFunctions <: ModelingToolkit.AbstractODESystem end
    abstract type Abstract_R3BPModel <: Abstract_DynamicalModel end
    abstract type Abstract_R4BPModel <: Abstract_DynamicalModel end

    # Reference frames
    abstract type Abstract_ReferenceFrame end

    # Differential correctors
    abstract type Abstract_DifferentialCorrector end

    #----------#
    # INCLUDES #
    #----------#
    include("orbital_trajectories.jl")
    include("reference_frames.jl")
    include("sensitivities.jl")
    include("dynamical_models.jl")
    include("differential_correction.jl")
    include("plotting.jl")
end