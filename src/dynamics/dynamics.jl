module Dynamics
    using OrbitalTrajectories.SpiceUtils

    using DiffEqBase
    using DiffEqSensitivity
    using DifferentialEquations
    using FiniteDiff
    using ForwardDiff
    using LinearAlgebra
    using LineSearches
    using LoopVectorization
    using Memoize
    using ModelingToolkit
    using NLsolve
    using PhysicalConstants.CODATA2014: NewtonianConstantOfGravitation
    using Plots
    using ProgressMeter
    using RecipesBase
    using SciMLBase
    using SimpleTraits
    using SPICE
    using StaticArrays
    using Tullio
    using Unitful

    # Required to ensure that we can precompile ODEFunctions (which needs to be
    # done in our own cache).
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)

    #----------------#
    # ABSTRACT TYPES #
    #----------------#

    # Dynamical models
    abstract type Abstract_AstrodynamicalODESystem <: ModelingToolkit.AbstractODESystem end
    abstract type Abstract_VariationalEquationsODESystem{Order} <: Abstract_AstrodynamicalODESystem end
    abstract type Abstract_AstrodynamicalModel end
    abstract type Abstract_ModelProperties end

    # Specific abstract models
    abstract type Abstract_R3BPModel <: Abstract_AstrodynamicalModel end
    abstract type Abstract_R4BPModel <: Abstract_AstrodynamicalModel end

    # Reference frames
    abstract type Abstract_ReferenceFrame end

    # Differential correctors
    abstract type Abstract_DifferentialCorrector end

    # Sensitivities
    abstract type Abstract_StateTransitionTensor{N} end

    #----------#
    # INCLUDES #
    #----------#
    include("orbital_trajectories.jl")
    include("reference_frames.jl")
    include("sensitivities.jl")
    include("variational_equations.jl")
    include("dynamical_models.jl")
    include("differential_correction.jl")
    include("plotting.jl")
end