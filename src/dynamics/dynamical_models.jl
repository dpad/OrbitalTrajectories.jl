#-----------------------#
# GENERIC MODEL METHODS #
#-----------------------#

export AstrodynamicalSystem

import ..cse

struct AstrodynamicalSystem{M<:Abstract_AstrodynamicalModel} <: Abstract_AstrodynamicalSystem{M}
    model :: M
    dynamics :: AstrodynamicalSystemDynamics
end

struct AstrodynamicalSystemDynamics
    transformations :: Vector{Abstract_ModelTransformation}
    ode_system :: ODESystem
    ode_f :: ODEFunction
end

# Build a default system with no model transformations (except default)
AstrodynamicalSystem{M}(args...; transform=[], kwargs...) where {M<:Abstract_AstrodynamicalModel} =
    AstrodynamicalSystem(M(args...; kwargs...), transform)

# wrap_code: perform common substring elimination to improve performance of the resulting models.
@doc """Generic constructor for an AstroDynamicalSystem and its underlying ODEFunction."""
@memoize function ModelDynamics(model::M, transformations::Vector{Abstract_ModelTransformation}; 
        lower_order=true, wrap_code=(cse, cse)) where {M<:Abstract_AstrodynamicalModel}

    # First, build the ODE system for the given model
    ode_system = ODESystem(model)

    # Apply the common order-lowering transform
    if lower_order
        transformations = [OrderLoweringTransform()(model, ode_system), transformations...]
    end

    # Now, apply all the transformations to the system in order
    for transform in transformations
        ode_system = transform(model, ode_system)
    end

    # Generate the functions
    ode_f = ODEFunction(ode_system; 
        # Don't need to differentiate further
        jac=false, 
        tgrad=false, 
        # Make sure the function is evaluated in current module (especially for package precompile)
        eval_expression=false, 
        eval_module=@__MODULE__,
        # Perform any further transformations on the code
        wrap_code
    )

    AstrodynamicalSystem{M}(model, transformations, ode_system, ode_f)
end


@doc """Generic model function."""
(model::Abstract_AstrodynamicalSystem)(args...; kwargs...) = model.ode_f(args...; kwargs...)

#--------------------------#
# ORDER LOWERING TRANSFORM #
#--------------------------#

struct OrderLoweringTransform <: Abstract_StateExtension end

@memoize function (::OrderLoweringTransform)(model::Abstract_AstrodynamicalModel, _ode::ODESystem)
    # Reduce high-order system to 1st-order (re-ordering to keep the original equations first)
    # NOTE: need to expand the RHS derivatives before calling order-lowering
    eqs = expand_derivatives.(equations(_ode))
    num_orig_eqs = length(eqs)
    eqs, dvs = ode_order_lowering(eqs, independent_variable(_ode), ModelingToolkit.states(_ode)) 
    eqs = [eqs[end - num_orig_eqs + 1:end]..., eqs[1:num_orig_eqs]...]
    dvs = [dvs[end - num_orig_eqs + 1:end]..., dvs[1:num_orig_eqs]...]
    ODESystem(simplify.(eqs), independent_variable(_ode), dvs, parameters(_ode))
end

#---------#
# DISPLAY #
#---------#

Base.show(io::IO, ::MIME"text/plain", ::Type{T}) where {T<:Abstract_AstrodynamicalModel} = 
    print(io, nameof(T))

function Base.show(io::IO, M::MIME"text/plain", x::Abstract_AstrodynamicalModel)
    print(io, string(SciMLBase.TYPE_COLOR, nameof(typeof(x)), SciMLBase.NO_COLOR))
    show(io, M, x.props)
    show(io, M, x.ode)
end

# Base.show(io::IO, ::MIME"text/plain", ::Type{T}) where {T<:Abstract_AstrodynamicalODESystem} = 
#     print(io, nameof(T))

# Base.show(io::IO, ::MIME"text/plain", x::Abstract_AstrodynamicalODESystem) = 
#     nothing

ModelingToolkit.varmap_to_vars(model::Abstract_AstrodynamicalModel, varmap) = ModelingToolkit.varmap_to_vars(varmap, parameters(model))
DiffEqBase.isinplace(f::Abstract_AstrodynamicalModel) = true
DiffEqBase.isinplace(f::Abstract_AstrodynamicalModel, _) = isinplace(f)

# XXX: Need these due to new ModelingToolkit interface.
# function Base.getproperty(sys::Abstract_AstrodynamicalODESystem, name::Symbol)
#     # XXX: This is very much needed to avoid the depwarn introduced in
#     # ModelingToolkit, especially in Pkg.test() environments!
#     return getfield(sys, name)
# end
# ModelingToolkit.get_systems(::Abstract_AstrodynamicalODESystem) = []
# ModelingToolkit.get_eqs(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_eqs(f.ode_system)
# ModelingToolkit.get_states(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_states(f.ode_system)
# ModelingToolkit.get_ps(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_ps(f.ode_system)
# Base.nameof(f::Abstract_AstrodynamicalODESystem) = nameof(typeof(f))

state_length(m::Abstract_AstrodynamicalModel) = length(states(m.ode))
state_length(s::State) = state_length(s.system)

#----------#
# INCLUDES #
#----------#
include("models/3bp.jl")
# include("models/4bp.jl")
# include("models/nbp.jl")