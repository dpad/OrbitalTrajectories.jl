#-----------------------#
# GENERIC MODEL METHODS #
#-----------------------#

import ..cse

# wrap_code: perform common substring elimination to improve performance of the resulting models.
@doc """Generic constructor for a DynamicalModel's underlying ODEFunctions."""
function (T::Type{<:Abstract_AstrodynamicalODESystem})(args...; wrap_code=(cse, cse), order=1, kwargs...)
    _ode = ODESystem(T, args...; kwargs...)

    # Reduce high-order system to 1st-order (re-ordering to keep the original equations first)
    # NOTE: need to expand the RHS derivatives before calling order-lowering
    eqs = expand_derivatives.(equations(_ode))
    num_orig_eqs = length(eqs)
    eqs, dvs = ode_order_lowering(eqs, independent_variable(_ode), ModelingToolkit.states(_ode)) 
    eqs = [eqs[end - num_orig_eqs + 1:end]..., eqs[1:num_orig_eqs]...]
    dvs = [dvs[end - num_orig_eqs + 1:end]..., dvs[1:num_orig_eqs]...]
    ode = ODESystem(simplify.(eqs), independent_variable(_ode), dvs, parameters(_ode); name=nameof())

    # Generate variational equations (if any)
    var_eq_systems = variational_equation_ODESystems(ode; order)

    # Generate the functions
    # TODO: Add support for tgrad (need to define derivative(get_pos) for EphemerisNBP)
    ode_f = ODEFunction(ode; 
        # Don't need to differentiate further
        jac=false, 
        tgrad=false, 
        # Make sure the function is evaluated in current module (especially for package precompile)
        eval_expression=false, 
        eval_module=@__MODULE__,
        # Perform any further transformations on the code
        wrap_code
    )

    T(ode, ode_f)
end

@doc """Generic model function."""
(model::Abstract_AstrodynamicalModel)(args...; kwargs...) = model.ode.ode_f(args...; kwargs...)

#---------#
# DISPLAY #
#---------#
Base.show(io::IO, x::Abstract_AstrodynamicalModel) = show(io, typeof(x))
Base.show(io::IO, x::Type{<:Abstract_AstrodynamicalModel}) = print(io, nameof(x))
ModelingToolkit.varmap_to_vars(model::Abstract_AstrodynamicalModel, varmap) = ModelingToolkit.varmap_to_vars(varmap, parameters(model))
DiffEqBase.isinplace(f::Abstract_AstrodynamicalModel) = true
DiffEqBase.isinplace(f::Abstract_AstrodynamicalModel, _) = isinplace(f)

# XXX: Need these due to new ModelingToolkit interface.
function Base.getproperty(sys::Abstract_AstrodynamicalODESystem, name::Symbol)
    # XXX: This is very much needed to avoid the depwarn introduced in
    # ModelingToolkit, especially in Pkg.test() environments!
    return getfield(sys, name)
end
ModelingToolkit.get_systems(::Abstract_AstrodynamicalODESystem) = []
ModelingToolkit.get_eqs(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_eqs(f.ode_system)
ModelingToolkit.get_states(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_states(f.ode_system)
ModelingToolkit.get_ps(f::Abstract_AstrodynamicalODESystem) = ModelingToolkit.get_ps(f.ode_system)
Base.nameof(f::Abstract_AstrodynamicalODESystem) = nameof(typeof(f))

state_length(m::Abstract_AstrodynamicalModel) = length(states(m.ode))
state_length(s::State) = state_length(s.model)

#----------#
# INCLUDES #
#----------#
include("models/3bp.jl")
include("models/4bp.jl")
include("models/nbp.jl")