#-----------------------#
# GENERIC MODEL METHODS #
#-----------------------#

import ..cse

# wrap_code: perform common substring elimination to improve performance of the resulting models.
@doc """Generic constructor for a DynamicalModel's underlying ODEFunctions."""
function (T::Type{<:Abstract_ModelODEFunctions})(args...; wrap_code=(cse, cse), kwargs...)
    _ode = ODESystem(T, args...; kwargs...)

    # Reduce high-order system to 1st-order (re-ordering to keep the original equations first)
    # NOTE: need to expand the RHS derivatives before calling order-lowering
    eqs = expand_derivatives.(equations(_ode))
    num_orig_eqs = length(eqs)
    eqs, dvs = ode_order_lowering(eqs, independent_variable(_ode), ModelingToolkit.states(_ode)) 
    eqs = [eqs[end - num_orig_eqs + 1:end]..., eqs[1:num_orig_eqs]...]
    dvs = [dvs[end - num_orig_eqs + 1:end]..., dvs[1:num_orig_eqs]...]
    ode = ODESystem(simplify.(eqs), independent_variable(_ode), dvs, parameters(_ode))

    # Generate the functions
    # TODO: Add support for tgrad (need to define derivative(get_pos) for EphemerisNBP)
    ode_f = ODEFunction(ode; jac=true, tgrad=false, eval_expression=false, eval_module=@__MODULE__, wrap_code)
    ode_stm_f = STM_ODEFunction(ode, ode_f; wrap_code)
    T(ode, ode_f, ode_stm_f)
end

@doc """Generic model function."""
(model::Abstract_DynamicalModel)(args...; kwargs...) = model.ode.ode_f(args...; kwargs...)

@traitdef HasJacobian{X}
@traitimpl HasJacobian{X} <- has_jacobian(X)
has_jacobian(X::Type{<:DiffEqBase.ODEFunction}) = !isnothing(fieldtype(X, :jac))

@doc """
    Function to generate an ODE Function that computes the State Transition
    Matrix simultaneously with the given function f. This requires f to have
    a computable Jacobian function.
"""
@traitfn function STM_ODEFunction(ode::ModelingToolkit.AbstractODESystem, f::::HasJacobian; kwargs...)
    # TODO: Output just the variational equations, join as a SplitODEProblem
    # TODO: Memoize this function! It's very slightly slow for EphemerisNBP

    # Convert Variable -> Operation
    iv  = independent_variable(ode)
    dvs = states(ode)
    params = parameters(ode)

    # NOTE: Depends on the Jacobian, corresponding to A(t) matrix (for the
    # State Transition Matrix) [Parker & Anderson 2014].
    @variables   ϕ[1:length(dvs),1:length(dvs)](iv)
    D = Differential(iv)

    # Get the Jacobian matrix (A(t))
    A = Base.invokelatest(f.jac, Num.(dvs), Num.(params), Num(iv))

    # The State Transition Matrix (STM) ODE function is defined as follows, including the N^2 Jacobian equations +
    # the N first-order equations of motion. [Koon 2011]
    # NOTE: the Differential is defined element-wise and flattened to a list.
    stm_eqs = collect(simplify.(D.(ϕ) .~ A * ϕ))

    # Create the ODE system and generate its functions
    stm_ode = ODESystem(
        [equations(ode)..., stm_eqs...], # Append the ODE equations.
        iv,
        [dvs..., ϕ...],  # Append the STM and motion state variables
        params)
    stm_f = ODEFunction(stm_ode; sparse=true, eval_expression=false, eval_module=@__MODULE__, kwargs...)
end

#---------#
# DISPLAY #
#---------#
Base.show(io::IO, x::Abstract_DynamicalModel) = show(io, typeof(x))
Base.show(io::IO, x::Type{<:Abstract_DynamicalModel}) = print(io, nameof(x))
ModelingToolkit.varmap_to_vars(model::Abstract_DynamicalModel, varmap) = ModelingToolkit.varmap_to_vars(varmap, parameters(model))
DiffEqBase.isinplace(f::Abstract_DynamicalModel, _) = true

# XXX: Need these due to new ModelingToolkit interface.
function Base.getproperty(sys::Abstract_ModelODEFunctions, name::Symbol)
    # XXX: This is very much needed to avoid the depwarn introduced in
    # ModelingToolkit, especially in Pkg.test() environments!
    return getfield(sys, name)
end
ModelingToolkit.get_systems(::Abstract_ModelODEFunctions) = []
ModelingToolkit.get_eqs(f::Abstract_ModelODEFunctions) = ModelingToolkit.get_eqs(f.ode_system)
ModelingToolkit.get_states(f::Abstract_ModelODEFunctions) = ModelingToolkit.get_states(f.ode_system)
ModelingToolkit.get_ps(f::Abstract_ModelODEFunctions) = ModelingToolkit.get_ps(f.ode_system)
Base.nameof(f::Abstract_ModelODEFunctions) = nameof(typeof(f))

#----------#
# INCLUDES #
#----------#
include("models/3bp.jl")
include("models/4bp.jl")
include("models/nbp.jl")