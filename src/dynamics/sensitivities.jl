export AD, FD, VE, sensitivity, sensitivity_trace, extract_STMs, extract_stability, has_variational_equations

#---------------------------#
# SENSITIVITIES (i.e. STMs) #
#---------------------------#

# TODO: Use Symbols instead of these constants?
const AD = ForwardDiff
const FD = FiniteDiff
const VE = Val(:VariationalEquations)

@traitdef HasVE{X}
@traitimpl HasVE{X} <- has_variational_equations(X)
has_variational_equations(X::Type{<:Abstract_ModelODEFunctions}) = hasfield(X, :ode_stm_f)
has_variational_equations(X::Type{<:Abstract_DynamicalModel}) = has_variational_equations(fieldtype(X, :ode))
has_variational_equations(X::Type{<:State}) = has_variational_equations(fieldtype(X, :model))

@doc """ Sensitivity of the final state with respect to initial state. """
sensitivity(m::Module, args...; kwargs...) = sensitivity(Val(first(fullname(m))), args...; kwargs...)
function sensitivity(::Val{:ForwardDiff}, state, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
    ForwardDiff.jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        return convert_to_frame(solve(new_state, alg; kwargs...), desired_frame).sol[end]
    end
end
function sensitivity(::Val{:FiniteDiff}, state, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
    FiniteDiff.finite_difference_jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        return convert_to_frame(solve(new_state, alg; kwargs...), desired_frame).sol[end]
    end
end
function sensitivity(v::Val{:VariationalEquations}, state, args...; kwargs...)
    STM_VE = sensitivity_trace(v, state, args...; kwargs...).sol[end]
    dim = length(state.u0)
    return reshape(STM_VE[dim+1:end], (dim, dim))
end

@doc """ Sensitivity of the full propagated trajectory with respect to initial state. """
sensitivity_trace(m::Module, args...; kwargs...) = sensitivity_trace(Val(first(fullname(m))), args...; kwargs...)
function sensitivity_trace(::Val{:ForwardDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; trace_time=false, kwargs...)
    values = trace_time ? [state.u0..., state.tspan[end] - state.tspan[begin]] : state.u0

    # Seed the values we want to trace with Dual numbers
    tag = typeof(state.model)
    duals = DiffEqSensitivity.seed_duals(values, tag)
    u0 = MVector{length(state.u0)}(duals[1:length(state.u0)])

    # Remake the state with the seeded values
    tspan = trace_time ? (state.tspan[begin], state.tspan[begin] + duals[end]) : state.tspan
    state_AD = remake(state; u0, tspan)

    # Solve and convert to the desired frame
    sol = solve(state_AD, alg; kwargs...)
    return convert_to_frame(sol, desired_frame)
end
@traitfn function sensitivity_trace(::Val{:VariationalEquations}, state::S, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...) where {S; HasVE{S}}
    dim = length(state.u0)
    I_flat = reshape(Matrix{Float64}(I, dim, dim), dim^2)
    u0_STM = MVector{dim^2 + dim}(vcat(state.u0, I_flat))
    prob_vareqns = State(state.model, state.frame, ODEProblem(state.model.ode.ode_stm_f, u0_STM, state.tspan, state.p))
    traj_vareqns_unconverted = solve(prob_vareqns, alg; kwargs...)

    # Convert to the desired frame
    return convert_to_frame(traj_vareqns_unconverted, desired_frame)
end

# TODO: Better functions for extracting STMs and stability index from Dual numbers
function extract_STMs(u)
    [hcat([Array(t.partials) for t in ut]...)' for ut in u]
end
function extract_stability(u)
    hcat(eigvals.(extract_STMs(u))...)
end