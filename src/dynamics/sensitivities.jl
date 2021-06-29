export AD, FD, VE, solve_sensitivity, get_sensitivity, stability_index, has_variational_equations

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

@doc """ Return the fully propagated trajectory including state sensitivities with respect to the initial state. """
solve_sensitivity(m::Module, args...; kwargs...) = solve_sensitivity(Val(first(fullname(m))), args...; kwargs...)
function solve_sensitivity(::Val{:ForwardDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; trace_time=false, kwargs...)
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
@traitfn function solve_sensitivity(::Val{:VariationalEquations}, state::S, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...) where {S; HasVE{S}}
    dim = length(state.u0)
    I_flat = reshape(Matrix{Float64}(I, dim, dim), dim^2)
    u0_STM = MVector{dim^2 + dim}(vcat(state.u0, I_flat))
    prob_vareqns = State(state.model, state.frame, ODEProblem(state.model.ode.ode_stm_f, u0_STM, state.tspan, state.p))
    traj_vareqns_unconverted = solve(prob_vareqns, alg; kwargs...)

    # Convert to the desired frame
    return convert_to_frame(traj_vareqns_unconverted, desired_frame)
end

@doc """ Solve and return the sensitivity of the final state (at t=state.tspan[2]) with respect to the given state. """
get_sensitivity(m::Module, args...; kwargs...) = get_sensitivity(Val(first(fullname(m))), args...; kwargs...)
function get_sensitivity(m::Val, state::State, args...; kwargs...)
    sol = solve_sensitivity(m, state, args...; kwargs...)
    get_sensitivity(sol, sol.t[end])
end
function get_sensitivity(::Val{:FiniteDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
    FiniteDiff.finite_difference_jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        return convert_to_frame(solve(new_state, alg; kwargs...), desired_frame).sol[end]
    end
end

@doc """ Sensitivity of the state at time t with respect to initial state. """
function get_sensitivity(sol::Trajectory, t)
    extract_STMs(sol, t)
end

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
function get_sensitivity(sol::Trajectory)
    extract_STMs(sol)
end

@doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
function get_sensitivity(sol::Trajectory, t1, t2)
    @assert t1 <= t2  "Expected t1 <= t2"
    if t2 == t1
        return Matrix(1.0 * I, 6, 6)  # TODO: Make this type and size-generic, static
    else
        return get_sensitivity(sol, t2) / get_sensitivity(sol, t1)
    end
end

# TODO: Better functions for extracting STMs and stability index from Dual numbers
function extract_STMs(sol::Trajectory, t)
    extract_STMs(sol(t))
end
function extract_STMs(sol::Trajectory)
    SVector{length(sol.t)}([extract_STMs(sol[i]) for i in 1:length(sol.t)])
end
function extract_STMs(state::State)
    T = eltype(state.u0).parameters[1]  # Get the tag type
    ForwardDiff.extract_jacobian(T, state.u0, values.(state.u0))
end
function stability_index(sol::Union{Trajectory,State})
    hcat(eigvals.(extract_STMs(sol))...)
end