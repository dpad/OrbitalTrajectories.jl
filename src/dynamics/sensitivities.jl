export AD, FD, VE, solve_sensitivity, stability_index, has_variational_equations
export StateTransitionMatrix, StateTransitionTensor, STM, STT

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
function solve_sensitivity(::Val{:ForwardDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; order=1, trace_time=false, kwargs...)
    trace_time && !isnothing(get(kwargs, :callback, nothing)) && @warn("trace_time sensitivity does not seem to play well with callbacks!")

    values = trace_time ? [state.u0..., state.tspan[end]] : state.u0

    # Seed the values we want to trace with Dual numbers
    tag = typeof(state.model)
    duals = copy(values)
    for _ in 1:order
        duals = DiffEqSensitivity.seed_duals(duals, tag)
    end
    u0 = MVector{length(state.u0)}(duals[1:length(state.u0)])

    # Remake the state with the seeded values
    tspan = trace_time ? (state.tspan[begin], duals[end]) : state.tspan
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

function stability_index(sol::Union{Trajectory,State})
    hcat(eigvals.(extract_STMs(sol))...)
end


# State Transition Tensors
struct StateTransitionTensor{N,V<:AbstractArray,T}
    tspan::Tuple{T,T}
    vals::V

    # Constructor
    function StateTransitionTensor(tspan, vals::AbstractArray, order=get_order(eltype(vals)))
        real_order = get_order(eltype(vals))
        order > 0 || throw("Expected order > 0, got $(order).")
        order <= real_order || throw("Values in given STT only support order <= $(real_order), got $(order).")
        return new{order, typeof(vals), eltype(tspan)}(tuple(tspan...), vals)
    end
end
const STT = StateTransitionTensor

# State Transition Matrix
const StateTransitionMatrix{V} = StateTransitionTensor{1, V}
const STM = StateTransitionMatrix
StateTransitionMatrix(other::StateTransitionTensor) = StateTransitionTensor(other, 1)
StateTransitionMatrix(args...; kwargs...) = StateTransitionTensor(args...; kwargs...)

# Properties
Base.axes(stm::StateTransitionMatrix) = (Base.OneTo(length(stm.vals)), Base.OneTo(length(stm.vals[1].partials)))
Base.axes(stm::StateTransitionMatrix, d) = axes(stm)[d]

function Base.getindex(stm::StateTransitionMatrix, idx1, idx2)
    # TODO: Make this generic, currently expects ForwardDiff.Duals
    @view ForwardDiff.value.(ForwardDiff.partials.(stm.vals, transpose(1:length(stm.vals[1].partials))))[idx1, idx2]
end

Base.show(io::IO, x::StateTransitionTensor{N,V}) where {N,V} = print(io, "STT{$(nameof(eltype(V))), order=$(N)}$(recursive_value.(x.vals))")
Base.show(io::IO, x::StateTransitionMatrix{V}) where {V} = print(io, "STM{$(nameof(eltype(V)))}$(recursive_value.(x.vals))")

function (Base.:+)(stm::StateTransitionMatrix, dx::AbstractArray{<:Number})
    @error "Addition not yet supported with new STT/STM types!"
#     T = eltype(state.u0).parameters[1]  # Get the tag type
    # STM = ForwardDiff.extract_jacobian(T, stm.vals, values.(dx))
    @info "mult!"
end

# Copy constructor
StateTransitionTensor(other::STT, new_order=N) where {N,STT<:StateTransitionTensor{N}} = StateTransitionTensor(other.tspan, other.vals, new_order)

@doc """ Extracts the State Transition Tensor from the state (solved with solve_sensitivity). """
StateTransitionTensor(state::State) = StateTransitionTensor(state.prob.tspan, state.u0)

@doc """ Solve and return the sensitivity of the final state (at t=state.tspan[2]) with respect to the given state. """
StateTransitionTensor(m::Module, args...; kwargs...) = StateTransitionTensor(Val(first(fullname(m))), args...; kwargs...)
function StateTransitionTensor(m::Val, state::State, args...; kwargs...)
    sol = solve_sensitivity(m, state, args...; kwargs...)
    StateTransitionTensor(sol, sol.t[end])
end
function StateTransitionMatrix(::Val{:FiniteDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
    FiniteDiff.finite_difference_jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        return convert_to_frame(solve(new_state, alg; kwargs...), desired_frame).sol[end]
    end
end

@doc """ Sensitivity of the interpolated states (at times t) with respect to initial state. """
function StateTransitionTensor(sol::Trajectory, t)
    StateTransitionTensor(sol(t))
end

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
function StateTransitionTensor(sol::Trajectory)
    SVector{length(sol.t)}([StateTransitionTensor(sol[i]) for i in 1:length(sol.t)])
end

@doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
function StateTransitionMatrix(sol::Trajectory, t1, t2)
    @assert t1 <= t2  "Expected t1 <= t2"
    if t2 == t1
        return Matrix(1.0 * I, 6, 6)  # TODO: Make this type and size-generic, static
    else
        return StateTransitionMatrix(sol, t2)[:,:] / StateTransitionMatrix(sol, t1)[:,:]
    end
end

recursive_value(val) = val
recursive_value(val::ForwardDiff.Dual) = recursive_value(ForwardDiff.value(val))

get_order(valtype::Type{<:ForwardDiff.Dual}) = 1 + get_order(valtype.parameters[2])
get_order(::Type{<:Number}) = 0