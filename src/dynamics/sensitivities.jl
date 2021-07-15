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
struct StateTransitionTensor{Order,Size<:Tuple,TensorsTuple<:Tuple} <: Abstract_StateTransitionTensor{Order}
    tensors::TensorsTuple

    function StateTransitionTensor(tensors::TensorsTuple) where {TensorsTuple<:Tuple}
        order = length(tensors)
        tensor_size = size(tensors[1])
        new{order, Tuple{tensor_size...}, TensorsTuple}(tensors)
    end
end

# Constructor

@doc """ Extracts the State Transition Tensor from the state (solved with solve_sensitivity). """
function StateTransitionTensor(state::State; order=get_order(eltype(state.u0)))
    real_order = get_order(eltype(state.u0))
    order > 0 || error("Expected order > 0, got $(order).")
    order <= real_order || error("Values in given STT only support order <= $(real_order), got $(order).")

    # TODO: Make this generic (assumes ForwardDiff partials)

    input_length = length(state.u0)
    output_length = ForwardDiff.npartials(eltype(state.u0))

    tensors = SArray[]
    for o in 1:order
        output_dims = ntuple(x -> output_length, o)

        coeff = 1 / factorial(o)
        tensor = SArray{Tuple{input_length,output_dims...}}(
            ForwardDiff.value.(
                ForwardDiff.partials.(
                    state.u0, 
                    ntuple(i -> reshape(1:output_dims[i], (ntuple(x -> 1, i)..., output_dims[i])), o)...
                )
            )
        )
        
        push!(tensors, coeff * tensor)
    end
    return StateTransitionTensor(tuple(tensors...))
end
const STT = StateTransitionTensor

# State Transition Matrix
const StateTransitionMatrix = StateTransitionTensor{1}
const STM = StateTransitionMatrix
StateTransitionMatrix(args...; kwargs...) = StateTransitionTensor(args...; kwargs..., order=1)

# Properties
# Base.axes(stm::StateTransitionTensor, args...) = axes(stm.tensor, args...)

# Base.getindex(stm::StateTransitionMatrix, idx1, idx2) = @view stm.tensor[idx1, idx2]

Base.show(io::IO, x::StateTransitionTensor{Order,Size}) where {Order,Size} = print(io, "STT{order=$(Order)}$(tuple(Size.parameters...))")
Base.show(io::IO, x::StateTransitionMatrix{Size}) where {Size} = print(io, "STM$(tuple(Size.parameters...))")

function (Base.:*)(stm::StateTransitionTensor{1}, dx)
    A = stm.tensors[1]
    @tullio(DX[i] := A[i,j] * dx[j])
    DX
end

function (Base.:*)(stm::StateTransitionTensor{2}, dx)
    A, B = stm.tensors
    @tullio(DX[i] := A[i,j] * dx[j])
    @tullio(DX[i] += B[i,j,k] * dx[j] * dx[k])
    DX
end

function (Base.:*)(stt1::StateTransitionTensor, stt2::StateTransitionTensor)
    StateTransitionTensor(stt1.tensor * stt2.tensor)
end

@doc """ Sensitivity of the interpolated states (at times t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory, t; kwargs...) = StateTransitionTensor(sol(t); kwargs...)

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory; kwargs...) = SVector{length(sol.t)}([StateTransitionTensor(sol[i]; kwargs...) for i in 1:length(sol.t)])

@doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
function StateTransitionMatrix(sol::Trajectory, t1, t2; kwargs...)
    @assert t1 <= t2  "Expected t1 <= t2"
    if t2 == t1
        return Matrix(1.0 * I, 6, 6)  # TODO: Make this type and size-generic, static
    else
        return StateTransitionMatrix(sol, t2; kwargs...)[:,:] / StateTransitionMatrix(sol, t1; kwargs...)[:,:]
    end
end

recursive_value(val) = val
recursive_value(val::ForwardDiff.Dual) = recursive_value(ForwardDiff.value(val))

get_order(valtype::Type{<:ForwardDiff.Dual}) = 1 + get_order(valtype.parameters[2])
get_order(::Type{<:Number}) = 0