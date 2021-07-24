export AD, FD, VE
export solve_sensitivity, supports_sensitivity_variational_equations
export StateTransitionMatrix, StateTransitionTensor, STM, STT

#---------------------------#
# SOLVING WITH SENSITIVITY  #
#---------------------------#

# TODO: Use Symbols instead of these constants?
const AD = ForwardDiff
const FD = FiniteDiff
const VE = Val(:VariationalEquations)

@traitdef HasVE{X}
@traitimpl HasVE{X} <- supports_sensitivity_variational_equations(X)
supports_sensitivity_variational_equations(X::Type{<:Abstract_ModelODEFunctions}) = hasfield(X, :ode_stm_f)
supports_sensitivity_variational_equations(X::Type{<:Abstract_DynamicalModel}) = supports_sensitivity_variational_equations(fieldtype(X, :ode))
supports_sensitivity_variational_equations(X::Type{<:State}) = supports_sensitivity_variational_equations(fieldtype(X, :model))

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
@traitfn function solve_sensitivity(::Val{:VariationalEquations}, state::S, desired_frame=state.frame, alg=DEFAULT_ALG; order=1, kwargs...) where {S; HasVE{S}}
    order == 1 || error("Variational equations only support order=1 for now.")

    dim = length(state.u0)
    I_flat = reshape(Matrix{Float64}(I, dim, dim), dim^2)
    u0_STM = MVector{dim^2 + dim}(vcat(state.u0, I_flat))
    prob_vareqns = State(state.model, state.frame, ODEProblem(state.model.ode.ode_stm_f, u0_STM, state.tspan, state.p))
    traj_vareqns_unconverted = solve(prob_vareqns, alg; kwargs...)

    # Convert to the desired frame
    return convert_to_frame(traj_vareqns_unconverted, desired_frame)
end
solve_sensitivity(::Val{:FiniteDiff}, args...; kwargs...) = error("solve_sensitivity(state) is not defined for FiniteDiff. Call STM(FD, state) instead.")

#--------------------------------#
# State Transition Tensors (STT) #
#--------------------------------#
# An STT represents the sensitivities of a set of output functions with respect
# to some input variables (around some reference point).
struct StateTransitionTensor{Order,TensorsTuple<:Tuple{Vararg{<:TensorMap,Order}}} <: Abstract_StateTransitionTensor{Order}
    t::Float64             # The time represented by the STT -- NOTE: we do not necessarily know the reference time.
    tensors::TensorsTuple  # tensors[i] holds the ith-order tensor.
end
const STT = StateTransitionTensor

# State Transition Matrix
const StateTransitionMatrix = StateTransitionTensor{1}
const STM = StateTransitionMatrix

# Constructors
@doc """ Extracts the State Transition Tensor from the state (solved with solve_sensitivity). """
function StateTransitionTensor(state::State, t; order=get_order(eltype(state.u0)))
    tensors = extract_sensitivity_tensors(state.u0, state_length(state), order)
    return StateTransitionTensor(recursive_value(t), tensors)
end
StateTransitionTensor(state::State; kwargs...) = StateTransitionTensor(state, state.tspan[begin]; kwargs...)
StateTransitionTensor(m::Module, args...; kwargs...) = StateTransitionTensor(Val(first(fullname(m))), args...; kwargs...)
StateTransitionTensor(v::Val, args...; kwargs...) = StateTransitionTensor(solve_sensitivity(v, args...; kwargs...))
StateTransitionMatrix(args...; kwargs...) = StateTransitionTensor(args...; kwargs..., order=1)

function extract_sensitivity_tensors(u0::AbstractArray, codomain_length, order)
    order == 0 || order == 1 || error("Given state only contains Variational Equation sensitivities of up to order 1, but got $(order).")

    error("Not implemented")

    # Variational Equations
    expected_length = codomain_length^2 + codomain_length
    length(u0) == expected_length || error("Unknown sensitivity type, expected $(expected_length) variables but found $(length(u0)).")

    tensor = SMatrix{codomain_length,codomain_length}(reshape(u0[codomain_length+1:end], codomain_length, codomain_length))
    return (tensor,)
end

function extract_sensitivity_tensors(u0::AbstractArray{DualType}, codomain_length, order) where {DualType<:ForwardDiff.Dual}
    real_order = get_order(DualType)

    order > 0 || error("Expected order > 0, got $(order).")
    order <= real_order || error("Given state only contains sensitivities of up to order $(real_order), but got $(order).")

    domain_length = numin(DualType)

    tensors = TensorMap[]
    for o in 1:order
        input_dims = ntuple(x -> domain_length, o)

        tensor = TensorMap(
            ForwardDiff.value.(
                ForwardDiff.partials.(
                    u0, 
                    ntuple(i -> reshape(1:input_dims[i], (ntuple(x -> 1, i)..., input_dims[i])), o)...
                )
            ), ⊗([ℝ^domain_length for _ in 1:o]...) → ℝ^codomain_length
        )
        
        push!(tensors, tensor)
    end
    return tuple(tensors...)
end

TensorKit.numin(dualType::Type{<:ForwardDiff.Dual}) = ForwardDiff.npartials(dualType)

tensor_order(::Type{<:StateTransitionTensor{Order}}) where {Order} = Order
tensor_order(stt::StateTransitionTensor) = tensor_order(typeof(stt))

# Define STT axes as (Out, In) dimensions.
# NOTE: We use only the 1st-order tensor, because we don't want to subset on higher-orders.
Base.axes(stt::StateTransitionTensor) = tuple(SOneTo.(dim.([codomain(stt.tensors[1])..., domain(stt.tensors[1])...]))...)
Base.axes(stt::StateTransitionTensor, i) = axes(stt)[i]

# Make STMs (1st-order only) broadcastable
Broadcast.broadcastable(stm::StateTransitionMatrix) = convert(Array, stm.tensors[1])

# Base.length(stm::StateTransitionMatrix) = length(stm.tensors[1])
# Base.iterate(stm::StateTransitionMatrix, args...) = iterate(stm.tensors[1], args...)

# Define subsets of STTs.
# We only subset on the (Out, In) dimensions (i.e. up to 1st-order)
Base.getindex(stm::StateTransitionTensor, domain_to_codomain::Pair) = stm[domain_to_codomain[2], domain_to_codomain[1]]
function Base.getindex(stm::StateTransitionTensor, codomain_idx, domain_idx)
    # Work out size of the newly-indexed codomain and domain. Need this to avoid:
    # 1) scalar indices, which will cause a dimension to be dropped
    # 2) Colon() indices, which have unknown size without context
    tensor_axes = axes(stm)  # Original tensor axes (to give context)
    codomain_size, domain_size = @. length(getindex(tensor_axes, (codomain_idx, domain_idx)))

    # Create the vector space based on the above sizes
    codomain_dim, domain_dim = (ℝ^codomain_size, ℝ^domain_size)

    new_tensors = []
    for tensor in stm.tensors
        tensor_array = convert(Array, tensor)[codomain_idx, domain_idx]
        tensor_map = TensorMap(
            collect(tensor_array), # XXX: Need the collect() to avoid scalars when the dims are (1, 1)
            ⊗([domain_dim for _ in 1:numin(tensor)]...) → codomain_dim
        )
        push!(new_tensors, tensor_map)
    end
    StateTransitionTensor(stm.t, tuple(new_tensors...))
end

# Display
Base.show(io::IO, x::StateTransitionTensor{Order}) where {Order} = print(io, string(SciMLBase.TYPE_COLOR, "STT", SciMLBase.NO_COLOR, "($(space(x.tensors[end])), t=$(x.t))"))
Base.show(io::IO, x::StateTransitionMatrix) = print(io, string(SciMLBase.TYPE_COLOR, "STM", SciMLBase.NO_COLOR, "($(space(x.tensors[end])), t=$(x.t))"))

# Comparison
# Base.isapprox(stt1::StateTransitionTensor{N}, stt2::StateTransitionTensor{N}; kwargs...) where {N} = all(isapprox.(stt1.tensors, stt2.tensors; kwargs...))
# Base.isapprox(stt1::StateTransitionTensor{N}, stt2::StateTransitionTensor{M}; kwargs...) where {N,M} = false

# Generate multiplication functions for STTs up to MAX_STT_ORDER
#----------------------------------------------------------------
# Each multiplication function uses Tullio.@tullio to perform a tensor contraction using einstein summation notation,
# contracting each tensor in STT.tensors with the appropriate number of vector values.
@generated function (Base.:*)(STT::StateTransitionTensor{Order}, DX::Tensor) where {Order}
    # The ith-order STT should be contracted with i variables.
    # For example, an STM (order=1) is DX[i] = STM[i,j] * dx[j]
    # But an STT(order=2) is DX[i] = STM[i,j,k] * dx[j] * dx[k]
    # Here, we generate an automatic list of indices (e.g. (j, k)), and from
    # that a list of [dx[j], dx[k]...].
    tensor_indices = [gensym() for _ in 1:Order]
    dxs = [:(DX[$(idx)]) for idx in tensor_indices]

    # TODO: Update this to use TensorKit.permute instead of @tensor, which should speed things up.
    # Example:
    # function contract(tensor_order1, tensor_order2, dx)
    #     DX = tensor_order1 * dx
    #     DX += 1/factorial(2) * permute(tensor_order2 * dx, (1,), (2,)) * dx
    #     return DX
    # end
    # (Or alternatively, go back to using @tullio).

    exprs = []
    for order in 1:Order
        push!(exprs, quote
            tensor = STT.tensors[$(order)]
            coeff = $(1 / factorial(order))
            @tensor RES[i] += *(coeff, tensor[i,$(tensor_indices[1:order]...)], $(dxs[1:order]...))
        end)
    end

    return quote
        RES = Tensor(zeros, domain(STT.tensors[1]))
        $(exprs...)
        return convert(Array, RES)
    end
end

# (Base.:*)(coeff::Number, stt::StateTransitionTensor) = StateTransitionTensor(stt.)
(Base.:*)(stt::StateTransitionTensor, dx) = stt * Tensor(collect(dx), domain(stt.tensors[1]))

# Inverses
(Base.inv)(stm::StateTransitionMatrix) = StateTransitionTensor(stm.t, (inv(stm.tensors[1]),))

# Multiplication of STMs.
# TODO: Need a much better explanation of this.
# NOTE: Assumes that the timespan represented by each STM is non-overlapping but
# continuous over a full time interval of [t_0, t_f], such that:
# -- STM1 represents a time interval [t_k, t_f]
# -- STM2 represents a time interval [t_0, t_k]
# -- (STM1 * STM2) results in the STM for the time interval [t_0, t_f]
# function (Base.:*)(stm1::StateTransitionMatrix, stm2::StateTransitionMatrix)
#     tensor1 = stm1.tensors[1]
#     tensor2 = stm2.tensors[1]
#     new_tensor = similar(tensor1)
#     @tullio new_tensor[i,a] = tensor1[i,alpha] * tensor2[alpha,a]
#     StateTransitionTensor((new_tensor,))
# end

# function (Base.:*)(stt1::S, stt2::S) where {S<:StateTransitionTensor{2}}
#     # 1st order
#     tensor1 = stt1.tensors[1]
#     tensor2 = stt2.tensors[1]
#     new_tensor1 = similar(tensor1)
#     @tullio new_tensor1[i,a] = tensor1[i,alpha] * tensor2[alpha,a]

#     # 2nd order
#     tensor1_2 = stt1.tensors[2]
#     tensor2_2 = stt2.tensors[2]
#     new_tensor2 = similar(tensor1_2)
#     @tullio new_tensor2[i,a,b] = tensor1[i,alpha] * tensor2_2[alpha,a,b] + tensor1_2[i,alpha,beta] * tensor2[alpha,a] * tensor2[beta,b]

#     StateTransitionTensor((new_tensor1, new_tensor2))
# end

# STTs of orders not generated above are not supported automatically. Users
# should define their contractions manually. (I wanted to do this completely
# automatically but haven't been able to yet.)
# (Base.:*)(::StateTransitionTensor{Order}, _) where {Order} = error("Multiplication for STTs is only defined up to order $(MAX_STT_ORDER), got order $(Order). Please define the multiplication manually.")

@doc """ Sensitivity of the interpolated states (at times t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory, t; kwargs...) = StateTransitionTensor(sol(t); kwargs...)
StateTransitionTensor(sol::Trajectory, t::Number; kwargs...) = StateTransitionTensor(sol(t), t; kwargs...)

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory; kwargs...) = [StateTransitionTensor(sol[i], sol.t[i]; kwargs...) for i in 1:length(sol.t)]

# @doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
# function StateTransitionMatrix(sol::Trajectory, t1, t2; kwargs...)
#     @assert t1 <= t2  "Expected t1 <= t2"
#     if t1 == t2
#         sol0 = sol[begin]
#         input_length = state_length(sol0)
#         output_length = ForwardDiff.npartials(eltype(sol0.u0))  # TODO: Make generic
#         return StateTransitionTensor((SMatrix{input_length, output_length}(1.0 * I),))
#     elseif t1 == sol.t[begin]
#         return StateTransitionMatrix(sol, t2; kwargs...)
#     else
#         return StateTransitionMatrix(sol, t2; kwargs...) * inv(StateTransitionMatrix(sol, t1; kwargs...))
#     end
# end

# Compute STMs by solving the given state first
# function StateTransitionTensor(::Val{:FiniteDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
#     # NOTE: Only supports order=1
#     tensors = (FiniteDiff.finite_difference_jacobian(state.u0) do u0
#         new_state = remake(state, u0=u0)
#         trajectory = solve(new_state, alg; kwargs...)
#         trajectory_converted = convert_to_frame(trajectory, desired_frame)
#         end_state_u = trajectory_converted.sol[end]
#     end,)
#     StateTransitionTensor(tensors)
# end

recursive_value(val) = val
recursive_value(val::ForwardDiff.Dual) = recursive_value(ForwardDiff.value(val))

get_order(valtype::Type{<:ForwardDiff.Dual}) = 1 + get_order(valtype.parameters[2])
get_order(::Type{<:Number}) = 0

# LinearAlgebra.eigvals(stm::StateTransitionMatrix) = eigvals(Array(stm.tensors[1]))