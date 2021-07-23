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
#
# An STT represents the sensitivities of a set of output functions with respect
# to some input variables (around some reference point).
# The Size tuple is (Out, In), where Out is the number of output functions, and
# In is the number of input variables.
#
# For example, a dynamical model with variables [x, y, z, ẋ, ẏ, ż] has Out=6,
# and usually these are defined with respect to u0 = [x(0), y(0), z(0), ẋ(0),
# ẏ(0), ż(0)], so In=6 (or it may also include t, in which case In=7).
struct StateTransitionTensor{Order,Size<:Tuple,TensorsTuple<:Tuple} <: Abstract_StateTransitionTensor{Order}
    # tensors[i] holds the ith-order tensor, already multiplied by its 1/factorial(i) coefficient.
    #
    # For example, the 1st-order tensor is simply the State Transition Matrix of
    # size (Out, In). Note that ndims=2 and hence it is called a "Matrix".
    #
    # The 2nd-order tensor is 1/(2!) * (2nd-order derivatives arrays of size
    # (Out, In, In)). Note that ndims=3.
    #
    # Similarly, each ith-order tensor has a size of (Out, In...) where ndims=i+1.
    tensors::TensorsTuple

    function StateTransitionTensor(tensors::TensorsTuple) where {TensorsTuple<:Tuple}
        order = length(tensors)
        tensor_size = size(tensors[1])
        new{order, Tuple{tensor_size...}, TensorsTuple}(tensors)
    end
end
const STT = StateTransitionTensor

# State Transition Matrix
const StateTransitionMatrix = StateTransitionTensor{1}
const STM = StateTransitionMatrix

# Constructors
StateTransitionTensor(m::Module, args...; kwargs...) = StateTransitionTensor(Val(first(fullname(m))), args...; kwargs...)
StateTransitionTensor(v::Val, args...; kwargs...) = StateTransitionTensor(solve_sensitivity(v, args...; kwargs...))

StateTransitionMatrix(args...; kwargs...) = StateTransitionTensor(args...; kwargs..., order=1)

@doc """ Extracts the State Transition Tensor from the state (solved with solve_sensitivity). """
function StateTransitionTensor(state::S; order=get_order(eltype(state.u0))) where {S<:State}
    tensors = extract_sensitivity_tensors(state.u0, state_length(state), order)
    return StateTransitionTensor(tensors)
end

function extract_sensitivity_tensors(u0::AbstractArray, input_length, order)
    order == 0 || order == 1 || error("Given state only contains Variational Equation sensitivities of up to order 1, but got $(order).")

    # Variational Equations
    expected_length = input_length^2 + input_length
    length(u0) == expected_length || error("Unknown sensitivity type, expected $(expected_length) variables but found $(length(u0)).")

    tensor = SMatrix{input_length,input_length}(reshape(u0[input_length+1:end], input_length, input_length))
    return (tensor,)
end

function extract_sensitivity_tensors(u0::AbstractArray{<:ForwardDiff.Dual}, input_length, order)
    real_order = get_order(eltype(u0))
    order > 0 || error("Expected order > 0, got $(order).")
    order <= real_order || error("Given state only contains sensitivities of up to order $(real_order), but got $(order).")

    output_length = ForwardDiff.npartials(eltype(u0))

    tensors = SArray[]
    for o in 1:order
        output_dims = ntuple(x -> output_length, o)

        tensor = SArray{Tuple{input_length,output_dims...}}(
            ForwardDiff.value.(
                ForwardDiff.partials.(
                    u0, 
                    ntuple(i -> reshape(1:output_dims[i], (ntuple(x -> 1, i)..., output_dims[i])), o)...
                )
            )
        )
        
        coeff = 1 / factorial(order)
        push!(tensors, coeff .* tensor)
    end
    return tuple(tensors...)
end

tensor_order(::Type{<:StateTransitionTensor{Order}}) where {Order} = Order
tensor_order(stt::StateTransitionTensor) = tensor_order(typeof(stt))
tensor_size(::Type{<:StateTransitionTensor{Order,Size}}) where {Order,Size} = tuple(Size.parameters...)
tensor_size(stt::StateTransitionTensor{Order,Size}) where {Order,Size} = tensor_size(typeof(stt))

# Define subsets of STTs
Base.axes(::StateTransitionTensor{O,Size}) where {O,Size} = tuple(SOneTo.(Size.parameters)...)
Base.axes(::StateTransitionTensor{O,Size}, i) where {O,Size} = SOneTo.(Size.parameters[i])

# Define subsets of STMs
Broadcast.broadcastable(stm::StateTransitionMatrix) = stm.tensors[1]
Base.length(stm::StateTransitionMatrix) = length(stm.tensors[1])
Base.iterate(stm::StateTransitionMatrix, args...) = iterate(stm.tensors[1], args...)
function Base.getindex(stm::StateTransitionMatrix, idx1, idx2)
    # TODO: Ensure this is a matrix or array of size = length(idx1), length(idx2), ...
    tensor_views = [view(tensor, idx1, idx2) for tensor in stm.tensors]
    StateTransitionTensor(tuple(tensor_views...))
end

# Display
Base.show(io::IO, x::StateTransitionTensor{Order,Size}) where {Order,Size} = print(io, "STT{order=$(Order)}$(tuple(Size.parameters...))")
Base.show(io::IO, x::StateTransitionMatrix{Size}) where {Size} = print(io, "STM$(tuple(Size.parameters...))")

# Comparison
Base.isapprox(stt1::StateTransitionTensor{N}, stt2::StateTransitionTensor{N}; kwargs...) where {N} = all(isapprox.(stt1.tensors, stt2.tensors; kwargs...))
Base.isapprox(stt1::StateTransitionTensor{N}, stt2::StateTransitionTensor{M}; kwargs...) where {N,M} = false

# Generate multiplication functions for STTs up to MAX_STT_ORDER
#----------------------------------------------------------------
# Each multiplication function uses Tullio.@tullio to perform a tensor contraction using einstein summation notation,
# contracting each tensor in STT.tensors with the appropriate number of vector values.
# XXX: I tried doing this with @generated functions, but it seems that the @tullio macro creates closures, meaning
# it can't be used in the output of an @generated function.
const MAX_STT_ORDER = 4
let exprs = []
    for STT_ORDER in 1:MAX_STT_ORDER
        # The ith-order STT should be contracted with i variables.
        # For example, an STM (order=1) is DX[i] = STM[i,j] * dx[j]
        # But an STT(order=2) is DX[i] = STM[i,j,k] * dx[j] * dx[k]
        #
        # Here, we generate an automatic list of indices (e.g. (j, k)), and from
        # that a list of [dx[j], dx[k]...].
        tensor_indices = [gensym() for _ in 2:(STT_ORDER+1)]
        dxs = [:(dx[$(tensor_indices[i])]) for i in 1:STT_ORDER]

        push!(exprs, quote
            tensor = stt.tensors[$(STT_ORDER)]
            @tullio DX[i] += tensor[i,$(tensor_indices...)] * *($(dxs...))
        end)
        eval(quote
            # Create a function specialised to the specific STT order.
            function (Base.:*)(stt::StateTransitionTensor{$(STT_ORDER)}, dx::AbstractVector)
                DX = zeros(eltype(stt.tensors[1]), tensor_size(stt)[1])
                $(exprs...)
                DX
            end
        end)
    end
end

function (Base.inv)(stm::StateTransitionMatrix)
    StateTransitionTensor((inv(stm.tensors[1]),))
end

function (Base.:*)(stm1::StateTransitionMatrix, stm2::StateTransitionMatrix)
    tensor1 = stm1.tensors[1]
    tensor2 = stm2.tensors[1]
    new_tensor = similar(tensor1)
    @tullio new_tensor[i,a] = tensor1[i,alpha] * tensor2[alpha,a]
    StateTransitionTensor((new_tensor,))
end

# STTs of orders not generated above are not supported automatically. Users
# should define their contractions manually. (I wanted to do this completely
# automatically but haven't been able to yet.)
(Base.:*)(::StateTransitionTensor{Order}, _) where {Order} = error("Multiplication for STTs is only defined up to order $(MAX_STT_ORDER), got order $(Order). Please define the multiplication manually.")

@doc """ Sensitivity of the interpolated states (at times t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory, t; kwargs...) = StateTransitionTensor(sol(t); kwargs...)

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory; kwargs...) = [StateTransitionTensor(sol[i]; kwargs...) for i in 1:length(sol.t)]

@doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
function StateTransitionMatrix(sol::Trajectory, t1, t2; kwargs...)
    @assert t1 <= t2  "Expected t1 <= t2"
    if t1 == t2
        sol0 = sol[begin]
        input_length = state_length(sol0)
        output_length = ForwardDiff.npartials(eltype(sol0.u0))  # TODO: Make generic
        return StateTransitionTensor((SMatrix{input_length, output_length}(1.0 * I),))
    elseif t1 == sol.t[begin]
        return StateTransitionMatrix(sol, t2; kwargs...)
    else
        return StateTransitionMatrix(sol, t2; kwargs...) * inv(StateTransitionMatrix(sol, t1; kwargs...))
    end
end

# Compute STMs by solving the given state first
function StateTransitionTensor(::Val{:FiniteDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; kwargs...)
    # NOTE: Only supports order=1
    tensors = (FiniteDiff.finite_difference_jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        trajectory = solve(new_state, alg; kwargs...)
        trajectory_converted = convert_to_frame(trajectory, desired_frame)
        end_state_u = trajectory_converted.sol[end]
    end,)
    StateTransitionTensor(tensors)
end

recursive_value(val) = val
recursive_value(val::ForwardDiff.Dual) = recursive_value(ForwardDiff.value(val))

get_order(valtype::Type{<:ForwardDiff.Dual}) = 1 + get_order(valtype.parameters[2])
get_order(::Type{<:Number}) = 0

LinearAlgebra.eigvals(stm::StateTransitionMatrix) = eigvals(Array(stm.tensors[1]))