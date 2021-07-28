export AD, FD, VE
export solve_sensitivity
export StateTransitionMatrix, StateTransitionTensor, STM, STT

#---------------------------#
# SOLVING WITH SENSITIVITY  #
#---------------------------#

# ForwardDiff Dual utilities
recursive_value(val) = val
recursive_value(val::ForwardDiff.Dual) = recursive_value(ForwardDiff.value(val))
get_order(valtype::Type{<:ForwardDiff.Dual}) = 1 + get_order(valtype.parameters[2])
get_order(::Type{<:Number}) = 0

const AD = ForwardDiff
const FD = FiniteDiff
const VE = Val(:VariationalEquations)

@traitdef ModelSupports{X}
@traitimpl ModelSupports{X} <- model_supports_sensitivity_variational_equations(X)
# model_supports_sensitivity_variational_equations(X::Type{<:Abstract_AstrodynamicalODESystem}) = hasfield(X, :ode_stm_f)
model_supports_sensitivity_variational_equations(X::Type{<:Abstract_AstrodynamicalModel}) = model_supports_sensitivity_variational_equations(fieldtype(X, :ode))
model_supports_sensitivity_variational_equations(X::Type{<:State}) = model_supports_sensitivity_variational_equations(fieldtype(X, :model))

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
# @traitfn function solve_sensitivity(::Val{:VariationalEquations}, state::S, desired_frame=state.frame, alg=DEFAULT_ALG; order=1, kwargs...) where {S; ModelSupports{S}}
#     order == 1 || error("Variational equations only support order=1 for now.")

#     dim = length(state.u0)
#     I_flat = reshape(Matrix{Float64}(I, dim, dim), dim^2)
#     u0_STM = MVector{dim^2 + dim}(vcat(state.u0, I_flat))
#     prob_vareqns = State(state.model, state.frame, ODEProblem(state.model.ode.ode_stm_f, u0_STM, state.tspan, state.p))
#     traj_vareqns_unconverted = solve(prob_vareqns, alg; kwargs...)

#     # Convert to the desired frame
#     return convert_to_frame(traj_vareqns_unconverted, desired_frame)
# end
solve_sensitivity(::Val{:FiniteDiff}, args...; kwargs...) = error("solve_sensitivity(state) is not defined for FiniteDiff. Call STM(FD, state) instead.")

#--------------------------------#
# State Transition Tensors (STT) #
#--------------------------------#
# An STT represents the sensitivities of a set of {Out} output functions with
# respect to some {In} input variables (around some reference point).
# 
# You should be able to multiply an STT with a vector of length {In} to get a
# vector of length {Out}, as follows:
# 
# > dx::SVector{7}  # Some small error in the initial 7-dimensional inputs.
# > (some_STT::STT{N,6,7} * dx) => ::SVector{6}
# > # Result is an N-th order estimate of 6-dimensional output state.
#
# In practice, it is a Taylor series of tensors representing linear maps.
# Each ith-order tensor maps a domain=(ℝ^In) → codomain=(ℝ^Out⊗ℝ^Out...).
# For example, the "1st-order tensor" (the State Transition Matrix) is simply a
# linear map from a vector to a vector -- hence, it is technically a 2nd-order
# tensor (a matrix), but in the context of the Taylor series, it contains the
# 1st-order derivatives.
#
# NOTE: We should try to use a tensor library like TensorKit.jl. This was
# implemented previously (see commits 10d2b21 through to 94a05f9), however there
# was some performance issues and annoyances, so I reverted back to SArrays for now.
struct StateTransitionTensor{Order,Out,In,TensorsTuple<:Tuple{Vararg{SArray{<:Tuple{Out,In,Vararg}},Order}}} <: Abstract_StateTransitionTensor{Order}
    tspan::NTuple{2,Float64}  # The timespan represented by the STT.
    tensors::TensorsTuple     # tensors[i] holds the ith-order tensor.
end
const STT = StateTransitionTensor

# State Transition Matrix
const StateTransitionMatrix = StateTransitionTensor{1}
const STM = StateTransitionMatrix

# Constructors
@doc """ Extracts the State Transition Tensor from the state (solved with solve_sensitivity). """
function StateTransitionTensor(state::State, tspan=state.tspan; order=get_order(eltype(state.u0)))
    tensors = extract_sensitivity_tensors(state.u0, state_length(state), order)
    return StateTransitionTensor(recursive_value.(tspan), tensors)
end
StateTransitionTensor(m::Module, args...; kwargs...) = StateTransitionTensor(Val(first(fullname(m))), args...; kwargs...)
StateTransitionTensor(v::Val, args...; kwargs...) = StateTransitionTensor(solve_sensitivity(v, args...; kwargs...))
StateTransitionMatrix(args...; kwargs...) = StateTransitionTensor(args...; kwargs..., order=1)

function extract_sensitivity_tensors(u0::AbstractArray, codomain_length, order)
    order == 0 || order == 1 || error("Given state only contains Variational Equation sensitivities of up to order 1, but got $(order).")

    # Variational Equations
    expected_length = codomain_length^2 + codomain_length
    length(u0) == expected_length || error("Unknown sensitivity type, expected $(expected_length) variables but found $(length(u0)).")

    tensor = SMatrix{codomain_length,codomain_length}(reshape(u0[codomain_length+1:end], codomain_length, codomain_length))
    return (tensor,)
end

function extract_sensitivity_tensors(u0::AbstractArray{DualType}, codomain_dim, order) where {DualType<:ForwardDiff.Dual}
    real_order = get_order(DualType)

    order > 0 || error("Expected order > 0, got $(order).")
    order <= real_order || error("Given state only contains sensitivities of up to order $(real_order), but got $(order).")

    domain_dim = ForwardDiff.npartials(DualType)

    # Extract up to order tensors from the Duals
    tensors = SArray[]
    for N in 1:order
        # Input dimensions of the tensor.
        # For example, if the domain dim (i.e. input size) is 6, and order is 3, then domain=(6,6,6).
        domain = ntuple(_->domain_dim, N)

        # Extract the Partial values directly from the multi-level partials.
        tensor = SArray{Tuple{codomain_dim,domain...}}(
            ForwardDiff.value.(
                ForwardDiff.partials.(
                    u0,
                    # To extract the N-th order partials as an N-dimensional
                    # array, we need to pass in N tuples that look like:
                    # (6,), (1,6), (1,1,6), and so on...
                    ntuple(i -> reshape(
                        1:domain[i],
                        (ntuple(_->1, i)..., domain[i])), 
                    N)...
                )
            )
        )
        
        push!(tensors, tensor)
    end
    return tuple(tensors...)
end

# Build identity STMs (not STTs)
Base.one(::S) where {Out,In,S<:StateTransitionMatrix{Out,In}} = convert(S, I)

# Define subsets of STTs, with axes as (Out, In) dimensions.
# NOTE: We assume all the higher-order tensors just repeat the In dimension, and
# we only subset on thee 1st-order.
Base.size(stt::StateTransitionTensor{Order,Out,In}) where {Order,Out,In} = (Out, In)
Base.axes(stt::StateTransitionTensor) = SOneTo.(size(stt))
Base.axes(stt::StateTransitionTensor, i) = axes(stt)[i]
Base.getindex(stt::StateTransitionTensor, domain_to_codomain::Pair) = stt[domain_to_codomain.second, domain_to_codomain.first]
function Base.getindex(stt::StateTransitionTensor{Order}, codomain, domain) where {Order}
    views = [begin
            # Preserve all original singleton dimensions by getting the subset size...
            codomain_dim, domain_dim = length.(to_indices(tensor, (codomain, domain)))
            new_dims = (codomain_dim, ntuple(_->domain_dim, Order)...)
            SArray{Tuple{new_dims...}}(view(tensor, codomain, domain))
        end for tensor in stt.tensors
    ]
    StateTransitionTensor(stt.tspan, tuple(views...))
end

# Make STMs (1st-order only) broadcastable
Broadcast.broadcastable(stm::StateTransitionMatrix) = stm.tensors[1]

# Display STTs/STMs.
Base.show(io::IO, ::MIME"text/plain", x::StateTransitionTensor{Order,Out,In}) where {Order,Out,In} = print(io, string(SciMLBase.TYPE_COLOR, "STT{$(Order)}", SciMLBase.NO_COLOR, "($(In)=>$(Out), t=$(x.tspan))"))
Base.show(io::IO, ::MIME"text/plain", x::StateTransitionMatrix{Out,In}) where {Out,In} = print(io, string(SciMLBase.TYPE_COLOR, "STM", SciMLBase.NO_COLOR, "($(In)=>$(Out), t=$(x.tspan))"))

# Comparison
Base.isapprox(stt1::S, stt2::S; kwargs...) where {N,Out,In,S<:StateTransitionTensor{N,Out,In}} = all(isapprox.(stt1.tensors, stt2.tensors; kwargs...))
Base.isapprox(stt1::S, stt2::S; kwargs...) where {N,S<:StateTransitionTensor{N}} = false
Base.isapprox(stt1::StateTransitionTensor{N}, stt2::StateTransitionTensor{M}; kwargs...) where {N,M} = false

# STT Multiplication & Inverse
#----------------------------------------------------------------
# Each of these generated functions uses Tullio.@tullio to perform tensor
# contraction on the STT.tensors using einstein summation notation.
#
# XXX: @generated functions for these would be ideal, but the @tullio macro creates
# closures so we can't do it in @generated. Can use TensorOperations.@tensor
# in an @generated function but it's up to ~10x slower, so we just generate up
# to MAX_STT_ORDER at package pre-compile time here.
const MAX_STT_ORDER = 3
const exprs = Expr[]
const inv_exprs = Expr[]
for ORDER in 1:MAX_STT_ORDER
    # Here, we generate an automatic list of indices (e.g. (j, k)), and from
    # that a list of [dx[j], dx[k]...].
    tensor_indices = [gensym() for _ in 1:ORDER]
    dxs = [:(dx[$(idx)]) for idx in tensor_indices]

    # [Refer to Park 2007, "Nonlinear trajectory navigation", eq.2.22]
    # The ith-order STT tensor should be contracted with i variables.
    # For example order=2 gives DX[i] = tensor[i,j,k] * dx[j] * dx[k].
    push!(exprs, quote
        tensor = stt.tensors[$(ORDER)]
        coeff = 1 / factorial($(ORDER))
        @tullio DX[i] += *(coeff, tensor[i,$(tensor_indices[1:ORDER]...)], $(dxs[1:ORDER]...))
    end)

    if ORDER == 2
        push!(inv_exprs, quote
            @tullio INV_TENSOR[i,a,b] := -inverse_tensors[1][i,alpha] * tensors[2][alpha,j1,j2] * inverse_tensors[1][j1,a] * inverse_tensors[1][j2,b]
            push!(inverse_tensors, INV_TENSOR)
        end)
    elseif ORDER == 3
        push!(inv_exprs, quote
            @tullio INV_TENSOR[i,a,b,c] := begin
                -(inverse_tensors[1][i,alpha] * tensors[2][alpha,j1,j2,j3] + inverse_tensors[2][i,alpha,beta] * (tensors[1][alpha,j1] * tensors[2][beta,j2,j3] + tensors[2][alpha,j1,j2] * tensors[1][beta,j3] + tensors[2][alpha,j1,j3] * tensors[1][beta,j2])) * (inverse_tensors[1][j1,a] * inverse_tensors[1][j2,b] * inverse_tensors[1][j3,c])
            end
            push!(inverse_tensors, INV_TENSOR)
        end)
    end

    # Create a function specialised to this specific STT order.
    eval(quote
        function (Base.:*)(stt::StateTransitionTensor{$(ORDER)}, dx::AbstractVector)
            DX = zeros(eltype(stt.tensors[1]), size(stt)[1])
            $(exprs[1:ORDER]...)
            DX
        end

        # # [Boone.2021b Appendix B | Park.2007 | Park.2007b]
        # exprs = []

        # if Order == 2
        #     push!(exprs, quote
        #             # The second order tensor...
        #             @tensor TENSOR2[i,a,b] := stt1.tensors[1][i,alpha] * stt2.tensors[2][alpha,a,b] + stt1.tensors[2][i,alpha,beta] * stt2.tensors[1][alpha,a] * stt2.tensors[1][beta,b]
        #             push!(multiplied_tensors, permute(TENSOR2, (1,), (2, 3)))
        #         end
        #     )
        # elseif Order > 2
        #     error("Multiplication for STTs of Order > 2 not yet implemented!")
        # end

        #     exprs = []


        # function (Base.:*)(stt1::StateTransitionTensor{Order}, stt2::StateTransitionTensor{Order}) where {Order}
        #     # Resulting timespan.
        #     # NOTE: Expects that stt1.tspan[1] == stt2.tspan[2] (i.e. the timespans line up like (a,b)*(c,a)).
        #     # Behaviour is unknown otherwise.
        #     multiplied_timespan = (stt2.tspan[1], stt1.tspan[2])

        #     # The first order tensors simply multiply together
        #     multiplied_tensors = TensorMap[]
        #     push!(multiplied_tensors, permute(stt1.tensors[1] * stt2.tensors[1], (1,), (2,)))

        #     $(exprs...)

        #     return StateTransitionTensor(multiplied_timespan, tuple(multiplied_tensors...))
        # end

        # Inverses and adjoints
        function (Base.inv)(stt::StateTransitionTensor{$(ORDER)})
            # [Refer to Park 2007, "Nonlinear trajectory navigation", Definition 2.2.12 (Inverse STTs)]
            inverse_tspan = (stt.tspan[end], stt.tspan[begin])
            inverse_tensors = SArray[]
            tensors = stt.tensors

            # The first-order tensor simply equals its inverse [Park.2007 eq.2.36]
            push!(inverse_tensors, inv(tensors[1]))

            $(inv_exprs[1:ORDER-1]...)

            StateTransitionTensor(inverse_tspan, tuple(inverse_tensors...))
        end

    end)
end

# Type promotion
Base.promote_rule(::Type{S}, ::Type{<:Union{AbstractMatrix,UniformScaling}}) where {Out,In,S<:StateTransitionMatrix{Out,In}} = SMatrix{Out,In,Float64}
Base.convert(::Type{S}, stm::StateTransitionMatrix) where {S<:AbstractMatrix} = convert(S, stm.tensors[1])
Base.convert(::Type{S}, uniform::UniformScaling) where {A,S<:SMatrix{A,A}} = S(uniform)

# Other forms of multiplication
(Base.:*)(stm1::StateTransitionMatrix{Out,In}, stm2::StateTransitionMatrix{Out,In}) where {Out,In} = StateTransitionTensor((stm2.tspan[1], stm1.tspan[2]), stm1.tensors .* stm2.tensors)
(Base.:*)(stt::StateTransitionTensor, other) = *(promote(stt, other)...)
(Base.:*)(other, stt::StateTransitionTensor) = *(promote(other, stt)...)

# Elementary operations on STMs
(Base.:+)(stm1::StateTransitionMatrix, other) = +(promote(stm1, other)...)
(Base.:+)(other, stm1::StateTransitionMatrix) = +(promote(other, stm1)...)

# Solving linear equation
(Base.:\)(stt::StateTransitionTensor, dx) = inv(stt) * dx

# Adjoint
(Base.adjoint)(stm::StateTransitionMatrix) = StateTransitionTensor(stm.tspan, adjoint.(stm.tensors))

@doc """ Sensitivity of the interpolated states (at times t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory, t; kwargs...) = StateTransitionTensor(sol(t); kwargs...)
StateTransitionTensor(sol::Trajectory, t::Number; kwargs...) = StateTransitionTensor(sol(t), (sol.t[begin], t); kwargs...)

@doc """ Sensitivity trace of each state (at times t=sol.t) with respect to initial state. """
StateTransitionTensor(sol::Trajectory; kwargs...) = [StateTransitionTensor(sol[i], (sol.t[begin], sol.t[i]); kwargs...) for i in 1:length(sol.t)]

@doc """ Sensitivity of the state at time t2 with respect to the state at time t1 <= t2. """
function StateTransitionTensor(sol::Trajectory, t1, t2; kwargs...)
    @assert t1 <= t2  "Expected t1 <= t2"
    if t1 == t2
        return one(STM(sol, t1))
    elseif t1 == sol.t[begin]
        return StateTransitionTensor(sol, t2; kwargs...)
    else
        # [Park.2007 eq.2.40]
        return StateTransitionTensor(sol, t2; kwargs...) * inv(StateTransitionTensor(sol, t1; kwargs...))
    end
end

# Compute STMs by solving the given state first
function StateTransitionTensor(::Val{:FiniteDiff}, state::State, desired_frame=state.frame, alg=DEFAULT_ALG; order=1, kwargs...)
    order == 1 || error("FiniteDiff STT only supports order 1 (STM), got $(order).")
    tensor = FiniteDiff.finite_difference_jacobian(state.u0) do u0
        new_state = remake(state, u0=u0)
        trajectory = solve(new_state, alg; kwargs...)
        trajectory_converted = convert_to_frame(trajectory, desired_frame)
        end_state_u = trajectory_converted.sol[end]
    end
    StateTransitionTensor(state.tspan, (SMatrix{size(tensor)...}(tensor),))
end

# Compute eigenvalues of an STT -- only uses the 1st-order tensor (STM)!
LinearAlgebra.eigvals(stt::StateTransitionTensor) = eigvals(Array(stt.tensors[1]))