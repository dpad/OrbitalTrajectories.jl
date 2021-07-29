#-------------------------------------------------#
# VARIATIONAL EQUATIONS FOR ASTRODYNAMICAL MODELS #
#-------------------------------------------------#

export with_var_eqs

const MAX_VE_ORDER = 2

struct VarEqModel_ODESystem{Order, F} <: Abstract_VariationalEquationsODESystem{Order}
    ode_system  :: ODESystem   # The full combined ODE system
    ode_f       :: F           # The full ODE system function
    VE_system   :: ODESystem   # Only the differential variational equations for this order
    jacobian    :: Array{Num}  # Stores the symbolic jacobian expressions
end

struct VarEqModel{Order, O<:Abstract_VariationalEquationsODESystem{Order}, M<:Abstract_AstrodynamicalModel} <: Abstract_AstrodynamicalModel
    ode         :: O  # The new ODE system
    model       :: M  # The original model
end

function Base.getproperty(x::T, b::Symbol) where {T<:VarEqModel}
    if hasfield(T, b)
        return getfield(x, b)
    else
        return getproperty(x.model, b)
    end
end
default_reference_frame(model::VarEqModel) = default_reference_frame(model.model)
ModelingToolkit.parameters(model::VarEqModel) = parameters(model.model)
state_length(m::VarEqModel) = state_length(m.model)
get_order(::State{M}) where {Order,M<:VarEqModel{Order}} = Order
primary_body(model::VarEqModel) = primary_body(model.model)
secondary_body(model::VarEqModel) = secondary_body(model.model)

Base.show(io::IO, ::MIME"text/plain", x::Abstract_VariationalEquationsODESystem{Order}) where {Order} = print(io, string(" with ", SciMLBase.TYPE_COLOR, "Order=$(Order) var.eqs.", SciMLBase.NO_COLOR))
function Base.show(io::IO, M::MIME"text/plain", x::VarEqModel{Order}) where {Order}
    print(io, string(SciMLBase.TYPE_COLOR, "VarEqModel", SciMLBase.NO_COLOR, "{"))
    show(io, M, x.model) 
    show(io, M, x.ode)
    print(io, SciMLBase.NO_COLOR, "}")
end

@doc """Generic State with default initial values for VarEq systems."""
function State(model::VarEqModel{Order}, reference_frame::Abstract_ReferenceFrame, u0::AbstractArray, tspan) where {Order}
    # Get the problem u0
    u0 = deepcopy(convert(Array, u0))

    # Get the model dimensions
    dim = state_length(model.model)

    if length(u0) == dim && Order > 0
        # User only provided initial state vector.
        # We need to enlarge it to include the variational equation states!

        # 1st-order sensitivities start as the Identity matrix
        uType = eltype(u0)
        append!(u0, vec(Matrix{eltype(uType)}(I, dim, dim)))

        # 2nd-order and above are just zeros
        for order in 2:Order
            append!(u0, zeros(eltype(uType), dim^(order+1)))
        end
    end

    # The final state size should match
    expected_dim = sum([dim^N for N in 1:(Order+1)])
    length(u0) == expected_dim || error("Expected u0 vector of length $(expected_dim), got $(length(u0)).")

    # Build the state
    State(model, reference_frame, ODEProblem(model, u0, tspan, parameters(model.model)))
end

function with_var_eqs(model::M, order=1; force=false, kwargs...) where {M<:Abstract_AstrodynamicalModel}
    # Create the original model and get the system ODE
    if model.ode isa Abstract_VariationalEquationsODESystem && !force
        error("Cannot build variational equations for an already variational system (unless you set force=true)!")
    end
    VE_system = VarEqModel_ODESystem(model.ode.ode_system, order; kwargs...)
    VarEqModel{order, typeof(VE_system), typeof(model)}(VE_system, model)
end

@doc """
    Generate ODE Systems that compute the State Transition
    Tensors simultaneously with the given system. This requires the system
    to have a computable Jacobian function.
"""
@memoize function VarEqModel_ODESystem(system::ODESystem, order; wrap_code=(cse, cse))
    if order < 1
        error("Expected order >= 1 for variational equations (got $(order)).")
    elseif order > MAX_VE_ORDER
        error("Variational equations not currently implemented beyond order-$(MAX_VE_ORDER) derivatives (got $(order)).")
    end

    # Get properties of the system
    iv  = independent_variable(system)
    dvs = states(system)
    dim = length(dvs)

    # Build all models for lower orders (should be memoized so should not be expensive)
    vareqs = [VarEqModel_ODESystem(system, N) for N in  1:order-1]
    systems = [system, [sys.VE_system for sys in vareqs]...]
    jacs = Array{Num}[reshape(sys.jacobian, ntuple(_->dim, N+1)) for (N, sys) in enumerate(vareqs)]
    ϕ = Array{Num}[reshape(states(sys), ntuple(_->dim, N+1)) for (N, sys) in enumerate(systems[2:end])]

    @info "Building VarEqModel of order $(order) for state dimension $(dim)"

    # Create the STT array for this order
    dims = ntuple(_->1:dim, order+1)
    ϕ_name = Symbol("ϕ_$(order)")
    new_ϕ = (@variables $(ϕ_name)[dims...](iv))[1]
    push!(ϕ, new_ϕ)

    # Get the Jacobian matrix (A(t)) of the previous system, relative to the original state.
    push!(jacs, Array(reshape(compute_jacobian(systems[end], dvs; simplify=true), dims...)))

    # Contract the tensors as needed
    # NOTE: Depends on the Jacobian, corresponding to A(t) matrix (for the
    # State Transition Matrix) [Parker & Anderson 2014].
    RESULT = STT_diff_contraction(Val(order), jacs, ϕ)

    # NOTE: the Differential is defined element-wise and flattened to a list.
    D   = Differential(iv)
    m_eqs = D.(ϕ[end]) .~ RESULT
    eqs = collect(m_eqs)

    # Create the differentiated ODE system
    push!(systems, ODESystem(vec(eqs), iv, vec(new_ϕ), []))

    # Create the full system combining 
    full_eqs = [eq for eqs in equations.(systems) for eq in eqs]
    full_dvs = [state for sys in systems for state in states(sys)]
    full_system = ODESystem(full_eqs, iv, full_dvs, parameters(system))

    # Generate the full system function
    ode_f = ODEFunction(full_system; 
        # Don't need to differentiate further
        jac=false, 
        tgrad=false, 
        # Make sure the function is evaluated in current module (especially for package precompile)
        eval_expression=false, 
        eval_module=@__MODULE__,
        # Perform any further transformations on the code
        wrap_code
    )

    return VarEqModel_ODESystem{order, typeof(ode_f)}(full_system, ode_f, systems[end], jacs[end])
end

function STT_diff_contraction(::Val{1}, jacs, ϕ)
    # Order-1 (STM -- simple matrix multiplication)
    jac_1 = jacs[1]
    ϕ_1 = ϕ[1]
    @tullio RESULT[i,a] := jac_1[i,α] * ϕ_1[α,a]
end

function STT_diff_contraction(::Val{2}, jacs, ϕ)
    # Order-2 STT
    jac_1, jac_2 = jacs
    ϕ_1, ϕ_2 = ϕ
    @tullio RESULT[i,a,b] := jac_1[i,α] * ϕ_2[α,a,b] + jac_2[i,α,β] * ϕ_1[α,a] * ϕ_1[β,b]
end

function compute_jacobian(system::ModelingToolkit.AbstractODESystem, dvs; kwargs...)
    rhs = [eq.rhs for eq ∈ equations(system)]
    return Symbolics.jacobian(rhs, dvs; kwargs...)
end