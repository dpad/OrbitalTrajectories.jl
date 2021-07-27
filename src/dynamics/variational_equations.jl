#-------------------------------------------------#
# VARIATIONAL EQUATIONS FOR ASTRODYNAMICAL MODELS #
#-------------------------------------------------#

@doc """
    Generate ODE Systems that compute the State Transition
    Tensors simultaneously with the given system. This requires the system
    to have a computable Jacobian function.
"""

const MAX_VE_ORDER = 2

struct VarEqODESystem{Order} <: Abstract_AstrodynamicalODESystem
    ode_system  :: ODESystem  # The full combined system
    VE_system   :: ODESystem  # Only the differential variational equations for this order
    jacobian    :: Array{Num}  # Stores the symbolic jacobian expressions
end

@memoize function VarEqODESystem(system::ODESystem, order)
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
    vareqs = [VarEqODESystem(system, N) for N in  1:order-1]
    systems = (system, [sys.VE_system for sys in vareqs]...)
    jacs = Array{Num}[reshape(sys.jacobian, ntuple(_->dim, N+1)) for (N, sys) in enumerate(vareqs)]
    ϕ = Array{Num}[reshape(states(sys), ntuple(_->dim, N+1)) for (N, sys) in enumerate(systems[2:end])]

    # Create the STT array for this order
    dims = ntuple(_->1:dim, order+1)
    ϕ_name = Symbol("ϕ_$(order)")
    new_ϕ = (@variables $(ϕ_name)[dims...](iv))[1]
    push!(ϕ, new_ϕ)

    # Get the Jacobian matrix (A(t)) of the previous system, relative to the original state.
    push!(jacs, Array(reshape(compute_jacobian(systems[end], dvs; simplify=false), dims...)))

    # Contract the tensors as needed
    # NOTE: Depends on the Jacobian, corresponding to A(t) matrix (for the
    # State Transition Matrix) [Parker & Anderson 2014].
    RESULT = STT_diff_contraction(Val(order), jacs, ϕ)

    # NOTE: the Differential is defined element-wise and flattened to a list.
    D   = Differential(iv)
    m_eqs = D.(ϕ[end]) .~ RESULT
    eqs = collect(m_eqs)

    # Create the differentiated ODE system
    VE_system = ODESystem(vec(eqs), iv, vec(new_ϕ), [])
    return VarEqODESystem{order}(system, VE_system, jacs[end])
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