#-------------------------------------------------#
# VARIATIONAL EQUATIONS FOR ASTRODYNAMICAL MODELS #
#-------------------------------------------------#

@doc """
    Function to generate ODE Systems that compute the State Transition
    Tensors simultaneously with the given system. This requires the system
    to have a computable Jacobian function.
"""
function variational_equation_ODESystems(system::ModelingToolkit.AbstractODESystem; order=0, kwargs...)
    # TODO: Memoize this function! It's very slightly slow for EphemerisNBP
    if order <= 0
        error("Expected order >= 1 for variational equations (got $(order)).")
    elseif order > MAX_STT_ORDER
        # In theory, higher orders should be supported, but we want to make sure
        # the user isn't accidentally shooting themselves in the foot (e.g.
        # accidentally supply order=1000). So for now, we just error out.
        error("Variational equations not currently implemented beyond order-$(MAX_STT_ORDER) derivatives (got $(order)).")
    end

    # Keep track of all the systems (including the original one)
    systems = ODESystem[system]

    for N in 1:order
        # The system we want to differentiate is the previous one
        ode = systems[N-1]

        # Get the existing ODE system's information.
        iv  = independent_variable(ode)
        dvs = states(ode)

        # The total dimension of the differentiated system depends on the length
        # of the states
        dim = length(dvs)

        # NOTE: Depends on the Jacobian, corresponding to A(t) matrix (for the
        # State Transition Matrix) [Parker & Anderson 2014].
        @variables   ϕ[1:dim,1:dim](iv)
        D = Differential(iv)

        # Get the Jacobian matrix (A(t))
        A = calculate_jacobian(ode)

        # The State Transition Matrix (STM) ODE function is defined as follows, including the N^2 Jacobian equations +
        # the N first-order equations of motion. [Koon 2011]
        # NOTE: the Differential is defined element-wise and flattened to a list.
        stm_eqs = collect(simplify.(D.(ϕ) .~ A * ϕ))

        # Create the differentiated ODE system
        stm_ode = ODESystem(stm_eqs, iv, ϕ, [])

        # Add this differentiated systems to our list
        push!(stm_ode, systems)
    end

    # Return just the differentiated systems (i.e. not including the original one)
    return systems[2:end]
end