using OrbitalTrajectories

using Test
using LinearAlgebra
using DiffEqSensitivity
using StaticArrays
using ForwardDiff

@testset "Pellegrini2016" begin
    # Values taken from Table 1 of [Pellegrini 2016, On the Computation and Accuracy of Trajectory State Transition Matrices]
    test_cases = [(
        system = (:Jupiter, :Europa),
        u0     = [-0.1017472008677258, 0., 0., 0., -0.01806028472285857, 0.],
        μ      = 2.528009215182033e-5,
        t_p    = 25.13898226959327,
        λ_max  = 2.468621114047195,
        STT2_norm = 1.881880182286279e4
    ), (
        system = (:Jupiter, :Europa),
        u0     = [0.04867586089512202, 0., 0., 0., -0.09354853663949217, 0.],
        μ      = 2.528009215182033e-5,
        t_p    = 70.53945041512506,
        λ_max  = 2.804814246519340e7,
        STT2_norm = 8.290465853749079e17
    ), (
        system = (:Earth, :Moon),
        u0     = [-0.013059020050628, 0., 0.07129515195874, 0., -0.526306975588415, 0.],
        μ      = 0.01215509906405700,
        t_p    = 2.517727406553485,
        λ_max  = 17.632688124231755,
        STT2_norm = 5.182997297997671e6
    )]

    case = test_cases[1]

    for case in test_cases
        cr3bp = CR3BP(case.system...; μ=case.μ)
        cr3bp_handcoded = HandcodedCR3BP(case.system...; μ=case.μ)
        
        # The u0 from the paper is given relative to secondary body, so convert here to primary body
        u0 = deepcopy(case.u0)
        u0[1] += (1 - cr3bp.props.μ)

        # Initial state
        state_cr3bp = State(cr3bp, SynodicFrame(), u0, (0., case.t_p))
        state_cr3bp_handcoded = State(cr3bp_handcoded, SynodicFrame(), u0, (0., case.t_p))

        # Compute STMs
        STM_AD = STM(AD, state_cr3bp)
        STM_VE = STM(VE, state_cr3bp)
        STM_VE_handcoded = STM(VE, state_cr3bp_handcoded)

        @test STM_AD[end] ≈ STM_VE[end] rtol=1e-5
        @test STM_VE[end] ≈ STM_VE_handcoded[end] rtol=1e-5

        # Compute the maximum eigenvalues
        λ_max_AD = maximum(norm.(eigvals(STM_AD[end])))
        @test λ_max_AD ≈ case.λ_max rtol=1e-4
        λ_max_VE = maximum(norm.(eigvals(STM_VE[end])))
        @test λ_max_VE ≈ case.λ_max rtol=1e-4
        λ_max_VE_handcoded = maximum(norm.(eigvals(STM_VE_handcoded[end])))
        @test λ_max_VE_handcoded ≈ case.λ_max rtol=1e-4

        # Try to get 2nd-order STT and norm of its 2nd-order Tensor
        STT_AD = STT(AD, state_cr3bp; order=2)
        STT_AD_2_norm = norm(STT_AD[end].tensors[2])
        @test STT_AD_2_norm ≈ case.STT2_norm rtol=1e-2
    end
end

@testset "EphemerisNBP STM computation" begin
    u0 = [0.8574053516112442, 0., 0., 0., 0.47, 0.]
    dim = length(u0)
    system_nbp = EphemerisNBP(:earth, :moon)
    prob = State(system_nbp, SynodicFrame(), u0, (0., 3600.0*24))

    # Compute final STM
    STM_AD = STM(AD, prob)[end]
    STM_VE = STM(VE, prob)[end]
    STM_FD = STM(FD, prob)
    @test STM_AD ≈ STM_FD rtol=1e-5
    @test STM_AD ≈ STM_VE rtol=1e-5

    # Compute STM trace
    STM_AD_trace = solve_sensitivity(AD, prob)
    STM_VE_trace = solve_sensitivity(VE, prob)
    @test STM(STM_AD_trace[end]) ≈ STM_AD rtol=1e-5
    @test STM(STM_AD_trace, STM_AD_trace.t[end]) ≈ STM_AD rtol=1e-5
end

@testset "Higher-order STTs" begin
    system = CR3BP(:moon, :earth)

    ode_system = system.ode.ode_system
    ode_f = system.ode.ode_f

    # using ModelingToolkit

    test = OrbitalTrajectories.Dynamics.VarEqODESystem(ode_system, 1)

    system = CR3BP(:moon, :earth; VE_order=1)

    OrbitalTrajectories.Dynamics.STM_ODEFunction(ode_system, ode_f)
end