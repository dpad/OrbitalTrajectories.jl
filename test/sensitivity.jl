using OrbitalTrajectories

using Test
using LinearAlgebra

@testset "Pellegrini2016" begin
    # Values taken from Table 1 of [Pellegrini 2016, On the Computation and Accuracy of Trajectory State Transition Matrices]
    test_cases = [(
        system = (:Jupiter, :Europa),
        u0     = [-0.1017472008677258, 0., 0., 0., -0.01806028472285857, 0.],
        μ      = 2.528009215182033e-5,
        t_p    = 25.13898226959327,
        λ_max  = 2.468621114047195,
    ), (
        system = (:Jupiter, :Europa),
        u0     = [0.04867586089512202, 0., 0., 0., -0.09354853663949217, 0.],
        μ      = 2.528009215182033e-5,
        t_p    = 70.53945041512506,
        λ_max  = 2.804814246519340e7,
    ), (
        system = (:Earth, :Moon),
        u0     = [-0.013059020050628, 0., 0.07129515195874, 0., -0.526306975588415, 0.],
        μ      = 0.01215509906405700,
        t_p    = 2.517727406553485,
        λ_max  = 17.632688124231755,
    )]

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
        STM_AD = get_sensitivity(AD, state_cr3bp)
        STM_VE = get_sensitivity(VE, state_cr3bp)
        STM_VE_handcoded = get_sensitivity(VE, state_cr3bp_handcoded)

        # Compute the maximum eigenvalues
        λ_max_AD = maximum(norm.(eigvals(Matrix(STM_AD))))
        @test λ_max_AD ≈ case.λ_max rtol=1e-4
        λ_max_VE = maximum(norm.(eigvals(STM_VE)))
        @test λ_max_VE ≈ case.λ_max rtol=1e-4
        λ_max_VE_handcoded = maximum(norm.(eigvals(STM_VE_handcoded)))
        @test λ_max_VE_handcoded ≈ case.λ_max rtol=1e-4
    end
end

@testset "EphemerisNBP STM computation" begin
    u0 = [0.8574053516112442, 0., 0., 0., 0.47, 0.]
    dim = length(u0)
    system_nbp = EphemerisNBP(:earth, :moon)
    prob = State(system_nbp, SynodicFrame(), u0, (0., 3600.0*24))

    # Compute final STM
    STM_AD = get_sensitivity(AD, prob)
    STM_FD = get_sensitivity(FD, prob)
    @test STM_AD ≈ STM_FD rtol=1e-5

    # Compute STM trace
    STM_AD_trace = solve_sensitivity(AD, prob)
    @test get_sensitivity(STM_AD_trace, STM_AD_trace.t[end]) ≈ STM_AD rtol=1e-5
end