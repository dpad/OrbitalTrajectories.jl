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

@testset "Higher-order STTs" begin

    # Now a custom multi-value function
    function f(vals)
        x, y = vals
        [x^2, 0.5x^2*y + y^3, cos(y)]
    end

    vals = SVector{2}([1.5, 9])

    # 1st-order Duals
    # f(vals_dual)[1] corresponds to x^2 --> the partials are ∂(x²)/∂xᵢ
    # f(vals_dual)[2] corresponds to 0.5x^2 + y^3 --> the partials are ∂(0.5x²y + y³)/∂xᵢ
    # f(vals_dual)[3] corresponds to cos(y) --> the partials are ∂(cos(y))/∂xᵢ
    struct MyTag end
    vals_dual = DiffEqSensitivity.seed_duals(vals, MyTag)
    T = eltype(vals_dual[1]).parameters[1]  # Get the tag type
    f_vals_dual = f(vals_dual)
    @test f_vals_dual[1].value ≈ vals[1]^2 rtol=1e-5
    @test f_vals_dual[1].partials[1] ≈ 2*vals[1] rtol=1e-5
    @test f_vals_dual[1].partials[2] ≈ 0 rtol=1e-5
    @test f_vals_dual[2].value ≈ 0.5*vals[1]^2*vals[2] + vals[2]^3 rtol=1e-5
    @test f_vals_dual[2].partials[1] ≈ 0.5*2*vals[1]*vals[2] rtol=1e-5
    @test f_vals_dual[2].partials[2] ≈ 0.5*vals[1]^2 + 3*vals[2]^2 rtol=1e-5
    @test f_vals_dual[3].value ≈ cos(vals[2]) rtol=1e-5
    @test f_vals_dual[3].partials[1] ≈ 0 rtol=1e-5
    @test f_vals_dual[3].partials[2] ≈ -sin(vals[2]) rtol=1e-5

    STM1 = StateTransitionTensor(f_vals_dual)
    STM2 = StateTransitionTensor(STM1)

    # So the Jacobian is aligned with the first axis (rows) is the output, and
    # second axis (cols) is the input.
    J = ForwardDiff.extract_jacobian(T, f_vals_dual, vals)

    # 2nd-order Duals
    # f(vals_dual_2ndorder)[1] corresponds to x^2 --> the partials are ∂(x²)/∂xᵢ
    # f(vals_dual_2ndorder)[2] corresponds to 0.5x^2 + y^3 --> the partials are ∂(0.5x²y + y³)/∂xᵢ
    # f(vals_dual_2ndorder)[3] corresponds to cos(y) --> the partials are ∂(cos(y))/∂xᵢ
    vals_dual_2ndorder = DiffEqSensitivity.seed_duals(DiffEqSensitivity.seed_duals(vals, MyTag), MyTag)
    T2 = eltype(vals_dual_2ndorder[1]).parameters[1]  # Get the tag type
    f_vals_dual_2ndorder = f(vals_dual_2ndorder)
    @test f_vals_dual_2ndorder[1].value.value ≈ vals[1]^2 rtol=1e-5
    @test f_vals_dual_2ndorder[1].value.partials[1] ≈ 2*vals[1] rtol=1e-5
    @test f_vals_dual_2ndorder[1].value.partials[2] ≈ 0 rtol=1e-5
    @test f_vals_dual_2ndorder[1].partials[1].value ≈ 2*vals[1] rtol=1e-5  # ∂(x²)/∂x = 2x
    @test f_vals_dual_2ndorder[1].partials[1].partials[1] ≈ 2 rtol=1e-5  # ∂²(x²)/∂x∂x = 2
    @test f_vals_dual_2ndorder[1].partials[1].partials[2] ≈ 0 rtol=1e-5  # ∂²(x²)/∂y∂x = 0
    @test f_vals_dual_2ndorder[1].partials[2].value ≈ 0 rtol=1e-5  # ∂(x²)/∂y = 0
    @test f_vals_dual_2ndorder[1].partials[2].partials[1] ≈ 0 rtol=1e-5  # ∂²(x²)/∂x∂y = 0
    @test f_vals_dual_2ndorder[1].partials[2].partials[2] ≈ 0 rtol=1e-5  # ∂²(x²)/∂y∂y = 0

    @test f_vals_dual_2ndorder[2].value.value ≈ 0.5*vals[1]^2*vals[2] + vals[2]^3 rtol=1e-5
    @test f_vals_dual_2ndorder[2].value.partials[1] ≈ vals[1]*vals[2] rtol=1e-5
    @test f_vals_dual_2ndorder[2].value.partials[2] ≈ 0.5*vals[1]^2 + 3*vals[2]^2 rtol=1e-5
    @test f_vals_dual_2ndorder[2].partials[1].value ≈ vals[1]*vals[2] rtol=1e-5  # ∂()/∂x = x*y
    @test f_vals_dual_2ndorder[2].partials[1].partials[1] ≈ vals[2] rtol=1e-5  # ∂²()/∂x∂x = y
    @test f_vals_dual_2ndorder[2].partials[1].partials[2] ≈ vals[1] rtol=1e-5  # ∂²()/∂y∂x = x
    @test f_vals_dual_2ndorder[2].partials[2].value ≈ 0.5*vals[1]^2 + 3*vals[2]^2 rtol=1e-5  # ∂()/∂y = 0.5x^2 + 3y^2
    @test f_vals_dual_2ndorder[2].partials[2].partials[1] ≈ vals[1] rtol=1e-5  # ∂²()/∂x∂y = x
    @test f_vals_dual_2ndorder[2].partials[2].partials[2] ≈ 6*vals[2] rtol=1e-5  # ∂²()/∂y∂y = 6*y

    @test f_vals_dual_2ndorder[3].value.value ≈ cos(vals[2]) rtol=1e-5
    @test f_vals_dual_2ndorder[3].value.partials[1] ≈ 0  rtol=1e-5
    @test f_vals_dual_2ndorder[3].value.partials[2] ≈ -sin(vals[2]) rtol=1e-5
    @test f_vals_dual_2ndorder[3].partials[1].value ≈ 0 rtol=1e-5  # ∂()/∂x = 0
    @test f_vals_dual_2ndorder[3].partials[1].partials[1] ≈ 0 rtol=1e-5  # ∂²()/∂x∂x = 0
    @test f_vals_dual_2ndorder[3].partials[1].partials[2] ≈ 0 rtol=1e-5  # ∂²()/∂y∂x = 0
    @test f_vals_dual_2ndorder[3].partials[2].value ≈ -sin(vals[2]) rtol=1e-5  # ∂()/∂y = -sin(y)
    @test f_vals_dual_2ndorder[3].partials[2].partials[1] ≈ 0 rtol=1e-5  # ∂²()/∂x∂y = 0
    @test f_vals_dual_2ndorder[3].partials[2].partials[2] ≈ -cos(vals[2]) rtol=1e-5  # ∂²()/∂y∂y = -cos(y)

    J2 = StateTransitionTensor(f_vals_dual_2ndorder)
    ForwardDiff.extract_jacobian(eltype(J2).parameters[1], J2[:,1], vals)

    STM1 = STT1{Val(1)};

    function extract_STT(vals)
    end

    # T2 = eltype(vals_dual_2ndorder[1]).parameters[1]  # Get the tag type

    u0 = [0.8574053516112442, 0., 0., 0., 0.47, 0.]
    dim = length(u0)
    system = CR3BP(:earth, :moon)
    prob = State(system, SynodicFrame(), u0, (0., 10.))
    sol_AD = solve_sensitivity(AD, prob; order=1)

    plot(sol_AD)

    # 1st-order STMs
    STM_AD = StateTransitionMatrix(sol_AD)


    STM_AD[1]


    using DifferentialEquations
    half_orbit = solve(prob, DiffCorrectAxisymmetric(); dc_tolerance=1e-10, verbose=true)







    STM_FD = StateTransitionMatrix(FD, prob)
    @test STM_AD[end] ≈ STM_FD rtol=1e-5

    # Compare to the ForwardDiff Jacobian
    STM_AD_jac = ForwardDiff.jacobian(prob.u0) do u0
        state = remake(prob; u0)
        solve(state).sol[end]
    end
    @test STM_AD[end] ≈ STM_AD_jac rtol=1e-5

    # 2nd-order STMs
    sol_AD_2 = solve_sensitivity(AD, prob; order=2)
    STM_AD_2 = StateTransitionMatrix(sol_AD_2)

    # Compare to the ForwardDiff Hessian
    STM_AD_hes = ForwardDiff.hessian(prob.u0) do u0
        state = remake(prob; u0)
        solve(state).sol[end]
    end
    # @test STM_AD[end] ≈ STM_AD_jac rtol=1e-5
end