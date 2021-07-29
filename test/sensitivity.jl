using OrbitalTrajectories

using Test
using LinearAlgebra
using DifferentialEquations
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

        # Try to get 2nd-order STT from variational equations model
        # XXX: This seemingly works, but causes instabilities, so the solution
        # stops short.
        sol_AD = solve_sensitivity(AD, state_cr3bp; order=2)
        sol_VE_2 = solve_sensitivity(VE, state_cr3bp; order=2)

        STT(sol_AD(0.5)).tensors[1]
        STT(sol_VE_2(0.5)).tensors[1]
        STT(sol_AD(0.5)).tensors[2][:,:,5:6]
        STT(sol_VE_2(0.5)).tensors[2][:,:,5:6]
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


@testset "Example of VarEqs on simple model" begin
    using ModelingToolkit

    struct MySimpleModel4{S,F} <: OrbitalTrajectories.Dynamics.Abstract_AstrodynamicalODESystem
        ode_system :: S
        ode_f      :: F
    end

    struct MyMODEL{O} <: OrbitalTrajectories.Dynamics.Abstract_AstrodynamicalModel
        ode :: O
        props :: Nothing
    end
    OrbitalTrajectories.Dynamics.default_reference_frame(model::MyMODEL) = SynodicFrame()
    ModelingToolkit.parameters(model::MyMODEL) = []

    function ModelingToolkit.ODESystem(::Type{MySimpleModel4})
        @parameters  t  # time
        @variables   x(t) y(t)
        D, D2 = Differential(t), Differential(t)^2

        return ODESystem([
                D2(x) ~ 5x^2 + y^2,
                D2(y) ~ x + x*y^3,
            ], 
            t,
            [x, y],
            []
        )
    end

    model = MySimpleModel4()
    model.ode_f.f.f_oop
    model_VE_1 = OrbitalTrajectories.Dynamics.VarEqModel_ODESystem(model.ode_system, 1)
    model_VE_2 = OrbitalTrajectories.Dynamics.VarEqModel_ODESystem(model.ode_system, 2)

    model_VE_1.ode_f.f.f_oop
    model_VE_2.ode_f.f.f_oop

    system = MyMODEL(model, nothing)
    system_1 = with_var_eqs(system)
    system_2 = with_var_eqs(system, 2)

    # x(t)'   = x_t(t)    
    # x_t(t)' = 5 * x(t)^2
    #
    # Jacobian should be:
    # [ 0   1] [1,1  1,2]  = [2,1            2,2]
    # [10x  0] [2,1  2,2]  = [10x(1,1)  10x(1,2)]
    model_VE_1.jacobian
    system_1.ode.jacobian
    # phi[]  = 
    # 2nd-order Jacobian (Hessian?)
    model_VE_2.jacobian
    system_2.ode.jacobian

    # Compute a contraction
    jacs = [system_1.ode.jacobian, system_2.ode.jacobian]
    ϕ = [reshape(states(system_1.ode.VE_system), length(u0), length(u0)), reshape(states(system_2.ode.VE_system), length(u0), length(u0), length(u0))]
    RESULT = OrbitalTrajectories.Dynamics.STT_diff_contraction(Val(2), jacs, ϕ)
    RESULT[4,4,4]

    i, a, b = (4, 4, 4)
    x = Num(0)
    for α in 1:4
        x += jacs[1][i,α] * ϕ[2][α,a,b]
        for β in 1:4
            x += jacs[2][i,α,β] * ϕ[1][α,a] * ϕ[1][β,b]
        end
    end
    x

    # Use ForwardDiff to propagate the duals
    using ForwardDiff
    using DiffEqSensitivity
    u0 = [0.5, 0.3, 3.9, 4.5]
    val = model.ode_f.f.f_oop(u0, [], 0.)
    val_VE_1 = model_VE_1.ode_f.f.f_oop([u0..., vec(Array(I, length(u0), length(u0)))...], [], 0.)
    tensor_1 = OrbitalTrajectories.Dynamics.extract_sensitivity_tensors(val_VE_1, length(u0), 1)[1]
    val_VE_2 = model_VE_2.ode_f.f.f_oop([u0..., vec(Array(I, length(u0), length(u0)))..., zeros(length(u0)^3)...], [], 0.)
    tensor_2 = OrbitalTrajectories.Dynamics.extract_sensitivity_tensors(val_VE_2, length(u0), 2)[2]

    # Now with states
    state = State(system, u0, (0., 0.5))
    sol = solve(state)
    sol_1 = solve_sensitivity(VE, state; order=1)
    sol_2 = solve_sensitivity(VE, state; order=2)
    sol_AD_1 = solve_sensitivity(AD, state; order=1)
    sol_AD_2 = solve_sensitivity(AD, state; order=2)

    # STTs
    STT_VE = STT(sol_2)
    STT_AD = STT(sol_AD_2)

    @test STT_VE[end].tensors[1] ≈ STT_AD[end].tensors[1]
    STT_VE[end].tensors[2] 
    STT_AD[end].tensors[2]

    # Plot to check!
    using Plots
    plot([sol.t, sol.t, sol.t, sol.t], [sol.sol[1,:], sol.sol[2,:], sol.sol[3,:], sol.sol[4,:]])
    plot([sol_1.t, sol_1.t, sol_1.t, sol_1.t], [sol_1.sol[1,:], sol_1.sol[2,:], sol_1.sol[3,:], sol_1.sol[4,:]])
    plot([sol_2.t, sol_2.t, sol_2.t, sol_2.t], [sol_2.sol[1,:], sol_2.sol[2,:], sol_2.sol[3,:], sol_2.sol[4,:]])
    plot([sol_AD_1.t, sol_AD_1.t, sol_AD_1.t, sol_AD_1.t], [ForwardDiff.value.(sol_AD_1.sol[1,:]), ForwardDiff.value.(sol_AD_1.sol[2,:]), ForwardDiff.value.(sol_AD_1.sol[3,:]), ForwardDiff.value.(sol_AD_1.sol[4,:])])
    plot([sol_AD_2.t, sol_AD_2.t, sol_AD_2.t, sol_AD_2.t], [sol_AD_2.sol[1,:], sol_AD_2.sol[2,:], sol_AD_2.sol[3,:], sol_AD_2.sol[4,:]])

    dual_u0 = DiffEqSensitivity.seed_duals(copy(u0), typeof(model))
    dual_val = model.ode_f.f.f_oop(dual_u0, [], 0.)
    dual_tensor_1 = OrbitalTrajectories.Dynamics.extract_sensitivity_tensors(dual_val, length(u0), 1)[1]
    dual_u0_2 = DiffEqSensitivity.seed_duals(copy(dual_u0), typeof(model))
    dual_val_2 = model.ode_f.f.f_oop(dual_u0_2, 0., [])
    dual_tensor_2 = OrbitalTrajectories.Dynamics.extract_sensitivity_tensors(dual_val_2, length(u0), 2)[2]
end