using OrbitalTrajectories

using Test
using DifferentialEquations
using BenchmarkTools

@testset "convert_to_frame" begin
    u0 = [0.8574053516112442, 0., 0., 0., 0.47, 0.]
    system_nbp = EphemerisNBP(:earth, :moon)
    prob = State(system_nbp, SynodicFrame(), u0, (0., 3600.0*24*30))

    # Test that converting to and from inertial frames works
    inert_prob = convert_to_frame(prob, InertialFrame())
    synod_prob = convert_to_frame(inert_prob, SynodicFrame())
    @test synod_prob.prob.u0 ≈ u0

    # Test that converting to/from non-normalised frames also works
    synod_nonorm_prob = convert_to_frame(prob, SynodicFrame(false))
    synod_norm_prob = convert_to_frame(synod_nonorm_prob, SynodicFrame(true))
    @test synod_norm_prob.prob.u0 ≈ prob.prob.u0
    inert_prob2 = convert_to_frame(synod_nonorm_prob, InertialFrame())
    synod_nonorm_prob2 = convert_to_frame(inert_prob2, SynodicFrame(false))
    @test synod_nonorm_prob2.prob.u0 ≈ synod_nonorm_prob.prob.u0

    # Roughly in decreasing order of allocations (30μs~40μs, allocations 8 to 38)
    # @btime convert_to_frame(prob, InertialFrame())  # norm -> inert
    # @btime convert_to_frame(inert_prob, SynodicFrame())  # inert -> norm
    # @btime convert_to_frame(prob, SynodicFrame(false))  # norm -> unorm
    # @btime convert_to_frame(synod_nonorm_prob, SynodicFrame(true))  # unorm -> norm
    # @btime convert_to_frame(synod_nonorm_prob, InertialFrame())  # unorm -> inert
    # @btime convert_to_frame(inert_prob2, SynodicFrame(false))  # inert -> unorm

    # Solve and test that converting the solution to synodic also matches
    sol = solve(prob)
    sol2 = convert_to_frame(sol, SynodicFrame())
    @test sol2.sol.u[begin] ≈ u0
end

@testset "EphemerisNBP performance with/without CSE" begin
    u0 = [0.8574053516112442, 0., 0., 0., 0.47, 0.]
    bodies = (:earth, :moon)

    # NOTE: This test only works for EphemerisNBP (other models remain equally performant)
    model = EphemerisNBP
    sys = model(bodies...)

    tspan = (0., 3600.0*24*30) 
    prob = State(sys, SynodicFrame(), u0, tspan)

    # Test default performance of prob (with CSE)
    time_new = @benchmark solve($prob)

    # Test the old performance of prob (turn off Common Subexpression Elimination)
    sys_old = model(bodies...; wrap_code=(nothing, nothing))
    prob_old = State(sys_old, SynodicFrame(), u0, tspan)
    time_old = @benchmark solve($prob_old)

    @test 10 * median(time_new.times) < median(time_old.times)  # At least a 10x improvement
end