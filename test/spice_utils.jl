using OrbitalTrajectories

using Test

@testset "SpiceUtils" begin
    @test OrbitalTrajectories.SpiceUtils.mass_fraction("jupiter", "europa") ≈ 2.528e-5 atol=1e-5
    @test OrbitalTrajectories.SpiceUtils.mass_fraction("mars", "phobos") ≈ 1.0659e16 / (1.0659e16 + 6.4171e23) atol=1e-5
end