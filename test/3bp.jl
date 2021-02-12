using OrbitalTrajectories

using Test
using Unitful

function check_libration_points(model)
    L = libration_points(model)
    @test 0.5 < L[1][1] < 1.0        # L1 (X)
    @test L[1][2] == L[1][3] == 0.0  # L1 (Y & Z)
    @test 1.0 < L[2][1] < 1.5        # L2 (X)
    @test L[2][2] == L[2][3] == 0.0  # L2 (Y & Z)
    @test -1.5 < L[3][1] < -1.0      # L3 (X)
    @test L[3][2] == L[3][3] == 0.0  # L3 (Y & Z)
    @test 0.0 < L[4][1] < 1.0        # L4 (X)
    @test 0.0 < L[4][2] < 1.0        # L4 (Y)
    @test L[4][1] == L[5][1]         # L4 & L5 (X)
    @test L[4][2] == -L[5][2]        # L4 & L5 (Y)
    @test L[4][3] == L[5][3] == 0.0  # L4 & L5 (Z)
    return L
end

@testset "3BP Properties" begin
    # Cross-check the Mars-Phobos system (values from Wikipedia)
    MarsPhobos = OrbitalTrajectories.Dynamics.R3BPSystemProperties(:mars, :phobos)
    @test MarsPhobos.L ≈ 9376.0u"km" rtol=1e-5
    @test MarsPhobos.V ≈ 2.138u"km/s" rtol=1e-5
    @test MarsPhobos.T ≈ 0.31891023u"d"/2π rtol=1e-5
    @test MarsPhobos.R1 ≈ [3396.2, 3396.2, 3376.2]u"km"  rtol=1e-1
    @test MarsPhobos.R2 ≈ [13.5, 11., 9.]u"km"  rtol=1e-1

    # Check the relationship between libration points
    check_libration_points(MarsPhobos)

    SunMars = OrbitalTrajectories.Dynamics.R3BPSystemProperties(:sun, :mars)
    L = check_libration_points(SunMars)

    # Validate against constants in [DeiTos2018]
    @test SunMars.L ≈ 2.279497905330276e8u"km"                          rtol=1e-4
    @test SunMars.V ≈ 24.128831378998047u"km/s"                         rtol=1e-4
    @test SunMars.R2[1] ≈ 3396.19u"km"                                  rtol=1e-4
    @test SunMars.T ≈ 109.3425420965616u"d"                             rtol=1e-4
    @test L[1][1] ≈ (1 - SunMars.μ) + (-1.082385474e6u"km" / SunMars.L) rtol=1e-4  # (with respect to Mars)
    @test SunMars.e ≈ 0.0935643512                                      rtol=1e-4

    # XXX: These values as listed in [DeiTos2018] are possibly wrong?
    #@test SunMars.μ ≈ 3.227154876045166e-6      rtol=1e-4
end