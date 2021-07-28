export ER3BP
export elliptical_potential, centrifugal_potential

#-------------#
# ER3BP MODEL #
#-------------#

struct ER3BP{P<:R3BPSystemProperties} <: Abstract_R3BPModel
    props :: P
end

# Convenience helper
# const ER3BP = AstrodynamicalSystem{ER3BP_Model}
ER3BP(args...; kwargs...) = ER3BP(R3BPSystemProperties(args...; kwargs...))

#---------------------------#
# ER3BP EQUATIONS OF MOTION #
#---------------------------#

function ModelingToolkit.ODESystem(::Type{ER3BP})
    @parameters  μ  # Mass fraction
    @parameters  e  # Eccentricity
    @parameters  f  # True anomaly
    @variables   x(f) y(f) z(f)
    D, D2 = Differential(f), Differential(f)^2
    Dx, Dy, Dz = Differential.((x, y, z))

    # [DeiTos2017, Eqs. 14]
    # NOTE: need to expand the RHS derivatives before calling order-lowering
    ω = elliptical_potential(μ, (x, y, z), f, e)
    return ODESystem([
            D2(x) ~ +2D(y) + Dx(ω),
            D2(y) ~ -2D(x) + Dy(ω),
            D2(z) ~        + Dz(ω)
        ], 
        f,
        [x, y, z],
        [μ, e]
    )
end

#------------------#
# HELPER FUNCTIONS #
#------------------#

@doc "Centrifugal potential [DeiTos2017, Eq.13; Ichinomiya 2018, Eq. 2.2]"
function centrifugal_potential(μ, (x, y, z))
    r = distance_vectors(μ, (x, y, z))
    0.5*(x^2 + y^2) + (1 - μ)/r[1] + μ/r[2] + 0.5*(1 - μ)*μ
end

@doc "Elliptical potential in barycentric frame [DeiTos2017, Eq. 15]"
function elliptical_potential(μ, (x, y, z), f, e)
    Ω_3 = centrifugal_potential(μ, (x, y, z))
    return (Ω_3 - 0.5*z^2*e*cos(f)) / (1 + e*cos(f))
end

# XXX: The parameters overload below gives us a performance boost.
ModelingToolkit.parameters(model::ER3BP) = SVector(model.props.μ, model.props.e)

#-------------#
# COMPILATION #
#-------------#

# Build the equations at pre-compile time
const ER3BP_ODESystem = ODESystem(ER3BP)

# Retrieve the ER3BP ODESystem at runtime
ModelingToolkit.ODESystem(::ER3BP) = ER3BP_ODESystem