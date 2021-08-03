export BC4BP

#---------------------------#
# BC4BP EQUATIONS OF MOTION #
#---------------------------#
struct BC4BP_ODESystem{S,F} <: Abstract_AstrodynamicalODESystem
    ode_system :: S
    ode_f      :: F
end

function ModelingToolkit.ODESystem(::Type{BC4BP_ODESystem})
    @parameters  μ   # Mass fraction (smaller 2 bodies)
    @parameters  μ2  # Mass fraction (note: inverse, as per [DeiTos2018])
    @parameters  a3  # Distance to 3rd body
    @parameters  α0  # Initial phase angle of 3rd body
    @parameters  ω3  # Angular velocity of 3rd body
    @parameters  t   # Time
    @variables   x(t) y(t) z(t)
    D, D2 = Differential(t), Differential(t)^2

    # [Dei Tos 2018]
    p   = [x, y, z]                  # Spacecraft position
    α   = ω3 * t
    p_S = a3 .* [cos(α + α0), sin(α + α0), 0]  # Position of 3rd body
    x̂   = [1, 0, 0]                  # X unit vector

    vectors = @. [
        # (Coefficient,     vector to be normalised) 
        (-(μ2 - 1),         p - p_S),
        (-(μ2 - 1)*(1 - μ), μ*x̂ + p_S),
        (+(μ2 - 1)*μ,       (1 - μ)*x̂ - p_S),
        (-(1 - μ),          p + μ*x̂),
        (-μ,                p - (1 - μ)*x̂)
    ]

    # Add up the equations [DeiTos2018 Eq. 21]
    forces = sum([coeff .* (v / norm(v)^3) for (coeff, v) in vectors])
    eqs    = @. D2(p) ~ [+2D(y) + x, -2D(x) + y, 0] + forces

    # Simplify and return the system
    return ODESystem(simplify.(eqs), t, [x, y, z], [μ, μ2, a3, α0, a3, ω3])
end

# Build the equations at pre-compile time
const BC4BP_ODEFunctions = BC4BP_ODESystem()

#-------------#
# BC4BP MODEL #
#-------------#
struct BC4BP{O<:BC4BP_ODESystem,P<:R4BPSystemProperties} <: Abstract_R4BPModel
    ode   :: O
    props :: P
end
BC4BP(args...; kwargs...) = BC4BP(BC4BP_ODEFunctions, R4BPSystemProperties(args...; kwargs...))