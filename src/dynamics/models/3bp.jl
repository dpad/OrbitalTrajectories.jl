export libration_points, pos_to_inertial, state_to_inertial
export R3BPSystemProperties

#---------------------------#
# RESTRICTED 3-BODY SYSTEMS #
#---------------------------#

default_reference_frame(::Abstract_R3BPModel) = SynodicFrame()

Base.show(io::IO, x::Abstract_R3BPModel) = print(io, "$(nameof(typeof(x)))$(x.props)")

#-------------------#
# SYSTEM PROPERTIES # 
#-------------------#
struct R3BPSystemProperties{L <: Unitful.Length, V <: Unitful.Velocity, T <: Unitful.Time}
    b1 :: Symbol        # Identifier of body 1
    b2 :: Symbol        # Identifier of body 1
    μ  :: Float64       # Mass ratio
    L  :: L             # Distance between centers
    V  :: V             # Orbital velocity of larger body
    T  :: T             # Time unit: orbital period of smaller body (per radian)
    e  :: Float64       # Eccentricity
    R1 :: SVector{3,L}  # Radius of central body
    R2 :: SVector{3,L}  # Radius of smaller body
end

@memoize function R3BPSystemProperties(a::Symbol, b::Symbol; kwargs...)
    # Extract the static data for this system
    a, b = @. Symbol(lowercase(String((a, b))))
    data = get(_3BPSystemPropertyData, (a, b),
           get(_3BPSystemPropertyData, (b, a), nothing))
    isnothing(data) && error("No system data available for $a, $b")

    # Get user-provided values
    data = merge(data, kwargs)

    # Run-time compute remaining values
    μ = get(data, :μ, SpiceUtils.mass_fraction(a, b))
    R1 = get(data, :R1, SVector{3}(bodvrd(String(a), "RADII")u"km"))
    R2 = get(data, :R2, SVector{3}(bodvrd(String(b), "RADII")u"km"))
    # Ensure the Length attribute matches the Radius type
    L = convert(eltype(R1), data[:L])

    R3BPSystemProperties(a, b, μ, L, data[:V], data[:T], data[:e], R1, R2)
end
R3BPSystemProperties(system::Abstract_DynamicalModel) = R3BPSystemProperties(primary_body(system), secondary_body(system))

# TODO: Compute these constants at run-time?
const _3BPSystemPropertyData = Dict(
    (:mars, :phobos) => Dict(
        :L => 9.376e3u"km",
        :V => 2.138u"km/s",
        :T => 0.31891023u"d" / (2π),
        :e => 0.0151
    ),
    (:earth, :moon) => Dict(
        :L => 384399u"km",
        :V => 1.022u"km/s",
        :T => 27.321661u"d" / (2π),
        :e => 0.0549
    ),
    (:sun, :mars) => Dict(
        :L => 2.279497905330276e8u"km",
        :V => 24.128831378998047u"km/s",
        :T => 686.99u"d" / (2π),
        :e => 0.0935643512
    ),
    (:sun, :earth) => Dict(
        # Values from Wikipedia
        :L => 149598023u"km",
        :V => 29.78u"km/s",
        :T => 365.256363004u"d" / (2π),
        :e => 0.0167086
    ),
    (:jupiter, :europa) => Dict(
        # Values from Wikipedia, [Pellegrini.2006]
        :L => 670900u"km",
        :V => 13.740u"km/s",
        :T => 3.551181u"d" / (2π),
        :e => 0.009,
        :μ => 2.528009215182033e-5
    )
)

Base.show(io::IO, x::R3BPSystemProperties) = print(io, (x.b1, x.b2))
ModelingToolkit.parameters(model::Abstract_R3BPModel) = [getfield(model.props, i.name) for i in ModelingToolkit.parameters(model.ode.ode_system)]

#---------#
# METHODS #
#---------#
@memoize function libration_points(props::R3BPSystemProperties)
    # Begin with some guesses for the libration points
    μ = props.μ
    L_guess = [
        [(1 - μ)^2, 0., 0.],         # L1
        [(1 - μ)^-1, 0., 0.],        # L2
        [-1., 0., 0.],               # L3
        [(1 - μ)^2, (1 - μ^2), 0.],  # L4
        [(1 - μ)^2, -(1 - μ^2), 0.]  # L5
    ]
    
    # Solve for gradient(libration point guess) == 0 using forward AD
    libration_points = nlsolve.(
        (pos) -> ForwardDiff.gradient(p -> centrifugal_potential(μ, p), pos), 
        L_guess, autodiff=:forward)

    # Return the zero solution (i.e. the position) for each libration point
    map(l -> l.zero, libration_points)
end
libration_points(model::Abstract_R3BPModel) = libration_points(model.props)

function distance_vectors(μ, pos)
    x, y, z = pos
    ((μ + x)^2     + y^2 + z^2)^(1/2), # Distance to the central body
    ((1 - μ - x)^2 + y^2 + z^2)^(1/2)  # Distance to the smaller body
end

@doc "Compute the Hill sphere radius for two bodies"
function hill_sphere_radius(system::Abstract_R3BPModel)
    mass_1 = bodvrd(system.props.b1, "GM")[1]
    mass_2 = bodvrd(system.props.b2, "GM")[1]
    return  system.props.L * (1 - system.props.e) * cbrt(mass_2 / (3 * mass_1))
end


@doc """ Returns the position rotation matrix from a rotating 3BP frame to an inertial frame. """
pos_to_inertial(t) = [
    # [Pavlak 2013, eq. 2.47]
    cos(t)  -sin(t)  0;
    sin(t)   cos(t)  0;
    0        0       1
]
pos_from_inertial(t) = inv(pos_to_inertial(t))

@doc """ Returns the state rotation matrix from a rotating 3BP frame to an inertial frame. """
state_to_inertial(t) = [
    # [Pavlak 2013, eq. 2.50]
     cos(t)  -sin(t)  0  0        0       0;
     sin(t)   cos(t)  0  0        0       0;
     0        0       1  0        0       0;
    -sin(t)  -cos(t)  0  cos(t)  -sin(t)  0;
     cos(t)  -sin(t)  0  sin(t)   cos(t)  0;
     0        0       0  0        0       1
]
state_from_inertial(t) = inv(state_to_inertial(t))

primary_body(m::Abstract_R3BPModel) = m.props.b1
secondary_body(m::Abstract_R3BPModel) = m.props.b2

function collision(system::Abstract_R3BPModel, body, dist=bodvrd(String(body), "RADII")[1]; radii=1.)
    circ_props = R3BPSystemProperties(primary_body(system), secondary_body(system))
    diam = (radii * dist) ./ ustrip(u"km", circ_props.L)
    return ContinuousCallback((integrator) -> terminate!(integrator, :Crashed); interp_points=10) do u, t, integrator
        check_distance(u, t, system, body, diam)
    end
end

function check_distance(u, t, system::Abstract_R3BPModel, body::Symbol, distance)
    body_pos = (body == primary_body(system)) ? SVector{3}(-system.props.μ, 0., 0.) : SVector{3}(1 - system.props.μ, 0., 0.)
    return norm(@view(u[1:3]) - body_pos) - distance
end

function convert_to_frame(state::State{<:Abstract_R3BPModel,<:SynodicFrame}, frame::InertialFrame)
    to_inertial = state_to_inertial(state.tspan[1])
    converted_u0 = to_inertial * state.u0
    new_state = State(state.model, frame, converted_u0, state.tspan)
    return new_state
end

#----------#
# INCLUDES #
#----------#
include("er3bp.jl")
include("cr3bp.jl")  # Depends on ER3BP