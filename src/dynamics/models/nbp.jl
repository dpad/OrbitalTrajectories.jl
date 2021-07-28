export EphemerisNBP
export get_state

#----------------#
# N-BODY SYSTEMS #
#----------------#
struct NBPSystemProperties{N}
    center :: Symbol
    center_id :: Integer
    bodies :: SVector{N, Symbol}
    body_ids :: SVector{N, <: Integer}
    μ :: SVector{N, <: Real}
end

function NBPSystemProperties(center::Symbol, bodies::Vararg{Symbol})
    # Load ephemerides for all bodies...
    map(SpiceUtils.load_ephemerides, bodies)
    bodies = SVector(bodies...)

    center_id = bodn2c(String(center))
    body_ids = @. bodn2c(String(bodies))
    μ = SpiceUtils.get_GM.(body_ids)
    NBPSystemProperties(center, center_id, bodies, body_ids, μ)
end

Base.show(io::IO, ::MIME"text/plain", x::NBPSystemProperties) = print(io, "($(x.center), acc=$(x.bodies))")

#-----------------------------------#
# EPHEMERIS-NBP EQUATIONS OF MOTION #
#-----------------------------------#
struct NBP_ODESystem{S,F} <: Abstract_AstrodynamicalODESystem
    ode_system :: S
    ode_f      :: F
end

# XXX: Need the "T" in place for @memoize to work
@memoize function ModelingToolkit.ODESystem(T::Type{NBP_ODESystem}, props::NBPSystemProperties)
    @parameters t  # Time in J2000 epoch
    @variables x(t) y(t) z(t)
    D2 = Differential(t)^2

    accelerations = zeros(Num, 3)
    pos = [x, y, z]

    for (body, μ) in zip(props.bodies, props.μ)
        # See also [DeiTos2017, Eqs. 14], [JTOP implementation], [Ozaki2017, Eq.2]

        # Get the position of the body
        if body == props.center
            body_pos = [0., 0., 0.]
        else
            # Add the non-inertial center's acceleration (i.e. center is also pulled by this body) [?]
            # NOTE: the body and center need to be converted to Num{String} for the @registered function to work.
            #       Num{Symbol} won't work because they convert to variable names instead of Symbols.
            body_pos_term = SpiceUtils.get_pos(t, String(body), String(props.center))

            # XXX: Need to manually build the vector since the compiler won't know the size of body_pos otherwise.
            # TODO: Fix the allocation cost of this. view() and reshape() are not defined for Nums...
            body_pos = [body_pos_term[1], body_pos_term[2], body_pos_term[3]]
            accelerations += -μ * body_pos / norm(body_pos)^3
        end

        # Compute the acceleration of the spacecraft towards this body
        distance_vector = body_pos - pos
        accelerations += μ * distance_vector / norm(distance_vector)^3
    end

    # Equations of motion: sum of all accelerations
    eqs = @. D2(pos) ~ sum(accelerations)

    # Build the 2nd-order ODE props
    return ODESystem(eqs, t, pos, [])
end

#---------------------#
# EPHEMERIS-NBP MODEL #
#---------------------#
struct EphemerisNBP{O<:NBP_ODESystem,P<:NBPSystemProperties} <: Abstract_AstrodynamicalModel{O}
    ode   :: O
    props :: P
end
function EphemerisNBP(bodies::Vararg{Symbol}; center=nothing, kwargs...)
    (!isnothing(center) && isnothing(bodn2c(center)) || center ∈ bodies) && error("Unsupported center (must not be specified in bodies)")
    all_bodies = isnothing(center) ? (bodies[1], bodies...) : (center, bodies...)
    bodies_symbols = @. Symbol(lowercase(String(all_bodies)))
    props = NBPSystemProperties(bodies_symbols...)
    EphemerisNBP(NBP_ODESystem(props; kwargs...), props)
end

ModelingToolkit.parameters(model::EphemerisNBP) = SVector{0,Float64}()

# HELPERS

# TODO: Check that target is in the system?
get_state(system::EphemerisNBP, target::Symbol, t) =
    State(system, SpiceUtils.get_state(t, target, system.props.center), (t, t))

# TRAITS
default_reference_frame(::EphemerisNBP) = InertialFrame()
default_synodic_reference_frame(::EphemerisNBP) = SynodicFrame(false)

primary_body(m::EphemerisNBP) = m.props.center
secondary_body(m::EphemerisNBP) = m.props.bodies[2]

convert_to_frame(state::State{<:EphemerisNBP,F}, frame::F) where {F<:Abstract_ReferenceFrame} = state
function convert_to_frame(state::State{<:EphemerisNBP,<:Abstract_ReferenceFrame}, frame::Abstract_ReferenceFrame)
    # Get the synodic initial state of the secondary body
    # TODO: Fix up the SpiceUtils interface so the arguments are more consistently ordered.
    to_synodic = SpiceUtils.state_to_synodic(primary_body(state), secondary_body(state), state.prob.tspan[1])
    inv_synodic = inv(to_synodic)

    @assert length(state.prob.u0) >= STATE_DIMS
    # XXX: @inbounds @views is sometimes causing crashes below!!!
    # @inbounds @views begin
    converted_u0 = state_to_frame(state, frame, to_synodic, inv_synodic)
    # end

    prob1 = remake(state.prob; u0=converted_u0)
    state = State(state.model, frame, prob1)
    return state
end

function state_to_frame(state::State{<:EphemerisNBP,InertialFrame}, frame::SynodicFrame{true}, to_synodic, inv_synodic)
    u1 = state_to_frame(state, SynodicFrame(false), to_synodic, inv_synodic)
    return state_to_frame(State(state.model, SynodicFrame(false), u1, state.prob.tspan), frame, to_synodic, inv_synodic)
end

function state_to_frame(state::State{<:EphemerisNBP,SynodicFrame{false}}, ::SynodicFrame{true}, to_synodic, inv_synodic)
    secondary_u0 = SpiceUtils.get_state(state.prob.tspan[1], secondary_body(state), primary_body(state))
    synod_secondary_u0 = to_synodic * secondary_u0
    # Normalise the state
    circ_props = R3BPSystemProperties(primary_body(state), secondary_body(state))
    u1 = copy(state.prob.u0)
    P, V = norm(secondary_u0[1:3]), norm(secondary_u0[4:6])
    u1[1:3] .= (u1[1:3] - synod_secondary_u0[1:3]) ./ P + [1 - circ_props.μ, 0., 0.]
    u1[4:6] .= (u1[4:6] - synod_secondary_u0[4:6]) ./ V

    if length(u1) == (STATE_DIMS + STATE_DIMS*STATE_DIMS)
        # XXX: If there are vareqn. variables, cross-multiply them appropriately
        vareqns = reshape(@view(u1[STATE_DIMS+1:end]), (STATE_DIMS, STATE_DIMS))
        pos_vareqns = @view vareqns[1:3, :]
        vel_vareqns = @view vareqns[4:6, :]
        pos_vareqns .= pos_vareqns ./ P
        vel_vareqns .= vel_vareqns ./ V
    end

    return u1
end

function state_to_frame(state::State{<:EphemerisNBP,<:SynodicFrame{true}}, ::SynodicFrame{false}, to_synodic, inv_synodic)
    secondary_u0 = SpiceUtils.get_state(state.prob.tspan[1], secondary_body(state), primary_body(state))
    synod_secondary_u0 = to_synodic * secondary_u0
    # De-normalise
    circ_props = R3BPSystemProperties(primary_body(state), secondary_body(state))
    u1 = copy(state.prob.u0)
    P, V = norm(secondary_u0[1:3]), norm(secondary_u0[4:6])
    u1[1:3] .= (u1[1:3] - [1 - circ_props.μ, 0., 0.]) .* P + synod_secondary_u0[1:3]
    u1[4:6] .= u1[4:6] .* V + synod_secondary_u0[4:6]

    if length(u1) == (STATE_DIMS + STATE_DIMS*STATE_DIMS)
        # XXX: If there are vareqn. variables, cross-multiply them appropriately
        vareqns = reshape(@view(u1[STATE_DIMS+1:end]), (STATE_DIMS, STATE_DIMS))
        pos_vareqns = @view vareqns[1:3, :]
        vel_vareqns = @view vareqns[4:6, :]
        pos_vareqns .= pos_vareqns .* P
        vel_vareqns .= vel_vareqns .* V
    end

    return u1
end

function state_to_frame(state::State{<:EphemerisNBP,<:SynodicFrame{true}}, frame, to_synodic, inv_synodic)
    u1 = state_to_frame(state, SynodicFrame(false), to_synodic, inv_synodic)
    return state_to_frame(State(state.model, SynodicFrame(false), u1, state.prob.tspan), frame, to_synodic, inv_synodic)
end

state_to_frame(state::State{<:EphemerisNBP,<:InertialFrame}, ::SynodicFrame{false}, to_synodic, inv_synodic) =
    u0_to_frame(state.prob.u0, to_synodic)
state_to_frame(state::State{<:EphemerisNBP,<:SynodicFrame{false}}, ::InertialFrame, to_synodic, inv_synodic) =
    u0_to_frame(state.prob.u0, inv_synodic)

# Dispatch on MVectors
const STATE_DIMS = 6
u0_to_frame(u0::StaticVector{STATE_DIMS}, transform::AbstractMatrix) = transform * u0
function u0_to_frame(u0, transform::AbstractMatrix)
    # TODO: What happens if u0 is a StaticVector? Things break :(

    # XXX: Only transform the state vector, not any var.eqns.
    new_u0 = copy(u0)
    new_u0[begin:STATE_DIMS] .= transform * @view(new_u0[1:STATE_DIMS])

    if length(new_u0) == (STATE_DIMS + STATE_DIMS*STATE_DIMS)
        # XXX: If there are vareqn. variables, cross-multiply them appropriately
        vareqns = reshape(@view(new_u0[STATE_DIMS+1:end]), (STATE_DIMS, STATE_DIMS))
        vareqns .= transform * vareqns
    end

    return new_u0
end


function check_distance(u, t, system::EphemerisNBP, body)
    diam = bodvrd(String(body), "RADII")[1]  # km
    return check_distance(u, t, system, body, diam)
end

function check_distance(u, t, system::EphemerisNBP, body, distance)
    if body == primary_body(system)
        body_pos = SVector{3}(0., 0., 0.)
    else
        body_pos = SpiceUtils.get_pos(t, body, primary_body(system))
    end
    return @inbounds(norm(@view(u[1:3]) - body_pos) - distance)
end