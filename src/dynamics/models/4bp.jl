#---------------------------#
# RESTRICTED 4-BODY SYSTEMS #
#---------------------------#

default_reference_frame(::Abstract_R4BPModel) = SynodicFrame()

#-------------------#
# SYSTEM PROPERTIES # 
#-------------------#
struct R4BPSystemProperties{N<:Number} <: Abstract_ModelProperties
    b1 :: Symbol  # Identifier of body 1 (central body)
    b2 :: Symbol  # Identifier of body 2 (smaller body)
    b3 :: Symbol  # Identifier of body 3 (large distance body)
    m3 :: N       # Normalised mass of body 3 (relative to mass of R3BP)
    a3 :: N       # Normalised distance between body 3 and R3BP barycenter
    ω3 :: N       # Angular velocity of body 3 in synodic coordinates
    μ  :: N       # Mass ratio of the r3bp system
    μ2 :: N       # Mass ratio of the whole system
    α0 :: N       # Initial phase angle of body 3
end

# @memoize function R4BPSystemProperties(a::Symbol, b::Symbol, c::Symbol)
function R4BPSystemProperties(a::Symbol, b::Symbol, c::Symbol, α0::Number=0.; kwargs...)
    # NOTE: c is the 3rd body (e.g. the Sun), a and b are the primaries
    # Get the 3-body problem properties
    r3bp = R3BPSystemProperties(a, b)

    # Run-time compute remaining values
    GM = @. SpiceUtils.get_GM((a, b, c))u"km^3/s^2"
    m3 = GM[3] / (GM[1] + GM[2])
    map(SpiceUtils.load_ephemerides, (a, b, c))
    state_3 = Array(SpiceUtils.get_state(0., a, c))
    elements_3 = oscltx(state_3, 0., ustrip(u"km^3/s^2", GM[3]))
    a3 = get(kwargs, :a3, elements_3[10]u"km" / r3bp.L)
    ω3 = get(kwargs, :ω3, √((1 + m3) / a3^3) - 1)  # [DeiTos2017 eq.18]
    μ2 = get(kwargs, :μ2, sum(GM) / (GM[1] + GM[2]))
    μ = get(kwargs, :μ, r3bp.μ)

    R4BPSystemProperties(a, b, c, m3, a3, ω3, μ, μ2, α0)
end

Base.show(io::IO, ::MIME"text/plain", x::R4BPSystemProperties) = print(io, parameters(x))#(a=x.b1, b=x.b2, c=x.b3, α0=x.α0))
ModelingToolkit.parameters(props::R4BPSystemProperties) = NamedTuple([i => getfield(props, i) for i in propertynames(props)])
ModelingToolkit.parameters(model::Abstract_R4BPModel) = [getfield(model.props, i.name) for i in ModelingToolkit.parameters(model.ode.ode_system)]

#---------#
# METHODS #
#---------#

primary_body(m::Abstract_R4BPModel) = m.props.b1
secondary_body(m::Abstract_R4BPModel) = m.props.b2

#----------#
# INCLUDES #
#----------#
include("bc4bp.jl")