export CR3BP, HandcodedCR3BP
export jacobi

#---------------------------#
# CR3BP EQUATIONS OF MOTION #
#---------------------------#
struct CR3BP_ODESystem{S,F} <: Abstract_AstrodynamicalODESystem
    ode_system :: S
    ode_f      :: F
end

function ModelingToolkit.ODESystem(::Type{CR3BP_ODESystem})
    # Build from the ER3BP equations (with eccentricity = 0 for circular)
    eqs_er3bp = ODESystem(ER3BP_ODESystem)
    (μ, e) = parameters(eqs_er3bp)
    eqs = [eq.lhs ~ simplify(substitute(eq.rhs, e => 0)) for eq in equations(eqs_er3bp)]
    return ODESystem(eqs, independent_variable(eqs_er3bp), states(eqs_er3bp), [μ]; name=:CR3BP)
end

# Build the equations at pre-compile time
const CR3BP_ODEFunctions = CR3BP_ODESystem()

#-------------#
# CR3BP MODEL #
#-------------#

struct CR3BP{O<:Abstract_AstrodynamicalODESystem,P<:R3BPSystemProperties} <: Abstract_R3BPModel
    ode   :: O
    props :: P
end

CR3BP(args...; kwargs...) = CR3BP(R3BPSystemProperties(args...; kwargs...))

# TODO: Move hand-coded functions into :analytical of CR3BP
function CR3BP(props::R3BPSystemProperties; kwargs...)
    if length(kwargs) == 0
        ode = CR3BP_ODEFunctions
    else
        ode = CR3BP_ODESystem(; kwargs...)
    end
    CR3BP{typeof(ode), typeof(props)}(ode, props)
end

# XXX: The parameters overload below gives us a performance boost.
ModelingToolkit.parameters(model::CR3BP) = SVector(model.props.μ)

#---------#
# METHODS #
#---------#

@doc "Jacobi integral = 2U - V^2 = -2E [Oshima 2019, Eq. 3; Koon, Eq. 2.3.13; Wakker, Eq. 3.54]"
jacobi(μ, vel, pos) = 2*centrifugal_potential(μ, pos) - norm(vel)^2


#---------------------------#
# HANDCODED CR3BP FUNCTIONS #
#---------------------------#

# Hand-coded CR3BP
struct HandcodedCR3BP_ODESystem{F} <: Abstract_AstrodynamicalODESystem
    ode_f      :: F
end

struct HandcodedCR3BP_VarEqODESystem{Order,F} <: Abstract_VariationalEquationsODESystem{Order}
    ode_f      :: F
end

const CR3BP_with_Handcoded = CR3BP{<:Union{HandcodedCR3BP_ODESystem,HandcodedCR3BP_VarEqODESystem}}

function Base.show(io::IO, M::MIME"text/plain", x::CR3BP_with_Handcoded)
    print(io, string(SciMLBase.TYPE_COLOR, nameof(typeof(x)), SciMLBase.NO_COLOR))
    show(io, M, x.props)
    print(io, " with hand-coded equations")
end

HandcodedCR3BP(args...; kwargs...) = HandcodedCR3BP(R3BPSystemProperties(args...; kwargs...))
function HandcodedCR3BP(props::R3BPSystemProperties; kwargs...)
    ode = HandcodedCR3BP_ODESystem(handcodedCR3BP_f)
    CR3BP{typeof(ode), typeof(props)}(ode, props)
end

ModelingToolkit.parameters(model::CR3BP_with_Handcoded) = SVector(model.props.μ)
state_length(::CR3BP_with_Handcoded) = 6

function with_var_eqs(model::CR3BP{<:HandcodedCR3BP_ODESystem}, order=1)
    order == 1 || error("HandcodedCR3BP only supports 1st-order variational equations!")
    ode = HandcodedCR3BP_VarEqODESystem{1, typeof(handcodedCR3BP_withSTM_f)}(handcodedCR3BP_withSTM_f)
    VarEqModel{order, typeof(ode), typeof(model)}(ode, model)
end

function handcodedCR3BP_f(du, u, p, t)
    @inbounds begin
        μ = p[1]
        x, y, z, x̂, ŷ, ẑ = u
        r1, r2 = (1 - μ, μ) ./ (distance_vectors(μ, (x, y, z)).^3)
        du[1] = x̂
        du[2] = ŷ
        du[3] = ẑ
        du[4] =  2ŷ  + x - r1*(μ + x) + r2*(1 - μ - x)
        du[5] = -2x̂  + y - r1*y       - r2*y
        du[6] =          - r1*z       - r2*z
    end
end

function handcodedCR3BP_withSTM_f(du, u, p, t)
    @inbounds begin
        μ = p[1]
        x, y, z, x̂, ŷ, ẑ = u
        (d, r) = distance_vectors(μ, (x, y, z))
        d3, d5 = d^3, d^5
        r3, r5 = r^3, r^5
        r1, r2 = (1 - μ, μ) ./ ((d, r).^3)
        A, B = 3*(1-μ)/d5, 3*μ/r5

        du[1] = x̂
        du[2] = ŷ
        du[3] = ẑ
        du[4] =  2ŷ  + x - r1*(μ + x) + r2*(1 - μ - x)
        du[5] = -2x̂  + y - r1*y       - r2*y
        du[6] =          - r1*z       - r2*z

        XX = 1 - (1-μ)/d3 - μ/r3 + A*(x+μ)^2 + B*(x - 1 + μ)^2
        XY = A*(x+μ)*y + B*(x-1+μ)*y
        XZ = A*(x+μ)*z + B*(x-1+μ)*z
        YY = 1 - (1-μ)/d3 - μ/r3 + A*y^2 + B*y^2
        YZ = A*y*z + B*y*z
        ZZ = -(1-μ)/d3 - μ/r3 + A*z^2 + B*z^2

        du[7] = u[10]
        du[8] = u[11]
        du[9] = u[12]
        du[10] = XX*u[7] + XY*u[8] + XZ*u[9] + 2*u[11]
        du[11] = XY*u[7] + YY*u[8] + YZ*u[9] - 2*u[10]
        du[12] = XZ*u[7] + YZ*u[8] + ZZ*u[9]

        du[13] = u[16]
        du[14] = u[17]
        du[15] = u[18]
        du[16] = XX*u[13] + XY*u[14] + XZ*u[15] + 2*u[17]
        du[17] = XY*u[13] + YY*u[14] + YZ*u[15] - 2*u[16]
        du[18] = XZ*u[13] + YZ*u[14] + ZZ*u[15]

        du[19] = u[22]
        du[20] = u[23]
        du[21] = u[24]
        du[22] = XX*u[19] + XY*u[20] + XZ*u[21] + 2*u[23]
        du[23] = XY*u[19] + YY*u[20] + YZ*u[21] - 2*u[22]
        du[24] = XZ*u[19] + YZ*u[20] + ZZ*u[21]

        du[25] = u[28]
        du[26] = u[29]
        du[27] = u[30]
        du[28] = XX*u[25] + XY*u[26] + XZ*u[27] + 2*u[29]
        du[29] = XY*u[25] + YY*u[26] + YZ*u[27] - 2*u[28]
        du[30] = XZ*u[25] + YZ*u[26] + ZZ*u[27]

        du[31] = u[34]
        du[32] = u[35]
        du[33] = u[36]
        du[34] = XX*u[31] + XY*u[32] + XZ*u[33] + 2*u[35]
        du[35] = XY*u[31] + YY*u[32] + YZ*u[33] - 2*u[34]
        du[36] = XZ*u[31] + YZ*u[32] + ZZ*u[33]

        du[37] = u[40]
        du[38] = u[41]
        du[39] = u[42]
        du[40] = XX*u[37] + XY*u[38] + XZ*u[39] + 2*u[41]
        du[41] = XY*u[37] + YY*u[38] + YZ*u[39] - 2*u[40]
        du[42] = XZ*u[37] + YZ*u[38] + ZZ*u[39]
    end
end