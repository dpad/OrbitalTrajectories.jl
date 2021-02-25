export State, Trajectory
export primary_body, secondary_body
export collision, check_distance, crashed

#------------------#
# ORBITAL PROBLEMS #
#------------------#

struct State{M<:Abstract_DynamicalModel,F<:Abstract_ReferenceFrame,O<:DiffEqBase.DEProblem} <: DiffEqBase.DEProblem
    model :: M
    frame :: F  # Reference frame that the problem's u0 is defined in
    prob :: O
end
State(model::Abstract_DynamicalModel, reference_frame::Abstract_ReferenceFrame, u0::AbstractArray, tspan) =
    State(model, reference_frame, MArray{Tuple{size(u0)...}}(u0), tspan)
State(model::Abstract_DynamicalModel, reference_frame::Abstract_ReferenceFrame, u0::StaticArray, tspan) =
    State(model, reference_frame, ODEProblem(model, u0, tspan, parameters(model)))
State(model, u0, tspan) = State(model, default_reference_frame(model), u0, tspan)

struct Trajectory{M<:Abstract_DynamicalModel,F<:Abstract_ReferenceFrame,T,N,A,O<:DiffEqBase.AbstractTimeseriesSolution{T,N,A},} <: DiffEqBase.AbstractTimeseriesSolution{T,N,A}
    model :: M
    frame :: F  # Reference frame that the solution is defined in
    sol :: O
end

State(traj::Trajectory) = State(traj.model, traj.frame, traj.sol.prob)

primary_body(state::State) = primary_body(state.model)
primary_body(traj::Trajectory) = primary_body(traj.model)
secondary_body(state::State) = secondary_body(state.model)
secondary_body(traj::Trajectory) = secondary_body(traj.model)

DiffEqBase.remake(state::State; kwargs...) = State(state.model, state.frame, remake(state.prob; kwargs...))
ModelingToolkit.parameters(state::State) = ModelingToolkit.parameters(state.model)

#---------------#
# INTERPOLATION #
#---------------#

# Interpolation
(traj::Trajectory)(t::Number) = State(traj.model, traj.frame, traj.sol(t), (t, t))

# Indexing
Base.getindex(state::State, idx...) = getindex(state.prob.u0, idx...)
Base.getindex(traj::Trajectory, idx...) = State(traj.model, traj.frame, getindex(traj.sol, idx...), (traj.sol.t[idx...], traj.sol.t[idx...]))
Base.getindex(traj::Trajectory, idx::Int) = State(traj.model, traj.frame, getindex(traj.sol, idx), (traj.sol.t[idx], traj.sol.t[idx]))
Base.getindex(traj::Trajectory, idx::AbstractArray{Int}) = State(traj.model, traj.frame, getindex(traj.sol, idx), (traj.sol.t[idx], traj.sol.t[idx]))
Base.axes(state::State, idx...) = axes(state.prob, idx...)
Base.axes(traj::Trajectory, idx...) = axes(traj.sol, idx...)
Base.firstindex(traj::Trajectory) = firstindex(traj.sol)
Base.firstindex(traj::Trajectory, idx) = firstindex(traj.sol, idx)
Base.lastindex(traj::Trajectory) = lastindex(traj.sol)
Base.lastindex(traj::Trajectory, idx) = lastindex(traj.sol, idx)
Base.size(traj::Trajectory) = size(traj.sol)

function Base.getproperty(x::T, b::Symbol) where {T<:Union{State,Trajectory}}
    if hasfield(T, b)
        return getfield(x, b)
    else
        return getproperty(isa(x, State) ? x.prob : x.sol, b)
    end
end

#---------#
# DISPLAY #
#---------#

Base.show(io::IO, A::State) = println(io, summary(A))
Base.summary(state::State) = string(
    SciMLBase.TYPE_COLOR, nameof(typeof(state)), SciMLBase.NO_COLOR, " in ",
    SciMLBase.TYPE_COLOR, state.model, SciMLBase.NO_COLOR, " in ",
    SciMLBase.TYPE_COLOR, state.frame, "\n",
    "    ", SciMLBase.NO_COLOR, "Underlying ", summary(state.prob), "\n",
    "    ", SciMLBase.NO_COLOR, "tspan = ", state.prob.tspan, "\n",
    "    u0    = ", state.prob.u0)
Base.show(io::IO, m::MIME"text/plain", A::Trajectory) = show(io, A) # XXX: Required because also defined in DiffEqBase
function Base.show(io::IO, A::Trajectory)
    println(io, string(
        SciMLBase.TYPE_COLOR, nameof(typeof(A)), SciMLBase.NO_COLOR, " in ",
        SciMLBase.TYPE_COLOR, A.model, SciMLBase.NO_COLOR, " in ",
        SciMLBase.TYPE_COLOR, A.frame, SciMLBase.NO_COLOR))
    println(io, string(
        "    retcode  = $(A.retcode) [$(length(A.t)) timesteps]\n",
        "    t        = ($(A.t[begin]), $(A.t[end]))\n",
        "    u[begin] = $(A.u[begin])\n",
        "    u[end]   = $(A.u[end])\n"
    ))
end

#-------#
# SOLVE #
#-------#

const DEFAULT_ALG = Vern7();

# XXX: Required to support solving a State problem.
DiffEqBase.solve(state::State, args...; reltol=1e-10, abstol=1e-10, kwargs...) =
    DiffEqBase.__solve(state, args...; reltol, abstol, kwargs...)

DiffEqBase.__solve(state::State; kwargs...) = DiffEqBase.__solve(state, DEFAULT_ALG; kwargs...)

# The __solve() method does the actual heavy lifting, including converting to a Trajectory.
function DiffEqBase.__solve(state::State, alg::DiffEqBase.DEAlgorithm; userdata=Dict(), kwargs...)
    default_frame = default_reference_frame(state.model)
    real_state = convert_to_frame(state, default_frame)

    # Pass the default state into the underlying solver
    # TODO: Remove the need for this in DiffCorrectAxisymmetric
    userdata_new = deepcopy(userdata)
    userdata_new[:real_state] = real_state

    # Call the underlying solver
    raw_sol = solve(real_state.prob, alg; userdata=userdata_new, kwargs...)
    return Trajectory(state.model, default_frame, raw_sol)
end