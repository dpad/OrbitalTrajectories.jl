export State, Trajectory
export primary_body, secondary_body
export collision, check_distance, crashed

#----------------#
# ORBITAL STATES #
#----------------#

struct State{S<:Abstract_AstrodynamicalSystem,F<:Abstract_ReferenceFrame,uType,tType,isinplace,O<:SciMLBase.AbstractODEProblem{uType,tType,isinplace}} <: SciMLBase.AbstractODEProblem{uType,tType,isinplace}
    system :: S
    frame :: F  # Reference frame that the problem's u0 is defined in
    prob :: O
end

State(system::Abstract_AstrodynamicalSystem, reference_frame::Abstract_ReferenceFrame, u0::AbstractArray, tspan) =
    State(system, reference_frame, MArray{Tuple{size(u0)...}}(u0), tspan)

State(system::Abstract_AstrodynamicalSystem, reference_frame::Abstract_ReferenceFrame, u0::StaticArray, tspan) =
    State(system, reference_frame, ODEProblem(system, u0, tspan, parameters(system)))

State(system::Abstract_AstrodynamicalSystem, u0::AbstractArray, tspan) = 
    State(system, default_reference_frame(system), u0, tspan)

# Remake a State
DiffEqBase.remake(state::State; kwargs...) = State(state.system, state.frame, remake(state.prob; kwargs...))

# Indexing into a State
Base.getindex(state::State, idx...) = getindex(state.prob.u0, idx...)
Base.axes(state::State, idx...) = axes(state.prob, idx...)

#----------------------#
# ORBITAL TRAJECTORIES #
#----------------------#

struct Trajectory{S<:Abstract_AstrodynamicalSystem,F<:Abstract_ReferenceFrame,T,N,A,O<:DiffEqBase.AbstractTimeseriesSolution{T,N,A},} <: DiffEqBase.AbstractTimeseriesSolution{T,N,A}
    system :: S
    frame :: F  # Reference frame that the solution is defined in
    sol :: O
end

# Convert a Trajectory back to a State (allowing it to be propagated again)
State(traj::Trajectory) = State(traj.system, traj.frame, traj.u[begin], (traj.t[begin], traj.t[end]))

# Properties of a Trajectory
primary_body(traj::Trajectory) = primary_body(traj.system)
secondary_body(traj::Trajectory) = secondary_body(traj.system)

# Indexing into a Trajectory
Base.getindex(traj::Trajectory, idx...) = State(traj.system, traj.frame, getindex(traj.sol, idx...), (traj.sol.t[idx...], traj.sol.t[idx...]))
Base.getindex(traj::Trajectory, idx::Int) = State(traj.system, traj.frame, getindex(traj.sol, idx), (traj.sol.t[idx], traj.sol.t[idx]))
Base.getindex(traj::Trajectory, idx::AbstractArray{Int}) = State(traj.system, traj.frame, getindex(traj.sol, idx), (traj.sol.t[idx], traj.sol.t[idx]))

Base.firstindex(traj::Trajectory) = firstindex(traj.sol)
Base.firstindex(traj::Trajectory, idx) = firstindex(traj.sol, idx)

Base.lastindex(traj::Trajectory) = lastindex(traj.sol)
Base.lastindex(traj::Trajectory, idx) = lastindex(traj.sol, idx)

Base.axes(traj::Trajectory, idx...) = axes(traj.sol, idx...)
Base.size(traj::Trajectory) = size(traj.sol)

#----------------------#
# PROPERTIES           #
#----------------------#

function Base.getproperty(x::T, b::Symbol) where {T<:Union{State,Trajectory}}
    if hasfield(T, b)
        return getfield(x, b)
    else
        return getproperty(isa(x, State) ? x.prob : x.sol, b)
    end
end

#---------------#
# INTERPOLATION #
#---------------#

# Interpolation of a Trajectory at time t
(traj::Trajectory)(t::Number) = State(traj.system, traj.frame, traj.sol(t), (t, t))

# Interpolation of a Trajectory at any list of times t
function (traj::Trajectory)(t::AbstractArray{<:Number})
    new_sol = DiffEqBase.build_solution(traj.prob, traj.alg, t, traj.sol(t); interp=traj.interp, retcode=traj.retcode)
    return Trajectory(traj.system, traj.frame, new_sol)
end

#---------#
# DISPLAY #
#---------#

Base.show(io::IO, ::MIME"text/plain", ::Type{S}) where {M2,R,S<:State{<:Abstract_AstrodynamicalSystem{M2},R}} =
    print(io, "State{$(nameof(M2)),$(nameof(R))}")

function Base.show(io::IO, M::MIME"text/plain", state::State{S}) where {M2,S<:Abstract_AstrodynamicalSystem{M2}}
    if get(io, :compact, false)
        print(io, "State{$(nameof(M2))}(t=$(state.prob.tspan), u0=$(state.prob.u0))")
    else
        print(io, string(
            SciMLBase.TYPE_COLOR, nameof(typeof(state)), SciMLBase.NO_COLOR, " in "))
        show(io, M, state.system)
        print(io, string(SciMLBase.NO_COLOR, " in ", SciMLBase.TYPE_COLOR))
        show(io, M, state.frame)
        println("\n", string(
        "    ", SciMLBase.NO_COLOR, "Underlying ", summary(state.prob), "\n",
        "    ", SciMLBase.NO_COLOR, "tspan = ", state.prob.tspan, "\n",
        "    u0    = ", state.prob.u0)
        )
    end
end

Base.show(io::IO, ::MIME"text/plain", ::Type{T}) where {M,R,T<:Trajectory{M,R}} =
    print(io, "Trajectory{$(nameof(M)),$(nameof(R))}")

function Base.show(io::IO, ::MIME"text/plain", A::Trajectory)
    if get(io, :compact, false)
        print(io, "Trajectory{$(nameof(typeof(A.system)))}(t=$((A.t[begin], A.t[end])))")
    else
        println(io, string(
            SciMLBase.TYPE_COLOR, nameof(typeof(A)), SciMLBase.NO_COLOR, " in ",
            SciMLBase.TYPE_COLOR, A.system, SciMLBase.NO_COLOR, " in ",
            SciMLBase.TYPE_COLOR, A.frame, SciMLBase.NO_COLOR))
        println(io, string(
            "    retcode  = $(A.retcode) [$(length(A.t)) timesteps]\n",
            "    t        = ($(A.t[begin]), $(A.t[end]))\n",
            "    u[begin] = $(A.u[begin])\n",
            "    u[end]   = $(A.u[end])\n"
        ))
    end
end

#-----------------------------------#
# PROPAGATION (STATE -> TRAJECTORY) #
#-----------------------------------#

const DEFAULT_ALG = Vern7();

# XXX: Required to support solving a State problem.
DiffEqBase.solve(state::State, alg=DEFAULT_ALG, args...; reltol=1e-10, abstol=1e-10, kwargs...) =
    DiffEqBase.__solve(state.system, state, alg, args...; reltol, abstol, kwargs...)

# The __solve() method does the actual heavy lifting, including converting to a Trajectory.
function DiffEqBase.__solve(state::State, alg::OrdinaryDiffEqAlgorithm; 
        userdata=Dict(), callback=nothing, kwargs...)

    default_frame = default_reference_frame(state.system)
    real_state = convert_to_frame(state, default_frame)

    # Pass the default state into the underlying solver
    # TODO: Remove the need for this in DiffCorrectAxisymmetric
    userdata = deepcopy(userdata)
    userdata[:real_state] = real_state

    # Copy the callbacks (for thread-safety)
    callback = deepcopy(callback)

    # Call the underlying solver
    raw_sol = solve(real_state.prob, alg; userdata, callback, kwargs...)

    # Save as a Trajectory
    return Trajectory(state.system, default_frame, raw_sol)
end