export InertialFrame, SynodicFrame, convert_to_frame

#------------------#
# REFERENCE FRAMES #
#------------------#

struct InertialFrame <: Abstract_ReferenceFrame end
struct SynodicFrame{IsNormalised} <: Abstract_ReferenceFrame
    SynodicFrame(isnormalised::Bool=true) = new{isnormalised}()
end

default_reference_frame(model::Abstract_AstrodynamicalModel) = error("$(model) has no default reference frame assigned!")
default_synodic_reference_frame(::Abstract_AstrodynamicalModel) = SynodicFrame()

#----------------------------------#
# CONVERT BETWEEN REFERENCE FRAMES #
#----------------------------------#
convert_to_frame(state::State, frame::Abstract_ReferenceFrame) = error("Cannot convert $(nameof(typeof(state))){$(state.model),$(state.frame)} to $(frame).")
convert_to_frame(state::State{<:Abstract_AstrodynamicalModel,T}, ::T) where {T<:Abstract_ReferenceFrame} = state
convert_to_frame(traj::Trajectory{<:Abstract_AstrodynamicalModel,T}, ::T) where {T<:Abstract_ReferenceFrame} = traj
function convert_to_frame(traj::Trajectory, frame::Abstract_ReferenceFrame)
    prob0 = traj.prob
    times = traj.t

    u1 = deepcopy(traj.u)
    convert_u!(u1, times, traj, frame)

    # TODO: Convert this into a type so that we can make an interp_summary(::interpType) dispatch (to print out converted sols)
    function interp_and_convert(times,idxs,deriv::Type{Val{0}},p,continuity::Symbol=:left)
        u = deepcopy(traj.interp(times, idxs, deriv, p, continuity))
        convert_u!(u, times, traj, frame)
        return u
    end

    new_sol = DiffEqBase.build_solution(prob0, traj.alg, times, u1; interp=interp_and_convert, retcode=traj.retcode)

    return Trajectory(traj.model, frame, new_sol)
end
function convert_u!(u::AbstractArray, times::AbstractArray, traj, frame)
    for (i, t) in enumerate(times)
        convert_u!(u[i], t, traj, frame)
    end
end
function convert_u!(u, t, traj, frame)
    prob1 = remake(traj.sol.prob; u0=u, tspan=(t, t))
    new_state = convert_to_frame(State(traj.model, traj.frame, prob1), frame)
    u .= new_state.prob.u0
end

const CRASHED_RETCODE = :Crashed
function collision(system::Abstract_AstrodynamicalModel, body, dist=bodvrd(String(body), "RADII")[1]; radii=1., interp_points=10)
    diam = radii * dist
    ContinuousCallback((integrator) -> terminate!(integrator, CRASHED_RETCODE); interp_points) do u, t, integrator
        check_distance(u, t, system, body, diam)
    end
end
check_distance(u, t, system::Abstract_AstrodynamicalModel, body) = error("check_distance not defined for $(nameof(typeof(system)))")
crashed(sol::Trajectory) = crashed(sol.sol)
crashed(sol::ODESolution) = sol.retcode == CRASHED_RETCODE