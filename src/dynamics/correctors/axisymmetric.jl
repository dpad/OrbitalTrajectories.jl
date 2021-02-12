export DiffCorrectAxisymmetric

#-------#
# TYPES #
#-------#

abstract type Abstract_AxisymmetricCorrector <: Abstract_DifferentialCorrector end

struct DiffCorrectAxisymmetric{D<:DiffEqBase.DEAlgorithm,F<:Abstract_ReferenceFrame} <: Abstract_AxisymmetricCorrector
    u0_free :: SVector{2,Integer}  # Indices of the free variables in the initial condition
    u1_free :: SVector{3,Integer}  # Indices of the free variables in the final condition
    alg :: D
    frame :: F

    # Search perturbations for axisymmetric case [Russell 2006]
    # NOTE: Only valid if the reference sol is already near-axisymmetric
    DiffCorrectAxisymmetric(u0_free, u1_free, alg, frame) = new{typeof(alg),typeof(frame)}(u0_free, u1_free, alg, frame)
    DiffCorrectAxisymmetric() = DiffCorrectAxisymmetric((5,6), (2,3,4), DEFAULT_ALG, SynodicFrame(true))
end

#-----------#
# FUNCTIONS #
#-----------#

abstract type AxisymmetricException <: Exception end
struct DidNotCross <: AxisymmetricException end
struct Crashed <: AxisymmetricException end
struct BadTspan <: AxisymmetricException end

function orbit_xcrossing(u, t, integrator)
    state = get(integrator.opts.userdata, :real_state, nothing)
    isnothing(state) && error("orbit_xcrossing expects a `real_state` to be provided in userdata.")
    new_state = remake(state; u0=u, tspan=(t, t))
    synodic_state = convert_to_frame(new_state, default_synodic_reference_frame(state.model))

    # XXX: This cannot return the Dual, or the downstream Duals (i.e. sol.u) get very broken when this triggers!
    return @inbounds ForwardDiff.value.(synodic_state.prob.u0[2])
end

function terminate_after_N_crossings!(integrator)
    t = integrator.t
    tspan = integrator.sol.prob.tspan
    if (t - tspan[begin]) / (tspan[end] - tspan[begin]) >= 0.01
        N = get!(integrator.opts.userdata, :crossings, 1)
        integrator.opts.userdata[:crossings] = N - 1
        if integrator.opts.userdata[:crossings] == 0
            terminate!(integrator, :NCrossings)
        end
    end
end

crossed(sol::Trajectory) = sol.retcode == :NCrossings

function corrector_callback(::Abstract_AxisymmetricCorrector, system::EphemerisNBP)
    # TODO: Play with interp_points, interp_points=0 halves runtime/memory, but might cause issues due to oscillation around y-axis
    # TODO: idxs=[2] makes things very slow, even though it seems it should speed it up.
    return VectorContinuousCallback(
        (integrator, idx) -> begin
            if idx == 3 # Crossing!
                terminate_after_N_crossings!(integrator)
            else
                terminate!(integrator, :Crashed)
            end
        end, 3;
        interp_points=10, save_positions=(false, false)) do du, u, t, integrator
            diam1 = maximum(bodvrd(String(primary_body(system)), "RADII"))  # km
            diam2 = maximum(bodvrd(String(secondary_body(system)), "RADII"))  # km
            du[1] = check_distance(u, t, system, primary_body(system), diam1)
            du[2] = check_distance(u, t, system, secondary_body(system), diam2)
            du[3] = orbit_xcrossing(u, t, integrator)

            # XXX: I think that some of the underlying SPICE functions are slightly discontinuous,
            # maybe there are differences in number of bits being used between here and in the
            # SPICE library. Basically, this causes issues in event callback detection, because it
            # gets stuck in a value close, but not quite equal to zero (e.g. on a value like -2.0e-14),
            # despite incrementing the timestep -- see https://github.com/SciML/DiffEqBase.jl/issues/646.
            # This fix is very hacky, but I have no idea what might be a better option for now.
            @. du[isapprox(du, 0; atol=1e-4)] = 0.
        end
end

function corrector_callback(::Abstract_AxisymmetricCorrector, ::Abstract_DynamicalModel)
    # TODO: Play with interp_points, interp_points=0 halves runtime/memory, but might cause issues due to oscillation around y-axis
    # TODO: idxs=[2] makes things very slow, even though it seems it should speed it up.
    return ContinuousCallback(orbit_xcrossing, terminate_after_N_crossings!; interp_points=10, save_positions=(false, false))
end

function corrector_solve(corrector::Abstract_AxisymmetricCorrector, state::State; strict=false, verbose=false, crossings=1, kwargs...)
    callback = strict ? corrector_callback(corrector, state.model) : nothing
    userdata = Dict{Symbol,Any}(:crossings => crossings)  # XXX: Dict type specified for easy merging

    sol = sensitivity_trace(AD, state, corrector.frame, corrector.alg; trace_time=true, abstol=1e-12, reltol=1e-12, callback, userdata, kwargs...)
    sol_tspan = ForwardDiff.value.((sol.t[begin], sol.t[end]))

    # Error-checking
    if strict
        try
            crashed(sol) && throw(Crashed())
            !crossed(sol) && throw(DidNotCross())
        catch ex
            if verbose && isa(ex, AxisymmetricException)
                display(plot(solve(remake(state, tspan=sol_tspan)), corrector.frame; title="$(ex)"))
            end
            rethrow(ex)
        end
    end

    if verbose
        display(plot(sol, corrector.frame))
    end

    return convert_to_frame(sol, corrector.frame), sol_tspan
end


function (corrector::Abstract_AxisymmetricCorrector)(state::State, F, J, x; kwargs...)
    tspan = (state.tspan[begin], state.tspan[end] + x[end])
    if (tspan[end] - tspan[begin]) < 0 || (tspan[end] - tspan[begin]) > 1e8
        # Don't like this tspan, it's either negative time direction, or too long
        throw(BadTspan())
    end

    # Recreate the state with the parameters provided in x
    correction_state = convert_to_frame(state, corrector.frame)
    u0 = deepcopy(correction_state.u0)
    u0[corrector.u0_free] .= x[1:length(corrector.u0_free)]
    correction_state = remake(correction_state; u0, tspan)

    # Compute sensitivities!
    sol, sol_tspan = corrector_solve(corrector, correction_state; strict=false, kwargs...)

    if !isnothing(F)
        # Calculate how far away we are from the desirable final state (i.e. the residual).
        end_u0 = ForwardDiff.value.(sol.u[end])
        F .= [end_u0[corrector.u1_free]...]
    end

    if !isnothing(J)
        # As per [Russell.2006, Eq.(10) to (13)], we need to add the time variation term, 
        # which requires the values of [ẏ, ż, ẍ]. Here, these are provided in the STM
        # on the "end" column (since we're calling corrector_solve with trace_time=true).
        STM = extract_STMs([sol.u[end]])[1]
        J .= STM[[corrector.u1_free...], [corrector.u0_free..., end]]
    end
end