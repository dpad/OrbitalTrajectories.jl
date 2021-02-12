#-----------------------------------#
# DIFFERENTIAL CORRECTION FUNCTIONS #
#-----------------------------------#

function corrector_solve end
function corrector_callback end

function DiffEqBase.__solve(state::State, corrector::Abstract_DifferentialCorrector;
                            dc_maxiters=10, dc_tolerance=1e-6, dc_xtol=dc_tolerance^2,
                            dc_method=:newton, dc_linesearch=LineSearches.BackTracking(),
                            verbose=false, kwargs...)

    # Run once with strict=true to retrieve the initial guess.
    correction_state = convert_to_frame(state, corrector.frame)
    guess_sol, guess_tspan = corrector_solve(corrector, correction_state; strict=true, verbose, kwargs...)

    # The state we will correct begins with a tspan up to the x-crossing
    correction_state = remake(correction_state; tspan=guess_tspan)
    u0 = deepcopy(correction_state.u0)

    # Solve iteratively to find a corrected.
    sol = nlsolve(
        # Function for computing residuals (F) and Jacobian (J)
        only_fj!((F, J, x) -> corrector(correction_state, F, J, x; verbose, kwargs...)),  
        # Initial guess (x) -- includes the initial guess tspan's final time
        [u0[corrector.u0_free]..., 0.];
        # Options
        method=dc_method, linesearch=dc_linesearch, ftol=dc_tolerance, xtol=dc_xtol,
        show_trace=verbose, extended_trace=verbose, iterations=dc_maxiters
    )
    verbose && @info sol

    # Propagate out the last solution.
    u0[corrector.u0_free] .= sol.zero[1:length(corrector.u0_free)]
    tspan = (guess_tspan[begin], guess_tspan[end] + sol.zero[end])
    final_state = remake(correction_state; u0, tspan)
    traj = solve(final_state, corrector.alg; callback=corrector_callback(corrector, final_state.model), kwargs...)

    return traj
end

#----------#
# INCLUDES #
#----------#

include("correctors/axisymmetric.jl")
include("correctors/continuation.jl")