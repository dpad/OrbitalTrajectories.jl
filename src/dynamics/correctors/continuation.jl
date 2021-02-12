export continuation_simple

function continuation_simple(state::State; x_perturbation=1e-2, ops=(:zero, :+, :-), dc_tolerance=1e-3, kwargs...)
    # NOTE: Only designed for DiffCorrectAxisymmetric!
    corrector = DiffCorrectAxisymmetric()

    # Try to do a simple continuation from sols
    sols = Dict{Symbol,Vector{Trajectory}}(:zero => [], :+ => [], :- => [])
    orig_state = deepcopy(state)

    # Properties
    props = R3BPSystemProperties(state.model)

    # Limits
    min_x = uconvert(NoUnits, props.R1[1] / props.L) - props.μ
    max_x = (1 - props.μ) - uconvert(NoUnits, props.R2[1] / props.L)

    # Solutions to return
    new_sols::Array{Trajectory} = []

    # Threads.@threads for id in ops
    p = ProgressThresh(0., "Continuing orbits:")
    try
        for id in ops
            op = getfield(Base, id)
            p = ProgressThresh(id == :+ ? max_x : min_x, "Continuing orbits along $(op):")
            last_state = orig_state
            while true
                # Build a new solution by perturbing the last solution
                new_u0 = convert_to_frame(last_state, corrector.frame).u0
                new_u0[1] += op(x_perturbation)
                new_state = State(state.model, corrector.frame, new_u0, orig_state.tspan)

                # Check conditions
                state_check = convert_to_frame(new_state, SynodicFrame(true))
                if state_check.u0[1] >= max_x || state_check.u0[1] <= min_x
                    # Initial position as not between primary and secondary body...
                    break
                end

                try
                    new_sol = solve(new_state, corrector; dc_tolerance, kwargs...)
                    if new_sol.u[end][1] < (1 - props.μ)
                        # Solution did not cross x-axis on the right side of secondary body...
                        break
                    end

                    # Append the new solution
                    push!(sols[id], new_sol)
                    ProgressMeter.update!(p, state_check.u0[1])
                    last_state = State(new_sol)
                catch ex
                    if isa(ex, AxisymmetricException)
                        break
                    else
                        rethrow(ex)
                    end
                end

                id == :zero && break
            end
            ProgressMeter.finish!(p)
        end
    catch ex
        if isa(ex, InterruptException)
            @warn "Interrupting continuation! Found $(length(sols)) solutions."
        else
            rethrow()
        end
    end

    new_sols = vcat(values(sols)...)  # Concatenate all solutions
    sort!(new_sols, by=x -> x.prob.u0[1])  # Sort solutions by X value

    return new_sols
end