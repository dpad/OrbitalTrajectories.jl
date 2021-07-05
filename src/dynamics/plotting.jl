# Define a circle by radius r
circle(x, y, r) = ellipse_by_eccentricity(x, y, a = r, e = 0)

# Define an ellipse with: semi-major axis a, eccentricity e
ellipse_by_eccentricity(x, y; a, e) = ellipse_by_axis(x, y, a = a, b = √(a^2 * (1 - e^2)))

function ellipse_by_axis(x, y; a, b, num=500)
    # Define an ellipse with:
    #   semi-major axis a
    #   semi-minor axis b
    θ = LinRange(0, 2π, num)
    return @. (x + a * cos(θ), y + b * sin(θ))
end

@recipe function f(sols::AbstractArray{<:Trajectory})
    sols, sols[1].frame
end

@recipe function f(sols::AbstractArray{<:Trajectory}, frame::Abstract_ReferenceFrame)
    values = get(plotattributes, :values, nothing)
    plot_cbar = !isnothing(values) && length(values) == length(sols)

    user_xlim = pop!(plotattributes, :xlims, nothing)
    user_ylim = pop!(plotattributes, :ylims, nothing)

    legend --> false # XXX: Legend seems to mess up the plot badly!

    # Work out the maximum extent of the orbit
    xlim = (Inf, -Inf)
    ylim = (Inf, -Inf)

    for (i, sol) in enumerate(sols)
        current_xlim, current_ylim = get_margin_lims(convert_to_frame(sol, frame), plotattributes)
        if isnothing(user_xlim)
            xlim = (min(xlim[1], current_xlim[1]), max(xlim[2], current_xlim[2]))
        else
            xlim = user_xlim
        end
        if isnothing(user_ylim)
            ylim = (min(ylim[1], current_ylim[1]), max(ylim[2], current_ylim[2]))
        else
            ylim = user_ylim
        end

        @series begin
            label := false
            if plot_cbar
                seriescolor --> cgrad(:thermal)
                arrow --> false
                linewidth --> 3
                label --> false
                colorbar --> :left
                line_z --> values[i]
            else
                seriescolor --> :blue
                arrow --> true
                linewidth --> 1.5
                label --> "Trajectory"
            end
            xlims := xlim
            ylims := ylim
            if i < length(sols)
                nomodel := true
            end
            sol, frame
        end
    end
end

# ------------------
# NEW STYLE PLOTTING
# ------------------

@recipe function f(sol::Trajectory)
    trace_vars = get(plotattributes, :trace, get(plotattributes, :trace_stability, false))
    if trace_vars === false
        # Plot the frame
        if !get(plotattributes, :nomodel, false)
            @series begin
                seriesalpha := 1.0
                (sol.model, sol.frame)
            end
            framestyle --> :zerolines
        else
            framestyle --> :none
        end

        # Plot the trajectory
        vars --> (1, 2)  # (x, y)
        denseplot --> get(plotattributes, :denseplot, true)

        xlim, ylim = get_margin_lims(sol, plotattributes)
        xlims --> xlim
        ylims --> ylim
        
        # Formatting
        arrow --> true
        dpi --> 150
        size --> (400, 500)
        legend --> false
        xaxis --> (rotation=45)
        aspect_ratio --> 1
        ticks --> false
    end

    sol.sol
end

@recipe function f(sol::Trajectory{<:Abstract_DynamicalModel,<:Abstract_ReferenceFrame,T}) where {T <: ForwardDiff.Dual}
    trace_vars = get(plotattributes, :trace, false)
    trace_stability = get(plotattributes, :trace_stability, false)
    denseplot = get(plotattributes, :denseplot, true)
    plotdensity = get(plotattributes, :plotdensity, 1000)

    tspan = ForwardDiff.value.(denseplot ? range(sol.t[begin], sol.t[end], length=plotdensity) : sol.t)
    tspan_norm = @. (tspan - tspan[begin]) / (tspan[end] - tspan[begin])
    sol_interpolated = sol(tspan)

    if trace_vars
        STMs = get_sensitivity(sol_interpolated)
        @series begin
            label --> ""
            legend --> false
            tspan_norm, STMs
        end
    elseif trace_stability
        stability_indices = norm.(stability_index(sol_interpolated))'
        sorted_indices = hcat(map(sort, eachslice(stability_indices, dims=1))...)[4:6,:]'
        max_eigenvalues = map(maximum, eachslice(sorted_indices, dims=1))
        @series begin
            label --> ""
            legend --> false
            tspan_norm, max_eigenvalues
        end
    else
        vars = get(plotattributes, :vars, (1,2))
        ([[u[v].value for u in sol_interpolated.u] for v in vars]...,)
    end
end

@recipe function f(model::Abstract_DynamicalModel, frame::Abstract_ReferenceFrame)
    nothing
end

@recipe function f(model::Abstract_DynamicalModel, frame::SynodicFrame)
    # User arguments
    nolabels = get(plotattributes, :nolabels, get(plotattributes, :nomodel, false))
    plot_libration = get(plotattributes, :libration_points, true)
    circ_props = R3BPSystemProperties(primary_body(model), secondary_body(model))

    vars = get(plotattributes, :vars, (1,2))

    @series begin
        seriestype := :shape
        seriescolor := get(plotattributes, :primary_color, :blue)
        linecolor := :black
        linewidth := 0.5
        label := nolabels ? "" : titlecase(String(primary_body(model)))
        line_z := nothing

        primary_pos = (-circ_props.μ, 0., 0.)
        ellipse_by_axis(primary_pos[vars[1]], primary_pos[vars[2]];
                        a = circ_props.R1[vars[1]] / circ_props.L,
                        b = circ_props.R1[vars[2]] / circ_props.L)
    end

    @series begin
        seriestype := :shape
        seriescolor := get(plotattributes, :secondary_color, :brown)
        linecolor := :black
        linewidth := 0.5
        label := nolabels ? "" : titlecase(String(secondary_body(model)))
        line_z := nothing

        secondary_pos = (1 - circ_props.μ, 0., 0.)
        ellipse_by_axis(secondary_pos[vars[1]], secondary_pos[vars[2]];
                        a = circ_props.R2[vars[1]] / circ_props.L,
                        b = circ_props.R2[vars[2]] / circ_props.L)
    end

    if get(plotattributes, :origin_primary, true)
        @series begin
            seriestype := :hline
            seriescolor := :black
            linestyle := :dash
            linewidth := 0.75
            line_z := nothing
            label := ""
            [0]
        end
    end

    if get(plotattributes, :origin_secondary, true)
        @series begin
            seriestype := :vline
            seriescolor := :black
            linestyle := :dash
            linewidth := 0.75
            label := nolabels ? "" : "Origin of $(titlecase(String(secondary_body(model))))"
            line_z := nothing
            if vars[1] == 1
                [1 - circ_props.μ]
            else
                [0.]
            end
        end
    end

    @series begin
        # Find the libration points
        L = libration_points(circ_props)
        if !isnothing(L) && plot_libration
            seriestype := :scatter
            markercolor := :yellow
            markerstrokewidth := 1
            markerstrokecolor := :red
            markershape := :diamond
            markersize := 3
            label := nolabels ? "" : "Libration points"
            line_z := nothing
            [l[vars[1]] for l in L], [l[vars[2]] for l in L]
        else
            [], []
        end
    end

end

@recipe function f(sol::Trajectory, frame::Abstract_ReferenceFrame)
    convert_to_frame(sol, frame)
end

@recipe function f(state::State)
    seriestype := :scatter
    ([state.prob.u0[1]], [state.prob.u0[2]])
end

@recipe function f(state::State, frame::Abstract_ReferenceFrame)
    @series begin
        state.model, frame
    end
    convert_to_frame(state, frame)
end

function get_margin_lims(sol::Trajectory, plotattributes)
    margins = get(plotattributes, :padding, 0.10)
    a, b = get(plotattributes, :vars, (1, 2))

    # Work out the maximum extent of the orbit
    x, y = (ForwardDiff.value.(sol.sol[a,:]), ForwardDiff.value.(sol.sol[b,:]))
    xlim = (minimum(x), maximum(x))
    ylim = (minimum(y), maximum(y))

    # Add some margins to the figure
    xdiff = margins * (xlim[2] - xlim[1])
    ydiff = margins * (ylim[2] - ylim[1])

    xlims = (xlim[1] - xdiff, xlim[2] + xdiff)
    ylims = (ylim[1] - ydiff, ylim[2] + ydiff)

    return xlims, ylims
end