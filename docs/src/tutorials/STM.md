# [State Transition Matrices (STMs)](@id STMs)

This tutorial will briefly introduce computing and plotting State Transition Matrices (STMs) for
any of the astrodynamical models provided by OrbitalTrajectories.jl.

Here, we will reproduce Figure 4 from the OrbitalTrajectories.jl paper[^Padilha2021], which shows
example trajectories and traces out the maximum eigenvalue of their STMs.

[^Padilha2021]: Dan Padilha, Diogene A. Dei Tos, Nicola Baresi, Junichiro Kawaguchi, "Modern Numerical Programming with Julia for Astrodynamic Trajectory Design", 31st AAS/AIAA Space Flight Mechanics Meeting, 2021.


## Load packages

First, we load all necessary packages:

```@example 1
using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using LinearAlgebra          # To compute matrix eigenvalues
```


## Orbit data

We build a structure to store the orbit data necessary. This data comes from Table 1 of this paper.[^Pellegrini2016]
Note that we include the Sun as a (fourth) body for higher-order propagation in `BC4BP` and `EphemerisNBP`.

[^Pellegrini2016]: Etienne Pellegrini, Ryan P. Russell, "On the Computation and Accuracy of Trajectory State Transition Matrices", Journal of Guidance, Control, and Dynamics, Vol. 39 Issue 11, Pg. 2485-2499, 2016, DOI:10.2514/1.G001920.

```@example 1
test_cases = [(
    system = (:Jupiter, :Europa, :Sun),
    u0     = [-0.1017472008677258, 0., 0., 0., -0.01806028472285857, 0.],
    μ      = 2.528009215182033e-5,
    t_p    = 25.13898226959327,
    λ_max  = 2.468621114047195,
    plot_args = (label="7:4 resonant\n(Jupiter-Europa)", xlim=(-1.15, 1.07))
), (
    system = (:Jupiter, :Europa, :Sun),
    u0     = [0.04867586089512202, 0., 0., 0., -0.09354853663949217, 0.],
    μ      = 2.528009215182033e-5,
    t_p    = 70.53945041512506,
    λ_max  = 2.804814246519340e7,
    plot_args = (label="8:11 + flybys\n(Jupiter-Europa)", xlim=(-1.37, 1.35))
), (
    system = (:Earth, :Moon, :Sun),
    u0     = [-0.013059020050628, 0., 0.07129515195874, 0., -0.526306975588415, 0.],
    μ      = 0.01215509906405700,
    t_p    = 2.517727406553485,
    λ_max  = 17.632688124231755,
    plot_args = (label="L1 halo orbit\n(Earth-Moon)", xlim=(0.8, 1.015))
)]
nothing # hide
```

## Propagate & plot trajectory

For demonstration, we first propagate and plot one of the test cases above in the **Circular Restricted 3-Body Problem** (`CR3BP`).

```@example 1
case = test_cases[1]  # Pick the first test case

# Build the system
# NOTE: overriding the default mass parameter μ with the one provided by [Pellegrini2016]
system = CR3BP(case.system[1:2]...; μ=case.μ)

# Propagation timespan
t_init = π/4  # Arbitrary initial time, does not matter for CR3BP
tspan  = (t_init, t_init + case.t_p)

# Initial state vector in [Pellegrini2016] is relative to secondary body
u0 = copy(case.u0)
u0[1] += 1 - case.μ

# Build the initial state and propagate it
state = State(system, SynodicFrame(), u0, tspan)
trajectory = solve(state)
plot(trajectory; legend=:bottomleft, case.plot_args...)
```

## Compute State Transition Matrix (STM)

We can compute the STM using one of several methods:

* **Automatic Differentiation** (`AD`) -- using the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package.
* **Finite Differencing** (`FD`) -- using the [FiniteDiff](https://github.com/JuliaDiff/FiniteDiff.jl) package.
* **Symbolic Variational Equations** (`VE`) -- using the [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) package.
* **Hand-coded Variational Equations** (`VE`) -- when using the `HandcodedCR3BP` model.

Each of the above methods is implemented in OrbitalTrajectories.jl. Below, we demonstrate
how to use the Automatic Differentiation (AD) method to compute and plot the STM:

```@example 1
# Compute the STM at the final point of the trajectory (with respect to initial state)
STM = sensitivity(AD, state)
```

We can check that the stability index (maximum eigenvalue) matches the value from the paper.[^Pellegrini2016]

```@example 1
# Check if STM eigenvalues matches what [Pellegrini2016] gives
λ_max = maximum(norm.(eigvals(Matrix(STM))))
@assert isapprox(λ_max, case.λ_max; rtol=1e-5)
```

We can also compute and plot the STM along the full trajectory. Note that this trace is also
interpolatable at any point.

```@example 1
# Compute the STM trace (i.e. the STM at every point along the trajectory)
STM_trace = sensitivity_trace(AD, state)

# Plot the STM elements
plot(STM_trace; trace=true)
```

We can also plot the stability index directly:

```@example 1
# Plot the STM stability index (i.e. maximum eigenvalue)
plot(STM_trace; trace_stability=true)
```

## Plot all the trajectories and STM stability traces

The models we will test are the following:

```@example 1
models = (EphemerisNBP, BC4BP, ER3BP, CR3BP, HandcodedCR3BP)
nothing # hide
```

The test cases are plotted and traced for each of the above models as follows.

```@example 1
using Unitful                # Units (u"km", u"d" (days)) and unit conversion
using LaTeXStrings           # So we can use Latex formatting in strings
using Plots.PlotMeasures     # To change plot margins using "mm" units

frame = SynodicFrame()
plot_attrs = (legendfontsize=8, titlefontsize=10, bg_color_legend=RGBA(1, 1, 1, 0.5), fg_color_legend=nothing)

p_orbits = []
p_traces = []

for (j, case) in enumerate(test_cases)
    push!(p_orbits, plot())
    push!(p_traces, plot())

    for (i, model) in enumerate(models)
        if model == HandcodedCR3BP || model <: Union{CR3BP, ER3BP}
            bodies = case.system[1:2]
            system = model(bodies...; μ=case.μ)
        else
            system = model(case.system...)
        end

        # Plotting options for this model
        sys_name = String(nameof(typeof(system))) * (isa(system, CR3BP{true}) ? " (hand-coded)" : "")
        alpha = isa(system, CR3BP) ? 1. : 0.7
        VE_label = " only"

        # Set an appropriate timespan
        t_init = π/4  # arbitrary initial epoch
        tspan = (t_init, t_init + case.t_p)
        if isa(system, EphemerisNBP)
            # tspan[1] is effectively an arbitrarily-chosen initial epoch time
            props = R3BPSystemProperties(primary_body(system), secondary_body(system))
            tspan = (tspan[1], case.t_p * ustrip(u"s", props.T) + tspan[1])
            alpha = 1.0
        end

        # Convert initial state relative to secondary body --> relative to primary body
        u0 = copy(case.u0)
        u0[1] += 1 - case.μ

        # Generate initial state and trajectory
        state = State(system, frame, u0, tspan)
        traj = solve(state)

        # Plot the trajectory
        color = [:orange, :green, :red, :blue, :black][i]
        !isa(system, CR3BP{true}) && plot!(p_orbits[end], traj, frame; arrow=false, title=case.plot_args.label,
            linewidth=isa(system, CR3BP{false}) ? 1.5 : 1, color, label="", alpha, nolabels=true,
            origin_secondary=false, case.plot_args..., plot_attrs...)

        # Build trace using Variational Equations (VE)
        if has_variational_equations(typeof(state))
            traj_VE = sensitivity_trace(VE, state, frame)
            VE_label = " & VE"

            # Extract the sensitivities
            tspan = range(traj_VE.t[begin], traj_VE.t[end], length=1000)
            tnorm = range(0., 1., length=length(tspan))
            u_vals = traj_VE.sol.(tspan)
            dim = length(state.u0)
            eigenvalues = eigvals.([reshape(ut[dim+1:end], (dim, dim)) for ut in u_vals])
            stability_indices = map(x -> maximum(norm.(x)), eigenvalues)

            # Plot!
            !isa(system, CR3BP{true}) && plot!(p_traces[end], tnorm, stability_indices; linestyle=:dash, color, label="", alpha)
        end
        VE_label = isa(system, CR3BP{false}) ? "$(VE_label) & hand-coded" : VE_label

        # Build trace using Automatic Differentiation (AD)
        traj_AD = sensitivity_trace(AD, state, frame)
        !isa(system, CR3BP{true}) && plot!(p_traces[end], traj_AD, frame; trace_stability=true, label="$(sys_name) (AD$(VE_label))", color,
            link=:y, yscale=:log10, ylabel=(j == 1) ? L"$\lambda_{\max}(\textrm{STM})$" : "", ytickfontcolor=(j == 1) ? RGBA(0,0,0,1) : RGBA(0,0,0,0),
            ytickfontsize=(j == 1) ? 8 : 0, ygrid=true, legend=(j==1) ? :topleft : false, alpha, plot_attrs...)
        j == 2 && plot!(p_traces[end]; xlabel="Time (normalised to 1 period)")
        
        # Build trace using Finite Differencing (FD)
        STM_FD = Matrix(sensitivity(FD, state, frame))
        λ_max_FD = maximum(norm.(eigvals(STM_FD)))
        isa(system, CR3BP{false}) && scatter!(p_traces[end], [1.], [λ_max_FD]; color=:black, markersize=6, markerstrokewidth=0, markershape=:star5, 
            label=L"$\lambda_{\max,\textrm{end}}\;\textrm{(FD)}$")
        !isa(system, CR3BP{true}) && scatter!(p_traces[end], [1.], [λ_max_FD]; color, markersize=6, markerstrokewidth=0, markershape=:star5, label="")
    end
end

# Plot: trajectories and STM traces
plot(p_orbits..., p_traces...; layout=(2,3), thickness_scaling=0.75, left_margin=2mm, right_margin=4mm)
```