# [Astrodynamical Models](@id dynamical_models)

This tutorial will briefly introduce propagating and plotting orbits in the astrodynamical models
provided by OrbitalTrajectories.jl

Here, we will reproduce Figure 3 from the OrbitalTrajectories.jl paper[^Padilha2021], which shows
example trajectories in the Earth-Moon(-Sun) system.

[^Padilha2021]: Dan Padilha, Diogene A. Dei Tos, Nicola Baresi, Junichiro Kawaguchi, "Modern Numerical Programming with Julia for Astrodynamic Trajectory Design", 31st AAS/AIAA Space Flight Mechanics Meeting, 2021.


## Load packages

First, we load all necessary packages:

```@example 1
using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using Unitful                # Units (u"km", u"d" (days)) and unit conversion
```


## Build astrodynamical models

We build a set of astrodynamical models (`systems`) in the Earth-Moon(-Sun) system. This includes:

* **Circular Restricted 3-Body Problem** (`CR3BP`) for Earth-Moon
* **Elliptic Restricted 3-Body Problem** (`ER3BP`) for Earth-Moon
* **Bi-Circular Restricted 4-Body Problem** (`BC4BP`) for Earth-Moon-Sun
* **Ephemeris Restricted N-Body Problem** (`EphemerisNBP`) for Earth-Moon-Sun

```@example 1
a, b, c = (:Earth, :Moon, :Sun)

# The EphemerisNBP model requires a timespan in seconds, rather than radians.
tspan       = (0.88π, 2.83π)
tspan_epoch = ustrip.(u"s", tspan .* R3BPSystemProperties(a, b).T)

# Define the systems to test.
# xoffset/yoffset are used to offset the model labels in the resulting plot.
systems = Dict(
  CR3BP(a, b)           => (color=:blue,   tspan=(0.88π, 2.88π), xoffset=-0.06, yoffset=-0.04),
  ER3BP(a, b)           => (color=:red,    tspan=(0.88π, 3.00π), xoffset=-0.06, yoffset=-0.04),
  BC4BP(a, b, c)        => (color=:green,  tspan=(0.88π, 2.85π), xoffset=-0.05, yoffset=-0.06),
  EphemerisNBP(a, b, c) => (color=:orange, tspan=tspan_epoch,    xoffset=-0.20, yoffset=-0.02)
)
nothing # hide
```

## Initial state and reference frame

We specify an initial state simply as a vector of ``[x_0, y_0, z_0, \dot{x}_0, \dot{y}_0, \dot{z}_0]``.
The state we provide is defined in a normalised synodic rotating frame, so we specify the initial frame as
`SynodicFrame()`.

```@example 1
u0       = [0.76710535, 0., 0., 0., 0.47262724, 0.]
u0_frame = SynodicFrame()
nothing # hide
```

## Simple propagation and plotting

We can propagate a trajectory by creating an initial state object `State`, providing the desired
astrodynamical model (`sys`), initial state vector (`u0`), initial state reference frame (`u0_frame`),
and propagation timespan (`tspan`).

To propagate, we simply call `solve(state)`.

For example, to plot in the **Circular Restricted 3-Body Problem** (`CR3BP`):

```@example 1
# Create the initial state
state = State(CR3BP(:Earth, :Moon), u0_frame, u0, tspan)

# Propagate the state to compute a trajectory
trajectory = solve(state)

# Plot the trajectory
plot(trajectory; legend=:bottomright)
```

## Propagate and plot in all the models

Here we propagate and plot trajectories in all of the aforementioned models.

First, set up some plotting options as follows.

```@example 1
plot_attrs = (thickness_scaling=0.75, legendfontsize=9, tickfontsize=9, guidefontsize=10,
              titlefontsize=10, bg_color_legend=RGBA(1, 1, 1, 0.8), fg_color_legend=nothing,
              legend=:bottomright)
p = plot();

# Keep track of the maximum x and y limits
xlim = (Inf, -Inf)
ylim = (Inf, -Inf)
nothing # hide
```

To plot in the original reference frame (i.e. normalised synodic rotating
frame), we simply call `plot(trajectory, u0_frame)`.

```@example 1
for (i, (sys, opts)) in enumerate(systems)
    # Create the initial state
    state = State(sys, u0_frame, u0, opts.tspan)

    # Propagate it to compute the trajectory
    trajectory = solve(state)

    # Plot the trajectory
    plot!(p, trajectory, u0_frame; label="", color=opts.color, nolabels=(i!=1))

    # Add a label for the model name
    traj_converted = convert_to_frame(trajectory, u0_frame)
    annotate!(p, [(traj_converted.sol[end][1] + opts.xoffset,
                   traj_converted.sol[end][2] - opts.yoffset,
                   Plots.text("$(nameof(typeof(sys)))" * 
                    ((isa(sys, EphemerisNBP) || isa(sys, BC4BP)) ? "\n(+ $(String(c)))" : ""),
                   opts.color, :left, 9, 0.))]);

    # Add a label for barycenter
    (i == 1) && hline!(p, [0.]; linecolor=:black, linestyle=:solid,
                       label="$(titlecase(String(a)))-$(titlecase(String(b))) Barycenter")
    
    # Expand the plot limits to the biggest size
    global xlim, ylim  # we want to override the global variables outside this scope
    xlim = (min(xlim[1], xlims()[1]), max(xlim[2], xlims()[2]))
    ylim = (min(ylim[1], ylims()[1]), max(ylim[2], ylims()[2]))
    plot!(p; xlims=xlim, ylims=ylim, plot_attrs...)
end
p  # Display the final plot
```