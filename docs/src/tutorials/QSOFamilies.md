# [Families of Quasi-Satellite Orbits (QSOs)](@id QSOs)

This tutorial will briefly introduce computing parametrised families of Quasi-Satellite
Orbits (QSOs) in any of the astrodynamical models provided by OrbitalTrajectories.jl.

## Load packages

First, we load all necessary packages:

```@example 1
using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using LinearAlgebra          # To compute stability index (eigenvalue)
```

## Differential Correction

We can correct an initial guess trajectory to find an axisymmetric orbit.[^Russell2006] We begin
with an initial guess in a CR3BP (Mars-Phobos) system:

[^Russell2006]: R. P. Russell, “Global search for planar and three-dimensional periodic orbits near Europa,” The Journal of the Astronautical Sciences, Vol. 54, No. 2, 2006, pp. 199–226, 10.1007/BF03256483. 

```@example 1
# Initial guess state
u0 = [0.994, 0., 0., 0., 0.0117, 0.]
state = State(CR3BP(:mars, :phobos), u0, (0., 2π))
trajectory = solve(state)
p = plot(trajectory; color=:blue, alpha=0.25)
```

We can correct this initial guess using a `DiffCorrectAxisymmetric` corrector:

```@example 1
trajectory = solve(state, DiffCorrectAxisymmetric())
plot!(p, trajectory; color=:blue)
```

## Parameter Continuation for DRO/QSO Orbit Families

Simple parameter continuation works by beginning with a known orbit, then perturbing it very slightly and correcting it to find
a similar orbit. Here, we use `continuation_simple` to perturb our orbit along its initial ``x``-axis position, returning a set
of many similar orbits along the axis.

!!! warning "Current limitations of `continuation_simple`"
    Note that the `continuation_simple` method is not currently designed to be robust. It will only look for orbits that satisfy
    several conditions, such as ``x_0 \in (R_1 - \mu, 1 - \mu)`` (given the primary body's normalised ``x`` radius of ``R_1``) and
    ``x_f \in (1 - \mu, \infty)`` (for the final ``x``-crossing position). Thus, it will essentially only look for and find
    orbits like DROs/QSOs, perturbing only along the initial ``x`` position.

```@example 1
trajectories = continuation_simple(state; x_perturbation=3e-4, dc_tolerance=1e-4)
length(trajectories)
```

The example above returned a large number of orbits, so to make plotting faster, we only plot every 10th orbit using the `1:10:end`
selector. In addition, we also provide a `values` array to provide a colour to each individual orbit depending on its stability
index.

```@example 1
selected_trajectories = trajectories[1:10:end]

# Compute stability index for each orbit
stability_index = map(selected_trajectories) do traj
    STM = Matrix(sensitivity(AD, State(traj)))
    eigenvalues = norm.(eigvals(STM))
    return maximum(eigenvalues)
end

# Propagate each trajectory to a full orbit
full_orbits = [
    solve(remake(State(traj); tspan=(traj.t[begin], traj.t[end] * 2.)))
    for traj in selected_trajectories
]

# Plot the orbit family, coloured by its stability index
plot(full_orbits; values=stability_index, primary_color=:red, legend=:bottomleft)
```