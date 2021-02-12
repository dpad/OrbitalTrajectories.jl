# [Sun-Mars CR3BP periodic orbits](@id SunMars_cr3bp_orbits)

This tutorial will briefly introduce propagating and plotting orbits in the Circular Restricted 3-Body Problem (CR3BP).

Here, we will reproduce Figure 4 from a Mars ballistic capture paper[^DeiTos2018], which plots a few different types
of periodic orbits in the Sun-Mars 3-body system.

[^DeiTos2018]: Diogene A. Dei Tos, Ryan P. Russell, and Francesco Topputo, "Survey of Mars Ballistic Capture Trajectories Using Periodic Orbits as Generating Mechanisms", Journal of Guidance, Control, and Dynamics 2018 41:6, 1227-1242

!!! note "Coding style"
    Note that we've split this example into individual steps, and purposefully introduced the
    use of Julia `Dictionaries`, `Arrays`, `NamedTuples`, and so on. In actual code,
    you should be able to simplify this code as desired.

## Load packages

First, we load all necessary packages:

```@example 1
using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using Unitful                # Units (u"km", u"d" (days)) and unit conversion
```

## Compute initial ``y``-velocity given a Jacobi energy

The orbit parameters given in the paper are the initial ``x``-axis crossing position (in ``10^6``km from Mars)
and Jacobi energy ``J`` (non-dimensional).

First, we compute the initial state vector by normalising the initial ``x`` position (``x0``), then
computing the initial ``y``-velocity (``v0``) from the Jacobi energy.

We return a new spacecraft `State`(@ref) with the initial state (given ``x0`` and ``v0``) and expected
propagation time (which we set to some default ``tmax`` value, noting that it is normalised to radians
for the CR3BP model).

```@example 1
function compute_initial_state(system; x0, J, tmax=365u"d")
    x0 = 1 - system.props.μ + (x0*(10^6)u"km" / system.props.L)
    v0 = √(2*centrifugal_potential(system.props.μ, [x0, 0., 0.]) - J)
    State(system, [x0, 0., 0., 0., v0, 0.], (0, tmax / system.props.T))
end
```

## Create initial states

Now we create the initial states by mapping our function above to the set of orbit parameters given in the paper. Note that our initial states are guesses; they are not completely numerically accurate, and therefore will not lead to periodic orbits.

```@example 1
SunMars = CR3BP(:sun, :mars)

all_params = Dict(
    "g2" => [
        (x0 = 0.172, J = 3.000034),  # updated to get a similar trajectory
        (x0 = 0.432, J = 3.000118),
        (x0 = 0.556, J = 3.000203),
        (x0 = 0.775, J = 3.000203),
        (x0 = 0.841, J = 3.0001804), # highly sensitive, updated to similar
    ],
    "g1" => [
        (x0 = 0.006, J = 3.011712, tmax=1u"d"),
        (x0 = 0.198, J = 3.000410),
        (x0 = 0.161, J = 3.000203),
        (x0 = 0.081, J = 3.000203),
        (x0 = 0.027, J = 3.000196), # highly sensitive, updated to similar
    ],
    "DRO" => [
        (x0 = -0.008, J = 3.009513, tmax=1u"d"),
        (x0 = -0.256, J = 3.000255),
        (x0 = -0.944, J = 3.000015),
    ]
)

# Map the compute_initial_state function to all the parameters in each family above
all_states = Dict(family => begin
    [compute_initial_state(SunMars; p...) for p in params]
end for (family, params) in all_params)
nothing # hide
```

## Propagate and plot a trajectory

We can see an example trajectory arising from one of our initial state guesses by
propagating it with the default propagator options, and then plotting it. Let's
try to plot the first orbit in the "g2" family:

```@example 1
trajectory = solve(all_states["g2"][1])
p1 = plot(trajectory)
```

## Correct the orbits into periodic orbits

We will now find the exact periodic orbits corresponding to our initial state guesses
by using an axisymmetric single-shooting differential corrector scheme. This is done
simply by solving our states using the `DiffCorrectAxisymmetric`(@ref) corrector.

```@example 1
all_half_orbits = Dict(family => begin
    [solve(s, DiffCorrectAxisymmetric()) for s in states]
end for (family, states) in all_states)
```

!!! note "Differential corrector orbits"
    Note that `DiffCorrectAxisymmetric` differential corrector orbits are returned propagated
    only for half their period (i.e. until the ``x``-axis crossing).

Now, we can super-impose this corrected trajectory on our previous plot:

```@example 1
plot!(p1, all_half_orbits["g2"][1]; color=:blue)
```

## Propagate and plot full periodic orbits

Finally, we end by propagating our `all_half_orbits` above for a full orbital period each,
and plot them.

```@example 1
plots = Dict(family => begin
    # half_orbits are given in the (0., t[end]) timespan.
    # We convert these to full orbits by propagating them in the (0., t[end]*2) timespan.
    orbits = [solve(remake(State(o), tspan=(0., o.t[end] * 2))) for o in half_orbits]

    # Plot each of the orbits individually.
    plots = [plot(o; color=:blue) for o in orbits]

    # Plot the whole family as a single row.
    plot(plots...; layout=(1, length(plots)))
end for (family, half_orbits) in all_half_orbits)
nothing # hide
```

Now we can plot the full figure by laying out our plots appropriately:

```@example 1
plot(plots["g2"], plots["g1"], plots["DRO"]; layout=(3, 1),
     legend=false, grid=false, minorgrid=false, showaxis=false, size=(500,500))
```