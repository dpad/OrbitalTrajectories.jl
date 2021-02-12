# [Koopman Expectation under Uncertainty](@id Koopman)

This tutorial will briefly introduce computing expectation values from orbits under uncertainty, using the Koopman 
expectation introduced by [DiffEqUncertainty.jl](https://github.com/SciML/DiffEqUncertainty.jl).[^Gerlach2020]

[^Gerlach2020]: A. R. Gerlach, A. Leonard, J. Rogers, and C. Rackauckas, “The Koopman Expectation: An Operator Theoretic Method for Efficient Analysis and Optimization of Uncertain Hybrid Dynamical Systems,” arXiv: 2008.08737, 2020.

## Load packages

First, we load all necessary packages:

```@example 1
using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using Unitful                # Units (u"km", u"d" (days)) and unit conversion
using DiffEqUncertainty      # To compute the Koopman expectation
using Distributions          # To define uncertainty distributions
nothing # hide
```

## Differential Correction

We begin with a Europa fly-by test case[^Pellegrini2016] in the CR3BP (Jupiter-Europa) system:

[^Pellegrini2016]: E. Pellegrini and R. P. Russell, “On the Computation and Accuracy of Trajectory State Transition Matrices,” Journal of Guidance, Control, and Dynamics, Vol. 39, No. 11, 2016, pp. 2485–2499, 10.2514/1.G001920.

```@example 1
system = CR3BP(:jupiter, :europa)

# Periodic orbit around Jupiter with flybys of Europa
u0 = [1.0486505808029702, 0., 0., 0., -0.09354853663949217, 0.]
full_tspan = (0., 37.)
full_state = State(system, u0, full_tspan)
full_trajectory = solve(full_state)

# Truncate the trajectory to just a single flyby section
tspan = (33., 37.)
state = State(system, full_trajectory.sol(tspan[1]), tspan)
trajectory = solve(state)
plot(trajectory)
```

## Uncertainty in Initial State

We can introduce some uncertainty in the initial state by specifying some distributions which define
the uncertainty present in the initial values. Then, we can draw random samples from these to plot
a statistical ensemble of trajectories (this is a form of Monte Carlo sampling).

```@example 1
# Define some Normal distributions for the initial values of x, y, and ẋ, ẏ
pos_stdev = uconvert(NoUnits, 100u"km" / system.props.L)
vel_stdev = uconvert(NoUnits, 30u"m/s" / system.props.V)
N_pos = [truncated(Normal(u, pos_stdev), u-3*pos_stdev, u+3*pos_stdev) for u in state.u0[1:2]]
N_vel = [truncated(Normal(u, vel_stdev), u-3*vel_stdev, u+3*vel_stdev) for u in state.u0[4:5]]

# New problems will be defined by randomly sampling from the above distributions
prob_func(prob, i, r) = remake(prob, u0 = [rand.(N_pos)..., state.u0[3], rand.(N_vel)..., state.u0[6]])

# Check for collisions with Europa
callback = collision(system, secondary_body(system))

# Define and solve a probabilistic ensemble of these trajectories
ensemble = EnsembleProblem(state; prob_func)
trajectories = solve(ensemble, OrbitalTrajectories.Dynamics.DEFAULT_ALG, EnsembleThreads();
                        trajectories=200, callback)

# Plot the ensemble and original trajectory
p = plot(trajectories; alpha=0.1, arrow=true, color_palette=:linear_blue_5_95_c73_n256, nomodel=true)
plot!(p, trajectory; c=:black, lw=3)
```

## Computing Expectations (Monte Carlo)

Note that some of the random trajectories crash into the secondary body (Europa). We would like to
compute the expectation (i.e. average) for the likelihood that our trajectory will crash into
the secondary body, given the uncertainties defined for its initial state.

The first way is to compute this from our Monte Carlo samples:

```@example 1
# Compute the running mean probability of crashing from the Monte Carlo samples
# NOTE: This uses the `crashed()` function from OrbitalTrajectories.
collision_prob_MC = cumsum(map(crashed, trajectories)) ./ (1:length(trajectories))
p = plot(collision_prob_MC; label="Monte Carlo expectation",
            xlabel="Monte Carlo simulations", ylabel="Collision probability (%)")
```

## Computing Expectations (Koopman)

Alternatively, we can use the `Koopman()` operator from the `DiffEqUncertainty.jl` library to
compute the expectation directly, without needing to sample from the statistical distribution of
trajectories above.

```@example 1
# Define the initial uncertain state, which just contains Normally-distributed and constant values
u0_uncertain = [N_pos..., state.u0[3], N_vel..., state.u0[6]]

# Compute the true expectation value for `crashed()` function using the Koopman operator
collision_prob_Koopman = expectation(crashed, state.prob, u0_uncertain, state.p, Koopman(),
                            OrbitalTrajectories.Dynamics.DEFAULT_ALG; callback)

# Plot the Koopman collision probability
hline!(p, [collision_prob_Koopman.u]; label="Koopman expectation")
```