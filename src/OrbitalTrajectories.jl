module OrbitalTrajectories

    using Reexport

    #----------#
    # INCLUDES #
    #----------#
    include("piracy.jl")

    include("SPICE/spice.jl")
    include("dynamics/dynamics.jl")

    #----------------#
    # IMPORT MODULES #
    #----------------#
    @reexport using OrbitalTrajectories.Dynamics
    @reexport using OrbitalTrajectories.SpiceUtils

end # module