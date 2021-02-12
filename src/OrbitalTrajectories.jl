module OrbitalTrajectories

    using Reexport
    using GeneralizedGenerated: NGG

    function __init__()
        # XXX: Ensure that NGG is initialized appropriately, so that
        # GeneralizedGenerated functions work properly. See the NGG_modules constant.
        map(NGG.module_index, NGG_modules)
        nothing
    end

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

    #-----------#
    # CONSTANTS #
    #-----------#
    # Output of this gives the order that we should touch module_index() in __init__
    const NGG_modules = NGG._modules

end # module