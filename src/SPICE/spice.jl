module SpiceUtils

    using Downloads
    using Dates
    using SPICE
    using ModelingToolkit
    using ProgressMeter
    using Pkg.Artifacts
    using ForwardDiff
    using StaticArrays
    using Memoize

    export mass_fraction
    export state_to_synodic, state_from_synodic
    export pos_to_synodic, pos_from_synodic

    const artifacts_toml = find_artifacts_toml(@__DIR__)
    const kernel_bodies = Dict(
        # List of bodies included in each default kernel file (as specified in build.jl)
        "default" => tuple(1:10..., 199, 299, 301, 399),
        "mars"    => tuple(401:402..., 499),
        "jupiter" => tuple(501:505..., 514:516..., 599),
        "saturn"  => tuple(601:609..., 612:614..., 632, 634, 699),
        "uranus"  => tuple(701:705..., 799),
        "neptune" => tuple(801:808..., 814, 899),
        "pluto"   => tuple(901:905..., 999)
    )
    const body_to_kernel = Dict(idx => kernel for (kernel, bodies) in kernel_bodies for idx in bodies)

    function __init__()
        load_kernel_artifact("meta_kernels")
        load_kernel_artifact("default_kernels")
        nothing
    end

    @memoize function load_ephemerides(body_name)
        body_name = lowercase(String(body_name))
        body_ID   = bodn2c(body_name)

        if body_ID ∉ keys(body_to_kernel)
            @warn """No default kernel specified for "$(titlecase(body_name))".
            You may need to manually download relevant kernels and load them with `SPICE.furnsh(path_to_kernel)`."""
            return false
        end

        artifact_name = "$(body_to_kernel[body_ID])_kernels"
        return load_kernel_artifact(artifact_name)
    end

    @memoize function load_kernel_artifact(artifact_name)
        # XXX: This is a bit of a hacky workaround of the existing Pkg.Artifacts system,
        # which supports lazy artifact loading but assumes the files it downloads are
        # gzipped. However, our kernels are just files straight from the NAIF site, so we
        # instead load them lazily here. artifact"(body)_kernels" will only work after this.

        # Get artifact details.
        meta = artifact_meta(artifact_name, artifacts_toml)
        isnothing(meta) && error("Expected to find $(artifact_name) in Artifacts.toml, but could not find it!")
        
        # Get kernel details from this artifact.
        meta_hash     = Base.SHA1(meta["git-tree-sha1"])
        download_URLs = [d["url"] for d in meta["download"]]
        kernel_dir    = artifact_path(meta_hash)
        kernel_paths  = [basename(url) for url in download_URLs]

        # Check if the kernel has already been downloaded before, and if not, download it to a temporary directory.
        temp_dir = artifact_exists(meta_hash) ? nothing : mktempdir(first(Artifacts.artifacts_dirs()))
        if !isnothing(temp_dir)
            # Create a progress bar.
            progress = (name) -> begin
                max_n = 10000
                bar = Progress(max_n; desc="Downloading NAIF kernel: $(name)", color=Base.info_color())
                (total, now) -> update!(bar, total > 0 ? round(Int, (now / total) * max_n) : 0)
            end

            for (kernel, url) in zip(kernel_paths, download_URLs)
                try
                    Downloads.download(url, joinpath(temp_dir, kernel); progress=progress(kernel))
                finally
                    finish!(bar)
                end
            end

            # All kernels are now downloaded; move them into the final kernel folder.
            mv(temp_dir, kernel_dir)
        end

        # Load each kernel into SPICE.
        for (kernel, url) in zip(kernel_paths, download_URLs)
            kernel_path = joinpath(kernel_dir, kernel)
            if !isfile(kernel_path)
                error("Could not find kernel file '$(kernel)'. Delete the $(kernel_dir) directory and try again.")
            end
            furnsh(kernel_path)
        end

        return true
    end

    @doc "Get a target body position from SPICE kernels."
    get_pos(t, target, reference; relative="ECLIPJ2000") = SVector{3,Float64}(spkpos(String(target), t, relative, "none", String(reference))[1])
    get_pos(t::ForwardDiff.Dual, target, reference; kwargs...) = get_pos(ForwardDiff.value(t), String(target), String(reference); kwargs...)
    get_pos(t::Num, tgt, ref) = get_pos(t, Num(tgt), Num(ref))
    ModelingToolkit.@register get_pos(t, target, reference)

    @doc "Get a target body state from SPICE kernels"
    get_state(t, target, reference; relative="ECLIPJ2000") = SVector{6,Float64}(spkezr(String(target), t, relative, "none", String(reference))[1])
    get_state(t::ForwardDiff.Dual, target, reference; kwargs...) = get_state(ForwardDiff.value(t), target, reference; kwargs...)

    @doc "Get body gravitational constant (GM) from SPICE kernels"
    get_GM(body::String) = bodvrd(body, "GM")[1]
    get_GM(body::Symbol) = get_GM(String(body))
    get_GM(body::Number) = get_GM(bodc2n(body))

    @doc "Get the mass fraction of two bodies (smaller mass divided by total mass)"
    mass_fraction(b1::Symbol, b2::Symbol) = mass_fraction(String(b1), String(b2))
    function mass_fraction(body_1::String, body_2::String)
        mass_1 = bodvrd(body_1, "GM")[1]
        mass_2 = bodvrd(body_2, "GM")[1]
        min(mass_1, mass_2) ./ (mass_1 .+ mass_2)
    end

    # Dynamic frame ID should be between 1400000 to 2000000 (added a bit of margin for rand) as per
    # https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#Selecting%20a%20Frame%20ID
    const INIT_DYNAMIC_FRAME_IDX = rand(1400000:1900000)
    const DYNAMIC_FRAMES_TYPE = Dict{Tuple{Symbol,Symbol},Tuple{Int,String}}
    DYNAMIC_FRAMES = DYNAMIC_FRAMES_TYPE()

    @doc """ 
    Automatically build and return a dynamic frame for a synodic
    reference (b1 is the observer, b2 is the target).
    """
    function dynamic_synodic_frame(b1::Symbol, b2::Symbol; reference="ECLIPJ2000")
        if (b1, b2) ∉ keys(DYNAMIC_FRAMES::DYNAMIC_FRAMES_TYPE)
            # XXX: Maximum size of the frame name is limited so we use the body numbers.
            bodies = bodn2c.((b1, b2))
            frame_name = "$(bodies[1])_$(bodies[2])_SYNODIC"
            path, io = mktemp()
            ID = INIT_DYNAMIC_FRAME_IDX + length(DYNAMIC_FRAMES::DYNAMIC_FRAMES_TYPE)
            write(io, 
            """ \\begindata

                FRAME_$(frame_name)        = $(ID)
                FRAME_$(ID)_NAME           = '$(frame_name)'
                FRAME_$(ID)_CLASS          = 5
                FRAME_$(ID)_CLASS_ID       = $(ID)
                FRAME_$(ID)_CENTER         = $(bodies[1])
                FRAME_$(ID)_RELATIVE       = '$(reference)'
                FRAME_$(ID)_DEF_STYLE      = 'PARAMETERIZED'
                FRAME_$(ID)_ROTATION_STATE = 'ROTATING'
                FRAME_$(ID)_FAMILY         = 'TWO-VECTOR'
                FRAME_$(ID)_PRI_AXIS       = 'X'
                FRAME_$(ID)_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
                FRAME_$(ID)_PRI_OBSERVER   = $(bodies[1])
                FRAME_$(ID)_PRI_TARGET     = $(bodies[2])
                FRAME_$(ID)_PRI_ABCORR     = 'NONE'
                FRAME_$(ID)_SEC_AXIS       = 'Y'
                FRAME_$(ID)_SEC_VECTOR_DEF = 'OBSERVER_TARGET_VELOCITY'
                FRAME_$(ID)_SEC_OBSERVER   = $(bodies[1])
                FRAME_$(ID)_SEC_TARGET     = $(bodies[2])
                FRAME_$(ID)_SEC_FRAME      = '$(reference)'
                FRAME_$(ID)_SEC_ABCORR     = 'NONE'
            """)
            close(io)

            (DYNAMIC_FRAMES::DYNAMIC_FRAMES_TYPE)[(b1, b2)] = (ID, frame_name)
            @info "Importing new $frame_name kernel from $(path)"
            furnsh(path)
        else
            ID, frame_name = (DYNAMIC_FRAMES::DYNAMIC_FRAMES_TYPE)[(b1, b2)]
        end
        return frame_name
    end

    # TODO: Fix up the SpiceUtils interface so the arguments are more consistently ordered.
    state_to_synodic(b1, b2, t; reference="ECLIPJ2000") = sxform(reference, dynamic_synodic_frame(b1, b2; reference=reference), t)
    state_to_synodic(b1, b2, t::ForwardDiff.Dual; kwargs...) = state_to_synodic(b1, b2, ForwardDiff.value(t); kwargs...)
    state_from_synodic(args...; kwargs...) = inv(state_to_synodic(args...; kwargs...))

    pos_to_synodic(b1, b2, t; reference="ECLIPJ2000") = pxform(reference, dynamic_synodic_frame(b1, b2; reference=reference), t)
    pos_to_synodic(b1, b2, t::ForwardDiff.Dual; kwargs...) = pos_to_synodic(b1, b2, ForwardDiff.value(t); kwargs...)
    pos_from_synodic(args...; kwargs...) = inv(pos_to_synodic(args...; kwargs...))

end