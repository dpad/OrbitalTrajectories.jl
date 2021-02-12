using Pkg.Artifacts

artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
spice_kernels = artifact_hash("spice_kernels", artifact_toml)

if isnothing(spice_kernels) || !artifact_exists(spice_kernels)
    spice_kernels = create_artifact() do artifact_dir
        url_base = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
        kernels = [
            "pck/pck00010.tpc",
            "pck/gm_de431.tpc",
            "lsk/latest_leapseconds.tls",
            "spk/planets/de430.bsp",
            "spk/satellites/mar097.bsp",
            "spk/satellites/jup310.bsp",
        ]
        for kernel in kernels 
            file_path = joinpath(artifact_dir, basename(kernel))
            if !isfile(file_path)
                download(joinpath(url_base, kernel), joinpath(artifact_dir, basename(kernel)))
            end
        end
    end

    bind_artifact!(artifact_toml, "spice_kernels", spice_kernels, force=true)
end

spice_path = artifact_path(spice_kernels)
@info "SPICE kernels in $(spice_path)"