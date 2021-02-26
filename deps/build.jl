using Pkg.Artifacts

artifact_toml = find_artifacts_toml(@__DIR__)
SPICE_url_base = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"

# General planetary ephemerides are downloaded by default.
default_kernels = artifact_hash("spice_kernels", artifact_toml)
if isnothing(default_kernels) || !artifact_exists(default_kernels)
    spice_kernels = create_artifact() do artifact_dir
        kernels = [
            "pck/pck00010.tpc",
            "pck/gm_de431.tpc",
            "lsk/latest_leapseconds.tls",
            "spk/planets/de430.bsp",
        ]
        for kernel in kernels 
            file_path = joinpath(artifact_dir, basename(kernel))
            if !isfile(file_path)
                download(joinpath(SPICE_url_base, kernel), joinpath(artifact_dir, basename(kernel)))
            end
        end
    end
    bind_artifact!(artifact_toml, "spice_kernels", spice_kernels; force=true)
end

# System-specific kernels are downloaded lazily
kernels = Dict(
    "mars" => "spk/satellites/mar097.bsp",
    "jupiter" => "spk/satellites/jup310.bsp",
)
for (system, kernel) in kernels
    artifact_name = "$(system)_ephemerides"

    artifact = artifact_hash(artifact_name, artifact_toml)
    if isnothing(artifact)
        # Only download the kernel if it doesn't exist in Artifacts.toml. Otherwise, we want to
        # download it lazily in SpiceUtils.
        kernel_hash = ""
        artifact = create_artifact() do artifact_dir
            download(joinpath(SPICE_url_base, kernel), joinpath(artifact_dir, basename(kernel)))
        end
        bind_artifact!(artifact_toml, artifact_name, artifact; download_info=[(joinpath(SPICE_url_base, kernel), "")], lazy=true, force=true)
    end
end