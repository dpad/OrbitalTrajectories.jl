using Pkg.Artifacts

const artifact_toml = find_artifacts_toml(@__DIR__)
const SPICE_url_base = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"

# General planetary ephemerides are downloaded by default.
const default_kernels = artifact_hash("spice_kernels", artifact_toml)
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
                download(join((SPICE_url_base, kernel), "/"), joinpath(artifact_dir, basename(kernel)))
            end
        end
    end
    bind_artifact!(artifact_toml, "spice_kernels", spice_kernels; force=true)
end

# System-specific kernels are downloaded lazily
const kernels = Dict(
    "mars"    => "mar097.bsp",
    "jupiter" => "jup310.bsp",
    "saturn"  => "sat427.bsp",
    "uranus"  => "ura111.bsp",
    "neptune" => "nep095.bsp",
    "pluto"   => "plu055.bsp"
)
for (system, kernel) in kernels
    artifact_name = "$(system)_ephemerides"

    artifact = artifact_hash(artifact_name, artifact_toml)
    if isnothing(artifact)
        # Only download the kernel if it doesn't exist in Artifacts.toml. Otherwise, we want to
        # download it lazily in SpiceUtils.
        download_URL = join((SPICE_url_base, "spk", "satellites", kernel), "/")
        artifact = create_artifact() do artifact_dir
            download(download_URL, joinpath(artifact_dir, basename(kernel)))
        end
        bind_artifact!(artifact_toml, artifact_name, artifact; download_info=[(download_URL, "")], lazy=true, force=true)
    end
end