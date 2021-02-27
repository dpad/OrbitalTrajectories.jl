using Pkg.Artifacts
const artifact_toml = find_artifacts_toml(@__DIR__)

# Kernels to be downloaded lazily
const SPICE_url_base = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"
const all_kernels = Dict(
    "meta"    => ["pck/pck00010.tpc",
                  "pck/gm_de431.tpc",
                  "lsk/latest_leapseconds.tls"],  # Required for all other kernels
    "default" => ["spk/planets/de430.bsp"],
    "mars"    => ["spk/satellites/mar097.bsp"],
    "jupiter" => ["spk/satellites/jup310.bsp"],
    "saturn"  => ["spk/satellites/sat427.bsp"],
    "uranus"  => ["spk/satellites/ura111.bsp"],
    "neptune" => ["spk/satellites/nep095.bsp"],
    "pluto"   => ["spk/satellites/plu055.bsp"]
)
for (system, kernels) in all_kernels
    artifact_name = "$(system)_kernels"
    artifact = artifact_hash(artifact_name, artifact_toml)
    if isnothing(artifact)
        # Only download the kernel if it doesn't exist in Artifacts.toml. Otherwise, we want to
        # download it lazily in SpiceUtils.
        download_info = [(join((SPICE_url_base, kernel), "/"), "") for kernel in kernels]
        artifact = create_artifact() do artifact_dir
            for (download_URL, _) in download_info
                @info "Downloading $(basename(download_URL))..."
                download(download_URL, joinpath(artifact_dir, basename(download_URL)))
            end
        end
        if Sys.islinux()
            bind_artifact!(artifact_toml, artifact_name, artifact; download_info, lazy=true, force=true)
        else
            @warn "Could not bind artifact, need to run on Linux so the hashes are correct!"
        end
    end
end