using Documenter
using OrbitalTrajectories

get!(ENV, "JULIA_DEBUG", Documenter)  # Print @debug notices from Documenter
get!(ENV, "GKSwstype", "100")  # Avoid Plots GR issues (see Documenter.jl#1583)

makedocs(
    sitename = "OrbitalTrajectories",
    format = Documenter.HTML(),
    modules = [OrbitalTrajectories],
    pages = [
        "OrbitalTrajectories.jl" => "index.md",
        "Examples" => [
            "tutorials/DynamicModels.md",
            "tutorials/SunMarsCR3BP.md",
            "tutorials/STM.md",
            "tutorials/QSOFamilies.md",
            "tutorials/Koopman.md"
        ],
        "Package API" => "API.md"
    ],
    strict = true
)

REPO = get(ENV, "GITHUB_REPOSITORY", nothing)
@info "REPO: $(REPO)"
if !isnothing(REPO)
    deploydocs(
        repo = "github.com/$(REPO)",
        forcepush = true
    )
end