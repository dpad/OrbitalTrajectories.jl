using SafeTestsets

# --- TESTS ---
@safetestset "SPICE Utils" begin include("spice_utils.jl") end
@safetestset "NBP" begin include("nbp.jl") end
@safetestset "3BP" begin include("3bp.jl") end
@safetestset "Document tests" begin include("doctests.jl") end
@safetestset "Sensitivity" begin include("sensitivity.jl") end