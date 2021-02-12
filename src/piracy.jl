using ModelingToolkit
using Unitful
using LinearAlgebra
using ForwardDiff

# TODO: Check which of these are still in use and/or necessary.

# Avoid Plots compatibility issue with Unitful
Base.occursin(s::Any, q::Unitful.Quantity) = occursin(s, repr(q))
Base.pointer(q::Unitful.Quantity) = pointer(repr(q))

# Ensure that ForwardDiff values can get sent to C calls
Base.unsafe_convert(T::Type{<:Any}, x::ForwardDiff.Dual) = T(ForwardDiff.value.(x))

# Compute norms for arrays of symbolic variables (as needed by EphemerisNBP)
LinearAlgebra.norm(a::AbstractArray{<:ModelingToolkit.Num}) = sum(a .^ 2)^(1/2)