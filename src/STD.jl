


import Statistics: std
import Loess: loess, predict
import DataFrames: DataFrame

# function loess_smooth(f::AbstractVector{Float64})
#     n = length(f)
#     predict(loess(Float64.(1:n), f), Float64.(1:n))
# end

# function loess_smooth(f::AbstractVector{Float64}, x::AbstractVector{Float64})
#     predict(loess(x, f), x)
# end

function loess_smooth(f::AbstractVector{Float64}; tol=0.0)
    out = copy(f)
    idx = findall(v -> abs(v) > tol, f)
    length(idx) < 3 && return out
    xs = Float64.(idx)
    out[idx] = predict(loess(xs, f[idx]), xs)
    return out
end

function loess_smooth(f::AbstractVector{Float64}, x::AbstractVector{Float64}; tol=0.0)
    out = copy(f)
    idx = findall(v -> abs(v) > tol, f)
    length(idx) < 3 && return out
    xs = x[idx]
    out[idx] = predict(loess(xs, f[idx]), xs)
    return out
end

function loess_smooth(f::AbstractMatrix{Float64})
    mapslices(loess_smooth, f, dims=1)
end




function STD(sample::AbstractArray; kwargs...)
    err = std(sample; kwargs...)
    if ndims(sample) == 3
        err = err |> x -> dropdims(x, dims=3)
    end
    return loess_smooth(err)
end

function df_STD(sample::AbstractArray, name::AbstractVector{<:Union{Symbol, String}}; kwargs...)
    err = STD(sample; kwargs...)
    df = DataFrame(err, name)
    return df
end




