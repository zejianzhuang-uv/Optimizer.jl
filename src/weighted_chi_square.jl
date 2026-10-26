

function χsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64})
    dy = @. (yth - yexp) * inv_σ # broadcast
    chi2 = dot(dy, dy)
    return chi2
end

function χsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractMatrix{Float64})
    dy = yth - yexp 
    chi2 = dot(dy, inv_σ * dy) # broadcast
    return chi2
end

function weightedχsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64}, nk::Int64)
    chi2 = χsq(yth, yexp, inv_σ)
    return chi2 / nk
end

function weightedχsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractMatrix{Float64}, nk::Int64)
    chi2 = χsq(yth, yexp, inv_σ)
    return chi2 / nk
end

"""
Weighted-reduced χ square
- Nsp: The number of the scattering processes: K⁻p → {π⁰Λ, ...} and π⁻p → {ηn, ...}
- nₖ: The number of the observables in the scattering process k
- Ntot: Total number of observables in Nsp scattering processes
"""
function red_weightedχsq(all_weight_chi2::Float64, Ntot::Int64, nfit::Int64, Nsp::Int64)
    chi2 = (Ntot / (Ntot - nfit)) * (1/Nsp) * all_weight_chi2
    return chi2
end





"""
Hobson–Bridle–Lahav χ²: −2 ln L after marginalising the dataset weight α with prior e^(−α)
"""
function hyperχsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64}, n::Int64)
    chi2 = χsq(yth, yexp, inv_σ)
    return (n + 2) * log(chi2 + 2)
end

"""
Scheme B: hyperχsq divided by n, same normalisation as weightedχsq
"""
function whyperχsq(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64}, n::Int64)
    return hyperχsq(yth, yexp, inv_σ, n) / n
end

"""
Effective weight of the dataset, α = (n + 2) / (χ² + 2)
"""
function hyperα(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64}, n::Int64)
    chi2 = χsq(yth, yexp, inv_σ)
    return (n + 2) / (chi2 + 2)
end