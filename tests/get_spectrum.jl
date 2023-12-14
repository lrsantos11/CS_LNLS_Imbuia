using StatsBase

function get_s(methods:: Vector{Symbol}, sample:: Vector{Float64}, iDFT:: Array{Float64}:: Vector{Float64}, interferogram; λ=1.0)
    # Your method here
    for mtd ∈ methods

    end
end

function get_spectogram_from_sample(methods:: Vector{Symbol}, sampleₙ:: Vector{Float64}, iDFTₙ:: Array{Float64}, interferogramₙ:: Vector{Float64}, percent:: Number, n:: Inf64; λ=1.0)
    # m sub-samples
    rows_id        = sample(1:n, round(Int,percent*n))
    tₘ             = sampleₙ[rows_id]
    interferogramₘ = interferogramₙ[rows_id]
    iDFTₘ          = iDFTₙ[rows_id, :]

    return [tₘ, interferogramₘ, iDFTₘ, get_s(methods, tₘ, iDFTₘ, interferogramₘ; λ=λ)]
end