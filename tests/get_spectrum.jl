using StatsBase
include("../tests/tests_utils.jl")

function get_s(methods:: Array{Symbol}, iDFT:: Array{Number}, interferogram:: Vector{Number}; λ=1.0)
    # Your method here

    for mtd ∈ methods
        func     = eval(mtd)
        spectrum = func(iDFT, interferogram; λ=λ)
    end
end

function get_spectogram_from_sample(methods:: Array{Symbol}, sampleₙ:: Vector{Float64}, iDFTₙ:: Array{Number}, interferogramₙ:: Vector{Number}, percent:: Number; λ=1.0)
    # m sub-samples
    n              = length(interferogramₙ)
    rows_id        = sample(1:n, round(Int,percent*n))
    tₘ             = sampleₙ[rows_id]
    interferogramₘ = interferogramₙ[rows_id]
    iDFTₘ          = iDFTₙ[rows_id, :]

    return [tₘ, interferogramₘ, iDFTₘ, get_s(methods, interferogramₘ, iDFTₘ; λ=λ)]
end

spectogram_data    = # Spectogram here
interferogram_data = ifft(real_data.+im.*imag_data)

[sampleₙ, interferogramₙ, iDFTₙ]=get_normally_dist_data(t_data, interferogram_data, n)

get_spectogram_from_sample(methods, sampleₙ, interferogramₙ, iDFTₙ, percent)