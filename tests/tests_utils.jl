function iDFTx(x:: Vector{Number}; n=n)
    n₂=round(Int, n/2)
    ω₁=[exp(2*π*im*x[m])/n for m=1:n]

    return reduce(hcat, [ω₁.^(k/n) for k=-n₂:n₂])
end

function get_normally_dist_data(t_data:: Vector{Number}, interferogram_data:: Vector{Number}, n:: Int64; iDFTx=iDFTx)
    # n random samples randomly distributed
    sampleₙ        = sample(t_data, n)
    interferogramₙ = interferogram_data[sampleₙ] 
    iDFTₙ          = iDFTx(sampleₙ)

    return [sampleₙ, interferogramₙ, iDFTₙ]
end