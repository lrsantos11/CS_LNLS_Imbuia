function get_normally_dist_data(t_data:: Vector{Number}, interferogram_data:: Vector{Number}, n:: Int64; iDFTx=iDFTx)
    # n random samples randomly distributed
    tₙ             = sample(t_data, n)
    interferogramₙ = interferogram_data[tₙ] 
    iDFTₙ          = iDFTx(tₙ)

    return [tₙ, interferogramₙ, iDFTₙ]
end

function iDFTx(x:: Vector{Number}; n=n, n₂=n₂)
    ω₁=[exp(2*π*im*x[m])/n for m=1:n]

    return reduce(hcat, [ω₁.^(k/n) for k=-n₂:n₂])
end