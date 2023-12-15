function iDFTx(x:: Vector{<:Number}; n=1)
    n₂ = round(Int, n/2)
    ω(x) = exp(2*π*im*x/n)
    ω₁ = ω.(x)

    return vander(ω₁, -n₂:n₂)
end

function get_normally_dist_data(t_data:: Vector{<:Number}, interferogram_data:: Vector{<:Number}, n:: Int64)
    # n random samples randomly distributed
    sampleₙ = sample(t_data, n)
    interferogramₙ = interferogram_data[sampleₙ] 
    iDFTₙ          = iDFTx(sampleₙ; n=n)

    return sampleₙ, interferogramₙ, iDFTₙ
end


function vander!(V::AbstractMatrix, v::AbstractVector, range_vander::UnitRange{Int64})
    S = eltype(v)
    m = length(v)
    n = length(range_vander)
    (m, n) == size(V) || throw(DimensionMismatch())
    n₂ = div(n, 2)
    med = range_vander[n₂+1] + n₂ + 1
    V[:, med] .= ones(S,m)

    for j in med+1:n
        @inbounds V[:, j] = v .* V[:, j-1]
    end
    for j in med-1:-1:1
        @inbounds V[:, j] = V[:, j+1] ./ v
    end
    return V
end

vander(x::AbstractVector, range) = vander!(Array{eltype(v)}(undef, length(v), n), x, range)