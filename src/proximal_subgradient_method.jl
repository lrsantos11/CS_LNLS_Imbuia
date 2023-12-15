include("group_sparse_functions.jl")
include("proximal_subgradient_solver.jl")

function proximal_gradient_method(iDFT:: Matrix{<:Number}, interferogram:: Vector{<:Number}; λ=1.0)
    k_max=100

    f(x:: Vector{<:Number}; iDFT=iDFT, interferogram=interferogram) = 0.5*norm(iDFT*x.-interferogram, 2)^2
    ∇f(x:: Vector{<:Number}; iDFT=iDFT, interferogram=interferogram) = iDFT'*(iDFT*x.-interferogram)
    Lf=1/2

    n=size(iDFT)[2]
    x₀=zeros(Complex{Float64}, n)

    Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, fx:: Number, ∂fx:: Array{<:Number})=L, proxhL(L, x.-∂fx./L, n, n, λ)

    return proximal_subgradient(f, ∇f, Lₖ, x₀, Lf, n, k_max)
end



