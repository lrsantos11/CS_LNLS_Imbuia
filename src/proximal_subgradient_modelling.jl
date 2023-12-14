include("group_sparse_functions.jl")
include("proximal_subgradient_solver.jl")

function proximal_gradient_method(iDFT:: Matrix{Number}, interferogram:: Vector{Number}; h=h, λ=1.0)
    k_max=1000

    f(x:: Vector{<:Number}; iDFT=iDFT, interferogram=interferogram)= 0.5*norm(iDFT*x.-interferogram, 2)^2

    ∇f(x:: Vector{<:Number})=iDFT'(iDFT*x.-interferogram)

    Lf=1/2

    x₀=zeros(length(interferogram))

    Lₖ(L:: Number, k:: Int64, x:: Array{<:Number}, fx:: Number, ∂fx:: Array{<:Number})=L, proxhL(L, x.-∂fx./L; λ=λ)

    proximal_subgradient(f, h, ∇f, Lₖ, x₀, Lf, k_max)
end



