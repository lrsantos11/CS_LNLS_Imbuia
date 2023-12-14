using Plots
using LaTeXStrings

function proximal_subgradient(f:: Function, g:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, s:: Number, k_max:: Int64; ϵ=eps, p=Inf)
    x=x₀
    L=s
    k=0
    solved=false
    
    while !(solved || k>=k_max)
        ∂fx=∂f(x)
        fx=f(x)
        L, x=Lₖ(L, k, x, fx, ∂fx) #Backtracking mais atualização

        k+=1
        solved=norm(∂fx, p)<ϵ
    end 

    return x
end