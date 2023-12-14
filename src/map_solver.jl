using BenchmarkTools
using ProximalOperators


"""
Using the Method of Alternating Projections for the Basis Pursuit problem
    
min ||x||₁ s.t. Ax = b 
"""
function BP_MAP(Affine;
    itmax::Int=100,
    ε::Number=1e-6,
    verbose::Bool=true,
    x₀::Vector{Float64}=Float64[],
    kwargs...)
    m, n = size(Affine.A)
    ProjAffine(x) = ProjectIndicator(Affine, x)
    if isempty(x₀)
        x₀ = zeros(n)
    end
    xMAP = ProjAffine(x₀)
    distance = norm(xMAP, 2)
    radius = 0.0
    solved = false
    tired = false
    it = 0
    inner_it = 0
    status = :Tired
    verbose && @info "m = $m, n = $n"
    while !(solved || tired)
        radius += distance
        BallL1 = IndBallL1(radius)
        Proj_BallL1(x) = ProjectIndicator(BallL1, x)
        resultMAP = MAP(xMAP, ProjAffine, Proj_BallL1, itmax=itmax, gap_distance=false, kwargs...)
        inner_it += resultMAP.iter_total
        zMAP = resultMAP.xApprox
        xMAP = ProjAffine(zMAP)
        distance = norm(xMAP - zMAP, 2)
        it += 1
        solved = distance < ε
        if solved
            verbose && @info "solved"
            verbose && @info "it = $it"
            verbose && @info "inner_it = $inner_it"
            status = :Solved
        end
        tired = it >= itmax
    end
    return xMAP, it, inner_it, status
end

"""
    BP_MAP(A, b; kwargs...)
    Basis Pursuit using the Method of Alternating Projections
    min ||x||₁ s.t. Ax = b 
"""
BP_MAP(A::AbstractMatrix, b::Vector; kwargs...) = BP_MAP(IndAffine(A,b), kwargs...)


#### UTILS FUNCTIONS

"""
    proj = ProjectIndicator(indicator,x)
    Projection using Indicator Function from `ProximalOperators.jl`
    """
function ProjectIndicator(indicator, x)
    proj, _ = prox(indicator, x)
    return proj
end
    

"""
        compute_error( x, xold; xsol, normytpe)
"""

function compute_error(x::AbstractVector, xold::AbstractVector, xsol::AbstractVector;
    norm_p::Number = 2)
    if isempty(xsol)
        return norm(x - xold, norm_p)
    else
        return norm(x - xsol, norm_p)
    end
end

"""
    MAP(x₀, ProjectA, ProjectB)

    Method of Alternating Projections
"""
function MAP(x₀::Vector, 
            ProjectA::Function, 
            ProjectB::Function;
    ε::Float64 = 1e-5,
    itmax::Int = 100,
    xSol::Vector = [],
    gap_distance::Bool = true)
    k = 0
    tolMAP = 1.0
    xMAP = copy(x₀)
    ProjA = ProjectA(xMAP)
    while tolMAP > ε && k < itmax
        ProjA = ProjectA(xMAP)
        if gap_distance
            xMAP = ProjectB(ProjA)
            tolMAP = norm(ProjA - xMAP)
        else
            xMAPOld = copy(xMAP)
            xMAP = ProjectB(ProjA)
            tolMAP = compute_error(xMAP, xMAPOld, xSol)
        end
        k += 2
    end
    return k, tolMAP, xMAP
end

