using JuMP
using HiGHS



function cs_LP_model(M::AbstractMatrix{T}, y::Vector{T}) where {T<:Number}
    num_row, num_col = size(M)
    size_y = length(y)
    @assert num_row == size_y

    model = Model(HiGHS.Optimizer)

    # Define variables
    @variables(model, begin
        x⁺[1:ncol] >= 0
        x⁻[1:ncol] >= 0
    end)

    # Define objective function
    @objective(model, Min, sum(x⁺) + sum(x⁻))

    # Define constraints

    @constraints(model, begin
        M * x⁺ - M * x⁻ .== y
    end)

    optimize!(model)

    return value.(x⁺) - value.(x⁻)
end