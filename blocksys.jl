# Patryk Rygiel
module blocksys

include("./structs/matrice.jl")
using .Matrice: SparseMatrice

function gauss_elimination(A, b::Vector, main_element=nothing)

    n = A.size
    
    # If no main element chosen perform classic algorithm
    if isnothing(main_element)
        for k in 1:(n - 1)
            for i in (k + 1):n
                
                m = A.data[i, k] / A.data[k, k]

                for j in (k + 1):n
                    A.data[i, j] = A.data[i, j] - m * A.data[k, j]
                end

                b[i] = b[i] - m * b[k]
            end
        end
    end

    # Compute vector x, having an upper triangle matrice
    x = Vector{Float64}(undef, length(b))
    x[n] = b[n] / A.data[n, n]

    for k in reverse(1:(n-1))
        x[k] = b[k]

        for j in (k + 1):n
            x[k] -= A.data[k, j] * x[j]
        end

        x[k] /= A.data[k, k]
    end

    return x
end

end # End blocksys