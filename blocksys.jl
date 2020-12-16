# Patryk Rygiel
module blocksys

include("./structs/matrice.jl")
using .Matrice: SparseMatrice

function gauss_elimination(A, b::Vector, main_element=nothing)

    n = A.size
    l = A.sub_size

    # If no main element chosen perform classic algorithm
    if isnothing(main_element)
        
        rows = l

        # First loop goes through all the rows - O(n)
        for k in 1:(n - 1)
            
            # Update only up to l next rows
            if k >= rows
                rows += l
            end

            # Second loop needs to update only up to  rows - O(l)
            for i in (k + 1):rows
                
                m = A.data[i, k] / A.data[k, k]

                # Update only up to l next columns
                cols = min(rows + l, n)
                
                # Third loop needs to update only up to l columns - O(l)
                for j in (k + 1):cols
                    A.data[i, j] = A.data[i, j] - m * A.data[k, j]
                end

                b[i] = b[i] - m * b[k]
            end
        end
    end

    # Compute vector x, having an upper triangle matrice
    x = Vector{Float64}(undef, length(b))
    x[n] = b[n] / A.data[n, n]

    # Iterate from last row up - O(n)
    for k in reverse(1:(n-1))
        
        x[k] = b[k]
        end_col = min(k + l, n) 

        # Iterate over columns, only up to l - O(l)
        for j in (k + 1):end_col
            x[k] -= A.data[k, j] * x[j]
        end

        x[k] /= A.data[k, k]
    end

    return x
end

end # End blocksys