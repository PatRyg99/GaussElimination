# Patryk Rygiel
module blocksys

function gauss_elimination(A, b::Vector; partial_choice::Bool = false)
    """Basic Gauss elimination taking into account specific matrice structure"""

    n = A.size
    l = A.sub_size

    rows = l

    # First loop goes through all the rows - O(n)
    for k in 1:(n - 1)

        # Partial choice of main element
        if partial_choice

            idx = k
            max = abs(A.data[k][k - A.shift[k]])

            end_row = k + l - (k % l) - 1
            
            # Iterate over rows with non-zeros and check for column max
            for i in k+1:end_row

                v = abs(A.data[i][k - A.shift[i]])

                if max < v
                    idx = i
                    max = v  
                end
            end
            
            # Swap i-th row with k-th
            if idx != k
                
                # If swap between sections - extend row to have additional
                # l zeros at the end
                if k % 4 == 0
                    z = zeros(l)
                    A.data[k] = [A.data[k], z]
                end

                temp_A = A.data[k]
                A.data[k] = A.data[idx]
                A.data[idx] = temp_A

                temp_shift = A.shift[k]
                A.shift[k] = A.shift[idx]
                A.shift[idx] = temp_shift

                temp_b = b[k]
                b[k] = b[idx]
                b[idx] = temp_b
            end
        end
        
        if k >= rows
            rows += l
        end

        # Second loop needs to update only up to l rows - O(l)
        for i in (k + 1):rows
            
            m = A.data[i][k - A.shift[i]] / A.data[k][k - A.shift[k]]
            cols = min(A.shift[k] + length(A.data[k]), n)
            
            # Third loop needs to update only up to 2*l - O(l)
            # Longest row in array has 2*l columns due to possible swaps
            for j in (k + 1):cols
                A.data[i][j - A.shift[i]] -= m * A.data[k][j - A.shift[k]]
            end

            b[i] -= m * b[k]
        end
    end

    # Compute vector x, having an upper triangle matrice
    x = Vector{Float64}(undef, length(b))
    x[n] = b[n] / A.data[n][n - A.shift[n]]

    # Iterate from last row up - O(n)
    for k in reverse(1:(n-1))
        
        x[k] = b[k]
        end_col = min(A.shift[k] + length(A.data[k]), n)

        # Iterate over columns, only up to 2*l - O(l)
        # Longest row in array has 2*l columns due to possible swaps
        for j in (k + 1):end_col
            x[k] -= A.data[k][j - A.shift[k]] * x[j]
        end

        x[k] /= A.data[k][k - A.shift[k]]
    end

    return x
end


function LU_decomposition(A; partial_choice::Bool = false)
    """Decomposition of a Matrice into L and U matrices"""

    n = A.size
    l = A.sub_size

    rows = l

    # Permutation vector
    P = collect(1:n)

    # First loop goes through all the rows - O(n)
    for k in 1:(n - 1)

        # Partial choice of main element
        if partial_choice

            idx = k
            max = abs(A.data[k][k - A.shift[k]])

            end_row = k + l - (k % l) - 1
            
            # Iterate over rows with non-zeros and check for column max
            for i in k+1:end_row

                v = abs(A.data[i][k - A.shift[i]])

                if max < v
                    idx = i
                    max = v  
                end
            end
            
            # Swap i-th row with k-th
            if idx != k
                
                # If swap between sections - extend row to have additional
                # l zeros at the end
                if k % 4 == 0
                    z = zeros(l)
                    A.data[k] = [A.data[k], z]
                end

                temp_A = A.data[k]
                A.data[k] = A.data[idx]
                A.data[idx] = temp_A

                temp_shift = A.shift[k]
                A.shift[k] = A.shift[idx]
                A.shift[idx] = temp_shift

                temp_p = P[k]
                P[k] = P[idx]
                P[idx] = temp_p
            end
        end

        if k >= rows
            rows += l
        end
        
        # Elements in column below are being divided by the diagonal element
        d = A.data[k][k - A.shift[k]]

        # Second loop needs to update only up to l rows - O(l)
        for i in (k + 1):rows
            
            A.data[i][k - A.shift[i]] /= d 
            cols = min(A.shift[k] + length(A.data[k]), n)
            
            # Third loop needs to update only up to 2*l - O(l)
            # Longest row in array has 2*l columns due to possible swaps
            for j in (k + 1):cols

                # Apply Schur complement
                A.data[i][j - A.shift[i]] -= A.data[i][k - A.shift[i]] * A.data[k][j - A.shift[k]]
            end
        end
    end

    return A, P
end


function LU_gauss_elimination(A, P::Vector, b::Vector)
    """Computing Gauss elimination on LU matrice - saved in A matrice"""

    n = A.size
    x = Vector{Float64}(undef, length(b))

    # Solving Ly=b and applying permutation vector - O(n)
    for k in 1:n
        x[k] = b[P[k]]

        # Iterate over columns, only up to 2*l - O(l)
        # Longest row in array has 2*l columns due to possible swaps
        for j in (A.shift[k] + 1):(k-1)
            x[k] -= A.data[k][j - A.shift[k]] * x[j]
        end 
    end

    # Solving Ux=y - O(n)
    for k in reverse(1:n)
    
        end_col = min(A.shift[k] + length(A.data[k]), n)

        # Iterate over columns, only up to 2*l - O(l)
        # Longest row in array has 2*l columns due to possible swaps
        for j in (k + 1):end_col
            x[k] -= A.data[k][j - A.shift[k]] * x[j]
        end

        x[k] /= A.data[k][k - A.shift[k]]
    end

    return x
end

end # End blocksys