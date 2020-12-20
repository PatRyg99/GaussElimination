# Patryk Rygiel

module Utils

function compute_side_vector(A)
    """
    Function computing b fomr equation: Ax = b.
    Having x = (1, ..., 1)
    """
    b = Vector{Float64}(undef, A.size)
    for row in 1:A.size
        b[row] = sum(A.data[row])
    end

    return b
end

end # End Utils