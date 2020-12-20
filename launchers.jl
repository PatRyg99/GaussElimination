# Patryk Rygiel

module Launchers

using TickTock
using LinearAlgebra: norm

include("./data_structs.jl")
include("./blocksys.jl")
include("./utils.jl")

using .DataStructs: load_sparse_matrice, load_vector, load_row_matrice
using .blocksys: gauss_elimination, LU_decomposition, LU_gauss_elimination
using .Utils: compute_side_vector

function launch_gauss(A_file::String; b_file::String = "", partial_choice::Bool = false)
    """
    Launching Gauss elimination on given A and b files.
    If no b vector given than compute it as A * (1, ..., 1) = b

    Inputs:
	    A_file: file with stored A matrice
        b_file: file with stored b vector
        partial_choice: whether to use algorithm with it or without
    
    Outputs:
        x': output vector
        timer: algorithm time_1
        err: (only when no b_file) error between x and x'

    """

    # Loading A matrice from file
    A = load_row_matrice(A_file)

    # Loading b vector from file if file given
    if sizeof(b_file) == 0
        b = compute_side_vector(A)
    else
        b = load_vector(b_file)
    end

    tick()
    x = gauss_elimination(A, b, partial_choice=partial_choice)
    timer = tok()
    
    x_file_name = dirname(A_file) * "/x.txt"

    # Saving output into file
    if sizeof(b_file) == 0
        err = norm(x - ones(size(x))) / norm(ones(size(x)))

        open(x_file_name, "w") do out
            write(out, string(err) * "\n")
            write(out, join(string.(x), "\n"))
        end

        return x, timer, err

    else
        open(x_file_name, "w") do out
            write(out, join(string.(x), "\n"))
        end

        return x, timer
    end
end


function launch_LU_gauss(A_file::String; b_file::String = "", partial_choice::Bool = false)
    """
    Launching LU Gauss elimination on given A and b files.
    If no b vector given than compute it as A * (1, ..., 1) = b

    Inputs:
        A_file: file with stored A matrice
        b_file: file with stored b vector
        partial_choice: whether to use algorithm with it or without

    Outputs:
        x': output vector
        timer: algorithm time_1
        err: (only when no b_file) error between x and x'

    """

    # Loading A matrice from file
    A = load_row_matrice(A_file)

    # Loading b vector from file if file given
    if sizeof(b_file) == 0
        b = compute_side_vector(A)
    else
        b = load_vector(b_file)
    end

    tick()
    A, P = LU_decomposition(A, partial_choice=partial_choice)
    x = LU_gauss_elimination(A, P, b)
    timer = tok()

    x_file_name = dirname(A_file) * "/x.txt"

    # Saving output into file
    if sizeof(b_file) == 0
        err = norm(x - ones(size(x))) / norm(ones(size(x)))

        open(x_file_name, "w") do out
            write(out, string(err) * "\n")
            write(out, join(string.(x), "\n"))
        end

        return x, timer, err
    else
        open(x_file_name, "w") do out
            write(out, join(string.(x), "\n"))
        end

        return x, timer
    end
end

end # End Launcher