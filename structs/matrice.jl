# Patryk Rygiel
module Matrice

using SparseArrays

export load_matrice, load_vector

struct SideVector
    size::Integer
    data::Vector
end

struct SparseMatrice
    size::Integer
    sub_size::Integer
    data::SparseMatrixCSC
end

function load_vector(file_path::String)
    open(file_path) do file
        
        size = 0
        data = []

        for (idx, line) in enumerate(eachline(file))
            
            # First row in file holds size
            if idx == 1
                size = parse(Int, line)
            
            # Other rows hold values
            else
                push!(data, parse(Float64, line))
            end
        end
        
        # Create custom SideVector type
        vector = SideVector(size, Float64.(vec(data)))
        println("Loaded vector with size $size.")

        return vector
    end
end

function load_matrice(file_path::String)
    open(file_path) do file
        
        size = 0
        sub_size = 0
        
        # Vectors to build a sparse array from
        I = []
        J = []
        V = []
        
        for (idx, line) in enumerate(eachline(file))
            line = split(line)
            
            # First row in file holds size and sub_size
            if idx == 1
                size = parse(Int, line[1])
                sub_size = parse(Int ,line[2])
            
            # Other rows hold [i, j, value] tuples
            else
                push!(I, parse(Int, line[1]))
                push!(J, parse(Int, line[2]))
                push!(V, parse(Float64, line[3]))
            end
        end
        
        # Create custom SparseMatrice type
        mat = SparseMatrice(size, sub_size, sparse(I, J, V))
        println("Loaded matrice with size $size and sub_size $(sub_size).")

        return mat
    end
end

end # end matrice