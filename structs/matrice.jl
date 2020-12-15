# Patryk Rygiel
module Matrice

using SparseArrays

export load_matrice, load_vector

struct SparseMatrice
    size::Integer
    sub_size::Integer
    data::SparseMatrixCSC
end

function load_vector(file_path::String)
    open(file_path) do file
        
        size = 0
        data = Float64[]

        for (idx, line) in enumerate(eachline(file))
            
            # First row in file holds size - omit
            if idx != 1
                push!(data, parse(Float64, line))
            end
        end
        
        # Create custom SideVector type
        vector = vec(data)
        println("Loaded vector with size $size.")

        return vector
    end
end

function load_matrice(file_path::String)
    open(file_path) do file
        
        size = 0
        sub_size = 0
        
        # Vectors to build a sparse array from
        I = Int64[]
        J = Int64[]
        V = Float64[]
        
        for (idx, line) in enumerate(eachline(file))
            line = split(line)
            
            # First row in file holds size and sub_size
            if idx == 1
                size = parse(Int64, line[1])
                sub_size = parse(Int64,line[2])
            
            # Other rows hold [i, j, value] tuples
            else
                push!(I, parse(Int64, line[1]))
                push!(J, parse(Int64, line[2]))
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