# Patryk Rygiel
module Matrice

using SparseArrays

export load_matrice, load_vector

struct RowMatrice
    size::Integer
    sub_size::Integer
    shift::Vector
    data::Array
end

function load_vector(file_path::String)
    """
    Loads side vector from given file
    """
    open(file_path) do file
    
        data = Float64[]

        for (idx, line) in enumerate(eachline(file))
            
            # First row in file holds size - omit
            if idx != 1
                push!(data, parse(Float64, line))
            end
        end
        
        # Create custom SideVector type
        vector = vec(data)
        println("Loaded vector with size $(length(vector)).")

        return vector
    end
end

function load_sparse_matrice(file_path::String)
    """
    Loads matrice from file as Julia SparseMatrixCSC
    """
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
        mat = sparse(I, J, V)
        println("Loaded matrice with size $size and sub_size $(sub_size).")

        return mat
    end
end


function load_row_matrice(file_path::String)
    """
    Loads matrice from file as custom data structure
    """
    open(file_path) do file
        
        # First row in file holds size and sub_size
        size, sub_size = split(readline(file))
        
        size = parse(Int64, size)
        sub_size = parse(Int64, sub_size)
        v = size ÷ sub_size

        # Construct row Vectors
        data = []
        shift = Vector{Int64}(undef, size)

        for i in 1:size
            if i <= sub_size
                dim = 3 * sub_size
            elseif i > size - sub_size
                dim = sub_size + 1
            else
                dim = sub_size * 3 + 1
            end
            
            vector = vec(zeros(dim))
            push!(data, vector)

            shift[i] = max(0, sub_size * ((i-1) ÷ sub_size) - 1)
        end

        # Other rows hold [i, j, value] tuples
        for (idx, line) in enumerate(eachline(file))
            
            if idx != 0
                line = split(line)

                i = parse(Int64, line[1])
                j = parse(Int64, line[2])
                v = parse(Float64, line[3])
                
                data[i][j - shift[i]] = v
            end
        end

        # Create custom Matrice type
        mat = RowMatrice(size, sub_size, shift, data)
        println("Loaded matrice with size $size and sub_size $(sub_size).")

        return mat
    end
end

end # end matrice