# Patryk Rygiel
module Matrice

using SparseArrays

export load_matrice, load_vector

struct CustomMatrice
    size::Integer
    sub_size::Integer
    A::Array
    B::Array
    C::Array
end
struct SparseMatrice
    size::Integer
    sub_size::Integer
    data::SparseMatrixCSC
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

function load_matrice(file_path::String)
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

        # Insert zeros where elements will be inserted
        # C matrices under diagonal
        
        # c1
        # 0 c2
        # 0 0 c3
        # 0 0 0 c4
        
        for row in 1:size-4
            if row % sub_size != 1

                start_col = sub_size * (row ÷ (sub_size + 1) + 1) + 1
                end_col = sub_size + row - 1

                for col in start_col:end_col
                    push!(I, row)
                    push!(J, col)
                    push!(V, 0.0)
                end
            end
        end

        # Create custom SparseMatrice type
        mat = SparseMatrice(size, sub_size, sparse(I, J, V))
        println("Loaded matrice with size $size and sub_size $(sub_size).")

        return mat
    end
end


function load_custom_matrice(file_path::String)
    """
    Loads matrice from file as custom data structure
    """
    open(file_path) do file
        
        # First row in file holds size and sub_size
        size, sub_size = split(readline(file))
        
        size = parse(Int64, size)
        sub_size = parse(Int64, sub_size)
        v = size ÷ sub_size

        # Allocate matrice
        A = zeros((v, sub_size, sub_size))
        B = zeros((v-1, sub_size))
        C = zeros((v-1, sub_size, sub_size))

        # Other rows hold [i, j, value] tuples
        for (idx, line) in enumerate(eachline(file))
            
            if idx != 0
                line = split(line)

                i = parse(Int64, line[1])
                j = parse(Int64, line[2])
                v = parse(Float64, line[3])

                indice = ((i-1) ÷ sub_size) + 1
                
                # Filling matrice A when value on diagonal (sub_size, sub_size)
                if abs(i - j) < sub_size
                    A[indice, (i-1) % sub_size + 1, (j-1) % sub_size + 1] = v
                
                # Filling matrice B when value under diagonal
                elseif i >= j
                    B[indice - 1, (j-1) % sub_size + 1] = v

                # Filling matrice C when value above diagonal
                else
                    C[indice, (i-1) % sub_size + 1, (i-1) % sub_size + 1] = v
                end

            end
        end

        # Create custom Matrice type
        mat = CustomMatrice(size, sub_size, A, B, C)
        println("Loaded matrice with size $size and sub_size $(sub_size).")

        return mat
    end
end

end # end matrice