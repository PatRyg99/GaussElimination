push!(LOAD_PATH, pwd())

include("./structs/matrice.jl")
using .Matrice: load_matrice, load_vector

matrice = load_matrice("data/Dane16/A.txt")
vector = load_vector("data/Dane16/b.txt")

println(matrice)
println(vector)