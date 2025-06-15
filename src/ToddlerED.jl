module ToddlerED

# Import necessary libraries
using LinearAlgebra
using SparseArrays
using Arpack # gives eigs

export dirac, spinop, spinx, spiny, spinz, spin0, anticomm, comm
export load_accelerated_libraries
export fermion_creation, fermion_annihilation

"""
@load_accelerated_libraries()
Return code to check the CPU model and load appropriate libraries for optimized performance.
"""
macro load_accelerated_libraries()
    quote
        let cpu_info_list = Sys.cpu_info()
            for cpu in cpu_info_list
                if occursin("Intel", cpu.model)
                    using MKL, MKLSparse
                    break
                elseif occursin("Apple M", cpu.model)
                    using AppleAccelerate
                    break
                end
            end
        end
    end
end

# Define Pauli matrices as sparse matrices
const σ₀ = sparse(ComplexF64[1 0; 0 1])
const σ₁ = sparse(ComplexF64[0 1; 1 0])
const σ₂ = sparse(ComplexF64[0 -im; im 0])
const σ₃ = sparse(ComplexF64[1 0; 0 -1])

"""
    dirac(i::Int, n::Int, nspin::Int=0)

Returns the Dirac matrix representing `i`-th Majorana in a space of with a total of `n`-Majoranas. `n` needs to be even, unless `n`=1.
{dirac(i,n), dirac(j,n)} = 2 δᵢⱼ
    For odd i, dirac(i,n) = σ₃ ⊗ σ₃ ⊗ ... ⊗ σ₃ ⊗ σ₁ ⊗ σ₀ ⊗ ... ⊗ σ₀
    For even i, dirac(i,n) = σ₃ ⊗ σ₃ ⊗ ... ⊗ σ₃ ⊗ σ₂ ⊗ σ₀ ⊗ ... ⊗ σ₀
They are all Hermitian and unitary, although not explictly typed so.
`nspin` is the number of spin-1/2's in the total hilbert space. Each spin has twice the quantum dimension of a Majorana.

This choice of basis makes number operators of complex fermions diagonal. Eg:
```julia
julia> c = [0.5 * dirac(2i - 1, 2lx) + 0.5im * dirac(2i, 2lx) for i in 1:lx]
julia> cdg = [0.5 * dirac(2i - 1, 2lx) - 0.5im * dirac(2i, 2lx) for i in 1:lx]
"""
function dirac(i::Int, nm::Int, nspin::Int=0)

    if i > nm
        error("Majorana index i=$i is greater than the total number of Majoranas nm=$nm.")
    end
    if nm == 1 && i == 1 && nspin == 0
        return σ₂
    elseif isodd(nm)
        error("Only even n or (n=1 & nspin=0) allowed.\nReceived nm=$(nm), nspin=$(nspin).")
    end

    matrices = [σ₃ for _ in 2:2:i-1]
    push!(matrices, iseven(i) ? σ₁ : σ₂)
    append!(matrices, [σ₀ for _ in i+1:2:nm-1+2nspin])
    return reduce(kron, matrices)
end

# Just to give an example of complex fermions
fermion_creation(i,nm) = 0.5 * dirac(2i - 1, nm) + 0.5im * dirac(2i, nm)
fermion_annihilation(i,nm) = 0.5 * dirac(2i - 1, nm) - 0.5im * dirac(2i, nm)

"""
    spinop(s::Int, i::Int, nspin::Int, nm::Int=0)
spinop(s,i) = σ₀ ⊗ σ₀ ⊗ ... ⊗ σ₀ ⊗ σₛ ⊗ σ₀ ⊗ ... ⊗ σ₀
s = 1,2,3 for σ₁,σ₂,σ₃, and 0 or 4 for σ₀
"""
function spinop(s::Int, i::Int, nspin::Int, nm::Int=0)
    if i > nspin
        error("Spin index i=$i is greater than the total number of spins nspin=$nm.")
    end

    ss = [σ₁, σ₂, σ₃, σ₀];

    matrices = [σ₀ for _ in 1:2:nm+2(i-1)]
    push!(matrices, ss[s == 0 ? 4 : s])
    append!(matrices, [σ₀ for _ in i+1:nspin])
    return reduce(kron, matrices)
end

"spinx(i,nspin,nm=0) is shorthand for spinop(1,i,nspin,nm)"
function spinx(args...)
    spinop(1, args...)
end
"spiny(i,nspin,nm=0) is shorthand for spinop(2,i,nspin,nm)"
function spiny(args...)
    spinop(2, args...)
end
"spinz(i,nspin,nm=0) is shorthand for spinop(3,i,nspin,nm)"
function spinz(args...)
    spinop(3, args...)
end
"spin0(i,nspin,nm=0) is shorthand for spinop(0,i,nspin,nm)"
function spin0(args...)
    spinop(4, args...)
end

"""
    anticomm(x, y)
Returns the anticommutator of two operators x and y, defined as {x,y}
"""
function anticomm(x, y)
    x * y + y * x
end

"""
    comm(x, y)
Returns the commutator of two operators x and y, defined as [x,y]
"""
function comm(x, y)
    x * y - y * x
end



end
