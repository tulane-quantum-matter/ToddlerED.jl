# ToddlerED

[![Build Status](https://github.com/tulane-quantum-matter/ToddlerED.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tulane-quantum-matter/ToddlerED.jl/actions/workflows/CI.yml?query=branch%3Amain)

Simple code for exact-diagonalization calculations with some number of (majorana) fermion and spin-1/2 operators.

Operators built up from tensor products of Pauli / Clifford algebra matrices.

Immediate predecessor is baby-ED by J,S-K.

## Installation

```julia
using Pkg;
Pkg.add(url="https://github.com/tulane-quantum-matter/ToddlerED.jl")
```

## Example

Solving the t-V model in 1D chain.
$$ H = -t \sum_i c_i^\dagger c_{i+1} + V (n_i-1/2) (n_{i+1}-1/2) $$

```julia
# t-V model in 1D with periodic boundary condition. Solve and compute RG invariant using staggered density-density correlator.

using LinearAlgebra
using SparseArrays
using ToddlerED
using Arpack # gives eigs
using FFTW

# load accelerator libraries
let cpu_info_list = Sys.cpu_info()
    for cpu in cpu_info_list
        if occursin("Intel", cpu.model)
            try
                using MKL, MKLSparse
                break
            catch
            end
        elseif occursin("Apple M", cpu.model)
            try
                using AppleAccelerate
                break
            catch
            end
        end
    end
end

lx = 10; # lattice size. 4n+2 is non-degenerate and admits well-defined "staggering"

# fermion operators
# Define creation and annihilation operators
c = [0.5 * dirac(2i - 1, 2lx) + 0.5im * dirac(2i, 2lx) for i in 1:lx]
cdg = [0.5 * dirac(2i - 1, 2lx) - 0.5im * dirac(2i, 2lx) for i in 1:lx]

# Define number operators
n = Diagonal.([real.(cdg[i] * c[i]) for i in 1:lx])
dn = Diagonal.([n[i] - 0.5I for i in 1:lx])
deln = Diagonal.([n[i] - n[mod1(i+1,lx)] for i in 1:lx])
## t-V model

# Define hopping term
t = 1.0
tp = sum([cdg[i] * c[i+1] for i in 1:lx-1]) + cdg[lx] * c[1]
tp = real.(tp + tp')

# Define interaction term
vp = sum([dn[i] * dn[i+1] for i in 1:lx-1]) + dn[lx] * dn[1]

# Define density operators
nn = [(-1)^(1 + i) * dn[1] * dn[i] for i in 1:lx]
delnn = [(-1)^(1 + i) * deln[1] * deln[i] for i in 1:lx]

# Define Hamiltonian as a function of interaction strength v
h(v) = Symmetric(real.(-tp + v * vp))

vs = 0:0.2:3;
result = [];
for v in vs # there 
    (ee, ev) = eigs(h(v), nev=2, which=:SR) # SR for smallest real part
    ee = ee[1:2]
    ev = ev[:, 1:2]
    nns = [dot(ev, nn[i], ev) for i in 1:lx]
    nnsk = fft(nns)
    push!(result, [v, real(nnsk[1] ./ nnsk[2]) - 1, nns, ee, ev])
end
result

```