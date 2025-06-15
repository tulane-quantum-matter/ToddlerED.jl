using ToddlerED
using LinearAlgebra, SparseArrays
using Test

# under construction
@testset "ToddlerED.jl" begin
    nspin = 3
    nm = 6
    rm = () -> mod1(rand(Int), nm)
    rs = () -> mod1(rand(Int), nspin)

    # spin operator clifford algebra
    for mm = 1:3
        @test 2I == let g = rs()
            anticomm(spinop(mm, g, nspin, nm), spinop(mm, g, nspin, nm))
        end
    end
    for mm = 1:3
        @test 0I == let g = rs()
            anticomm(spinop(mm, g, nspin, nm), spinop(mod1(mm + 1, 3), g, nspin, nm))
        end
    end

    # spin operator su2 algebra
    for mm = 1:3
        @test let g = rs()
            2im .* spinop(mm, g, nspin, nm) == comm(spinop(mod1(mm + 1, 3), g, nspin, nm), spinop(mod1(mm + 2, 3), g, nspin, nm))
        end
    end

    @test let i = rm()
        2I == anticomm(dirac(i, nm, nspin), dirac(i, nm, nspin)) &&
            2I == anticomm(dirac(mod1(i + 1, nm), nm, nspin), dirac(mod1(i + 1, nm), nm, nspin))
    end

    # commutator between spin operators
    for mm = 1:4
        @test let i = rs()
            0I == comm(spinop(1, i, nspin, nm), spinop(1, mod1(i + 1, nspin), nspin, nm))
        end
    end

    # anticommutator between majorana
    for mm = 1:4
        @test let i = rm(), j = rm()
            (i == j ? 2 : 0)I == anticomm(dirac(i, nm, nspin), dirac(j, nm, nspin))
        end
    end


    # commutator between spin and majorana
    for mm = 1:4
        @test let i = rm(), j = rs()
            0I == comm(dirac(i, nm, nspin), spinop(1, j, nspin, nm))
        end
    end
end
