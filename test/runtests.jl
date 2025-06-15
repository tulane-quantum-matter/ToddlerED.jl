using ToddlerED
using Test

@testset "ToddlerED.jl" begin
    function test_dirac()
        nm = 4
        nspin = 2
        for i in 1:nm
            for j in 1:nm
                println("$(i), $(j): ", anticomm(dirac(i, nm, nspin), dirac(j, nm, nspin)))
            end
        end
        for i in 1:nspin
            for j in 1:nspin
                println("$(i), $(j): ", anticomm(spinx(i, nm, nspin), spinx(j, nm, nspin)))
                println("$(i), $(j): ", anticomm(spiny(i, nm, nspin), spinx(j, nm, nspin)))
            end
        end
    end
end
