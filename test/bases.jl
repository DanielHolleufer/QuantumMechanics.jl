using QuantumMechanics
using Test

@testset "bases" begin
    @test length(GenericBasis(4)) == 4

    @testset "FockBasis" begin
        @test length(FockBasis(10)) == 11
        @test length(FockBasis(10, 3)) == 8
    end
end
