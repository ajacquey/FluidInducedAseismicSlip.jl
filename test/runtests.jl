using FluidInducedAseismicSlip
using Test

@testset "FluidInducedAseismicSlip.jl" begin
    # Test the calculation of λ for different values of T
    @test lambda_analytical_gs(0.1, 500)[1] ≈ 3.6403634658646724
    @test lambda_analytical_gs(0.2, 500)[1] ≈ 1.9124952153633206
    @test lambda_analytical_gs(0.3, 500)[1] ≈ 1.3474176336951285
    @test lambda_analytical_gs(0.4, 500)[1] ≈ 1.0252408283009076
    @test lambda_analytical_gs(0.5, 500)[1] ≈ 0.7916049887046924
    @test lambda_analytical_gs(0.6, 500)[1] ≈ 0.601228108399377
    @test lambda_analytical_gs(0.7, 500)[1] ≈ 0.4351313800638411
    @test lambda_analytical_gs(0.8, 500)[1] ≈ 0.2833769094180259
    @test lambda_analytical_gs(0.9, 500)[1] ≈ 0.1398124625560765
end
