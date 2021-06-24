using FluidInducedAseismicSlip
using Test

@testset "FluidInducedAseismicSlip.jl" begin
    # Write your tests here.
    FluidInducedAseismicSlip.injection_analytical(0.5, 100)
end
