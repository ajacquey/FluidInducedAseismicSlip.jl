using FluidInducedAseismicSlip

N = 500
L = 1.0
x = discretizeFaultSegment(L, N)
τ₀ = zeros(Real, N, 1)
σ₀ = zeros(Real, N, 1)
bem_problem = initializeProblem(x, τ₀, σ₀)
E = collocationBEMMatrix(x, 1.0, 2.0)