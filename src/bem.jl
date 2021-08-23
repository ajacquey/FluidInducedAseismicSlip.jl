"""
    discretizeFaultSegment(L::Real, N::Int)

Returns a discretized fault segment of length L with N elements
Inputs:
    L: length of the fault segment
    N: number of elements
Outputs:
    x: vector of elements centers
"""
function discretizeFaultSegment(L::Real, N::Int)
    # grid spacing
    dx = L / N
    # returns a vector centered around 0
    return -L / 2 .+ [1:N;] * dx .- dx / 2
end

"""
    collocationBEMMatrix(x::Array{T,1}, μ::Real, λ::Real) where T<:Real

Returns a global collocation matrix for boundary elements for a 2D problem
Inputs:
    x: vector of elements centers
    μ: effective shear modulus
    λ: effective normal modulus
Outputs:
    E: the global (2N x 2N) collocation BEM matrix
"""
function collocationBEMMatrix(x::Array{T,1}, μ::Real, λ::Real) where T<:Real
    # number of elements
    N = size(x)[1]
    # initialize global collocation BEM matrix
    E = zeros(Real, 2*N, 2*N)
    # grid spacing (not necessarily constant)
    Δx = x[2] - x[1]
    # shear and normal collocation matrices
    Eₛ = reshape([μ ./ (Δx * π) / ((i - j)^2 - 0.25) for i = 1:N for j = 1:N], N, N)
    Eₙ = reshape([λ ./ (Δx * π) / ((i - j)^2 - 0.25) for i = 1:N for j = 1:N], N, N)
    # global assembly
    E = globalBEMMatrixAssembly(Eₛ, Eₙ)
    return E
end

"""
    globalBEMMatrixAssembly(Eₛ::Matrix{T}, Eₙ::Matrix{T}) where T<:Real

Assembly of the global BEM matrix for 2D problems
Inputs:
    Eₛ: the global (N x N) collocation BEM matrix for shear components
    Eₙ: the global (N x N) collocation BEM matrix for normal components
Outputs:
    E: the global (2N x 2N) collocation BEM matrix
"""
function globalBEMMatrixAssembly(Eₛ::Matrix{T}, Eₙ::Matrix{T}) where T<:Real
    # number of elements
    N = size(Eₛ)[1]
    # initialize global collocation BEM matrix
    E = zeros(Real, 2*N, 2*N)
    # global assembly
    for i = 1:N
        for j = 1:N
            E[2*i-1, 2*j-1] = Eₛ[i,j]
            E[2*i, 2*j-1] = 0
            E[2*i-1, 2*j] = 0
            E[2*i, 2*j] = Eₙ[i,j]
        end
    end
    return E
end