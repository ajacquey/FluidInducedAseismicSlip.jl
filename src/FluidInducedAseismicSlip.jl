module FluidInducedAseismicSlip

using FastGaussQuadrature
using SpecialFunctions
using LinearAlgebra

export injection_analytical_gs
export lambda_analytical_gs

"""
    injection_analytical_gs(T::Float64, N::Int64 = 100)

Solves the slip distribution using Gauss-Chebyshev quadratures due to an
analytical pressure expression.
Inputs:
    T: the stress parameter
    N: number of quadrature points
Outputs:
    x: the Gauss-Checbyshev quadrature points
    δ: the slip distribution
    λ: the slip to fluid migration factor
"""
function injection_analytical_gs(T::Float64, N::Int64 = 100)
    # Solve for λ
    (λ, k) = lambda_analytical_gs(T, N)

    # Slip
    (x, δ) = slip_distribution_gs(λ, N)

    return (x, δ, λ)
end

"""
    lambda_analytical_gs(
        T::Float64,
        N::Int64,
        max_iters::Int64 = 200,
        abs_tol::Float64 = 1.0e-10,
        debug::Bool = False,
    )

Use Newton-Raphson iterations to solve for λ using Gauss-Chebyshev quadratures.
Inputs:
    T: the stress parameter
    N: number of quadrature points
    max_iters: the maximum number of Newton-Raphson iterations
    abs_tol: the absolute tolerance for the Newton solve
    debug: print convergence information
Outputs:
    λ: the slip to fluid migration factor
    k: the number of iterations used
"""
function lambda_analytical_gs(
    T::Float64,
    N::Int64,
    max_iters::Int64 = 200,
    abs_tol::Float64 = 1.0e-10,
    debug::Bool = false,
)
    # Initialization
    λ = 1.0
    k = 0
    # Create Gauss-Chebyshev quadrature points
    s, w = gausschebyshev(N, 1)

    while (k < max_iters)
        # Update residuals and jacobian
        Res = residual_analytical(T, s, w, λ)
        Jac = jacobian_analytical(T, s, w, λ)

        # Check convergence
        if (abs(Res) <= abs_tol)
            if (debug)
                println("Solve converged! λ = ", λ, " after ", k, " iterations.")
            end
            return (λ, k)
        end

        # Update λ
        dλ = -Res / Jac
        λ += dλ
        k += 1
    end

    throw(ErrorException("Maximm number of iterations reached in the Newton solve!!!\n"))
end

"""
    residual_analytical(T::Float64, s::Vector{Float64}, w::Vector{Float64}, λ::Float64)

The residual for calculating λ
"""
function residual_analytical(T::Float64, s::Vector{Float64}, w::Vector{Float64}, λ::Float64)
    f(s) = erfc(λ * abs(s))
    return dot(w, f.(s)) / pi - T
end

"""
    jacobian_analytical(T::Float64, s::Vector{Float64}, w::Vector{Float64}, λ::Float64)

The jacobian for calculating λ
"""
function jacobian_analytical(T::Float64, s::Vector{Float64}, w::Vector{Float64}, λ::Float64)
    f(s) = -2 * abs(s) / sqrt(pi) * exp(-λ^2 * abs(s)^2)
    return dot(w, f.(s)) / pi
end

"""
    dF(s::Float64, u::Float64, λ::Float64)

The function dF from Viesca and Garagash (2018)
"""
function dF(s::Float64, u::Float64, λ::Float64)
    return 1  / π * erfc(λ * abs(s)) / (u - s)
end

"""
    F(u::Float64, λ::Float64, N::Int64)

The function F from Viesca and Garagash (2018)
"""
function F(u::Float64, λ::Float64, N::Int64)
    # Gauss-Chebyshev quadrature points
    s, w = gausschebyshev(N, 1)
    return dot(w, dF.(s, u, λ))
end

"""
    slip_weigth(x::Float64, s::Vector{Float64}, N::Int64)

The slip weigth S_ij from Viesca and Garagash (2018)
"""
function slip_weigth(x::Float64, s::Vector{Float64}, N::Int64)
    θ = theta(x)
    Φ = [0.5 * (sin(k * θ) / k - sin((k + 2) * θ) / (k + 2)) for k = N:-1:1]
    B = [2 * sin(theta(sj)) * sin((k+1)*theta(sj)) / (N+1) for sj in s, k = N:-1:1]
    return B * Φ
end

"""
    slip_gs(x::Float64, N::Int64, λ::Float64)

Slip evaluated with Gauss-Chebyshev quadrature
"""
function slip_gs(x::Float64, N::Int64, λ::Float64)
    # Gauss-Chebyshev quadrature points
    (s, w) = gausschebyshev(N-1, 2)
    S = slip_weigth(x, s, N)

    return dot(S, F.(s, λ, N))
end

"""
    theta(x::Float64)

Theta function from Viesca and Garagash (2018)
"""
function theta(x::Float64)
    return acos(x::Float64)
end

"""
    slip_distribution_gs(λ::Float64, N::Int64)

Computes the slip distribution based on Gauss-Chebyshev quadratures
Inputs:
    λ: the slip to fluid migration factor
    N: number of quadrature points
Outputs
    x: the Gauss-Checbyshev quadrature points
    δ: the slip distribution
    
"""
function slip_distribution_gs(λ::Float64, N::Int64)
    x = range(-1.0, 1.0, length = N)

    δ = slip_gs.(x, N, λ)
    return (x, δ)
end

end
