module FluidInducedAseismicSlip

using Plots: length, size
using FastGaussQuadrature
using SpecialFunctions
using LinearAlgebra
using Plots
using Revise

export injection_analytical_gs
export lambda_analytical_gs
export discretizeFaultSegment, collocationBEMMatrix

include("bem.jl")

"""
    injection_analytical_gs(T::Float64, N::Int64 = 100)

Solves the slip distribution using Gauss-Chebyshev quadratures due to an
analytical pressure expression
"""
function injection_analytical_gs(T::Float64, N::Int64 = 100)
    # Solve for λ
    λ, k = lambda_analytical_gs(T, N)

    # Slip
    x, δ = slip_distribution_gs(λ, N)

    # # Plot slip gradient
    # plot(x, δ)
    # x_a = collect(0:0.05:1.0)
    # δ_a = [0.23366202, 0.23081724, 0.22426165, 0.21513089, 0.20404685, 0.19145812, 0.17772295, 0.16314401, 0.14798653, 0.13248893, 0.11686999, 0.101334019, 0.086075209, 0.071281873, 0.057141201, 0.043845559, 0.031602114, 0.020649782, 0.0112944788, 0.0040083095, 0.0]
    # plot!(x_a, δ_a)
    return x, δ, λ
end

"""
    lambda_analytical_gs(
        T::Float64,
        N::Int64,
        max_iters::Int64 = 200,
        abs_tol::Float64 = 1.0e-10,
        debug::Bool = False,
    )

Use Newton-Raphson iterations to solve for λ using Gauss-Chebyshev quadratures
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
                print("Solve converged! λ = ", λ, " after ", k, " iterations.\n")
            end
            return λ, k
        end

        # Update λ
        dλ = -Res / Jac
        λ += dλ
        k += 1
    end

    print("Maximm number of iterations reached in the Newton solve!!!\n")
end

"""
    residual_analytical(T, s, w, λ)

The residual for calculating λ
"""
function residual_analytical(T, s, w, λ)
    f(s) = erfc(λ * abs(s))
    return dot(w, f.(s)) / pi - T
end

"""
    jacobian_analytical(T, s, w, λ)

The jacobian for calculating λ
"""
function jacobian_analytical(T, s, w, λ)
    f(s) = -2 * abs(s) / sqrt(pi) * exp(-λ^2 * abs(s)^2)
    return dot(w, f.(s)) / pi
end

"""
    dF(s, u, λ)

The function dF from Viesca and Garagash (2018)
"""
function dF(s, u, λ)
    return 1  / π * erfc(λ * abs(s)) / (u - s)
end

"""
    F(u, λ, N)

The function F from Viesca and Garagash (2018)
"""
function F(u, λ, N)
    # Gauss-Chebyshev quadrature points
    s, w = gausschebyshev(N, 1)
    return dot(w, dF.(s, u, λ))
end

"""
    slip_weigth(x, s, N)

The slip weigth S_ij from Viesca and Garagash (2018)
"""
function slip_weigth(x, s, N)
    θ = theta(x)
    Φ = [0.5 * (sin(k * θ) / k - sin((k + 2) * θ) / (k + 2)) for k = N:-1:1]
    B = [2 * sin(theta(sj)) * sin((k+1)*theta(sj)) / (N+1) for sj in s, k = N:-1:1]
    return B * Φ
end

"""
    slip_gs(u, N, λ)

Slip evaluated with Gauss-Chebyshev quadrature
"""
function slip_gs(x, N, λ)
    # Gauss-Chebyshev quadrature points
    s, w = gausschebyshev(N-1, 2)
    S = slip_weigth(x, s, N)

    return dot(S, F.(s, λ, N))
end

"""
    theta(x)

Theta function from Viesca and Garagash (2018)
"""
function theta(x)
    return acos(x)
end

"""
    slip_gradient_distribution(λ, N)

Computes the slip gradient distribution based on Gauss-Chebyshev quadrature
"""
function slip_gradient_distribution(λ, N)
    dδ = zeros(N-1)
    # Gauss-Chebyshev quadrature points
    x, w = gausschebyshev(N-1, 2)

    for i in 1:N-1
        dδ[i] = slip_gradient(x[i], N, λ)
    end
    return x, dδ
end

"""
    slip_distribution_gs(λ, N)

Computes the slip distribution based on Gauss-Chebyshev quadratures
"""
function slip_distribution_gs(λ, N)
    δ = zeros(N)
    x = range(-1, 1, length = N)

    for i in 1:N
        δ[i] = slip_gs(x[i], N, λ)
    end
    return x, δ
end

end
