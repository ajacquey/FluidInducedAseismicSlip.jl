using FluidInducedAseismicSlip
using Plots

function test_slip()
    (x, δ, λ) = injection_analytical_gs(0.5, 500)

    # Analytical slip
    x_a = collect(0:0.05:1.0)
    δ_a = [0.23366202, 0.23081724, 0.22426165, 0.21513089, 0.20404685, 0.19145812, 0.17772295, 0.16314401, 0.14798653, 0.13248893, 0.11686999, 0.101334019, 0.086075209, 0.071281873, 0.057141201, 0.043845559, 0.031602114, 0.020649782, 0.0112944788, 0.0040083095, 0.0]

    # Plot
    plot(x, δ, label="This package")
    plot!(x_a, δ_a, label="Viesca (2021)")
end

test_slip()
