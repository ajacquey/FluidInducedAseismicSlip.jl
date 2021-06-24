import Pkg
Pkg.activate("/home/ajacquey/projects/FluidInducedAseismicSlip/")
using Revise
using FluidInducedAseismicSlip

injection_analytical(0.5, 500)
