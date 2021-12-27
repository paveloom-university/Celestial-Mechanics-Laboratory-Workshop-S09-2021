# This script finds a figure-8 solution to a
# three-body problem with arbitrary precision

println('\n', " "^4, "> Loading the packages...")

using LinearAlgebra # Norm
using Plots # Plotting

# Use the GR backend for plots
gr()

# Change the default font for plots
default(fontfamily = "Computer Modern", dpi = 300, legend = :topright)

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")

# Prepare the output directory
mkpath(PLOTS_DIR)

# Define the mass of each body (sum is 1)
const m = 1 / 3

# Define the value of ϰ = -Gm (G = 1)
const ϰ = -m

"""
Integrate the three-body problem using the
4th-order Runge-Kutta's method, return the
values of position and velocity on each step
"""
function rk4(
    r₀₁::Vector{F},
    r₀₂::Vector{F},
    r₀₃::Vector{F},
    v₀₁::Vector{F},
    v₀₂::Vector{F},
    v₀₃::Vector{F},
    h::F,
    n::I,
)::Tuple{Matrix{F},Matrix{F},Matrix{F},Matrix{F},Matrix{F},Matrix{F}} where
{F<:AbstractFloat,I<:Unsigned}
    # Determine the length of the input vectors
    N = length(r₀₁)
    # Prepare output matrixes
    r₁ = Matrix{F}(undef, n + 1, N)
    r₂ = Matrix{F}(undef, n + 1, N)
    r₃ = Matrix{F}(undef, n + 1, N)
    v₁ = Matrix{F}(undef, n + 1, N)
    v₂ = Matrix{F}(undef, n + 1, N)
    v₃ = Matrix{F}(undef, n + 1, N)
    r₁[1, :] = copy(r₀₁)
    r₂[1, :] = copy(r₀₂)
    r₃[1, :] = copy(r₀₃)
    v₁[1, :] = copy(v₀₁)
    v₂[1, :] = copy(v₀₂)
    v₃[1, :] = copy(v₀₃)
    # Compute the increments to the position
    function k(i, v₁₊, v₂₊, v₃₊)
        k₁ = v₁[i-1, :] + v₁₊
        k₂ = v₂[i-1, :] + v₂₊
        k₃ = v₃[i-1, :] + v₃₊
        return k₁, k₂, k₃
    end
    # Compute the increments to the velocity
    function a(i, r₁₊, r₂₊, r₃₊)
        r₁⁺, r₂⁺, r₃⁺ = r₁[i-1, :] + r₁₊, r₂[i-1, :] + r₂₊, r₃[i-1, :] + r₃₊
        r₁₂, r₁₃, r₂₃ = r₁⁺ - r₂⁺, r₁⁺ - r₃⁺, r₂⁺ - r₃⁺
        r₂₁, r₃₁, r₃₂ = -r₁₂, -r₁₃, -r₂₃
        ρr₁₂ = ρr₂₁ = norm(r₁₂)^3
        ρr₁₃ = ρr₃₁ = norm(r₁₃)^3
        ρr₂₃ = ρr₃₂ = norm(r₂₃)^3
        a₁ = ϰ * r₁₂ / ρr₁₂ + ϰ * r₁₃ / ρr₁₃
        a₂ = ϰ * r₂₃ / ρr₂₃ + ϰ * r₂₁ / ρr₂₁
        a₃ = ϰ * r₃₁ / ρr₃₁ + ϰ * r₃₂ / ρr₃₂
        return a₁, a₂, a₃
    end
    # Compute the solution
    for i in 2:(n+1)
        k₁₁, k₁₂, k₁₃ = k(i, ntuple(_ -> repeat([0], N), 3)...)
        a₁₁, a₁₂, a₁₃ = a(i, ntuple(_ -> repeat([0], N), 3)...)
        k₂₁, k₂₂, k₂₃ = k(i, h * a₁₁ / 2, h * a₁₂ / 2, h * a₁₃ / 2)
        a₂₁, a₂₂, a₂₃ = a(i, h * k₁₁ / 2, h * k₁₂ / 2, h * k₁₃ / 2)
        k₃₁, k₃₂, k₃₃ = k(i, h * a₂₁ / 2, h * a₂₂ / 2, h * a₂₃ / 2)
        a₃₁, a₃₂, a₃₃ = a(i, h * k₂₁ / 2, h * k₂₂ / 2, h * k₂₃ / 2)
        k₄₁, k₄₂, k₄₃ = k(i, h * a₃₁, h * a₃₂, h * a₃₃)
        a₄₁, a₄₂, a₄₃ = a(i, h * k₃₁, h * k₃₂, h * k₃₃)
        r₁[i, :] = r₁[i-1, :] + h / 6 * (k₁₁ + 2 * k₂₁ + 2 * k₃₁ + k₄₁)
        r₂[i, :] = r₂[i-1, :] + h / 6 * (k₁₂ + 2 * k₂₂ + 2 * k₃₂ + k₄₂)
        r₃[i, :] = r₃[i-1, :] + h / 6 * (k₁₃ + 2 * k₂₃ + 2 * k₃₃ + k₄₃)
        v₁[i, :] = v₁[i-1, :] + h / 6 * (a₁₁ + 2 * a₂₁ + 2 * a₃₁ + a₄₁)
        v₂[i, :] = v₂[i-1, :] + h / 6 * (a₁₂ + 2 * a₂₂ + 2 * a₃₂ + a₄₂)
        v₃[i, :] = v₃[i-1, :] + h / 6 * (a₁₃ + 2 * a₂₃ + 2 * a₃₃ + a₄₃)
    end
    return r₁, r₂, r₃, v₁, v₂, v₃
end

println(" "^4, "> Integrating the three-body problem...")

# Set the precision (in bits)
setprecision(BigFloat, 128)

# Define the initial values of the position and velocity
v₀₂ = BigFloat.(["0.7494421910777922289898659", "1.1501789857502275024030202"])
v₀₁ = [-0.5 * v₀₂[1], -0.5 * v₀₂[2]]
v₀₃ = [-0.5 * v₀₂[1], -0.5 * v₀₂[2]]
r₀₃ = [-5 * m^2 / (1 + m * (norm(v₀₁)^2 + norm(v₀₂)^2 + norm(v₀₃)^2)), 0.0]
r₀₁ = [-r₀₃[1], 0.0]
r₀₂ = BigFloat.(["0.0", "0.0"])

# Define the time step and number of iterations
h = BigFloat("0.0001")
n = ceil(UInt, BigFloat("1.676118923759279") / h / 3)

# Compute the energy constant
E = m / 2 * (norm(v₀₁)^2 + norm(v₀₂)^2 + norm(v₀₃)^2) -
    m^2 * (1 / norm(r₀₁ - r₀₂) + 1 / norm(r₀₁ - r₀₃) + 1 / norm(r₀₂ - r₀₃))

println(" "^4, "  > E: ", E)

# Integrate the problem
r₁, r₂, r₃, _, _, _ = rk4(r₀₁, r₀₂, r₀₃, v₀₁, v₀₂, v₀₃, h, n)

println(" "^4, "> Plotting the figure...")

# Plot the solution
p = plot(r₁[:, 1], r₁[:, 2]; label = "1", legend_position = :topleft)
plot!(p, r₂[:, 1], r₂[:, 2]; label = "2")
plot!(p, r₃[:, 1], r₃[:, 2]; label = "3")

# Save the figure
savefig(p, joinpath(PLOTS_DIR, "solution.png"))
savefig(p, joinpath(PLOTS_DIR, "solution.pdf"))

println()
