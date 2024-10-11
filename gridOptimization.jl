using JuMP
import Ipopt
using DataFrames
using CSV

# Create the model object
model = Model(Ipopt.Optimizer)

# Data
n = 9  # Number of generators
m = 11  # Number of nodes

# Generator data
max_capacity = [0.02, 0.15, 0.08, 0.07, 0.04, 0.17, 0.17, 0.26, 0.05]
cost_per_unit = [175.0, 100.0, 150.0, 150.0, 300.0, 350.0, 400.0, 300.0, 200.0]
locations = [2, 2, 2, 3, 4, 5, 7, 9, 9]

# Map generators to nodes
GeneratorsAtNode = Dict(k => [i for i in 1:n if locations[i] == k] for k in 1:m)

# Consumer data
consumer_nodes = [1, 4, 6, 8, 9, 10, 11]
active_power_demand = [0.10, 0.19, 0.11, 0.09, 0.21, 0.05, 0.04]

# Map demands to nodes
function safe_sum(iterator)
    s = zero(Float64)
    for x in iterator
        s += x
    end
    return s
end

DemandAtNode = Dict(k => safe_sum(active_power_demand[j] for j in 1:length(consumer_nodes) if consumer_nodes[j] == k) for k in 1:m)

# Voltage and phase angle limits
voltage_limits = (0.98, 1.02)
angle_limits = (-π, π)

# Edge data
edges = [(1,2), (1,11), (2,3), (2,11), (3,4), (3,9), (4,5), (5,6), (5,8), (6,7), (7,8), (7,9), (8,9), (9,10), (10,11)]
g_kl = [4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.71, 1.26, 1.11, 1.32, 2.01, 4.41, 2.14, 5.06]
b_kl = [-20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3, -9.2, -13.9, -8.7, -11.3, -7.7, -13.5, -26.7]

# Convert g_kl and b_kl to Float64
g_kl = Float64.(g_kl)
b_kl = Float64.(b_kl)

# Build admittance matrix
Y = zeros(ComplexF64, m, m)
for idx in 1:length(edges)
    (k_node, l_node) = edges[idx]
    y = Complex(g_kl[idx], b_kl[idx])
    Y[k_node,k_node] += y
    Y[l_node,l_node] += y
    Y[k_node,l_node] -= y
    Y[l_node,k_node] -= y
end

G = real(Y)
B = imag(Y)

# Variables without bounds
@variable(model, A[i=1:n])  # Active power produced by each generator
@variable(model, R[i=1:n])  # Reactive power produced

# Upper and lower bounds as constraints for A[i]
A_upper_bounds = @constraint(model, [i=1:n], A[i] <= max_capacity[i])
A_lower_bounds = @constraint(model, [i=1:n], A[i] >= 0)

# Upper and lower bounds as constraints for R[i]
R_upper_bounds = @constraint(model, [i=1:n], R[i] <= 0.03 * max_capacity[i])
R_lower_bounds = @constraint(model, [i=1:n], R[i] >= -0.03 * max_capacity[i])

# Voltage variables
@variable(model, voltage_limits[1] <= v[i=1:m] <= voltage_limits[2], start=1.0)  # Voltage magnitudes
@variable(model, angle_limits[1] <= θ[i=1:m] <= angle_limits[2], start=0.0)  # Voltage angles

# Objective: Minimize total production cost
@objective(model, Min, sum(cost_per_unit[i] * A[i] for i in 1:n))

# Expressions for active and reactive power injections at each node
@expression(model, P_k[k=1:m], 
    sum(v[k] * v[l] * (G[k,l] * cos(θ[k] - θ[l]) + B[k,l] * sin(θ[k] - θ[l])) for l in 1:m)
)

@expression(model, Q_k[k=1:m], 
    sum(v[k] * v[l] * (G[k,l] * sin(θ[k] - θ[l]) - B[k,l] * cos(θ[k] - θ[l])) for l in 1:m)
)

# Power balance constraints
@constraint(model, [k=1:m], P_k[k] == safe_sum(A[i] for i in GeneratorsAtNode[k]) - get(DemandAtNode, k, 0.0))
@constraint(model, [k=1:m], Q_k[k] == safe_sum(R[i] for i in GeneratorsAtNode[k]))

# Solve model
optimize!(model)

# Check solver status
if termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL
    # --- Generator Results ---
    generator_data = DataFrame(
        Generator = ["G$i" for i in 1:n],
        A_pu = [round(value(A[i]), digits=6) for i in 1:n],
        R_pu = [round(value(R[i]), digits=6) for i in 1:n],
        M_pu = max_capacity,
        h_SEK_pu = cost_per_unit
    )
    CSV.write("generator_results.csv", generator_data)
    
    # --- Node Voltage Results ---
    node_voltage_data = DataFrame(
        Node = 1:m,
        v_pu = [round(value(v[k]), digits=6) for k in 1:m],
        θ_rad = [round(value(θ[k]), digits=6) for k in 1:m]
    )
    CSV.write("node_voltage_results.csv", node_voltage_data)
    
    # --- Power Flow Results ---
    power_flow_rows = []
    for idx in 1:length(edges)
        (k, l) = edges[idx]
        p_kl = round(value(v[k]^2 * G[k,l] - v[k]*v[l]*(G[k,l]*cos(value(θ[k]-θ[l])) + B[k,l]*sin(value(θ[k]-θ[l])))), digits=6)
        p_lk = round(value(v[l]^2 * G[l,k] - v[l]*v[k]*(G[l,k]*cos(value(θ[l]-θ[k])) + B[l,k]*sin(value(θ[l]-θ[k])))), digits=6)
        q_kl = round(value(-v[k]^2 * B[k,l] - v[k]*v[l]*(G[k,l]*sin(value(θ[k]-θ[l])) - B[k,l]*cos(value(θ[k]-θ[l])))), digits=6)
        q_lk = round(value(-v[l]^2 * B[l,k] - v[l]*v[k]*(G[l,k]*sin(value(θ[l]-θ[k])) - B[l,k]*cos(value(θ[l]-θ[k])))), digits=6)
        push!(power_flow_rows, (Edge = "($k,$l)", p_kl_pu = p_kl, p_lk_pu = p_lk, q_kl_pu = q_kl, q_lk_pu = q_lk))
    end
    power_flow_data = DataFrame(power_flow_rows)
    CSV.write("power_flow_results.csv", power_flow_data)
    
    # --- Lagrange Multipliers for Generator Capacity Constraints ---
    lagrange_rows = []
    for i in 1:n
        mu_upper = dual(A_upper_bounds[i])
        mu_lower = dual(A_lower_bounds[i])
        push!(lagrange_rows, (Generator = "G$i", μ_upper_SEK_pu = round(mu_upper, digits=4), μ_lower_SEK_pu = round(mu_lower, digits=4)))
    end
    lagrange_data = DataFrame(lagrange_rows)
    CSV.write("lagrange_multipliers.csv", lagrange_data)
    
    println("Optimization results have been saved to CSV files:")
    println(" - generator_results.csv")
    println(" - node_voltage_results.csv")
    println(" - power_flow_results.csv")
    println(" - lagrange_multipliers.csv")
else
    println("Solver did not find an optimal solution. Termination status: ", termination_status(model))
end