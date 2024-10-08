using JuMP
import Ipopt


# Create the model object
model = Model(Ipopt.Optimizer)

# Data
n = 9  # Number of generators
m = 11  # Number of nodes

# Generator data
max_capacity = [0.02, 0.15, 0.08, 0.07, 0.04, 0.17, 0.17, 0.26, 0.05]
cost_per_unit = [175, 100, 150, 150, 300, 350, 400, 300, 200]
locations = [2, 2, 2, 3, 4, 5, 7, 9, 9]

# Consumer data
consumer_nodes = [1, 4, 6, 8, 9, 10, 11]
active_power_demand = [0.10, 0.19, 0.11, 0.09, 0.21, 0.05, 0.04]

# Voltage and phase angle limits
voltage_limits = (0.98, 1.02)
angle_limits = (-π, π)

# Table of edges and coefficients
edges = [(1,2), (1,11), (2,3), (2,11), (3,4), (3,9), (4,5), (5,6), (5,8), (6,7), (7,8), (7,9), (8,9), (9,10), (10,11)]
g_kl = [4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.71, 1.26, 1.11, 1.32, 2.01, 4.41, 2.14, 5.06]
b_kl = [-20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3, -9.2, -13.9, -8.7, -11.3, -7.7, -13.5, -26.7]

# Model
model = Model(Ipopt.Optimizer)

# Variables
@variable(model, 0 <= A[i=1:n] <= max_capacity[i])  # Active power produced by each generator
@variable(model, -0.03 * max_capacity[i] <= R[i=1:n] <= 0.03 * max_capacity[i])  # Reactive power produced
@variable(model, voltage_limits[1] <= v[i=1:m] <= voltage_limits[2])  # Voltage magnitudes at each node
@variable(model, angle_limits[1] <= θ[i=1:m] <= angle_limits[2])  # Voltage angles at each node

# Objective: Minimize total production cost
@objective(model, Min, sum(cost_per_unit[i] * A[i] for i in 1:n))

@variable(model, zero == 0)  # Define a variable that is non-negative
#@constraint(model, zero == 0)  # Add a constraint to fix it at zero


# Power balance constraints
for k in 1:m
    # Active power balance
    active_generation=[A[i] for i in 1:n if locations[i] == k]
    push!(active_generation,zero)
    active_generation = sum(active_generation)
    println("active_generation: ",active_generation)
    
    active_consumption=[active_power_demand[j] for j in 1:length(consumer_nodes) if consumer_nodes[j] == k]
    push!(active_consumption,0)
    active_consumption = sum(active_consumption)
    
    active_flow=[v[k]^2 * g_kl[idx] - v[k] * v[l] * g_kl[idx] * cos(θ[k] - θ[l]) - v[k] * v[l] * b_kl[idx] * sin(θ[k] - θ[l])
    for (idx, (k, l)) in enumerate(edges) if k == k || l == k]
    push!(active_flow,zero)
    active_flow = sum(active_flow)

    @constraint(model, active_generation - active_flow - active_consumption == 0)
    
    # Reactive power balance
    reactive_generation=[R[i] for i in 1:n if locations[i] == k]
    push!(reactive_generation,zero)
    reactive_generation = sum(reactive_generation)
    
    reactive_flow=[-v[k]^2 * b_kl[idx] + v[k] * v[l] * b_kl[idx] * cos(θ[k] - θ[l]) - v[k] * v[l] * g_kl[idx] * sin(θ[k] - θ[l])
    for (idx, (k, l)) in enumerate(edges) if k == k || l == k]
    push!(reactive_flow,zero)
    reactive_flow = sum(reactive_flow)
    @constraint(model, reactive_generation - reactive_flow == 0)
end

# Solve model
optimize!(model)

# Printing some of the results for further analysis
println("") # Printing white line after solver output, before printing
println("Termination statue: ", JuMP.termination_status(model))
println("Optimal(?) objective function value: ", JuMP.objective_value(model))
println("Optimal(?) A points: ", JuMP.value.(A))
println("Optimal(?) R points: ", JuMP.value.(R))
println(JuMP.dual.(JuMP.UpperBoundRef.(A)))

for i in 1:n
    println("Active power for generator $i: ", value(A[i]))
    println("Reactive power for generator $i: ", value(R[i]))
end

#=println("Dual variables/Lagrange multipliers corresponding to some of the constraints: ")
println(JuMP.dual.(SOS_constr))=#

