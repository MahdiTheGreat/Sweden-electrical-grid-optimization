using JuMP
using Ipopt
using DataFrames
using CSV

# Function to safely sum over a collection, returning 0.0 if empty
function safe_sum(collection)
    return isempty(collection) ? 0.0 : sum(collection)
end

# Function to check if a value is within bounds considering tolerance
function within_bounds(val, lower, upper, tol)
    return (val >= lower - tol) && (val <= upper + tol)
end

# Function to verify constraints
function verify_constraints(n, m, GeneratorsAtNode, DemandAtNode, max_capacity, G, B, v_values, θ_values, A_values, R_values, voltage_limits, angle_limits, tolerance=1e-6)
    println("\n--- Constraint verification ---")

    # Generator capacity constraints
    println("\nGenerator capacity constraints:")
    local all_gen_constraints = true
    for i in 1:n
        a = A_values[i]
        if !within_bounds(a, 0.0, max_capacity[i], tolerance)
            println(" - Generator G$i: A[$i] = $a exceeds bounds [0, $(max_capacity[i])]")
            all_gen_constraints = false
        end
    end
    if all_gen_constraints
        println("All generator capacity constraints are satisfied.")
    end

    # Reactive power constraints
    println("\nReactive power constraints:")
    local all_reactive_constraints = true
    for i in 1:n
        r = R_values[i]
        if !within_bounds(r, -0.03 * max_capacity[i], 0.03 * max_capacity[i], tolerance)
            println(" - Generator G$i: R[$i] = $r exceeds bounds [$( -0.03 * max_capacity[i]), $(0.03 * max_capacity[i])]")
            all_reactive_constraints = false
        end
    end
    if all_reactive_constraints
        println("All reactive power constraints are satisfied.")
    end

    # Voltage limits
    println("\nVoltage limits:")
    local all_voltage_constraints = true
    for k in 1:m
        v = v_values[k]
        if !within_bounds(v, voltage_limits[1], voltage_limits[2], tolerance)
            println(" - Node $k: v[$k] = $v pu exceeds bounds [$(voltage_limits[1]), $(voltage_limits[2])]")
            all_voltage_constraints = false
        end
    end
    if all_voltage_constraints
        println("All voltage constraints are satisfied.")
    end

    # Angle limits
    println("\nAngle limits:")
    local all_angle_constraints = true
    for k in 1:m
        θ_val = θ_values[k]
        if !within_bounds(θ_val, angle_limits[1], angle_limits[2], tolerance)
            println(" - Node $k: θ[$k] = $θ_val rad exceeds bounds [$(angle_limits[1]), $(angle_limits[2])]")
            all_angle_constraints = false
        end
    end
    if all_angle_constraints
        println("All angle constraints are satisfied.")
    end

    # Power balance constraints
    println("\nPower balance constraints:")
    local all_power_balance = true
    for k in 1:m
        # Calculate P_k and Q_k based on optimized values
        P_calc = sum(v_values[k] * v_values[l] * (G[k,l] * cos(θ_values[k] - θ_values[l]) + B[k,l] * sin(θ_values[k] - θ_values[l])) for l in 1:m)
        Q_calc = sum(v_values[k] * v_values[l] * (G[k,l] * sin(θ_values[k] - θ_values[l]) - B[k,l] * cos(θ_values[k] - θ_values[l])) for l in 1:m)

        # Expected P_k and Q_k from constraints using safe_sum
        P_expected = safe_sum([A_values[i] for i in GeneratorsAtNode[k]]) - get(DemandAtNode, k, 0.0)
        Q_expected = safe_sum([R_values[i] for i in GeneratorsAtNode[k]])

        # Check if the calculated values match the expected values within tolerance
        if abs(P_calc - P_expected) > tolerance
            println(" - Node $k: Active power balance violated. |Calculated: $P_calc, Expected: $P_expected|")
            all_power_balance = false
        end
        if abs(Q_calc - Q_expected) > tolerance
            println(" - Node $k: Reactive power balance violated. |Calculated: $Q_calc, Expected: $Q_expected|")
            all_power_balance = false
        end
    end
    if all_power_balance
        println("All power balance constraints are satisfied.")
    end

    # Overall verification status
    local overall_status = all_gen_constraints && all_reactive_constraints && all_voltage_constraints && all_angle_constraints && all_power_balance
    if overall_status
        println("\nAll constraints are successfully satisfied within the specified tolerance.")
    else
        println("\nSome constraints are violated. Please review the above messages for details.")
    end
end

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
DemandAtNode = Dict(k => safe_sum([active_power_demand[j] for j in 1:length(consumer_nodes) if consumer_nodes[j] == k]) for k in 1:m)

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
@constraint(model, [k=1:m], P_k[k] == safe_sum([A[i] for i in GeneratorsAtNode[k]]) - get(DemandAtNode, k, 0.0))
@constraint(model, [k=1:m], Q_k[k] == safe_sum([R[i] for i in GeneratorsAtNode[k]]))

# Solve the model
optimize!(model)

# Check solver status
status = termination_status(model)
if status == JuMP.OPTIMAL || status == JuMP.LOCALLY_SOLVED
    # Generator results
    generator_data = DataFrame(
        Generator = ["G$i" for i in 1:n],
        A_pu = [round(value(A[i]), digits=6) for i in 1:n],
        R_pu = [round(value(R[i]), digits=6) for i in 1:n],
        M_pu = max_capacity,
        h_SEK_pu = cost_per_unit
    )
    CSV.write("generator_results.csv", generator_data)
    
    # Node voltage results
    node_voltage_data = DataFrame(
        Node = 1:m,
        v_pu = [round(value(v[k]), digits=6) for k in 1:m],
        θ_rad = [round(value(θ[k]), digits=6) for k in 1:m]
    )
    CSV.write("node_voltage_results.csv", node_voltage_data)
    
    # Power Flow Results
    power_flow_rows = []
    for idx in 1:length(edges)
        (k, l) = edges[idx]
        g = g_kl[idx]
        b = b_kl[idx]
        θ_k = value(θ[k])
        θ_l = value(θ[l])
        v_k = value(v[k])
        v_l = value(v[l])
        # Compute the power flows
        p_kl = round(v_k^2 * g - v_k * v_l * (g * cos(θ_k - θ_l) + b * sin(θ_k - θ_l)), digits=6)
        p_lk = round(v_l^2 * g - v_l * v_k * (g * cos(θ_l - θ_k) + b * sin(θ_l - θ_k)), digits=6)
        q_kl = round(-v_k^2 * b - v_k * v_l * (g * sin(θ_k - θ_l) - b * cos(θ_k - θ_l)), digits=6)
        q_lk = round(-v_l^2 * b - v_l * v_k * (g * sin(θ_l - θ_k) - b * cos(θ_l - θ_k)), digits=6)
        push!(power_flow_rows, (Edge = "($k,$l)", p_kl_pu = p_kl, p_lk_pu = p_lk, q_kl_pu = q_kl, q_lk_pu = q_lk))
    end
    power_flow_data = DataFrame(power_flow_rows)
    CSV.write("power_flow_results.csv", power_flow_data)
    
    println("Optimization results have been saved to CSV files.")

    # Retrieve optimized variable values
    A_values = [value(A[i]) for i in 1:n]
    R_values = [value(R[i]) for i in 1:n]
    v_values = [value(v[k]) for k in 1:m]
    θ_values = [value(θ[k]) for k in 1:m]

    # Call the verification function
    verify_constraints(n, m, GeneratorsAtNode, DemandAtNode, max_capacity, G, B, v_values, θ_values, A_values, R_values, voltage_limits, angle_limits)

else
    println("Solver did not find an optimal solution. Termination status: ", status)
end
