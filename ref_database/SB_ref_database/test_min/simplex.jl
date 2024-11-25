using Combinatorics

# Function to generate a simplex in n dimensions with given step size
function generate_simplex(n, step)
    # Generate all possible points with the given step size
    points = collect(0:step:1)
    
    # Generate all combinations of points in n dimensions
    combinations = Iterators.product(ntuple(_ -> points, n)...)
    
    # Filter combinations to keep only those where the sum of coordinates is 1
    simplex = [collect(comb) for comb in combinations if sum(comb) ≈ 1.0]
    
    return simplex
end

# Function to adjust the simplex points
function adjust_simplex(simplex, eps)
    adjusted_simplex = []
    for point in simplex
        adjusted_point = copy(point)
        for i in 1:length(point)
            if point[i] == 0.0
                adjusted_point[i] = eps
            elseif point[i] == 1.0
                adjusted_point[i] = 1.0 - eps
            end
        end
        # Adjust the remaining coordinates to ensure the sum is still 1
        sum_adjusted = sum(adjusted_point)
        if sum_adjusted ≈ 1.0
            push!(adjusted_simplex, adjusted_point)
        else
            # Find the index of the largest coordinate
            max_index = argmax(adjusted_point)
            adjusted_point[max_index] += 1.0 - sum_adjusted
            push!(adjusted_simplex, adjusted_point)
        end
    end
    return adjusted_simplex
end

# Define the number of dimensions, step size, and epsilon
n       = 4
step    = 0.1
eps     = 1e-4

# Generate a simplex in n dimensions
simplex = generate_simplex(n, step)

# Adjust the simplex points
adjusted_simplex = adjust_simplex(simplex, eps)

# Print the vertices of the adjusted simplex
println("Vertices of the adjusted simplex:")
for vertex in adjusted_simplex
    println(vertex)
end
sum.(adjusted_simplex,dims=1)