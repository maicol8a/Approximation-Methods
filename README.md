# Universal Algorithms for Approximating Projection Depths
This repository provides implementations of three universal algorithms for approximating projection depths in multivariate statistics: Random Search, Refined Random Search, and Spherical Nelder-Mead. These algorithms are designed to work with any depth function that satisfies the projection property, making them versatile for various statistical analyses.

# Algorithms
## Random Search Universal
The Random Search Universal algorithm approximates projection depths by generating random directions on the unit sphere and computing the depth in each direction. The minimum depth found across all directions is returned as the approximate depth.

## Refined Random Search Universal
The Refined Random Search Universal algorithm improves upon the basic random search by iteratively refining the search within a shrinking neighborhood around the current best direction. This approach concentrates the search in regions with lower depth, leading to more efficient approximation.

## Spherical Nelder-Mead Universal
The Spherical Nelder-Mead Universal algorithm adapts the classic Nelder-Mead optimization method to the geometry of the unit sphere. It minimizes the depth along great circles, making it efficient for approximating projection depths in high-dimensional spaces.

# Usage
To use these algorithms, you need to adapt your depth functions to accept the following parameters:

data: The data matrix (n x d).
point: The point to evaluate (vector of length d).
direction: The direction vector on the unit sphere.
...: Any additional parameters required by the depth function.
The depth function should return the computed depth for the given data, point, and direction.

## Example
Here is an example of how to use the Random Search Universal algorithm with a custom depth function:

### Define a custom depth function
custom_depth <- function(data, point, direction, p = 0.1) {
  # Compute the depth using the provided direction
  # Your depth computation logic here
}

### Generate sample data
set.seed(123)
data <- matrix(rnorm(2000), ncol = 10)  # 200 x 10
point <- rnorm(10)  # Length 10

### Run the Random Search Universal algorithm
result <- random.search.universal(
  data = data,
  point = point,
  depth.func = custom_depth,
  n_it = 1000,
  p = 0.1
)

print(result)

# Acknowledgments
The code was created by I. Cascos (https://github.com/icascos) and M. Ochoa. The theoretical description of the procedures is available at: https://doi.org/10.1016/j.csda.2020.107166.
