import numpy as np
from scipy.spatial import Delaunay

def barycentric_coordinates(triangle, point):
    """
    Calculate the barycentric coordinates of a point relative to a triangle.
    Args:
        triangle: A 3x2 array with the coordinates of the triangle's vertices.
        point: A 1x2 array with the coordinates of the point.
    Returns:
        Barycentric coordinates (lambda1, lambda2, lambda3).
    """
    A = np.array([
        [triangle[0, 0], triangle[1, 0], triangle[2, 0]],
        [triangle[0, 1], triangle[1, 1], triangle[2, 1]],
        [1, 1, 1]
    ])
    b = np.array([point[0], point[1], 1])
    lambdas = np.linalg.solve(A, b)
    return lambdas

def interpolate_scalar(triangle, scalar_values, point):
    """
    Interpolate a scalar value at a given point using barycentric coordinates.
    Args:
        triangle: A 3x2 array with the coordinates of the triangle's vertices.
        scalar_values: A 1x3 array of scalar values at the triangle's vertices.
        point: A 1x2 array with the coordinates of the point.
    Returns:
        Interpolated scalar value.
    """
    lambdas = barycentric_coordinates(triangle, point)
    return np.dot(lambdas, scalar_values)

# Example source mesh (triangles and scalar values)
source_nodes = np.array([
    [0, 0],
    [1, 0],
    [0, 1],
    [1, 1]
])
source_triangles = np.array([
    [0, 1, 2],
    [1, 3, 2]
])
source_scalars = np.array([1.0, 2.0, 3.0, 4.0])  # Scalars at each source node

# Example target mesh
target_nodes = np.array([
    [0.25, 0.25],
    [0.75, 0.25],
    [0.25, 0.75]
])

# Find which source triangle each target node belongs to
tri = Delaunay(source_nodes)
target_scalars = np.zeros(len(target_nodes))

for i, point in enumerate(target_nodes):
    simplex = tri.find_simplex(point)
    if simplex == -1:
        print(f"Point {point} is outside the source mesh")
        continue
    vertices = source_triangles[simplex]
    triangle = source_nodes[vertices]
    scalar_values = source_scalars[vertices]
    target_scalars[i] = interpolate_scalar(triangle, scalar_values, point)

# Output interpolated scalar values at target nodes
print("Target scalars:", target_scalars)
