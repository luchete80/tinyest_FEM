import matplotlib.pyplot as plt

def plot_deformation(node_coords, displacements, scale=100):
    deformed_coords = node_coords + displacements.reshape(-1, 2) * scale
    plt.figure(figsize=(6,6))
    
    # Plot original shape
    plt.plot(*zip(*node_coords, node_coords[0]), 'bo-', label='Original')
    # Plot deformed shape
    plt.plot(*zip(*deformed_coords, deformed_coords[0]), 'ro-', label='Deformed')
    
    plt.legend()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Finite Element Deformation")
    plt.show()

plot_deformation(node_coords, X)

