import numpy as np

class ExplicitSolver:
    def __init__(self, num_nodes, dt, method="leapfrog"):
        self.num_nodes = num_nodes
        self.dt = dt
        self.method = method
        
        # State variables
        self.x = np.zeros((num_nodes, 3))  # Position
        self.v = np.zeros((num_nodes, 3))  # Velocity
        self.a = np.zeros((num_nodes, 3))  # Acceleration
        self.m = np.ones(num_nodes)  # Assume unit mass for simplicity
        
        # Previous step variables (for Verlet)
        self.x_prev = np.copy(self.x)
        
    def apply_forces(self):
        """Compute forces and update acceleration."""
        F = np.zeros((self.num_nodes, 3))  # External forces (e.g., gravity, springs, etc.)
        F[:, 1] -= 9.81  # Simple gravity in -y direction
        self.a = F / self.m[:, None]  # Acceleration = Force / Mass
        
    def step(self):
        if self.method == "leapfrog":
            self.leapfrog_step()
        elif self.method == "verlet":
            self.verlet_step()
        else:
            raise ValueError("Unsupported integration method")
        
    def leapfrog_step(self):
        """Perform a Leapfrog integration step."""
        self.v += 0.5 * self.dt * self.a  # Half-step velocity update
        self.x += self.dt * self.v  # Position update
        self.apply_forces()  # Compute new accelerations
        self.v += 0.5 * self.dt * self.a  # Complete velocity update
        
    def verlet_step(self):
        """Perform a Verlet integration step."""
        x_new = 2 * self.x - self.x_prev + (self.a * self.dt ** 2)
        self.x_prev = np.copy(self.x)  # Store previous position
        self.x = x_new  # Update position
        self.apply_forces()  # Compute new accelerations
        self.v = (self.x - self.x_prev) / (2 * self.dt)  # Approximate velocity
        
    def run(self, steps):
        """Run the solver for a given number of steps."""
        for _ in range(steps):
            self.step()
            print("Step completed. First node position:", self.x[0])

# Example usage
num_nodes = 10
dt = 0.01
method = "leapfrog"  # Change to "verlet" to use Verlet integration

solver = ExplicitSolver(num_nodes, dt, method)
solver.run(100)
