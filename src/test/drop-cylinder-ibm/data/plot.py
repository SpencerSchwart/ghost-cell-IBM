
# Fully corrected sessile drop analytical solution

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Given parameters
Rc = 1.0  # Cylinder radius
theta_deg = 60  # Contact angle in degrees
theta = np.radians(theta_deg)  # Convert to radians

# Function to solve for R_d (drop radius) and lc (drop center location)
def equations(vars):
    R_d, lc = vars
    # 1. Ensuring the drop's bottom aligns correctly with the cylinder
    eq1 = Rc - (lc - R_d * np.cos(theta))
    # 2. Correcting the intersection condition based on the correct geometry
    eq2 = (lc - Rc)**2 + (R_d * np.sin(theta))**2 - R_d**2  # Ensures curvature match
    return [eq1, eq2]

# Improved initial guess for R_d and lc
initial_guess = [Rc / np.cos(theta), Rc + 0.5]

# Solve for R_d and lc
solution = fsolve(equations, initial_guess)
R_d, lc = solution

# Generate cylinder coordinates
theta_cyl = np.linspace(0, 2 * np.pi, 300)
x_cyl = Rc * np.cos(theta_cyl)
y_cyl = Rc * np.sin(theta_cyl)

# Generate droplet coordinates (ensuring full arc)
theta_drop = np.linspace(-theta, theta, 300)
x_drop = lc - R_d * np.cos(theta_drop)
y_drop = R_d * np.sin(theta_drop)

# Compute the correct contact points where the drop touches the cylinder
contact_x1 = Rc * np.cos(theta)
contact_y1 = Rc * np.sin(theta)
contact_x2 = Rc * np.cos(-theta)
contact_y2 = Rc * np.sin(-theta)

# Plot
fig, ax = plt.subplots(figsize=(6,6))
ax.plot(x_cyl, y_cyl, 'k-', linewidth=2, label="Cylinder")
ax.plot(x_drop, y_drop, 'b-', linewidth=2, label="Sessile Drop")
ax.scatter([contact_x1, contact_x2], [contact_y1, contact_y2], color='r', zorder=3, label="Contact Points")

# Formatting
ax.set_aspect('equal')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.legend()
ax.set_title(f"Sessile Drop on a Cylinder (Contact Angle = {theta_deg}°)")

# Show plot
plt.show()

