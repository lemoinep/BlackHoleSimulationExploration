# Author(s): Dr. Patrick Lemoine

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch  # For legend handles

# Purpose of the simulation:
# This simulation aims to visually and simply illustrate the process of gravitational collapse of a matter distribution,
# leading to the formation of a black hole.
# We want to observe how, under its own gravity, a region of space can concentrate matter increasingly,
# forming an extremely dense object from which nothing, not even light, can escape.
# This simplified simulation does not claim to faithfully reproduce all the complexities of relativistic physics,
# but it allows us to visualize some key concepts.

# Physical constants (arbitrary units for simulation)
G = 1.0  # Gravitational constant (Newtonian) - determines the strength of gravitational attraction.
c = 1.0  # Speed of light (set to 1 in our units) - important for Schwarzschild radius calculation.
M_initial = 1.0  # Initial mass of the black hole seed - the black hole starts with a nonzero mass to initiate attraction.

# Simulation parameters
#n_particles = 1500  # Number of particles simulating matter - more particles means more realistic but more computationally expensive.
dt = 0.05  # Time step between simulation frames - affects simulation speed and numerical stability.
n_steps = 600  # Total number of simulation steps.
friction_coeff = 0.01  # Friction coefficient (energy dissipation) - simulates energy loss allowing matter to spiral inward.

# Initialize random particle positions inside a spherical shell
np.random.seed(0)  # Seed random number generator for reproducibility.
#radii = np.random.uniform(5, 20, n_particles)  # Random radii between 5 and 10 units.

n_particles_1 = 500
n_particles_2 = 500
n_particles_3 = 1000

radii_far_1 =  np.random.uniform(5, 10, n_particles_1)
radii_far_2 =  np.random.uniform(15, 20, n_particles_2)
radii_far_3 =  np.random.uniform(30, 40, n_particles_3)

n_particles = n_particles_1+n_particles_2+n_particles_3 # Number of particles simulating matter - more particles means more realistic but more computationally expensive.
radii = np.concatenate([radii_far_1,radii_far_2,radii_far_3])

theta = np.arccos(np.random.uniform(-1, 1, n_particles))  # Polar angles (theta), uniform distribution on sphere.
phi = np.random.uniform(0, 2 * np.pi, n_particles)  # Azimuthal angles (phi), uniform between 0 and 2pi.

# Convert spherical coordinates (r, theta, phi) to Cartesian coordinates (x, y, z)
x = radii * np.sin(theta) * np.cos(phi)
y = radii * np.sin(theta) * np.sin(phi)
z = radii * np.cos(theta)
pos = np.vstack((x, y, z))  # Positions matrix shape (3, n_particles)

# Initialize approximate orbital velocities perpendicular to radius vector
mass = M_initial  # Initial mass for velocity calculation, updated during simulation
v_circ = np.sqrt(G * mass / radii)  # Circular orbital velocity approximation: v = sqrt(GM/r)
vx = -v_circ * np.sin(phi)  # Velocity component x (tangential)
vy = v_circ * np.cos(phi)   # Velocity component y (tangential)
vz = np.zeros(n_particles)  # Velocity component z (initially zero)
vel = np.vstack((vx, vy, vz))  # Velocities matrix shape (3, n_particles)

# Create matplotlib figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial scatter plot of particles (blue)
scat = ax.scatter(pos[0], pos[1], pos[2], s=5, color='blue', alpha=0.6)

# Parameters for plotting the event horizon (sphere)
u = np.linspace(0, 2 * np.pi, 30)  # Azimuthal angle grid
v = np.linspace(0, np.pi, 30)      # Polar angle grid

def plot_horizon(ax, mass):
    # Schwarzschild radius Rs = 2GM/c^2
    Rs = 2 * G * mass / c**2
    X = Rs * np.outer(np.cos(u), np.sin(v))
    Y = Rs * np.outer(np.sin(u), np.sin(v))
    Z = Rs * np.outer(np.ones(np.size(u)), np.cos(v))
    return ax.plot_surface(X, Y, Z, color='red', alpha=0.3)  # Semi-transparent red sphere

horizon = plot_horizon(ax, mass)  # Initial event horizon plot

# Set axis limits and labels
lim1 = 30
ax.set_xlim(-lim1, lim1)
ax.set_ylim(-lim1, lim1)
ax.set_zlim(-lim1, lim1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
title = ax.set_title(f"Gravitational Collapse - Mass: {mass:.2f}")

# Dynamic text to display elapsed simulation time
time_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

# Legend showing color meaning
legend_elements = [
    Patch(facecolor='red', edgecolor='r', label='Event Horizon (Black Hole)'),
    Patch(facecolor='blue', edgecolor='b', label='Dark Matter Particles')
]
ax.legend(handles=legend_elements, loc='upper right')

def update(frame):
    global pos, vel, mass, horizon

    # Compute distances of particles from center
    r = np.linalg.norm(pos, axis=0)
    r_safe = np.maximum(r, 1e-5)  # Avoid division by zero

    # Compute gravitational acceleration (Newtonian gravity)
    a_grav = -G * mass * pos / r_safe**3  # a = -GM/r^3 * r_vector

    # Compute friction acceleration (energy dissipation)
    a_friction = -friction_coeff * vel

    # Total acceleration
    a = a_grav + a_friction

    # Update velocities and positions (Euler integration)
    vel += a * dt
    pos += vel * dt

    # Detect particles inside event horizon (radius Rs)
    Rs = 2 * G * mass / c**2
    inside_horizon = r < Rs

    if np.any(inside_horizon):
        # Calculate mass accreted this step (proportional to number of particles inside horizon)
        accreted_mass = np.sum(inside_horizon) * 0.005
        mass += accreted_mass  # Increase black hole mass

        # Remove accreted particles from simulation
        mask = ~inside_horizon
        pos = pos[:, mask]
        vel = vel[:, mask]

    # Update scatter plot positions
    scat._offsets3d = (pos[0], pos[1], pos[2])

    # Remove old horizon and plot new horizon with updated mass
    horizon.remove()
    horizon = plot_horizon(ax, mass)

    # Update title and elapsed time text
    title.set_text(f"Gravitational Collapse - Mass: {mass:.2f}\nParticles remaining: {pos.shape[1]}")
    elapsed_time = frame * dt
    time_text.set_text(f"Elapsed time: {elapsed_time:.2f} units")

    # Keep axis limits fixed for consistent view
    lim1 = 30
    ax.set_xlim(-lim1, lim1)
    ax.set_ylim(-lim1, lim1)
    ax.set_zlim(-lim1, lim1)

    # Return updated artists for animation
    return scat, horizon, title, time_text

# Create animation
ani = FuncAnimation(fig, update, frames=n_steps, interval=30, blit=False)

# Show plot with animation
plt.show()
