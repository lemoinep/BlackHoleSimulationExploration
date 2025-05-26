# Author(s): Dr. Patrick Lemoine

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Patch

# --- Physical constants and simulation parameters ---
G = 1.0       # Gravitational constant (Newtonian approximation)
c = 1.0       # Speed of light (unit system)
dt = 0.5     # Time step for numerical integration
n_steps = 600 # Total number of simulation frames
friction_coeff = 0.01  # Artificial friction to mimic energy dissipation (e.g., gravitational wave emission)

# --- Initialize matter particles ---
# Particles represent diffuse matter (e.g., gas, dust, or stars) surrounding black holes.
np.random.seed(0)  # Seed for reproducibility

# Distribute particles in three spherical shells at different radial distances
n_particles_1, n_particles_2, n_particles_3, n_particles_4 = 500, 500, 1000,1000
radii_far_1 = np.random.uniform(5, 10, n_particles_1)
radii_far_2 = np.random.uniform(15, 20, n_particles_2)
radii_far_3 = np.random.uniform(30, 40, n_particles_3)
radii_far_4 = np.random.uniform(70, 100, n_particles_4)
n_particles = n_particles_1 + n_particles_2 + n_particles_3 + n_particles_4
radii = np.concatenate([radii_far_1, radii_far_2, radii_far_3, radii_far_4])

# Uniform angular distribution on sphere (theta: polar angle, phi: azimuthal angle)
theta = np.arccos(np.random.uniform(-1, 1, n_particles))
phi = np.random.uniform(0, 2 * np.pi, n_particles)

# Convert spherical coordinates to Cartesian coordinates for particle positions
x = radii * np.sin(theta) * np.cos(phi)
y = radii * np.sin(theta) * np.sin(phi)
z = radii * np.cos(theta)
pos = np.vstack((x, y, z))  # Shape: (3, n_particles)

# Initialize approximate circular orbital velocities for particles around combined black hole mass
mass_total_initial = 2.0  # Sum of initial black hole masses
v_circ = np.sqrt(G * mass_total_initial / radii)  # Circular velocity approximation v = sqrt(GM/r)
vx = -v_circ * np.sin(phi)  # Velocity components tangential to radial vector
vy = v_circ * np.cos(phi)
vz = np.zeros(n_particles)  # Initial z-velocity zero for simplicity
vel = np.vstack((vx, vy, vz))  # Shape: (3, n_particles)

# --- Initialize two black holes ---
# Each black hole is modeled as a point mass with position, velocity, and mass
mass1, mass2 = 1.0, 1.5  # Initial black hole masses

pos_bh1 = np.array([-20.0, 0.0, 10.0])  # Initial position BH1 (left side)
pos_bh2 = np.array([30.0, 0.0, 0.0])   # Initial position BH2 (right side)

vel_bh1 = np.array([0.0, 0.15, 0.0])   # Initial velocity BH1 (upward)
vel_bh2 = np.array([0.0, -0.15, 0.0])  # Initial velocity BH2 (downward)

fused = False  # Flag indicating whether the black holes have merged

# --- Visualization setup ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot for matter particles (blue dots)
scat = ax.scatter(pos[0], pos[1], pos[2], s=5, color='blue', alpha=0.6)

# Angular grids for plotting spherical event horizons
u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, np.pi, 30)
U, V = np.meshgrid(u, v)

def plot_deformed_horizon(ax, center, other_center, mass, U, V,
                          max_deformation_distance=20, max_elongation=1.5, color='red', alpha=0.3):
    """
    Plot a deformable black hole event horizon as a sphere with radius equal to the Schwarzschild radius,
    deformed (elongated) along the axis pointing toward the other black hole to simulate tidal distortion.

    Parameters:
    - ax: matplotlib 3D axis
    - center: 3D coordinates of the black hole center
    - other_center: 3D coordinates of the other black hole (for deformation direction)
    - mass: mass of the black hole
    - U, V: meshgrid arrays for spherical angles
    - max_deformation_distance: distance scale over which deformation is significant
    - max_elongation: maximum elongation factor along the axis
    - color: color of the horizon sphere
    - alpha: transparency level
    """
    Rs = 2 * G * mass / c**2  # Schwarzschild radius: Rs = 2GM/c^2

    # Vector from this black hole to the other
    r_vec = other_center - center
    dist = np.linalg.norm(r_vec)
    if dist > 0:
        r_hat = r_vec / dist
    else:
        r_hat = np.array([1, 0, 0])  # Default direction if coincident

    # Unit vectors on the sphere surface
    x = np.sin(V) * np.cos(U)
    y = np.sin(V) * np.sin(U)
    z = np.cos(V)
    pos_vec = np.array([x, y, z])

    # Project each surface point onto the axis vector to get elongation factor
    dot_prod = np.sum(pos_vec * r_hat[..., None, None], axis=0)

    # Calculate deformation factor: stronger elongation when closer
    deformation_factor = 1.0 + (max_elongation - 1.0) * np.exp(-dist / max_deformation_distance) * dot_prod**2

    # Compute deformed horizon coordinates
    X = center[0] + Rs * deformation_factor * np.sin(V) * np.cos(U)
    Y = center[1] + Rs * deformation_factor * np.sin(V) * np.sin(U)
    Z = center[2] + Rs * deformation_factor * np.cos(V)

    return ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)

# Set plot limits and labels
lim = 40
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-lim, lim)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
title = ax.set_title("Collision of Two Black Holes with Deformable Horizons")
time_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

# Legend for particles and horizons
legend_elements = [
    Patch(facecolor='red', edgecolor='r', label='Deformable Black Hole Horizon'),
    Patch(facecolor='blue', edgecolor='b', label='Matter Particles')
]
ax.legend(handles=legend_elements, loc='upper right')

# Initial horizons
horizon1 = plot_deformed_horizon(ax, pos_bh1, pos_bh2, mass1, U, V)
horizon2 = plot_deformed_horizon(ax, pos_bh2, pos_bh1, mass2, U, V)

def safe_remove(artist):
    """Safely remove a matplotlib artist if it exists and is not already removed."""
    if artist is not None:
        try:
            artist.remove()
        except ValueError:
            pass

def update(frame):
    global pos, vel, mass1, mass2, pos_bh1, pos_bh2, vel_bh1, vel_bh2
    global horizon1, horizon2, fused

    # --- Black hole dynamics ---
    if not fused:
        # Vector and distance between black holes
        r_bh = pos_bh2 - pos_bh1
        dist_bh = np.linalg.norm(r_bh)
        dist_bh_safe = max(dist_bh, 1e-5)  # Avoid division by zero

        # Newtonian gravitational force between black holes
        force_bh = G * mass1 * mass2 / dist_bh_safe**3 * r_bh

        # Update velocities (Newton's 3rd law)
        vel_bh1 += force_bh / mass1 * dt
        vel_bh2 -= force_bh / mass2 * dt

        # Update positions
        pos_bh1 += vel_bh1 * dt
        pos_bh2 += vel_bh2 * dt
    else:
        dist_bh = np.nan  # Not used after fusion

    # --- Gravitational acceleration on matter particles ---
    if fused:
        # Single black hole gravitational field after merger
        r1 = pos - pos_bh1[:, np.newaxis]
        dist1 = np.linalg.norm(r1, axis=0)
        dist1_safe = np.maximum(dist1, 1e-5)
        a_grav = -G * mass1 * r1 / dist1_safe**3
    else:
        # Sum of gravitational fields from both black holes before merger
        r1 = pos - pos_bh1[:, np.newaxis]
        r2 = pos - pos_bh2[:, np.newaxis]
        dist1 = np.linalg.norm(r1, axis=0)
        dist2 = np.linalg.norm(r2, axis=0)
        dist1_safe = np.maximum(dist1, 1e-5)
        dist2_safe = np.maximum(dist2, 1e-5)
        a_grav1 = -G * mass1 * r1 / dist1_safe**3
        a_grav2 = -G * mass2 * r2 / dist2_safe**3
        a_grav = a_grav1 + a_grav2

    # Artificial friction to mimic dissipative processes allowing particles to spiral inward
    a_friction = -friction_coeff * vel

    # Total acceleration
    a = a_grav + a_friction

    # Update particle velocities and positions using Euler integration
    vel += a * dt
    pos += vel * dt

    # --- Accretion of matter by black holes ---
    Rs1 = 2 * G * mass1 / c**2  # Schwarzschild radius BH1
    Rs2 = 2 * G * mass2 / c**2  # Schwarzschild radius BH2

    # Identify particles inside each horizon (accreted)
    inside_horizon1 = np.linalg.norm(pos - pos_bh1[:, np.newaxis], axis=0) < Rs1
    inside_horizon2 = np.linalg.norm(pos - pos_bh2[:, np.newaxis], axis=0) < Rs2

    # Calculate mass gain proportional to number of accreted particles
    accreted_mass_1 = np.sum(inside_horizon1) * 0.005
    accreted_mass_2 = np.sum(inside_horizon2) * 0.005

    mass1 += accreted_mass_1
    mass2 += accreted_mass_2

    # Remove accreted particles from simulation
    mask = ~(inside_horizon1 | inside_horizon2)
    pos = pos[:, mask]
    vel = vel[:, mask]

    # --- Detect and process black hole merger ---
    if not fused and dist_bh < Rs1 + Rs2:
        fused = True
        total_mass = mass1 + mass2

        # Compute center of mass position and velocity for merged black hole
        pos_bh1 = (mass1 * pos_bh1 + mass2 * pos_bh2) / total_mass
        vel_bh1 = (mass1 * vel_bh1 + mass2 * vel_bh2) / total_mass

        mass1 = total_mass
        mass2 = 0.0

        # Move second black hole far away to avoid rendering issues
        pos_bh2 = np.array([1e10, 1e10, 1e10])
        vel_bh2 = np.zeros(3)

    # --- Update visualization ---
    scat._offsets3d = (pos[0], pos[1], pos[2])

    safe_remove(horizon1)
    safe_remove(horizon2)

    if fused:
        horizon1 = plot_deformed_horizon(ax, pos_bh1, pos_bh2, mass1, U, V,
                                        max_deformation_distance=20, max_elongation=1.0,
                                        color='darkred', alpha=0.5)
        horizon2 = None
    else:
        horizon1 = plot_deformed_horizon(ax, pos_bh1, pos_bh2, mass1, U, V,
                                        max_deformation_distance=20, max_elongation=1.5,
                                        color='red', alpha=0.4)
        horizon2 = plot_deformed_horizon(ax, pos_bh2, pos_bh1, mass2, U, V,
                                        max_deformation_distance=20, max_elongation=1.5,
                                        color='red', alpha=0.4)

    # Update title with current masses and remaining particle count
    title_text = "Collision of Two Black Holes with Deformable Horizons\n"
    if fused:
        title_text += f"Merger completed - Total mass: {mass1:.2f}\n"
    else:
        title_text += f"Mass BH1: {mass1:.2f} | Mass BH2: {mass2:.2f}\n"
    title_text += f"Remaining particles: {pos.shape[1]}"
    title.set_text(title_text)

    # Update elapsed time display
    elapsed_time = frame * dt
    time_text.set_text(f"Elapsed time: {elapsed_time:.2f} units")

    # Fix axis limits for consistent viewing
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    return scat, horizon1, horizon2, title, time_text

# Create animation object
ani = FuncAnimation(fig, update, frames=n_steps, interval=30, blit=False)

# Show the animation
plt.show()
