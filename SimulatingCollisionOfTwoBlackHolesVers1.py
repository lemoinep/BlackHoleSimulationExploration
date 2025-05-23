# Author(s): Dr. Patrick Lemoine

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch

# --- Physical constants and simulation parameters ---
G = 1.0       # Gravitational constant (Newtonian approximation)
c = 1.0       # Speed of light (set to 1 in simulation units)
dt = 0.05     # Time step for numerical integration
n_steps = 600 # Total number of simulation frames
friction_coeff = 0.01  # Artificial friction to simulate energy dissipation (e.g., gravitational wave emission)

# --- Initialize matter particles ---
# Particles represent diffuse matter around black holes, modeled as collisionless points
np.random.seed(0)  # Seed for reproducibility

# Distribute particles in three spherical shells at different radial distances
n_particles_1 = 500
n_particles_2 = 500
n_particles_3 = 1000

radii_far_1 = np.random.uniform(5, 10, n_particles_1)
radii_far_2 = np.random.uniform(15, 20, n_particles_2)
radii_far_3 = np.random.uniform(30, 40, n_particles_3)

n_particles = n_particles_1 + n_particles_2 + n_particles_3
radii = np.concatenate([radii_far_1, radii_far_2, radii_far_3])

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
mass1 = 1.0
mass2 = 1.0

pos_bh1 = np.array([-10.0, 0.0, 0.0])  # Initial position BH1 (left side)
pos_bh2 = np.array([10.0, 0.0, 0.0])   # Initial position BH2 (right side)

vel_bh1 = np.array([0.0, 0.15, 0.0])   # Initial velocity BH1 (upward)
vel_bh2 = np.array([0.0, -0.15, 0.0])  # Initial velocity BH2 (downward)

# --- Fusion control variables ---
fused = False                  # Flag indicating if black holes have merged
fusion_in_progress = False     # Flag for ongoing smooth fusion animation
fusion_steps = 30              # Number of frames over which fusion is interpolated
fusion_step_count = 0          # Counter for fusion progress
pos_fusion_start_1 = None      # Starting position BH1 at fusion start
pos_fusion_start_2 = None      # Starting position BH2 at fusion start
mass_fusion_start_1 = None     # Starting mass BH1 at fusion start
mass_fusion_start_2 = None     # Starting mass BH2 at fusion start
vel_fusion_start_1 = None      # Starting velocity BH1 at fusion start
vel_fusion_start_2 = None      # Starting velocity BH2 at fusion start

# --- Visualization setup ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

scat = ax.scatter(pos[0], pos[1], pos[2], s=5, color='blue', alpha=0.6)  # Matter particles

# Angular grids for plotting spherical event horizons
u = np.linspace(0, 2 * np.pi, 30)
v = np.linspace(0, np.pi, 30)

def plot_horizon(ax, center, mass, color='red', alpha=0.3):
    """
    Plot the event horizon of a black hole as a sphere with Schwarzschild radius.
    Rs = 2GM/c^2 defines the radius of the event horizon in Schwarzschild metric.
    """
    Rs = 2 * G * mass / c**2
    X = center[0] + Rs * np.outer(np.cos(u), np.sin(v))
    Y = center[1] + Rs * np.outer(np.sin(u), np.sin(v))
    Z = center[2] + Rs * np.outer(np.ones(np.size(u)), np.cos(v))
    return ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)

horizon1 = plot_horizon(ax, pos_bh1, mass1)
horizon2 = plot_horizon(ax, pos_bh2, mass2)

lim = 40
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-lim, lim)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
title = ax.set_title("Collision of Two Black Holes")

time_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

legend_elements = [
    Patch(facecolor='red', edgecolor='r', label='Black Hole Event Horizon'),
    Patch(facecolor='blue', edgecolor='b', label='Matter Particles')
]
ax.legend(handles=legend_elements, loc='upper right')

def safe_remove(artist):
    """Safely remove matplotlib artist if it exists and not already removed."""
    if artist is not None:
        try:
            artist.remove()
        except ValueError:
            pass

def update(frame):
    """
    Update function called every frame by matplotlib animation.
    Simulates gravitational interactions, accretion, and black hole merger.
    """
    global pos, vel, mass1, mass2, pos_bh1, pos_bh2, vel_bh1, vel_bh2
    global horizon1, horizon2, fused, fusion_in_progress, fusion_step_count
    global pos_fusion_start_1, pos_fusion_start_2, mass_fusion_start_1, mass_fusion_start_2
    global vel_fusion_start_1, vel_fusion_start_2

    # --- Black hole dynamics ---
    if not fused and not fusion_in_progress:
        # Compute vector and distance between black holes
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

    elif fusion_in_progress:
        # Smooth interpolation of fusion over multiple frames
        fusion_step_count += 1
        t = fusion_step_count / fusion_steps  # Normalized fusion progress [0,1]

        # Interpolate position, mass, and velocity of merged black hole
        pos_bh1 = (1 - t) * pos_fusion_start_1 + t * pos_fusion_start_2
        mass1 = (1 - t) * mass_fusion_start_1 + t * (mass_fusion_start_1 + mass_fusion_start_2)
        vel_bh1 = (1 - t) * vel_fusion_start_1 + t * vel_fusion_start_2

        # Deactivate second black hole during fusion
        pos_bh2 = np.array([np.nan, np.nan, np.nan])
        vel_bh2 = np.zeros(3)

        if fusion_step_count >= fusion_steps:
            fusion_in_progress = False
            fused = True
            mass2 = 0.0

    else:
        dist_bh = np.nan  # Not used after fusion

    # --- Gravitational acceleration on matter particles ---
    if fused or fusion_in_progress:
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

    # Add friction to simulate energy loss (e.g., gravitational wave emission)
    a_friction = -friction_coeff * vel

    # Total acceleration
    a = a_grav + a_friction

    # Update particle velocities and positions (Euler integration)
    vel += a * dt
    pos += vel * dt

    # --- Accretion of matter by black holes ---
    if fused or fusion_in_progress:
        Rs1 = 2 * G * mass1 / c**2  # Schwarzschild radius of merged BH
        inside_horizon1 = np.linalg.norm(pos - pos_bh1[:, np.newaxis], axis=0) < Rs1
        accreted_mass_1 = np.sum(inside_horizon1) * 0.005  # Mass gain per accreted particle
        mass1 += accreted_mass_1
        # Remove accreted particles from simulation
        mask = ~inside_horizon1
        pos = pos[:, mask]
        vel = vel[:, mask]
    else:
        Rs1 = 2 * G * mass1 / c**2
        Rs2 = 2 * G * mass2 / c**2

        inside_horizon1 = dist1 < Rs1
        inside_horizon2 = dist2 < Rs2

        accreted_mass_1 = np.sum(inside_horizon1) * 0.005
        accreted_mass_2 = np.sum(inside_horizon2) * 0.005

        mass1 += accreted_mass_1
        mass2 += accreted_mass_2

        mask = ~(inside_horizon1 | inside_horizon2)
        pos = pos[:, mask]
        vel = vel[:, mask]

    # --- Detect and initiate black hole merger ---
    if not fused and not fusion_in_progress and dist_bh < Rs1 + Rs2:
        fusion_in_progress = True
        fusion_step_count = 0
        pos_fusion_start_1 = pos_bh1.copy()
        pos_fusion_start_2 = pos_bh2.copy()
        mass_fusion_start_1 = mass1
        mass_fusion_start_2 = mass2
        vel_fusion_start_1 = vel_bh1.copy()
        vel_fusion_start_2 = vel_bh2.copy()

    # --- Update visualization ---
    scat._offsets3d = (pos[0], pos[1], pos[2])

    safe_remove(horizon1)
    if fused or fusion_in_progress:
        horizon1 = plot_horizon(ax, pos_bh1, mass1, color='darkred', alpha=0.5)
        safe_remove(horizon2)
    else:
        horizon1 = plot_horizon(ax, pos_bh1, mass1)
        safe_remove(horizon2)
        horizon2 = plot_horizon(ax, pos_bh2, mass2)

    # Update title with current masses and particle count
    title_text = "Collision of Two Black Holes\n"
    if fused:
        title_text += f"Merger completed - Total mass: {mass1:.2f}\n"
    elif fusion_in_progress:
        title_text += f"Merger in progress...\n"
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

    return scat, horizon1, title, time_text

ani = FuncAnimation(fig, update, frames=n_steps, interval=30, blit=False)

plt.show()
