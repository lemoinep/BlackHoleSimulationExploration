# Author(s): Dr. Patrick Lemoine

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Patch

# --- Physical constants and simulation parameters ---
G = 1.0       # Gravitational constant (Newtonian approximation)
c = 1.0       # Speed of light (unit system)
dt = 0.05     # Time step for numerical integration
n_steps = 800 # Number of simulation frames
friction_coeff = 0.01  # Artificial friction to mimic energy dissipation

# --- Initialize matter particles ---
np.random.seed(0)
n_particles_1, n_particles_2, n_particles_3 = 500, 500, 1000
radii_far_1 = np.random.uniform(5, 10, n_particles_1)
radii_far_2 = np.random.uniform(15, 20, n_particles_2)
radii_far_3 = np.random.uniform(30, 40, n_particles_3)
n_particles = n_particles_1 + n_particles_2 + n_particles_3
radii = np.concatenate([radii_far_1, radii_far_2, radii_far_3])
theta = np.arccos(np.random.uniform(-1, 1, n_particles))
phi = np.random.uniform(0, 2 * np.pi, n_particles)
x = radii * np.sin(theta) * np.cos(phi)
y = radii * np.sin(theta) * np.sin(phi)
z = radii * np.cos(theta)
pos = np.vstack((x, y, z))
mass_total_initial = 3.0  # Sum of initial BH masses
v_circ = np.sqrt(G * mass_total_initial / radii)
vx = -v_circ * np.sin(phi)
vy = v_circ * np.cos(phi)
vz = np.zeros(n_particles)
vel = np.vstack((vx, vy, vz))

# --- Initialize three black holes with different masses ---
masses = np.array([1.0, 1.5, 0.8])  # Different masses
positions = np.array([
    [-10.0, 0.0, 0.0],
    [10.0, 0.0, 0.0],
    [0.0, 15.0, 0.0]
])
velocities = np.array([
    [0.0, 0.15, 0.0],
    [0.0, -0.15, 0.0],
    [-0.15, 0.0, 0.0]
])

# Fusion tracking: keep track of which BHs remain active
active = np.array([True, True, True])

# --- Visualization setup ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scat = ax.scatter(pos[0], pos[1], pos[2], s=5, color='blue', alpha=0.6)
u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, np.pi, 30)
U, V = np.meshgrid(u, v)

def plot_deformed_horizon(ax, center, other_centers, mass, U, V,
                          max_deformation_distance=20, max_elongation=1.5, color='red', alpha=0.3):
    """
    Plot a deformable black hole horizon elongated toward the nearest other black hole.
    """
    Rs = 2 * G * mass / c**2
    # Find nearest other BH for deformation axis
    vecs = other_centers - center
    dists = np.linalg.norm(vecs, axis=1)
    if len(dists) == 0:
        r_hat = np.array([1, 0, 0])
    else:
        nearest_idx = np.argmin(dists)
        dist = dists[nearest_idx]
        r_vec = vecs[nearest_idx]
        r_hat = r_vec / dist if dist > 0 else np.array([1, 0, 0])

    x = np.sin(V) * np.cos(U)
    y = np.sin(V) * np.sin(U)
    z = np.cos(V)
    pos_vec = np.array([x, y, z])
    dot_prod = np.sum(pos_vec * r_hat[..., None, None], axis=0)
    deformation_factor = 1.0 + (max_elongation - 1.0) * np.exp(-dist / max_deformation_distance) * dot_prod**2

    X = center[0] + Rs * deformation_factor * np.sin(V) * np.cos(U)
    Y = center[1] + Rs * deformation_factor * np.sin(V) * np.sin(U)
    Z = center[2] + Rs * deformation_factor * np.cos(V)
    return ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)

lim = 50
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-lim, lim)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
title = ax.set_title("Collision of Three Black Holes with Deformable Horizons")
time_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

legend_elements = [
    Patch(facecolor='red', edgecolor='r', label='Deformable Black Hole Horizon'),
    Patch(facecolor='blue', edgecolor='b', label='Matter Particles')
]
ax.legend(handles=legend_elements, loc='upper right')

horizons = []

def safe_remove(artist):
    if artist is not None:
        try:
            artist.remove()
        except Exception:
            pass

def update(frame):
    global pos, vel, masses, positions, velocities, active, horizons

    n_bh = len(masses)
    # Update gravitational forces between active black holes
    for i in range(n_bh):
        if not active[i]:
            continue
        force = np.zeros(3)
        for j in range(n_bh):
            if i == j or not active[j]:
                continue
            r_vec = positions[j] - positions[i]
            dist = np.linalg.norm(r_vec)
            if dist < 1e-5:
                continue
            force += G * masses[i] * masses[j] / dist**3 * r_vec
        velocities[i] += force / masses[i] * dt

    # Update BH positions
    for i in range(n_bh):
        if active[i]:
            positions[i] += velocities[i] * dt

    # Compute gravitational acceleration on particles from all active BHs
    a_grav = np.zeros_like(pos)
    for i in range(n_bh):
        if not active[i]:
            continue
        r = pos - positions[i][:, None]
        dist = np.linalg.norm(r, axis=0)
        dist_safe = np.maximum(dist, 1e-5)
        a_grav += -G * masses[i] * r / dist_safe**3

    # Add friction to simulate energy loss
    a_friction = -friction_coeff * vel
    a_total = a_grav + a_friction

    # Update particle velocities and positions
    vel += a_total * dt
    pos += vel * dt

    # Accretion: remove particles inside any BH horizon and increase BH mass
    mask = np.ones(pos.shape[1], dtype=bool)
    for i in range(n_bh):
        if not active[i]:
            continue
        Rs = 2 * G * masses[i] / c**2
        dist_particles = np.linalg.norm(pos - positions[i][:, None], axis=0)
        inside = dist_particles < Rs
        accreted_mass = np.sum(inside) * 0.005
        masses[i] += accreted_mass
        mask &= ~inside  # Remove accreted particles

    pos = pos[:, mask]
    vel = vel[:, mask]

    # Check for mergers between black holes
    for i in range(n_bh):
        if not active[i]:
            continue
        for j in range(i + 1, n_bh):
            if not active[j]:
                continue
            dist_bh = np.linalg.norm(positions[j] - positions[i])
            Rs_i = 2 * G * masses[i] / c**2
            Rs_j = 2 * G * masses[j] / c**2
            if dist_bh < Rs_i + Rs_j:
                # Merge BH j into BH i
                total_mass = masses[i] + masses[j]
                positions[i] = (masses[i] * positions[i] + masses[j] * positions[j]) / total_mass
                velocities[i] = (masses[i] * velocities[i] + masses[j] * velocities[j]) / total_mass
                masses[i] = total_mass
                active[j] = False
                # Move merged BH far away to avoid rendering
                positions[j] = np.array([1e10, 1e10, 1e10])
                velocities[j] = np.zeros(3)

    # Update visualization
    scat._offsets3d = (pos[0], pos[1], pos[2])

    # Remove old horizons
    for h in horizons:
        safe_remove(h)
    horizons = []

    # Plot horizons for active BHs
    active_indices = np.where(active)[0]
    for idx in active_indices:
        others = np.delete(positions[active_indices], np.where(active_indices == idx), axis=0)
        h = plot_deformed_horizon(ax, positions[idx], others, masses[idx], U, V)
        horizons.append(h)

    # Update title with info
    title_text = f"Collision of Three Black Holes\n"
    for idx in active_indices:
        title_text += f"BH{idx+1}: Mass={masses[idx]:.2f}  "
    title_text += f"\nParticles remaining: {pos.shape[1]}"
    title.set_text(title_text)

    elapsed_time = frame * dt
    time_text.set_text(f"Elapsed time: {elapsed_time:.2f} units")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    return [scat] + horizons + [title, time_text]

ani = FuncAnimation(fig, update, frames=n_steps, interval=30, blit=False)
plt.show()

