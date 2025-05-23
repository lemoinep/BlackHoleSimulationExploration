# Author(s): Dr. Patrick Lemoine

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

# === Fundamental Constants and Simulation Parameters ===
# Units are normalized (G = c = 1) for simplicity in relativistic calculations.
# In these geometric units, masses and distances are implicitly in meters via the Schwarzschild radius Rs = 2GM/c².
G = 1.0  # Gravitational constant (dimensionless) [See: Schutz, "A First Course in General Relativity"]
c = 1.0  # Speed of light (dimensionless)

M_initial = 1.0  # Initial black hole mass

n_particles = 2000  # Number of baryonic matter particles
dt = 0.05           # Integration time step

# Tuned parameters for physical plausibility
mass_particle = 3e-4     # Mass of each baryonic particle
friction_coeff = 0.035   # Friction coefficient (models baryonic dissipation)
emission_rate = 0.5      # Photon emission rate (Hawking radiation)

# === Initial Conditions: Two Baryonic Matter Populations ===
# The bimodal distribution (near/far) models:
# 1) A fragmented accretion disk (radii 5-10)
# 2) A cold matter reservoir (radii 15-25)
# Reference: α-viscosity accretion disk model [Shakura & Sunyaev 1973]
np.random.seed(0)
n_near = int(n_particles * 0.6)
n_far = n_particles - n_near

radii_near = np.random.uniform(5, 10, n_near)   # Rapid coalescence region
radii_far = np.random.uniform(15, 25, n_far)    # Outer cold matter reservoir
radii = np.concatenate([radii_near, radii_far])

# Spherical coordinates for isotropic distribution
theta = np.arccos(np.random.uniform(-1, 1, n_particles))
phi = np.random.uniform(0, 2 * np.pi, n_particles)

# Conversion to Cartesian coordinates
x = radii * np.sin(theta) * np.cos(phi)
y = radii * np.sin(theta) * np.sin(phi)
z = radii * np.cos(theta)
pos = np.vstack((x, y, z))          # Shape: (3, n_particles)

mass = M_initial

# === Initial Velocities: Keplerian Circular Orbits ===
# v_circ = sqrt(GM/r) (Newtonian approximation)
# Note: In general relativity, the innermost stable circular orbit is at r = 6GM/c²,
# with maximum orbital velocity v = c/√3 [Bardeen et al. 1972]
v_circ = np.sqrt(G * mass / radii)
vx = -v_circ * np.sin(phi)
vy = v_circ * np.cos(phi)
vz = np.zeros(n_particles)
vel = np.vstack((vx, vy, vz))       # Shape: (3, n_particles)

# === Background Stars (for visual context) ===
n_stars = 500
star_radii = np.random.uniform(100, 150, n_stars)
star_theta = np.arccos(np.random.uniform(-1, 1, n_stars))
star_phi = np.random.uniform(0, 2 * np.pi, n_stars)

star_x = star_radii * np.sin(star_theta) * np.cos(star_phi)
star_y = star_radii * np.sin(star_theta) * np.sin(star_phi)
star_z = star_radii * np.cos(star_theta)
stars_pos = np.vstack((star_x, star_y, star_z))  # Shape: (3, n_stars)

# === Emitted Particles (Hawking Radiation) ===
emitted_pos = np.empty((3, 0))   # Positions of emitted photons
emitted_vel = np.empty((3, 0))   # Velocities of emitted photons

# === Rendering Functions (OpenGL) ===

def draw_sphere(radius, slices=30, stacks=30, color=(1, 0, 0, 0.3)):
    """
    Render a semi-transparent sphere (used for black hole event horizon/accretion region).
    Note: In reality, strong gravitational lensing would distort the appearance of the horizon.
    """
    glPushMatrix()
    glColor4f(*color)
    quad = gluNewQuadric()
    gluQuadricNormals(quad, GLU_SMOOTH)
    gluQuadricDrawStyle(quad, GLU_FILL)
    gluSphere(quad, radius, slices, stacks)
    gluDeleteQuadric(quad)
    glPopMatrix()

def draw_particles(positions, velocities, size=3):
    """
    Render baryonic particles as colored points.
    Color encodes velocity magnitude:
    - Blue: cold matter (v << v_escape)
    - Red: relativistic matter (v ~ c/2)
    """
    glPointSize(size)
    glBegin(GL_POINTS)
    for p, v in zip(positions.T, velocities.T):
        speed = np.linalg.norm(v)
        t = min(speed / 2.0, 1.0)
        r = t
        b = 1 - t
        glColor3f(r, 0, b)
        glVertex3f(*p)
    glEnd()

def draw_emitted_particles(positions):
    """
    Render emitted photons (Hawking radiation) as yellow points.
    In reality, Hawking radiation is thermal (blackbody spectrum at T_H = ħc³/(8πkGM)).
    Here, emission is simplified as isotropic test particles.
    """
    glPointSize(4)
    glBegin(GL_POINTS)
    glColor3f(1, 1, 0)
    for p in positions.T:
        glVertex3f(*p)
    glEnd()

def draw_stars(positions):
    """
    Render distant background stars as white points.
    """
    glPointSize(1)
    glBegin(GL_POINTS)
    glColor3f(1, 1, 1)
    for p in positions.T:
        glVertex3f(*p)
    glEnd()

def init_opengl(width, height):
    """
    Initialize OpenGL context, perspective, and lighting.
    """
    glViewport(0, 0, width, height)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glClearColor(0, 0, 0, 1)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45, (width / height), 0.1, 1000.0)
    glMatrixMode(GL_MODELVIEW)

    # Lighting configuration
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    light_pos = [10.0, 10.0, 10.0, 1.0]
    light_ambient = [0.5, 0.5, 0.5, 1.0]
    light_diffuse = [0.7, 0.7, 0.7, 1.0]
    light_specular = [1.0, 1.0, 1.0, 1.0]
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos)
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)
    glEnable(GL_COLOR_MATERIAL)
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [1,1,1,1])
    glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 50)

# === Main Simulation Loop ===

def main():
    """
    Main loop: handles event processing, physics updates, and rendering.
    """
    global pos, vel, mass, emitted_pos, emitted_vel

    pygame.init()
    display = (800, 600)

    # Enable anti-aliasing for smoother rendering
    pygame.display.gl_set_attribute(pygame.GL_MULTISAMPLEBUFFERS, 1)
    pygame.display.gl_set_attribute(pygame.GL_MULTISAMPLESAMPLES, 4)

    pygame.display.set_mode(display, DOUBLEBUF | OPENGL)
    init_opengl(*display)

    glEnable(GL_MULTISAMPLE)

    cam_distance = 60
    rot_x, rot_y = 40, 0
    mouse_down = False

    clock = pygame.time.Clock()
    running = True
    frame = 0

    font = pygame.font.SysFont("Arial", 18)

    while running:
        clock.tick(60)  # Limit to 60 FPS
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:
                    mouse_down = True
                elif event.button == 4:
                    cam_distance = max(5, cam_distance - 2)
                elif event.button == 5:
                    cam_distance = min(100, cam_distance + 2)
            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:
                    mouse_down = False
            elif event.type == pygame.MOUSEMOTION and mouse_down:
                dx, dy = event.rel
                rot_x += dy
                rot_y += dx

        # === Physics Update ===

        r = np.linalg.norm(pos, axis=0)             # Radial distance from origin
        r_safe = np.maximum(r, 1e-5)                # Avoid division by zero

        a_grav = -G * mass * pos / r_safe**3        # Newtonian gravitational acceleration
        # Note: In full GR, geodesic equations would be required.
        a_friction = -friction_coeff * vel          # Phenomenological friction (models turbulent viscosity, shocks, etc.)
        a = a_grav + a_friction

        vel += a * dt                               # Velocity Verlet integration
        pos += vel * dt

        Rs = 2 * G * mass / c**2                    # Schwarzschild radius (event horizon)
        # Beyond r ≈ 1.5Rs lies the photon sphere (unstable photon orbits)
        accretion_zone = Rs * 1.1                   # Slightly outside the event horizon

        inside_horizon = r < Rs                     # Particles inside the event horizon
        in_accretion_zone = (r >= Rs) & (r < accretion_zone)

        vel[:, in_accretion_zone] *= 0.90           # Damping in the accretion region

        if np.any(inside_horizon):
            # Inside the horizon, all geodesics point toward the singularity.
            # Accreted mass ΔM includes internal energy via E=mc² (energy-momentum conservation).
            accreted_mass = np.sum(inside_horizon) * mass_particle
            mass += accreted_mass                   # Increase BH mass by accreted particles
            mask = ~inside_horizon
            pos = pos[:, mask]
            vel = vel[:, mask]

        # === Hawking Evaporation (Simplified) ===
        # Evaporation rate: dM/dt ∝ ħc⁴/(G²M²) [Hawking 1974]
        # Here: simplified numerical law for visualization
        evaporation_rate = 1e-5 / (mass**2 + 1e-8)
        mass_loss = evaporation_rate * dt
        mass = max(mass - mass_loss, 0.1)           # Prevent negative mass

        # === Photon Emission (Hawking Radiation) ===
        n_emit = int(emission_rate)
        if n_emit > 0:
            emit_radii = Rs * 1.2 + np.random.uniform(-0.1, 0.1, n_emit)
            emit_theta = np.arccos(np.random.uniform(-1, 1, n_emit))
            emit_phi = np.random.uniform(0, 2 * np.pi, n_emit)

            emit_x = emit_radii * np.sin(emit_theta) * np.cos(emit_phi)
            emit_y = emit_radii * np.sin(emit_theta) * np.sin(emit_phi)
            emit_z = emit_radii * np.cos(emit_theta)
            new_emit_pos = np.vstack((emit_x, emit_y, emit_z))

            speed_emit = 2.0 / mass  # Emission speed related to effective temperature
            emit_vel_dir = new_emit_pos / np.linalg.norm(new_emit_pos, axis=0)
            new_emit_vel = emit_vel_dir * speed_emit

            emitted_pos = np.hstack((emitted_pos, new_emit_pos))
            emitted_vel = np.hstack((emitted_vel, new_emit_vel))

        if emitted_pos.shape[1] > 0:
            emitted_pos += emitted_vel * dt
            dist_emit = np.linalg.norm(emitted_pos, axis=0)
            keep_emit = dist_emit < 200             # Remove photons that have escaped far
            emitted_pos = emitted_pos[:, keep_emit]
            emitted_vel = emitted_vel[:, keep_emit]

        # === Rendering ===

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        # Camera setup
        gluLookAt(
            cam_distance * np.sin(np.radians(rot_y)) * np.cos(np.radians(rot_x)),
            cam_distance * np.sin(np.radians(rot_x)),
            cam_distance * np.cos(np.radians(rot_y)) * np.cos(np.radians(rot_x)),
            0, 0, 0,
            0, 1, 0
        )

        draw_stars(stars_pos)

        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE)
        glColor4f(1, 0, 0, 0.1)
        draw_sphere(Rs * 1.3)                      # Visualize accretion region
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        draw_sphere(Rs, color=(1, 0, 0, 0.6))      # Visualize event horizon

        draw_particles(pos, vel)                    # Render baryonic particles

        draw_emitted_particles(emitted_pos)         # Render Hawking radiation

        # === HUD: Display Simulation Info ===
        elapsed_time = frame * dt
        info_text = f"BH Mass: {mass:.4f} | Particles: {pos.shape[1]} | Emitted: {emitted_pos.shape[1]} | Time: {elapsed_time:.2f}"
        text_surface = font.render(info_text, True, (255, 255, 255))
        text_data = pygame.image.tostring(text_surface, "RGBA", True)
        glWindowPos2d(10, display[1] - 30)
        glDrawPixels(text_surface.get_width(), text_surface.get_height(), GL_RGBA, GL_UNSIGNED_BYTE, text_data)

        pygame.display.flip()
        frame += 1

    pygame.quit()

if __name__ == "__main__":
    main()
