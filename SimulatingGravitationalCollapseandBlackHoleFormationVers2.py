# Author(s): Dr. Patrick Lemoine

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

# Constants and parameters
G = 1.0
c = 1.0
M_initial = 1.0

n_particles = 2000
dt = 0.05
friction_coeff = 0.01

# Initialize particle positions in two groups: near and far
np.random.seed(0)

n_near = int(n_particles * 0.7)
n_far = n_particles - n_near

radii_near = np.random.uniform(5, 10, n_near)
radii_far = np.random.uniform(15, 25, n_far)
radii = np.concatenate([radii_near, radii_far])

theta = np.arccos(np.random.uniform(-1, 1, n_particles))
phi = np.random.uniform(0, 2 * np.pi, n_particles)

x = radii * np.sin(theta) * np.cos(phi)
y = radii * np.sin(theta) * np.sin(phi)
z = radii * np.cos(theta)
pos = np.vstack((x, y, z))

mass = M_initial
v_circ = np.sqrt(G * mass / radii)
vx = -v_circ * np.sin(phi)
vy = v_circ * np.cos(phi)
vz = np.zeros(n_particles)
vel = np.vstack((vx, vy, vz))

# Starfield initialization (distant stars)
n_stars = 500
star_radii = np.random.uniform(100, 150, n_stars)  # Large radius for distant stars
star_theta = np.arccos(np.random.uniform(-1, 1, n_stars))
star_phi = np.random.uniform(0, 2 * np.pi, n_stars)

star_x = star_radii * np.sin(star_theta) * np.cos(star_phi)
star_y = star_radii * np.sin(star_theta) * np.sin(star_phi)
star_z = star_radii * np.cos(star_theta)
stars_pos = np.vstack((star_x, star_y, star_z))

def draw_sphere(radius, slices=30, stacks=30, color=(1, 0, 0, 0.3)):
    glPushMatrix()
    glColor4f(*color)
    quad = gluNewQuadric()
    gluQuadricNormals(quad, GLU_SMOOTH)
    gluQuadricDrawStyle(quad, GLU_FILL)
    gluSphere(quad, radius, slices, stacks)
    gluDeleteQuadric(quad)
    glPopMatrix()

def draw_particles(positions):
    glPointSize(3)
    glBegin(GL_POINTS)
    glColor3f(0, 0, 1)
    for p in positions.T:
        glVertex3f(*p)
    glEnd()

def draw_stars(positions):
    glPointSize(1)
    glBegin(GL_POINTS)
    glColor3f(1, 1, 1)  # White stars
    for p in positions.T:
        glVertex3f(*p)
    glEnd()

def init_opengl(width, height):
    glViewport(0, 0, width, height)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glClearColor(0, 0, 0, 1)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45, (width / height), 0.1, 1000.0)  # Large far clipping plane for stars
    glMatrixMode(GL_MODELVIEW)

def main():
    pygame.init()
    display = (800, 600)

    # Request multisampling buffers before creating window
    pygame.display.gl_set_attribute(pygame.GL_MULTISAMPLEBUFFERS, 1)
    pygame.display.gl_set_attribute(pygame.GL_MULTISAMPLESAMPLES, 4)  # 4x MSAA

    pygame.display.set_mode(display, DOUBLEBUF | OPENGL)

    init_opengl(*display)

    # Enable multisampling in OpenGL
    glEnable(GL_MULTISAMPLE)

    cam_distance = 40
    rot_x, rot_y = 0, 0
    mouse_down = False

    global pos, vel, mass

    clock = pygame.time.Clock()
    running = True
    frame = 0

    font = pygame.font.SysFont("Arial", 18)

    while running:
        clock.tick(60)
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

        # Physics update
        r = np.linalg.norm(pos, axis=0)
        r_safe = np.maximum(r, 1e-5)

        a_grav = -G * mass * pos / r_safe**3
        a_friction = -friction_coeff * vel
        a = a_grav + a_friction

        vel += a * dt
        pos += vel * dt

        Rs = 2 * G * mass / c**2
        accretion_zone = Rs * 1.1

        inside_horizon = r < Rs
        in_accretion_zone = (r >= Rs) & (r < accretion_zone)

        vel[:, in_accretion_zone] *= 0.90

        if np.any(inside_horizon):
            accreted_mass = np.sum(inside_horizon) * 0.0015
            mass += accreted_mass
            mask = ~inside_horizon
            pos = pos[:, mask]
            vel = vel[:, mask]

        # Rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        gluLookAt(
            cam_distance * np.sin(np.radians(rot_y)) * np.cos(np.radians(rot_x)),
            cam_distance * np.sin(np.radians(rot_x)),
            cam_distance * np.cos(np.radians(rot_y)) * np.cos(np.radians(rot_x)),
            0, 0, 0,
            0, 1, 0
        )

        # Draw starfield first (background)
        draw_stars(stars_pos)

        # Draw event horizon and particles
        draw_sphere(Rs)
        draw_particles(pos)

        elapsed_time = frame * dt
        info_text = f"Mass: {mass:.3f} | Particles: {pos.shape[1]} | Time: {elapsed_time:.2f}"
        text_surface = font.render(info_text, True, (255, 255, 255))
        text_data = pygame.image.tostring(text_surface, "RGBA", True)
        glWindowPos2d(10, display[1] - 30)
        glDrawPixels(text_surface.get_width(), text_surface.get_height(), GL_RGBA, GL_UNSIGNED_BYTE, text_data)

        pygame.display.flip()
        frame += 1

    pygame.quit()

if __name__ == "__main__":
    main()
