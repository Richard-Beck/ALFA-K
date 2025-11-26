import numpy as np
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

# ==========================================
# PART 1: SHARED SIMULATION DATA
# ==========================================
print("Generating shared simulation data (200x200)...")
n_per_group = 5
sigma_cluster = 0.1

# Simulation Grid
x_sim = np.linspace(-2, 2, 200)
y_sim = np.linspace(-2, 2, 200)
X_sim, Y_sim = np.meshgrid(x_sim, y_sim)
Z_sim = np.exp(-(X_sim ** 2 + Y_sim ** 2)) + 0.3 * np.sin(3 * X_sim) * np.cos(3 * Y_sim)

interp_func = RegularGridInterpolator((x_sim, y_sim), Z_sim.T)
surface_z = lambda coords: interp_func(coords)

np.random.seed(42)
sample_coords = np.stack((np.random.uniform(-2, 2, 10), np.random.uniform(-2, 2, 10)), axis=1)
start_center = sample_coords[np.argmin(surface_z(sample_coords))]
x0 = np.stack((
    np.random.normal(loc=start_center[0], scale=sigma_cluster, size=n_per_group),
    np.random.normal(loc=start_center[1], scale=sigma_cluster, size=n_per_group)
), axis=1)


def clip_to_bounds(coords):
    return np.clip(coords, -2, 2)


x0 = clip_to_bounds(x0)


def s(group):
    z_vals = surface_z(group)
    centroid = group[np.argmax(z_vals)]
    new_coords = np.stack((
        np.random.normal(loc=centroid[0], scale=sigma_cluster, size=n_per_group),
        np.random.normal(loc=centroid[1], scale=sigma_cluster, size=n_per_group)
    ), axis=1)
    return clip_to_bounds(new_coords)


x_groups = [x0]
for _ in range(14):
    x_groups.append(s(x_groups[-1]))

visible_indices = list(range(0, 15, 2))
visible_groups = [x_groups[i] for i in visible_indices]
# Recalculate Z exactly on simulation grid
z_visible_sim = [surface_z(g) for g in visible_groups]
highlight_positions = np.vstack(visible_groups)

# Camera Setup
cell_positions = np.vstack(visible_groups)
cell_center = cell_positions.mean(axis=0)
z_center = np.mean([z.mean() for z in z_visible_sim])
cell_focus = np.array([cell_center[0], cell_center[1], z_center])

azimuth = 315
elevation = 30
radius = 4.0
phi = np.radians(90 - elevation)
theta = np.radians(azimuth)
eye_x = cell_focus[0] + radius * np.sin(phi) * np.cos(theta)
eye_y = cell_focus[1] + radius * np.sin(phi) * np.sin(theta)
eye_z = cell_focus[2] + radius * np.cos(phi)
eye_location = np.array([eye_x, eye_y, eye_z])

# ==========================================
# PART 2: RENDER 1 - GLOWING SPHERES (Unchanged)
# ==========================================
print("Rendering Image 1: Glowing Spheres...")

pl1 = pv.Plotter(off_screen=True, window_size=[2000, 2000])
pl1.set_background("black")
pl1.enable_anti_aliasing("ssaa")
pl1.enable_depth_peeling(number_of_peels=10)
pl1.enable_eye_dome_lighting()

grid_sim = pv.StructuredGrid(X_sim, Y_sim, Z_sim)

# BASELINE MESH: Dark Grey, Specular 0.2
pl1.add_mesh(
    grid_sim,
    color="#1a1a1a",
    show_edges=True,
    edge_color="#333333",
    line_width=0.5,
    lighting=True,
    specular=0.2
)

cmap = plt.get_cmap("viridis")
n_steps = len(visible_groups)

for i, (group, z_vals) in enumerate(zip(visible_groups, z_visible_sim)):
    color = cmap(i / (n_steps - 1))[:3]
    points = np.column_stack((group[:, 0], group[:, 1], z_vals + 0.02))
    cloud = pv.PolyData(points)

    core = cloud.glyph(geom=pv.Sphere(radius=0.03), scale=False, orient=False)
    pl1.add_mesh(core, color=color, lighting=False, opacity=1.0)
    halo1 = cloud.glyph(geom=pv.Sphere(radius=0.045), scale=False, orient=False)
    pl1.add_mesh(halo1, color=color, lighting=False, opacity=0.1)
    halo2 = cloud.glyph(geom=pv.Sphere(radius=0.06), scale=False, orient=False)
    pl1.add_mesh(halo2, color=color, lighting=False, opacity=0.1)

pl1.camera_position = [eye_location, cell_focus, (0, 0, 1)]
pl1.screenshot("lscape_spheres.png")
pl1.close()
print("Saved lscape_spheres.png")

# ==========================================
# PART 3: RENDER 2 - VIBRANT PAINT (Fixes applied here)
# ==========================================
print("Rendering Image 2: Vibrant Surface Paint...")

pl2 = pv.Plotter(off_screen=True, window_size=[2000, 2000])
pl2.set_background("black")
pl2.enable_anti_aliasing("ssaa")
pl2.enable_eye_dome_lighting()

# Use the EXACT same mesh object
grid_ren = pv.StructuredGrid(X_sim, Y_sim, Z_sim)
grid_points = grid_ren.points
gx = grid_points[:, 0]
gy = grid_points[:, 1]

highlight_radius = 0.03
min_dists = np.full(gx.shape, np.inf)

for cx, cy in highlight_positions:
    d = np.sqrt((gx - cx) ** 2 + (gy - cy) ** 2)
    min_dists = np.minimum(min_dists, d)

scalar_field = np.zeros_like(min_dists, dtype=int)
scalar_field[min_dists <= (highlight_radius * 4)] = 1
scalar_field[min_dists <= (highlight_radius * 2)] = 2
scalar_field[min_dists <= highlight_radius] = 3

# 0: Wireframe Match (#1a1a1a)
# 1: Electric Blue
# 2: Neon Cyan
# 3: Pure White
vibrant_colors = ["#1a1a1a", "#1F51FF", "#00FFFF", "#ffffff"]
grid_ren.point_data["colors"] = scalar_field

pl2.add_mesh(
    grid_ren,
    scalars="colors",
    cmap=vibrant_colors,
    n_colors=4,
    show_scalar_bar=False,

    # --- FIX 1: KILL THE BLOTCHES ---
    # Set specular to 0.0 to prevent any shiny highlights.
    # Set roughness to 1.0 for a fully matte surface.
    lighting=True,
    specular=0.0,
    roughness=1.0,

    # --- FIX 2: RESTORE SMOOTH PAINT ---
    # Removed interpolate_before_map=False.
    # Default behavior (True) will smooth colors between vertices.

    show_edges=True,
    edge_color="#333333",
    line_width=0.5
)

pl2.camera_position = [eye_location, cell_focus, (0, 0, 1)]
pl2.screenshot("lscape_painted.png")
pl2.close()
print("Saved lscape_painted.png")
print("Done.")