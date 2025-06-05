import pyvista as pv
import numpy as np
import os
from PIL import Image, ImageDraw, ImageFont
import pytesseract

# === Funzione per rifilare spazi bianchi ===
def crop_white_space(image_path):
    with Image.open(image_path) as im:
        im = im.convert("RGBA")
        bbox = im.getbbox()
        if bbox:
            cropped = im.crop(bbox)
            cropped.save(image_path)
            print(f"Spazio bianco rimosso: {image_path}")

# === Parametri personalizzabili ===
output_dir = "."  # <--- Modifica questo percorso a piacere
input_dir = "../outputs"  # <--- Modifica questo percorso a piacere
os.makedirs(output_dir, exist_ok=True)
filename = os.path.join(input_dir, "AT_RT.vtk")
image_filename_AT = os.path.join(output_dir, "AT.png")
image_filename_RT = os.path.join(output_dir, "RT.png")


# === Configura stile globale ===
#pv.global_theme.background = 'black' | 'white' | '#RRGGBB' | (r, g, b[, a])
#pv.global_theme.font.family = 'arial' | 'times' | 'courier'  NON DISPONIBILI| 'helvetica' | 'cambria' | 'calibri' | ...
#pv.global_theme.font.size = int (es: 10, 20, 30...)
#pv.global_theme.font.color = 'red' | '#00ff00' | (1, 0, 0) | ...
#pv.global_theme.show_edges = True | False
#pv.global_theme.cmap = 'viridis' | 'plasma' | 'magma' | 'jet' | 'coolwarm' | ...

pv.global_theme.background = 'white'
pv.global_theme.font.family = 'arial'
pv.global_theme.font.size = 30
pv.global_theme.font.color = 'black'
pv.global_theme.show_edges = True
pv.global_theme.cmap = 'coolwarm'


# === 1. Leggi il file VTK ===
grid = pv.read(filename)

# === 2. Estrai le dimensioni e slice centrale lungo Z ===
dims = grid.dimensions[::-1]  # (z, y, x)
spacing = grid.spacing
origin = grid.origin

nx, ny, nz = dims[2], dims[1], dims[0]
dx, dy, dz = spacing
ox, oy, oz = origin

z_index = nz // 2

# === 3. Estrai campi volumetrici ===
at_3d = grid.point_data['AT'].reshape((nz, ny, nx))  # Assicurati che esista

at_slice = at_3d[z_index, :, :]

# === 4. Costruisci griglia strutturata 2D ===
x = np.linspace(ox, ox + dx * (nx - 1), nx)
y = np.linspace(oy, oy + dy * (ny - 1), ny)
X, Y = np.meshgrid(x, y)
Z = np.full_like(X, oz + z_index * dz)
sg = pv.StructuredGrid(X, Y, Z)
sg['AT'] = at_slice.flatten(order='C')

# === 5. Crea il plot ===
plotter = pv.Plotter(off_screen=True)

# Mappa colori AT
plotter.add_mesh(sg, scalars='AT', style='points', edge_color='black', point_size=7, opacity=0.8, show_scalar_bar=False)

# Contorni
plotter.show_bounds(
    grid='back',
    location='outer',
    all_edges=True,
    xtitle='',
    ytitle='',
    bold=False,
    fmt="%.2f",
    font_size=13,
    color='black'
)

# Scalar bar laterale
plotter.add_scalar_bar(
    title="",
    n_labels=5,
    title_font_size=35,
    label_font_size=25,
    #font_family='times',
    width=0.08,
    height=0.6,
    position_x=0.85,
    position_y=0.2,
    vertical=True,
    use_opacity=False
)


# Contorni AT
contours = sg.contour(isosurfaces=25, scalars='AT')
plotter.add_mesh(contours, color='black', line_width=3.0, show_scalar_bar=False)

# Configura vista XY
plotter.view_xy()
plotter.trasparent_background = True
#plotter.set_background("white")
plotter.camera.zoom(1.2)

# Salva immagine
plotter.screenshot(image_filename_AT)
plotter.close()
print(f"Immagine salvata: {image_filename_AT}")

crop_white_space(image_filename_AT)




# === Configura stile globale ===
#pv.global_theme.background = 'black' | 'white' | '#RRGGBB' | (r, g, b[, a])
#pv.global_theme.font.family = 'arial' | 'times' | 'courier'  NON DISPONIBILI| 'helvetica' | 'cambria' | 'calibri' | ...
#pv.global_theme.font.size = int (es: 10, 20, 30...)
#pv.global_theme.font.color = 'red' | '#00ff00' | (1, 0, 0) | ...
#pv.global_theme.show_edges = True | False
#pv.global_theme.cmap = 'viridis' | 'plasma' | 'magma' | 'jet' | 'coolwarm' | ...

pv.global_theme.background = 'white'
pv.global_theme.font.family = 'arial'
pv.global_theme.font.size = 30
pv.global_theme.font.color = 'black'
pv.global_theme.show_edges = True
pv.global_theme.cmap = 'magma'


# === 1. Leggi il file VTK ===
grid = pv.read(filename)

# === 2. Estrai le dimensioni e slice centrale lungo Z ===
dims = grid.dimensions[::-1]  # (z, y, x)
spacing = grid.spacing
origin = grid.origin

nx, ny, nz = dims[2], dims[1], dims[0]
dx, dy, dz = spacing
ox, oy, oz = origin

z_index = nz // 2

# === 3. Estrai campi volumetrici ===
rt_3d = grid.point_data['RT'].reshape((nz, ny, nx))  # Assicurati che esista
rt_slice = rt_3d[z_index, :, :]

# === 4. Costruisci griglia strutturata 2D ===
x = np.linspace(ox, ox + dx * (nx - 1), nx)
y = np.linspace(oy, oy + dy * (ny - 1), ny)
X, Y = np.meshgrid(x, y)
Z = np.full_like(X, oz + z_index * dz)
sg = pv.StructuredGrid(X, Y, Z)
sg['RT'] = rt_slice.flatten(order='C')

# === 5. Crea il plot ===
plotter = pv.Plotter(off_screen=True)

# Mappa colori RT
plotter.add_mesh(sg, scalars='RT', style='points', edge_color='black', point_size=7, opacity=0.8, show_scalar_bar=False)

# Contorni
plotter.show_bounds(
    grid='back',
    location='outer',
    all_edges=True,
    xtitle='',
    ytitle='',
    bold=False,
    fmt="%.2f",
    font_size=13,
    color='black'
)

# Scalar bar laterale
plotter.add_scalar_bar(
    title="",
    n_labels=5,
    title_font_size=35,
    label_font_size=25,
    #font_family='times',
    width=0.08,
    height=0.6,
    position_x=0.85,
    position_y=0.2,
    vertical=True,
    use_opacity=False
)


# Contorni RT
contours = sg.contour(isosurfaces=25, scalars='RT')
plotter.add_mesh(contours, color='black', line_width=3.0, show_scalar_bar=False)

# Configura vista XY
plotter.view_xy()
plotter.trasparent_background = True
#plotter.set_background("white")
plotter.camera.zoom(1.2)

# Salva immagine
plotter.screenshot(image_filename_RT)
plotter.close()
print(f"Immagine salvata: {image_filename_RT}")

crop_white_space(image_filename_RT)


