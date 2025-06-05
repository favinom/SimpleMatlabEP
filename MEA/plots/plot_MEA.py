import pyvista as pv
import numpy as np
import os
from PIL import Image

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
input_dir = "../outputs"
output_dir = "."
os.makedirs(output_dir, exist_ok=True)
filename = os.path.join(input_dir, "electrodes.vtk")
image_filename = os.path.join(output_dir, "electrodes_image.png")

# === Configura stile globale ===
pv.global_theme.background = 'white'
pv.global_theme.font.family = 'arial'
pv.global_theme.font.size = 30
pv.global_theme.font.color = 'black'
pv.global_theme.show_edges = True
pv.global_theme.cmap = 'plasma'

# === 1. Leggi il file VTK ===
grid = pv.read(filename)

# === 2. Estrai CELL_DATA ===
cell_scalars = grid.cell_data['Example_Cell_Scalar']

# === 3. Converti in una griglia uniforme (volumetrica) ===
# Per visualizzare in 2D: seleziona uno slice o proietta
dims = grid.dimensions
spacing = grid.spacing
origin = grid.origin

# Ricostruisci volume da cell data
volume = grid.cast_to_unstructured_grid()

# === 4. Crea plotter ===
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(volume, scalars='Example_Cell_Scalar', style='surface', show_scalar_bar=False, opacity=1.0)

# Barra scalare
#plotter.add_scalar_bar(
#    title="",
#    n_labels=5,
#    title_font_size=35,
#    label_font_size=30,
#    width=0.08,
#    height=0.6,
#    position_x=0.85,
#    position_y=0.2,
#    vertical=True,
#    use_opacity=False
#)

plotter.show_bounds(
    grid='back',  # Aggiunge linee di griglia sul piano posteriore
    location='outer',  # Posiziona gli assi all'esterno del plot
    all_edges=True,  # Mostra tutte le linee di bordo
    #xtitle='X Axis',
    ytitle='',
    bold=False,
    #ztitle='Z Axis',
    #n_xlabels=10,  # Numero di etichette sull'asse X
    #n_ylabels=10,  # Numero di etichette sull'asse Y
    fmt="%.2f",  # Formato delle etichette
    font_size=13,
    color='black'
)



# === 5. Configura camera 2D (XY) ===
plotter.view_xy()
plotter.camera.zoom(1.2)

# === 6. Salva immagine ===
plotter.screenshot(image_filename)
plotter.close()
print(f"Immagine salvata: {image_filename}")

crop_white_space(image_filename)
