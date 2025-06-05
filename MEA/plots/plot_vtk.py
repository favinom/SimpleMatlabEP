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
output_dir = "./frames"
vtk_dir = "../outputs"
os.makedirs(output_dir, exist_ok=True)



# === Loop su tutti i file .vtk ===
for i in range(0, 10):  # imposta il range massimo previsto
#i = 1

    # === Configura stile globale ===
    pv.global_theme.background = 'white'
    pv.global_theme.font.family = 'arial'
    pv.global_theme.font.size = 30
    pv.global_theme.font.color = 'black'
    pv.global_theme.show_edges = True
    pv.global_theme.cmap = 'plasma'
    
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"sol_{istr}.vtk")
    image_filename_v = os.path.join(output_dir, f"sol_{istr}_v.png")
    image_filename_u = os.path.join(output_dir, f"sol_{istr}_u.png")

    #if not os.path.isfile(filename):
    #    continue

    print(f"Processing: {filename}")

    try:
        # === 1. Leggi il file VTK ===
        grid = pv.read(filename)

        # === 2. Estrai dati scalari disponibili ===
        v = grid.point_data.keys()[0] if grid.point_data else list(grid.cell_data.keys())[0]
        volume = grid.cast_to_unstructured_grid()

        # === 3. Crea plotter ===
        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(volume, scalars=v, style='surface', show_scalar_bar=False, opacity=1.0)

        # === 4. Configura assi ===
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

        # === 5. Camera XY ===
        plotter.view_xy()
        plotter.camera.zoom(1.2)

        # === 6. Salva immagine ===
        plotter.screenshot(image_filename_v)
        plotter.close()

        print(f"Immagine salvata: {image_filename_v}")
        crop_white_space(image_filename_v)



    except Exception as e:
        print(f"Errore nel file {filename}: {e}")



    # === Configura stile globale ===
    pv.global_theme.background = 'white'
    pv.global_theme.font.family = 'arial'
    pv.global_theme.font.size = 30
    pv.global_theme.font.color = 'black'
    pv.global_theme.show_edges = True
    pv.global_theme.cmap = 'turbo'


    try:
        # === 1. Leggi il file VTK ===
        grid = pv.read(filename)

        # === 2. Estrai dati scalari disponibili ===
        u = grid.point_data.keys()[1] if grid.point_data else list(grid.cell_data.keys())[1]
        volume = grid.cast_to_unstructured_grid()

        # === 3. Crea plotter ===
        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(volume, scalars=u, style='surface', show_scalar_bar=False, opacity=1.0, clim=[-4, 4])

        # === 4. Configura assi ===
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

        # === 5. Camera XY ===
        plotter.view_xy()
        plotter.camera.zoom(1.2)

        # === 6. Salva immagine ===
        plotter.screenshot(image_filename_u)
        plotter.close()

        print(f"Immagine salvata: {image_filename_u}")
        crop_white_space(image_filename_u)

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")




