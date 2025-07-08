import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import pyvista as pv
import os

# === Impostazioni iniziali ===
output_dir = "."
os.makedirs(output_dir, exist_ok=True)
vtk_files = [
    ("../confronto/bido_AT_RT.vtk", "BIDO"),
    ("../confronto/mea_AT_RT.vtk", "MEA")
]
font_path = "./fonts/Helvetica.ttf"

# Limiti colormap (modifica qui se vuoi)
clim_AT1 = (5, 60)
clim_AT2 = (1855, 1875)
clim_RT1 = None #(5, 60)
clim_RT2 = None #(1855, 1875)

# Dizionario per etichette personalizzate
field_labels = {
    "AT_1": "AT1",
    "AT_2": "AT2",
    "RT_1": "RT1",
    "RT_2": "RT2"
}

# === Font personalizzato ===
custom_font = fm.FontProperties(fname=font_path)

def get_slice(grid, field_name):
    dims = grid.dimensions[::-1]  # (z, y, x)
    spacing = grid.spacing
    origin = grid.origin
    nx, ny, nz = dims[2], dims[1], dims[0]
    dx, dy, dz = spacing
    ox, oy, oz = origin
    z_index = nz // 2  # Slicing al centro lungo Z
    x = np.linspace(ox, ox + dx * (nx - 1), nx)
    y = np.linspace(oy, oy + dy * (ny - 1), ny)
    X, Y = np.meshgrid(x, y)
    data_3d = grid.point_data[field_name].reshape((nz, ny, nx))
    slice_2d = data_3d[z_index, :, :]
    return X, Y, x, y, slice_2d



fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fields = [("AT_1", clim_AT1), ("AT_2", clim_AT2)]

for col, (vtk_file, titolo_colonna) in enumerate(vtk_files):
    grid = pv.read(vtk_file)
    for row, (field_name, clim) in enumerate(fields):
        ax = axes[row, col]
        X, Y, x, y, slice_2d = get_slice(grid, field_name)
        im = ax.imshow(
            slice_2d,
            origin='lower',
            cmap='coolwarm',
            extent=(x[0], x[-1], y[0], y[-1]),
            interpolation='none',
            vmin=clim[0] if clim else None,
            vmax=clim[1] if clim else None
        )
        ax.contour(
            X, Y, slice_2d,
            levels=25,
            colors='black',
            linewidths=1.0
        )
        ax.set_xlabel("X (mm)", fontproperties=custom_font, fontsize=10)
        ax.set_ylabel("Y (mm)", fontproperties=custom_font, fontsize=10)
        ax.set_aspect('equal')
        # Titolo della riga (campo) nella prima colonna
        #if col == 0:
        #    ax.annotate(
        #    field_labels.get(field_name, field_name),
        #    xy=(-0.22, 0.5), xycoords='axes fraction',
        #    ha='center', va='center', rotation=90,
        #    fontsize=16, fontproperties=custom_font
        #    )
        # Titolo della colonna (file) nella prima riga
        if row == 0:
            ax.annotate(titolo_colonna, xy=(0.5, 1.08), xycoords='axes fraction',
                        ha='center', va='center', fontsize=17, fontproperties=custom_font)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=11)
        cbar.ax.set_ylabel("")  # Rimuove l'etichetta verticale
        # Sposta il titolo della colorbar in verticale (sul lato sinistro della colorbar)
        
        if row == 0:
            cbar.set_label(
            field_labels.get(field_name, field_name),
            fontproperties=custom_font,
            fontsize=13,
            labelpad=-25,
            rotation=270,
            loc='center'
            )
        elif row == 1:
            cbar.set_label(
            field_labels.get(field_name, field_name),
            fontproperties=custom_font,
            fontsize=13,
            labelpad=-35,
            rotation=270,
            loc='center'
            )

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "AT_confronto.png"), dpi=300)
plt.close()
print(f"Salvata immagine: {os.path.join(output_dir, 'AT_confronto.png')}")





fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fields = [("RT_1", clim_RT1), ("RT_2", clim_RT2)]

for col, (vtk_file, titolo_colonna) in enumerate(vtk_files):
    grid = pv.read(vtk_file)
    for row, (field_name, clim) in enumerate(fields):
        ax = axes[row, col]
        X, Y, x, y, slice_2d = get_slice(grid, field_name)
        im = ax.imshow(
            slice_2d,
            origin='lower',
            cmap='gnuplot',
            extent=(x[0], x[-1], y[0], y[-1]),
            interpolation='none',
            vmin=clim[0] if clim else None,
            vmax=clim[1] if clim else None
        )
        ax.contour(
            X, Y, slice_2d,
            levels=25,
            colors='black',
            linewidths=1.0
        )
        ax.set_xlabel("X (mm)", fontproperties=custom_font, fontsize=10)
        ax.set_ylabel("Y (mm)", fontproperties=custom_font, fontsize=10)
        ax.set_aspect('equal')
        # Titolo della riga (campo) nella prima colonna
        #if col == 0:
        #    ax.annotate(
        #    field_labels.get(field_name, field_name),
        #    xy=(-0.22, 0.5), xycoords='axes fraction',
        #    ha='center', va='center', rotation=90,
        #    fontsize=16, fontproperties=custom_font
        #    )
        # Titolo della colonna (file) nella prima riga
        if row == 0:
            ax.annotate(titolo_colonna, xy=(0.5, 1.08), xycoords='axes fraction',
                        ha='center', va='center', fontsize=17, fontproperties=custom_font)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=11)
        cbar.ax.set_ylabel("")  # Rimuove l'etichetta verticale
        # Sposta il titolo della colorbar in verticale (sul lato sinistro della colorbar)
        
        if row == 0:
            cbar.set_label(
            field_labels.get(field_name, field_name),
            fontproperties=custom_font,
            fontsize=13,
            labelpad=-30,
            rotation=270,
            loc='center'
            )
        elif row == 1:
            cbar.set_label(
            field_labels.get(field_name, field_name),
            fontproperties=custom_font,
            fontsize=13,
            labelpad=-40,
            rotation=270,
            loc='center'
            )

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "RT_confronto.png"), dpi=300)
plt.close()
print(f"Salvata immagine: {os.path.join(output_dir, 'RT_confronto.png')}")


