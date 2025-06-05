import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import pyvista as pv
import os
from PIL import Image

# === Impostazioni iniziali ===
output_dir = "." #/test_matplotlib_AT_RT"
os.makedirs(output_dir, exist_ok=True)
vtk_file = "../outputs/AT_RT.vtk"  # <-- percorso file VTK
font_path = "./fonts/Helvetica.ttf"     # <-- percorso al font .ttf

# === Font personalizzato ===
custom_font = fm.FontProperties(fname=font_path)

# === Lettura file VTK e slicing ===
grid = pv.read(vtk_file)
dims = grid.dimensions[::-1]  # (z, y, x)
spacing = grid.spacing
origin = grid.origin

nx, ny, nz = dims[2], dims[1], dims[0]
dx, dy, dz = spacing
ox, oy, oz = origin
z_index = nz // 2

x = np.linspace(ox, ox + dx * (nx - 1), nx)
y = np.linspace(oy, oy + dy * (ny - 1), ny)
X, Y = np.meshgrid(x, y)

def plot_field(field_name, cmap, filename_out):
    data_3d = grid.point_data[field_name].reshape((nz, ny, nx))
    slice_2d = data_3d[z_index, :, :]

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(
        slice_2d,
        origin='lower',
        cmap=cmap,
        extent=(x[0], x[-1], y[0], y[-1]),
        interpolation='none'
    )

    contours = ax.contour(
        X, Y, slice_2d,
        levels=25,
        colors='black',
        linewidths=1.0
    )

    # Etichette degli assi
    ax.set_xlabel("X (mm)", fontproperties=custom_font, fontsize=14)
    ax.set_ylabel("Y (mm)", fontproperties=custom_font, fontsize=14)

    # Barre dei colori
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(field_name, rotation=270, labelpad=15,
                       fontproperties=custom_font, fontsize=16)

    #ax.set_title(f"Slice Z - Campo {field_name}", fontproperties=custom_font, fontsize=18)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(filename_out, dpi=300)
    plt.close()
    print(f"Salvata immagine: {filename_out}")

# === Plotta AT ===
plot_field("AT", cmap='coolwarm', filename_out=os.path.join(output_dir, "AT_new.png"))

# === Plotta RT ===
plot_field("RT", cmap='magma', filename_out=os.path.join(output_dir, "RT_new.png"))
