import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import pyvista as pv
import os

# === Parametri ===
output_dir = "."
os.makedirs(output_dir, exist_ok=True)
vtk_file = "../confronto/electrodes.vtk"  # <-- file VTK dal secondo script
font_path = "./fonts/Helvetica.ttf"     # <-- percorso al font .ttf

# === Carica font ===
custom_font = fm.FontProperties(fname=font_path)

# === Carica griglia da file ===
grid = pv.read(vtk_file)

# === Estrai dimensioni/griglia ===
dims = grid.dimensions[::-1]  # (z, y, x)
spacing = grid.spacing
origin = grid.origin

nx, ny, nz = dims[2], dims[1], dims[0]
dx, dy, dz = spacing
ox, oy, oz = origin

x = np.linspace(ox, ox + dx * (nx - 1), nx) * 10
y = np.linspace(oy, oy + dy * (ny - 1), ny) * 10
X, Y = np.meshgrid(x, y)

# === Estrai CELL_DATA, converti a volume ===
field_name = "Example_Cell_Scalar"
cell_data = grid.cell_data[field_name]

# === Converti cell data a point data ===
grid = grid.cell_data_to_point_data()

# === Ricostruisci campo 3D ===
data_3d = grid.point_data[field_name].reshape((nz, ny, nx))
z_index = nz // 2
slice_2d = data_3d[z_index, :, :]

# === Plot ===
fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(
    slice_2d,
    origin='lower',
    cmap='plasma',
    extent=(x[0], x[-1], y[0], y[-1]),
    interpolation='none'
)

#contours = ax.contour(
#    X, Y, slice_2d,
#    levels=25,
#    colors='black',
#    linewidths=1.0
#)

# Etichette assi e titolo
ax.set_xlabel("X (mm)", fontproperties=custom_font, fontsize=16)
ax.set_ylabel("Y (mm)", fontproperties=custom_font, fontsize=16)
ax.tick_params(axis='x', labelsize=14)  
ax.tick_params(axis='y', labelsize=14)  
#ax.set_title(f"Slice Z - Campo {field_name}", fontproperties=custom_font, fontsize=18)

# Colorbar
#cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
#cbar.ax.tick_params(labelsize=12)
#cbar.ax.set_ylabel(field_name, rotation=270, labelpad=15,
#                   fontproperties=custom_font, fontsize=16)

ax.set_aspect('equal')
plt.tight_layout()
filename_out = os.path.join(output_dir, "electrodes.png")
plt.savefig(filename_out, dpi=300)
plt.close()

print(f"Salvata immagine: {filename_out}")
