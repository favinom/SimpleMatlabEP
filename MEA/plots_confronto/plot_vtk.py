from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os


# === Parametri ===
vtk_dir = "../confronto"
frame_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
nrows, ncols = 2, 5

# === Setup della figura ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 8))
axes = axes.flatten()
last_im = None  # Per salvare l'ultimo imshow valido

# === Loop sui frame ===
for idx, i in enumerate(frame_indices):
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"bido_{istr}.vtk")

    if not os.path.exists(filename):
        print(f"File non trovato: {filename}")
        continue

    try:
        grid = pv.read(filename)
        dims = grid.dimensions[::-1]  # (z, y, x)
        spacing = grid.spacing
        origin = grid.origin

        nx, ny, nz = dims[2], dims[1], dims[0]
        dx, dy, dz = spacing
        ox, oy, oz = origin

        x = np.linspace(ox, ox + dx * (nx - 1), nx)
        y = np.linspace(oy, oy + dy * (ny - 1), ny)
        z_index = nz // 2

        v = grid.point_data["Transmembrane_Potential"]
        data_3d = v.reshape((nz, ny, nx))
        slice_2d = data_3d[z_index, :, :]

        ax = axes[idx]
        im = ax.imshow(
            slice_2d,
            origin="lower",
            cmap="plasma",
            extent=(x[0], x[-1], y[0], y[-1]),
            vmin=-80,
            vmax=40,
            interpolation="none"
        )
        ax.set_title(f"t = {i*20} ms", fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])

        last_im = im  # salva l'oggetto imshow valido

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")

# === Aggiungi barra colori comune se disponibile ===
if last_im is not None:
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(last_im, cax=cbar_ax, label="V [mV]")


fig.suptitle("Transmembrane Potential - bido", fontsize=18)  # Titolo globale
plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('bido_v', dpi=600)
plt.close()





# === Setup della figura ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 8))
axes = axes.flatten()
last_im = None  # Per salvare l'ultimo imshow valido

# === Loop sui frame ===
for idx, i in enumerate(frame_indices):
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"bido_{istr}.vtk")

    if not os.path.exists(filename):
        print(f"File non trovato: {filename}")
        continue

    try:
        grid = pv.read(filename)
        dims = grid.dimensions[::-1]  # (z, y, x)
        spacing = grid.spacing
        origin = grid.origin

        nx, ny, nz = dims[2], dims[1], dims[0]
        dx, dy, dz = spacing
        ox, oy, oz = origin

        x = np.linspace(ox, ox + dx * (nx - 1), nx)
        y = np.linspace(oy, oy + dy * (ny - 1), ny)
        z_index = nz // 2

        v = grid.point_data["Extracellular_Potential"]
        data_3d = v.reshape((nz, ny, nx))
        slice_2d = data_3d[z_index, :, :]

        ax = axes[idx]
        im = ax.imshow(
            slice_2d,
            origin="lower",
            cmap="plasma",
            extent=(x[0], x[-1], y[0], y[-1]),
            vmin=-4,
            vmax=4,
            interpolation="none"
        )
        ax.set_title(f"t = {i*20} ms", fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])

        last_im = im  # salva l'oggetto imshow valido

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")

# === Aggiungi barra colori comune se disponibile ===
if last_im is not None:
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(last_im, cax=cbar_ax, label="ue [mV]")

fig.suptitle("Extracellular Potential - bido", fontsize=18)  # Titolo globale

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('bido_ue', dpi=600)
plt.close()





# === Setup della figura ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 8))
axes = axes.flatten()
last_im = None  # Per salvare l'ultimo imshow valido

# === Loop sui frame ===
for idx, i in enumerate(frame_indices):
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"mea_{istr}.vtk")

    if not os.path.exists(filename):
        print(f"File non trovato: {filename}")
        continue

    try:
        grid = pv.read(filename)
        dims = grid.dimensions[::-1]  # (z, y, x)
        spacing = grid.spacing
        origin = grid.origin

        nx, ny, nz = dims[2], dims[1], dims[0]
        dx, dy, dz = spacing
        ox, oy, oz = origin

        x = np.linspace(ox, ox + dx * (nx - 1), nx)
        y = np.linspace(oy, oy + dy * (ny - 1), ny)
        z_index = nz // 2

        v = grid.point_data["Transmembrane_Potential"]
        data_3d = v.reshape((nz, ny, nx))
        slice_2d = data_3d[z_index, :, :]

        ax = axes[idx]
        im = ax.imshow(
            slice_2d,
            origin="lower",
            cmap="plasma",
            extent=(x[0], x[-1], y[0], y[-1]),
            vmin=-80,
            vmax=40,
            interpolation="none"
        )
        ax.set_title(f"t = {i*20} ms", fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])

        last_im = im  # salva l'oggetto imshow valido

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")

# === Aggiungi barra colori comune se disponibile ===
if last_im is not None:
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(last_im, cax=cbar_ax, label="V [mV]")


fig.suptitle("Transmembrane Potential - MEA", fontsize=18)  # Titolo globale

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('mea_v', dpi=600)
plt.close()





# === Setup della figura ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 8))
axes = axes.flatten()
last_im = None  # Per salvare l'ultimo imshow valido

# === Loop sui frame ===
for idx, i in enumerate(frame_indices):
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"mea_{istr}.vtk")

    if not os.path.exists(filename):
        print(f"File non trovato: {filename}")
        continue

    try:
        grid = pv.read(filename)
        dims = grid.dimensions[::-1]  # (z, y, x)
        spacing = grid.spacing
        origin = grid.origin

        nx, ny, nz = dims[2], dims[1], dims[0]
        dx, dy, dz = spacing
        ox, oy, oz = origin

        x = np.linspace(ox, ox + dx * (nx - 1), nx)
        y = np.linspace(oy, oy + dy * (ny - 1), ny)
        z_index = nz // 2

        v = grid.point_data["Extracellular_Potential"]
        data_3d = v.reshape((nz, ny, nx))
        slice_2d = data_3d[z_index, :, :]

        ax = axes[idx]
        im = ax.imshow(
            slice_2d,
            origin="lower",
            cmap="plasma",
            extent=(x[0], x[-1], y[0], y[-1]),
            vmin=-4,
            vmax=4,
            interpolation="none"
        )
        ax.set_title(f"t = {i*20} ms", fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])

        last_im = im  # salva l'oggetto imshow valido

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")

# === Aggiungi barra colori comune se disponibile ===
if last_im is not None:
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(last_im, cax=cbar_ax, label="ue [mV]")


fig.suptitle("Extracellular Potential - MEA", fontsize=18)  # Titolo globale

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('mea_ue', dpi=600)
plt.close()



# === Collage verticale dei potenziali di membrana ===
img1 = Image.open('bido_v.png')
img2 = Image.open('mea_v.png')
width = max(img1.width, img2.width)
height = img1.height + img2.height
collage_v = Image.new('RGB', (width, height), (255, 255, 255))
collage_v.paste(img1, (0, 0))
collage_v.paste(img2, (0, img1.height))
collage_v.save('collage_v.png')

# === Collage verticale dei potenziali extracellulari ===
img3 = Image.open('bido_ue.png')
img4 = Image.open('mea_ue.png')
width_ue = max(img3.width, img4.width)
height_ue = img3.height + img4.height
collage_ue = Image.new('RGB', (width_ue, height_ue), (255, 255, 255))
collage_ue.paste(img3, (0, 0))
collage_ue.paste(img4, (0, img3.height))
collage_ue.save('collage_ue.png')