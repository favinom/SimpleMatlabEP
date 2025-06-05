import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.ticker as ticker  
import pyvista as pv
import os
from moviepy.editor import ImageSequenceClip
from moviepy.editor import ImageClip, concatenate_videoclips
from PIL import Image, ImageDraw, ImageFont
import glob


# === Parametri ===
output_dir = "./frames"
vtk_dir = "../outputs"
os.makedirs(output_dir, exist_ok=True)
font_path = "./fonts/Helvetica.ttf"  # <-- aggiorna se necessario
custom_font = fm.FontProperties(fname=font_path)

# === Funzione per plottare un campo ===
def plot_field_matplotlib(grid, field_name, field_label, x, y, z_index, filename_out, cmap='plasma', clim=None, show_grid=True):
    data_3d = grid.point_data[field_name].reshape((nz, ny, nx))
    slice_2d = data_3d[z_index, :, :]
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(
        slice_2d,
        origin='lower',
        cmap=cmap,
        extent=(x[0], x[-1], y[0], y[-1]),
        interpolation='none',
        vmin=clim[0] if clim else None,
        vmax=clim[1] if clim else None
    )

    # Contorni neri
    #contours = ax.contour(X, Y, slice_2d, levels=25, colors='black', linewidths=0.8)

    # Etichette e titoli
    ax.set_xlabel("X (cm)", fontproperties=custom_font, fontsize=14)
    ax.set_ylabel("Y (cm)", fontproperties=custom_font, fontsize=14)
    #ax.set_title(f"Slice Z - Campo {field_name}", fontproperties=custom_font, fontsize=18)
    ax.set_title(f"t =  {int(istr) * 10} ms", fontproperties=custom_font, fontsize=18)

    # Barra colori
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(field_label, rotation=270, labelpad=15,
                       fontproperties=custom_font, fontsize=16)

    # === Aggiungi griglia ===
    # === Aggiungi griglia e controlla densità tick ===
    if show_grid:
        nticks = 10  # <-- cambia questo numero per più o meno tick
        xticks = np.linspace(x[0], x[-1], nticks)
        yticks = np.linspace(y[0], y[-1], nticks)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        ax.grid(True, which='both', color='gray', linestyle='--', linewidth=0.5)
        
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(filename_out, dpi=600)
    plt.close()
    print(f"Salvata immagine: {filename_out}")



# === Loop sui file ===
for i in range(100):
    istr = f"{i:04d}"
    filename = os.path.join(vtk_dir, f"sol_{istr}.vtk")
    if not os.path.isfile(filename):
        print(f"File non trovato: {filename}")
        continue

    print(f"Processing: {filename}")

    try:
        grid = pv.read(filename)

        # Assicurati che sia strutturato e con dati scalari
        if not grid.point_data:
            print(f"Nessun dato scalare in {filename}")
            continue

        field_names = list(grid.point_data.keys())
        if len(field_names) < 2:
            print(f"Meno di 2 campi scalari in {filename}")
            continue

        # Estrai geometria
        dims = grid.dimensions[::-1]  # (z, y, x)
        spacing = grid.spacing
        origin = grid.origin

        nx, ny, nz = dims[2], dims[1], dims[0]
        dx, dy, dz = spacing
        ox, oy, oz = origin

        x = np.linspace(ox, ox + dx * (nx - 1), nx)
        y = np.linspace(oy, oy + dy * (ny - 1), ny)
        z_index = nz // 2

        # Plot primo campo
        field_v = field_names[0]
        field_label = "V [mV]"
        out_v = os.path.join(output_dir, f"v_{istr}_new.png")
        plot_field_matplotlib(grid, field_v, field_label, x, y, z_index, out_v, cmap='plasma', clim=[-80, 40])

        # Plot secondo campo
        field_u = field_names[1]
        field_label = "ue [mV]"
        out_u = os.path.join(output_dir, f"u_{istr}_new.png")
        plot_field_matplotlib(grid, field_u, field_label, x, y, z_index, out_u, cmap='turbo', clim=[-4, 4])

    except Exception as e:
        print(f"Errore nel file {filename}: {e}")


# === Crea un video dalle immagini ===
def crea_video_da_immagini(cartella_frames, nome_video, fps=5, pattern="v_*.png"):
    import glob
    file_paths = sorted(glob.glob(os.path.join(cartella_frames, pattern)))
    if not file_paths:
        print("Nessuna immagine trovata per creare il video.")
        return

    print(f"Generazione video da {len(file_paths)} immagini...")
    clip = ImageSequenceClip(file_paths, fps=fps)
    clip.write_videofile(nome_video, codec='libx264', audio=False)



def crea_video_accoppiato(cartella_frames, nome_video, fps=5, font_path="./fonts/Helvetica.ttf"):
    v_paths = sorted(glob.glob(os.path.join(cartella_frames, "v_*_new.png")))
    u_paths = sorted(glob.glob(os.path.join(cartella_frames, "u_*_new.png")))

    if not v_paths or not u_paths or len(v_paths) != len(u_paths):
        print("Numero di immagini v e u non corrispondente o mancante.")
        return

    font_size = 100
    font = ImageFont.truetype(font_path, font_size)
    clips = []

    for v_img, u_img in zip(v_paths, u_paths):
        # Carica immagini
        im_v = Image.open(v_img)
        im_u = Image.open(u_img)

        # Combina le immagini affiancandole
        width, height = im_v.size
        combined = Image.new("RGB", (width * 2, height + 60), "white")

        # Incolla immagini
        combined.paste(im_v, (0, 60))
        combined.paste(im_u, (width, 60))

        ## Estrai tempo da nome file
        #frame_number = int(v_img.split("_")[-2])
        #tempo = f"t = {frame_number * 10} ms"

        ## Aggiungi il tempo centrato in alto
        #draw = ImageDraw.Draw(combined)
        #bbox = font.getbbox(tempo)
        #text_width = bbox[2] - bbox[0]
        #text_height = bbox[3] - bbox[1]
        #draw.text(
        #    ((width * 2 - text_width) / 2, 10),
        #    tempo,
        #    font=font,
        #    fill="black"
        #)


        # Converti a MoviePy ImageClip
        img_clip = ImageClip(np.array(combined)).set_duration(1 / fps)
        clips.append(img_clip)

    # Crea e salva il video
    final_clip = concatenate_videoclips(clips, method="compose")
    final_clip.write_videofile(nome_video, fps=fps, codec="libx264", audio=False)



# === Esegui ===
#crea_video_da_immagini(output_dir, "campo_v_video.mp4", fps=5, pattern="v_*_new.png")
#crea_video_da_immagini(output_dir, "campo_u_video.mp4", fps=5, pattern="u_*_new.png")

crea_video_accoppiato(output_dir, "campo_composito.mp4", fps=5, font_path=font_path)
