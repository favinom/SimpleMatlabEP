import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parametri
mesh_sizes = [64, 128, 256]#, 512, 1024]
nsd_x = [2, 3, 4, 5, 6, 7, 8, 9, 10]
results_dir = 'results'

# Liste dati
h_list = []
H_list = []
iter_list = []

# Ciclo su tutti i parametri
for mesh_size in mesh_sizes:
    h_piccolo = 1.0 / mesh_size
    for n in nsd_x:
        H_grande = 1.0 / n
        filename = f'results_iter_{h_piccolo:.4f}_{H_grande:.4f}.csv'
        filepath = os.path.join(filename)
        
        if os.path.exists(filepath):
            df = pd.read_csv(filepath, header=None)
            if not df.empty:
                #iter_value = df.iloc[0, 1]  # secondo valore della prima riga
                iter_value = df.iloc[0, 3]  # quart valore della prima riga
                h_list.append(h_piccolo)
                H_list.append(H_grande)
                iter_list.append(iter_value)

# Conversione in array numpy
h_array = np.array(h_list)
H_array = np.array(H_list)
iter_array = np.array(iter_list)

# Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot
sc = ax.scatter(h_array, H_array, iter_array, c=iter_array, cmap='viridis', s=50, label='Iterazioni SP')

# Superficie teorica
h_vals = np.unique(h_array)
H_vals = np.unique(H_array)
hh, HH = np.meshgrid(h_vals, H_vals)
zz = 0.03 / (hh * HH) + 30

# Plot superficie
#ax.plot_surface(hh, HH, zz, alpha=0.4, color='red')

# Etichette
ax.set_xlabel('h (mesh size)')
ax.set_ylabel('H (subdomain size)')
ax.set_zlabel('Numero iterazioni SP')
ax.set_title('Iterazioni vs 1/(h·H)')

# Legenda manuale
scatter_proxy = plt.Line2D([0], [0], linestyle="none", marker='o', color='blue')
surface_proxy = plt.Rectangle((0,0),1,1,fc="red", alpha=0.4)
#ax.legend([scatter_proxy, surface_proxy], ['Iterazioni SP', '1/(h·H)+30'])

plt.tight_layout()
plt.show()
