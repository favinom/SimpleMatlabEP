import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Cartella dove hai salvato i file CSV
results_dir = 'results'

# Lista di tutti i file CSV
csv_files = glob.glob(os.path.join('results.csv'))

# Liste per i dati
h_list = []
H_list = []
iter_list = []

# Legge ogni file e accumula i dati
for file in csv_files:
    df = pd.read_csv(file)
    h_list.extend(df['h'])
    H_list.extend(df['H'])
    iter_list.extend(df['iter_dd'])

## Plot 3D
#fig = plt.figure(figsize=(10, 7))
#ax = fig.add_subplot(111, projection='3d')
#
#print(h_list)
#ax.scatter(h_list, H_list, iter_list, c=iter_list, cmap='viridis', s=100)
#
#ax.set_xlabel('h (mesh size)')
#ax.set_ylabel('H (subdomain size)')
#ax.set_zlabel('Numero iterazioni SP')
#ax.set_title('Convergenza metodo Steklov-Poincaré')
#
#plt.tight_layout()
#plt.show()

# Conversione in numpy array
h_array = np.array(h_list)
H_array = np.array(H_list)
iter_array = np.array(iter_list)

# Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot dei dati computazionali
sc = ax.scatter(h_array, H_array, iter_array, c=iter_array, cmap='viridis', s=50, label='Iterazioni SP')

# Meshgrid per la superficie teorica z = 1/(h*H)
h_vals = np.unique(h_array)
H_vals = np.unique(H_array)
hh, HH = np.meshgrid(h_vals, H_vals)
zz = 0.03 / (hh * HH) + 30

# Plot superficie teorica
ax.plot_surface(hh, HH, zz, alpha=0.4, color='red', label='1/(h·H)')

# Etichette e titolo
ax.set_xlabel('h (mesh size)')
ax.set_ylabel('H (subdomain size)')
ax.set_zlabel('Numero iterazioni SP')
ax.set_title('Iterazioni vs 1/(h·H)')

# Legenda manuale
scatter_proxy = plt.Line2D([0], [0], linestyle="none", marker='o', color='blue')
surface_proxy = plt.Rectangle((0,0),1,1,fc="red", alpha=0.4)
ax.legend([scatter_proxy, surface_proxy], ['Iterazioni SP', '1/(h·H)+30'])

plt.tight_layout()
plt.show()