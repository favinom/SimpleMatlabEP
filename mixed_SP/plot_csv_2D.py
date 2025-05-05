import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parametri
mesh_sizes = [64, 128, 256, 512]
nsd_x = [2, 3, 4, 5, 6, 7, 8, 9, 10]
results_dir = 'results'

# Organizza i dati per h
data_per_h = {}

for mesh_size in mesh_sizes:
    h = 1.0 / mesh_size
    H_over_h_list = []
    iter_list = []
    for n in nsd_x:
        H = 1.0 / n
        H_over_h = H #1 / (H * h)
        filename = f'results_iter_{h:.4f}_{H:.4f}.csv'
        filepath = os.path.join(filename)

        if os.path.exists(filepath):
            df = pd.read_csv(filepath, header=None)
            if not df.empty:
                iter_value = df.iloc[0, 3]  # quarta colonna
                H_over_h_list.append(H_over_h)
                iter_list.append(iter_value)

    if H_over_h_list:
        data_per_h[h] = (np.array(H_over_h_list), np.array(iter_list))

# Plot 2D
plt.figure(figsize=(10, 6))
colors = plt.cm.viridis(np.linspace(0, 1, len(data_per_h)))

for idx, (h, (H_over_h_vals, iters)) in enumerate(data_per_h.items()):
    sorted_indices = np.argsort(H_over_h_vals)
    x_sorted = H_over_h_vals[sorted_indices]
    y_sorted = iters[sorted_indices]
    plt.plot(x_sorted, y_sorted, '-o', color=colors[idx], label=f'h = {h:.4f}')
    
    # Curva teorica: C(1 + log^2(H/h))
    C = y_sorted[0] / (1 + np.log(x_sorted[0])**2)  # fit la costante
    theory_curve = C * (1 + np.log(x_sorted)**2)
    #plt.plot(x_sorted, theory_curve, '--', color=colors[idx], label=f'Teoria (h={h:.4f})')


#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$Hh$')
plt.ylabel('Numero iterazioni SP')
plt.title(r'Confronto tra iterazioni osservate e teoriche ($\sim 1 + \log^2(Hh)$)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
