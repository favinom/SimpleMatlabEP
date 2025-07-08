import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Carica i dati
bido = pd.read_csv('../confronto/bido_traces.csv')
mea = pd.read_csv('../confronto/mea_traces.csv')  # Cambia in .cdv se necessario

# Lista dei canali
channels_V = [f'V{i}_mV' for i in range(1, 10)]
channels_u = [f'u{i}_mV' for i in range(1, 10)]

# Crea la figura
fig, axs = plt.subplots(2, 2, figsize=(16, 10), gridspec_kw={'height_ratios': [1, 1]})

# --- Plot sopra: segnali ---
for ch in channels_V:
    axs[0,0].plot(bido['Time_ms'], bido[ch], label=ch)
    axs[0,1].plot(mea['Time_ms'], mea[ch], label=ch)

axs[0,0].set_title('Bido traces')
axs[0,1].set_title('Mea traces')
axs[0,0].set_ylabel('mV')
#axs[0,0].legend(loc='upper right', fontsize=8)
axs[0,1].legend(loc='upper right', fontsize=8)

# --- Plot sotto: differenze ---
for ch in channels_u:
    axs[1,0].plot(bido['Time_ms'], bido[ch], label=ch)
    axs[1,1].plot(mea['Time_ms'], mea[ch], label=ch)

axs[1,0].set_title('Bido traces')
axs[1,1].set_title('Mea traces')
axs[1,0].set_ylabel('mV')
#axs[1,0].legend(loc='upper right', fontsize=8)
#axs[1,1].legend(loc='upper right', fontsize=8)


# Nascondi subplot vuoto
#axs[1,1].axis('off')
#plt.grid(True)
plt.savefig(os.path.join('.', 'V_ue.png'))
plt.close()

