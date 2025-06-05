import os
import pandas as pd
import matplotlib.pyplot as plt

# Percorso del file CSV
csv_path = '../outputs/voltage_traces.csv'

# Leggi il file CSV
df = pd.read_csv(csv_path)

# Crea la cartella per i plot
plot_dir = '.' #'../outputs'
os.makedirs(plot_dir, exist_ok=True)

# Imposta dimensione testo generale
label_fontsize = 20
title_fontsize = 22
legend_fontsize = 20
tick_labelsize = 20 


# Plot 1: Voltage_mV (Vsaved)
plt.figure(figsize=(10, 12))
plt.plot(df['Time_ms'], df['Voltage_mV'], label='Vsaved')
plt.xlabel('Time (ms)', fontsize=label_fontsize)
plt.ylabel('Voltage (mV)', fontsize=label_fontsize)
plt.title('V', fontsize=title_fontsize)
plt.tick_params(axis='both', labelsize=tick_labelsize)

plt.legend(fontsize=legend_fontsize)
plt.grid(True)
plt.savefig(os.path.join(plot_dir, 'Vsaved_plot.png'))
plt.close()

# Plot 2: ue_mV (usaved)
plt.figure(figsize=(10, 6))
plt.plot(df['Time_ms'], df['ue_mV'], label='usaved', color='orange')
plt.xlabel('Time (ms)', fontsize=label_fontsize)
plt.ylabel('ue (mV)', fontsize=label_fontsize)
plt.title('ue', fontsize=title_fontsize)
plt.tick_params(axis='both', labelsize=tick_labelsize)
plt.legend(fontsize=legend_fontsize)
plt.grid(True)
plt.savefig(os.path.join(plot_dir, 'uesaved_plot.png'))
plt.close()

# Plot 3: Tutti i campi U{k} (U0_mV, U1_mV, ...)
u_columns = [col for col in df.columns if col.startswith('U') and col.endswith('_mV')]

# Dizionario per etichette leggibili
label_map = {
    'U1_mV': 'FP in e1',
    'U5_mV': 'FP in e5',
    'U9_mV': 'FP in e9',
}

plt.figure(figsize=(10, 11))
for col in u_columns:
    if col in ['U1_mV', 'U5_mV', 'U9_mV']:
        plt.plot(df['Time_ms'], df[col], label=label_map[col], linewidth=3.5)
plt.xlabel('Time (ms)', fontsize=label_fontsize)
plt.ylabel('FP (mV)', fontsize=label_fontsize)
plt.title('FP signals', fontsize=title_fontsize)
plt.tick_params(axis='both', labelsize=tick_labelsize)
plt.legend(loc='best', fontsize=legend_fontsize)
plt.grid(True)
plt.savefig(os.path.join(plot_dir, 'FP_signals_plot.png'))
plt.close()

print(f"Plot salvati nella cartella: {plot_dir}")
