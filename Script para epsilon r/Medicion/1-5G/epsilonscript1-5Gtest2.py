import skrf as rf
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# --- PARÁMETROS FÍSICOS DEL RESONADOR A 1GHZ PARA ROGERS4350B ---
R_rogers = 19.13e-3             # Radio para 1.5GHZ
L_rogers = 2 * np.pi * R_rogers # longitud del anillo
c = 3e8                         # Velocidad de la luz en vacio en m/s
H_rogers = 1.524e-3             # Altura del sustrato
w_rogers = 3.3e-3               # Ancho de las pistas en metros (ejemplo 3 mm)

# --- PARÁMETROS FÍSICOS DEL RESONADOR A 1GHZ PARA FR4 ---
R_fr4 = 17.657e-3                # Radio para 1.5GHZ
L_fr4 = 2 * np.pi * R_fr4       # longitud del anillo
# c = 3e8 (ya definido)
H_fr4 = 1.56e-3                 # Altura del sustrato
w_fr4 = 2.85e-3                  # Ancho de las pistas en metros (ejemplo 3 mm)

# Modos de resonancia --> hasta 4.5GHZ (usando los primeros 3 modos n=1, 2, 3)
n_modos = [1, 2, 3]

# Array para los epsilon_eff y epsilon_r
epsilon_eff_rogers = np.zeros(3)
epsilon_r_rogers = np.zeros(3)

epsilon_eff_fr4 = np.zeros(3)
epsilon_r_fr4 = np.zeros(3)

# Variables para ajuste de busqueda de picos
distance_val = 100
prominence_val = 10

# --- CARGAR ARCHIVOS S2P ---
try:
    ntwk_rogers = rf.Network('Rogers_1-5GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'Rogers_1-5GHz_puntos.s2p'.")
    exit()

try:
    ntwk_fr4 = rf.Network('FR4_1-5GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'FR4_1-5GHz_puntos.s2p'.")
    exit()

print("Archivos cargados exitosamente.")

# --- Graficar ---
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
ntwk_rogers.s21.plot_s_db(ax=ax[0], label='$|S_{21}|$ (S21 Rogers 1.5GHZ)')
ax[0].set_title('ROGERS 4350B $|S_{21}|$')
ax[0].set_xlabel('Frecuencia (Hz)')
ax[0].set_ylabel('Magnitud (dB)')
ax[0].grid(True, which="both", ls="--", alpha=0.6)
ax[0].legend()

ntwk_fr4.s21.plot_s_db(ax=ax[1],label='$|S_{21}|$ (S21 FR4 1.5GHZ)')
ax[1].set_title('FR4 $|S_{21}|$')
ax[1].set_xlabel('Frecuencia (Hz)')
ax[1].set_ylabel('Magnitud (dB)')
ax[1].grid(True, which="both", ls="--", alpha=0.6)
ax[1].legend()

# Mostrar el gráfico
plt.suptitle('Comparación de Parámetros $|S_{21}|$ (Transmisión)') # Título general de la figura
plt.tight_layout(rect=[0, 0, 1, 0.96]) # Ajusta el espacio para el título principal
plt.show()

# --- PROCESAMIENTO ROGERS 4350B ---

frecuencias_rogers = ntwk_rogers.f  # En Hz
frecuencias_rogers_GHz = frecuencias_rogers / 1e9
s21_rogers = ntwk_rogers.s[:, 1, 0]
s21_dB_rogers = 20 * np.log10(np.abs(s21_rogers))

# Filtrar hasta 3.5 GHz
mask = frecuencias_rogers_GHz <= 5
frecuencias_filtradas_rogers = frecuencias_rogers_GHz[mask]
s21_dB_filtrado_rogers = s21_dB_rogers[mask]

# Detectar picos
picos_resonancia_rogers, _ = find_peaks(s21_dB_filtrado_rogers, distance=distance_val, prominence=prominence_val)
picos_a_usar_rogers = picos_resonancia_rogers[:len(n_modos)]
frecuencias_resonancia_rogers = frecuencias_filtradas_rogers[picos_a_usar_rogers]
frecuencias_resonancia_Hz_rogers = frecuencias_resonancia_rogers * 1e9

# Calcular epsilon efectivo y relativo
for i in range(3):
	epsilon_eff_rogers[i] = (n_modos[i] * c / (2 * np.pi * R_rogers * frecuencias_resonancia_Hz_rogers[i])) ** 2

for i in range(3):
	epsilon_r_rogers[i] = (2 * epsilon_eff_rogers[i] - (1 - (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))) / (1 + (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))

epsilon_r_promedio_rogers = np.mean(epsilon_r_rogers)

print("\n--- Resultados ROGERS 4350B ---")
print(f"Frecuencias de Resonancia (GHz): {frecuencias_resonancia_rogers}")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_rogers}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_rogers}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_rogers:.3f}")

# --- PROCESAMIENTO FR4 ---

frecuencias_fr4 = ntwk_fr4.f  # En Hz
frecuencias_fr4_GHz = frecuencias_fr4 / 1e9
s21_fr4 = ntwk_fr4.s[:, 1, 0]
s21_dB_fr4 = 20 * np.log10(np.abs(s21_fr4))

# Filtrar hasta 3.5 GHz
mask = frecuencias_fr4_GHz <= 5
frecuencias_filtradas_fr4 = frecuencias_fr4_GHz[mask]
s21_dB_filtrado_fr4 = s21_dB_fr4[mask]

# Detectar picos
picos_resonancia_fr4, _ = find_peaks(s21_dB_filtrado_fr4, distance=distance_val, prominence=prominence_val)
picos_a_usar_fr4 = picos_resonancia_fr4[:len(n_modos)]
frecuencias_resonancia_fr4 = frecuencias_filtradas_fr4[picos_a_usar_fr4]
frecuencias_resonancia_Hz_fr4 = frecuencias_resonancia_fr4 * 1e9

# Calcular epsilon efectivo y relativo
for i in range(3):
	epsilon_eff_fr4[i] = (n_modos[i] * c / (2 * np.pi * R_fr4 * frecuencias_resonancia_Hz_fr4[i])) ** 2

for i in range(3):
	epsilon_r_fr4[i] = (2 * epsilon_eff_fr4[i] - (1 - (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))) / (1 + (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))

epsilon_r_promedio_fr4 = np.mean(epsilon_r_fr4)

print("\n--- Resultados FR4 ---")
print(f"Frecuencias de Resonancia (GHz): {frecuencias_resonancia_fr4}")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_fr4}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_fr4}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_fr4:.3f}")

# --- GRÁFICO DE DEBUG (S21 con Picos Detectados) ---

fig, axs = plt.subplots(1, 2, figsize=(14, 6))
plt.suptitle(f'Depuración de Detección de Picos (Distancia={distance_val}, Prominencia={prominence_val})', fontsize=16)

# --- Subplot 1: ROGERS Debug ---
axs[0].plot(frecuencias_filtradas_rogers, s21_dB_filtrado_rogers, label='$|S_{21}|$ ROGERS')
picos_y_rogers = s21_dB_filtrado_rogers[picos_a_usar_rogers]
axs[0].plot(frecuencias_resonancia_rogers, picos_y_rogers, "x", color='red', markersize=10, label='Picos Detectados')
axs[0].set_title('ROGERS 4350B')
axs[0].set_xlabel('Frecuencia (GHz)')
axs[0].set_ylabel('Magnitud ($|S_{21}|$ en dB)')
axs[0].grid(True, which="both", ls="--", alpha=0.6)
axs[0].legend()

# --- Subplot 2: FR4 Debug ---
axs[1].plot(frecuencias_filtradas_fr4, s21_dB_filtrado_fr4, label='$|S_{21}|$ FR4', color='C1')
picos_y_fr4 = s21_dB_filtrado_fr4[picos_a_usar_fr4]
axs[1].plot(frecuencias_resonancia_fr4, picos_y_fr4, "x", color='red', markersize=10, label='Picos Detectados')
axs[1].set_title('FR4')
axs[1].set_xlabel('Frecuencia (GHz)')
axs[1].set_ylabel('Magnitud ($|S_{21}|$ en dB)')
axs[1].grid(True, which="both", ls="--", alpha=0.6)
axs[1].legend()

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# --- GRÁFICO COMPARATIVO ---
plt.figure(figsize=(10, 6))

# Curva ROGERS 4350B
plt.plot(frecuencias_resonancia_rogers, epsilon_r_rogers, 'o', color='r', markersize=10, label='εr ROGERS')
plt.axhline(y=epsilon_r_promedio_rogers, color='r', linestyle='--', label=f'Promedio εr ROGERS = {epsilon_r_promedio_rogers:.3f}')

# Etiquetas en cada punto ROGERS
for f, eps in zip(frecuencias_resonancia_rogers, epsilon_r_rogers):
    plt.text(f, eps + 0.03, f"{eps:.3f}\n{f:.3f} GHz", color='r', ha='center', va='bottom', fontsize=12)

# Curva FR4
plt.plot(frecuencias_resonancia_fr4, epsilon_r_fr4, 'o', color='b',markersize=10, label='εr FR4')
plt.axhline(y=epsilon_r_promedio_fr4, color='b', linestyle='--', label=f'Promedio εr FR4 = {epsilon_r_promedio_fr4:.3f}')

for f, eps in zip(frecuencias_resonancia_fr4, epsilon_r_fr4):
    plt.text(f, eps - 0.05, f"{eps:.3f}\n{f:.3f} GHz", color='b', ha='center', va='top', fontsize=12)

# Etiquetas y formato
plt.xlabel('Frecuencias de resonancia (GHz)')
plt.ylabel('εr(f)')
plt.title('Cálculo de Permitividad Relativa (εr): ROGERS 4350B vs FR4')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.ylim(3, 5) # Rango visual ajustado
plt.legend(title='Materiales')
plt.tight_layout()
plt.show()