import skrf as rf
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# Parámetros físicos del resonador a 1GHZ para ROGERS4350B
R_rogers = 28.2795114977997e-3             # Radio para 1GHZ
L_rogers = 2 * np.pi * R_rogers # longitud del anillo
c = 3e8                         # Velocidad de la luz en vacio en m/s
H_rogers = 1.524e-3             # Altura del sustrato 
w_rogers = 3.309e-3             # Ancho de las pistas en metros (ejemplo 3 mm)

# Parámetros físicos del resonador a 1GHZ para FR4
R_fr4 = 26.4853775112596e-3                # Radio para 1GHZ
L_fr4 = 2 * np.pi * R_fr4       # longitud del anillo
c = 3e8                         # Velocidad de la luz en vacio en m/s
H_fr4 = 1.56e-3                 # Altura del sustrato 
w_fr4 = 2.85e-3                 # Ancho de las pistas en metros (ejemplo 3 mm)

# Modos de resonancia --> hasta 6GHZ
n_modos = [1,2,3]

# Array para los epsilon_eff y epsilon_r
epsilon_eff_rogers = np.zeros(3)
epsilon_r_rogers = np.zeros(3)

epsilon_eff_fr4 = np.zeros(3)
epsilon_r_fr4 = np.zeros(3)

# Cargar el archivo .s2p 
ntwk_rogers = rf.Network('S21-Rogers-1G.s2p')
ntwk_fr4 = rf.Network('S21-FR4-1G.s2p')

# Extraer el vector de frecuencias (en Hz)
frecuencias_rogers = ntwk_rogers.f  # En Hz
frecuencias_rogers_GHz = frecuencias_rogers / 1e9  # Opcional: convertir a GHz

frecuencias_fr4 = ntwk_fr4.f  # En Hz
frecuencias_fr4_GHz = frecuencias_fr4 / 1e9        # Opcional: convertir a GHz

# Obtener el parámetro S21 (índice [0,1]) en dB
s21_rogers = ntwk_rogers.s[:, 1, 0]                # Complejo
s21_dB_rogers = 20 * np.log10(np.abs(s21_rogers))  # Módulo en dB

s21_fr4 = ntwk_fr4.s[:, 1, 0]                      # Complejo
s21_dB_fr4 = 20 * np.log10(np.abs(s21_fr4))        # Módulo en dB

# Filtrar hasta 6 GHz (6 picos para fundamental de 1GHZ)
mask = frecuencias_rogers_GHz <= 3.5
frecuencias_filtradas_rogers = frecuencias_rogers_GHz[mask]
s21_dB_filtrado_rogers = s21_dB_rogers[mask]

mask = frecuencias_fr4_GHz <= 3.5
frecuencias_filtradas_fr4 = frecuencias_fr4_GHz[mask]
s21_dB_filtrado_fr4 = s21_dB_fr4[mask]

# Detectar picos en el rango filtrado
picos_resonancia_rogers, _ = find_peaks(s21_dB_filtrado_rogers, distance=30, prominence=3)
frecuencias_resonancia_rogers = frecuencias_filtradas_rogers[picos_resonancia_rogers]
frecuencias_resonancia_Hz_rogers = frecuencias_resonancia_rogers * 1e9  # Convertir GHz a Hz

picos_resonancia_fr4, _ = find_peaks(s21_dB_filtrado_fr4, distance=30, prominence=3)
frecuencias_resonancia_fr4 = frecuencias_filtradas_fr4[picos_resonancia_fr4]
frecuencias_resonancia_Hz_fr4 = frecuencias_resonancia_fr4 * 1e9  # Convertir GHz a Hz

# Calcular epsilon efectivo y relativo para cada frecuencia
for i in range(3): # i va de 0 a 5
	epsilon_eff_rogers[i] = (n_modos[i] * c / (2 * np.pi * R_rogers * frecuencias_resonancia_Hz_rogers[i])) ** 2
for i in range(3): # i va de 0 a 5
	epsilon_r_rogers[i] = (2 * epsilon_eff_rogers[i] - (1 - (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))) / (1 + (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))

for i in range(3): # i va de 0 a 5
	epsilon_eff_fr4[i] = (n_modos[i] * c / (2 * np.pi * R_fr4 * frecuencias_resonancia_Hz_fr4[i])) ** 2
for i in range(3): # i va de 0 a 5
	epsilon_r_fr4[i] = (2 * epsilon_eff_fr4[i] - (1 - (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))) / (1 + (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))

# Calcular el promedio de epsilon_r
epsilon_r_promedio_rogers = np.mean(epsilon_r_rogers)
epsilon_r_promedio_fr4 = np.mean(epsilon_r_fr4)

print("\n--- Resultados ROGERS 4350B ---")
print(f"Frecuencias de Resonancia (GHz): {frecuencias_resonancia_rogers}")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_rogers}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_rogers}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_rogers:.3f}")

print("\n--- Resultados FR4 ---")
print(f"Frecuencias de Resonancia (GHz): {frecuencias_resonancia_fr4}")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_fr4}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_fr4}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_fr4:.3f}")

#Grafico
plt.figure()
# Curva ROGERS 4350B
plt.plot(frecuencias_resonancia_rogers, epsilon_r_rogers, color='r', label='εr ROGERS')
plt.axhline(y=epsilon_r_promedio_rogers, color='r', linestyle='--', label=f'Promedio εr ROGERS = {epsilon_r_promedio_rogers:.3f}')
# Curva FR4
plt.plot(frecuencias_resonancia_fr4, epsilon_r_fr4, color='b', label='εr FR4')
plt.axhline(y=epsilon_r_promedio_fr4, color='b', linestyle='--', label=f'Promedio εr FR4 = {epsilon_r_promedio_fr4:.3f}')
# Etiquetas y formato
plt.xlabel('Frecuencias de resonancia(GHz)')
plt.ylabel('εr(f)')
plt.title('ROGERS 4350B vs FR4 a 1GHZ')
plt.grid(True)
plt.ylim(3, 5)
plt.legend(title='Materiales')
plt.tight_layout()
plt.show()