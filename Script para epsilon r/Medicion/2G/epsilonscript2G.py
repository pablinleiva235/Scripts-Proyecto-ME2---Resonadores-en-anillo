import skrf as rf
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# --- PARÁMETROS FÍSICOS DEL RESONADOR A 2GHZ PARA ROGERS4350B ---
R_rogers = 14.347e-3             # Radio para 2GHZ
L_rogers = 2 * np.pi * R_rogers # longitud del anillo
c = 3e8                         # Velocidad de la luz en vacio en m/s
H_rogers = 1.524e-3             # Altura del sustrato
w_rogers = 3.3e-3               # Ancho de las pistas en metros (ejemplo 3 mm)

# --- PARÁMETROS FÍSICOS DEL RESONADOR A 2GHZ PARA FR4 ---
R_fr4 = 13.291e-3                # Radio para 2GHZ
L_fr4 = 2 * np.pi * R_fr4       # longitud del anillo
# c = 3e8 (ya definido)
H_fr4 = 1.56e-3                 # Altura del sustrato
w_fr4 = 2.6e-3                  # Ancho de las pistas en metros (ejemplo 3 mm)

# Modos de resonancia --> hasta 6GHZ (usando los primeros 3 modos n=1, 2, 3)
n_modos = [1, 2, 3]

# Array para los epsilon_eff y epsilon_r
epsilon_eff_rogers = np.zeros(3)
epsilon_r_rogers = np.zeros(3)

epsilon_eff_fr4 = np.zeros(3)
epsilon_r_fr4 = np.zeros(3)

frecuencias_resonancia_rogers = [1.966432, 3.902389, 5.789194]
# CONVERSIÓN a array de NumPy para permitir la multiplicación escalar
frecuencias_resonancia_array_rogers = np.array(frecuencias_resonancia_rogers)

# Multiplicación elemento por elemento
frecuencias_resonancia_Hz_rogers = frecuencias_resonancia_array_rogers * 1e9

# --- Datos de FR4 ---
frecuencias_resonancia_fr4 = [2.036855, 4.045878, 6.002182]
# CONVERSIÓN a array de NumPy
frecuencias_resonancia_array_fr4 = np.array(frecuencias_resonancia_fr4)

# Multiplicación elemento por elemento
frecuencias_resonancia_Hz_fr4 = frecuencias_resonancia_array_fr4 * 1e9

# Calcular epsilon efectivo y relativo
for i in range(3):
	epsilon_eff_rogers[i] = (n_modos[i] * c / (2 * np.pi * R_rogers * frecuencias_resonancia_Hz_rogers[i])) ** 2

for i in range(3):
	epsilon_r_rogers[i] = (2 * epsilon_eff_rogers[i] - (1 - (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))) / (1 + (1 / np.sqrt(1 + 12 * H_rogers / w_rogers)))

epsilon_eff_promedio_rogers = np.mean(epsilon_eff_rogers)
epsilon_r_promedio_rogers = np.mean(epsilon_r_rogers)

# Calcular epsilon efectivo y relativo
for i in range(3):
	epsilon_eff_fr4[i] = (n_modos[i] * c / (2 * np.pi * R_fr4 * frecuencias_resonancia_Hz_fr4[i])) ** 2

for i in range(3):
	epsilon_r_fr4[i] = (2 * epsilon_eff_fr4[i] - (1 - (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))) / (1 + (1 / np.sqrt(1 + 12 * H_fr4 / w_fr4)))

epsilon_eff_promedio_fr4 = np.mean(epsilon_eff_fr4)
epsilon_r_promedio_fr4 = np.mean(epsilon_r_fr4)

print("\n--- Resultados ROGERS 4350B ---")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_rogers}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_rogers}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_rogers:.3f}")

print("\n--- Resultados FR4 ---")
print(f"Valores de Epsilon efectivo (εeff): {epsilon_eff_fr4}")
print(f"Valores de Epsilon Relativo (εr): {epsilon_r_fr4}")
print(f"Epsilon Relativo Promedio (εr promedio): {epsilon_r_promedio_fr4:.3f}")

# --- GRÁFICO COMPARATIVO EPSILON EFF ---
plt.figure(figsize=(10, 6))

# Curva ROGERS 4350B
plt.plot(frecuencias_resonancia_rogers, epsilon_eff_rogers, 'o', color='r', markersize=10, label='εeff ROGERS')
plt.axhline(y=epsilon_eff_promedio_rogers, color='r', linestyle='--', label=f'Promedio εeff ROGERS = {epsilon_eff_promedio_rogers:.3f}')

# Etiquetas en cada punto ROGERS
for f, eps in zip(frecuencias_resonancia_rogers, epsilon_eff_rogers):
    plt.text(f, eps - 0.20, f"{eps:.3f}\n{f:.3f} GHz", color='r', ha='center', va='bottom', fontsize=12)

# Curva FR4
plt.plot(frecuencias_resonancia_fr4, epsilon_eff_fr4, 'o', color='b',markersize=10, label='εeff FR4')
plt.axhline(y=epsilon_eff_promedio_fr4, color='b', linestyle='--', label=f'Promedio εeff FR4 = {epsilon_eff_promedio_fr4:.3f}')

for f, eps in zip(frecuencias_resonancia_fr4, epsilon_eff_fr4):
    plt.text(f, eps + 0.2, f"{eps:.3f}\n{f:.3f} GHz", color='b', ha='center', va='top', fontsize=12)

# Etiquetas y formato
plt.xlabel('Frecuencias de resonancia (GHz)')
plt.ylabel('εeff(f)')
plt.title('Cálculo de Permitividad Efectiva (εeff) a 2GHZ: ROGERS 4350B vs FR4')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.ylim(2, 4) # Rango visual ajustado
plt.legend(title='Materiales')
plt.tight_layout()
plt.show()

# --- GRÁFICO COMPARATIVO EPSILON R ---
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
plt.title('Cálculo de Permitividad Relativa (εr) a 2GHZ: ROGERS 4350B vs FR4')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.ylim(3, 5) # Rango visual ajustado
plt.legend(title='Materiales')
plt.tight_layout()
plt.show()