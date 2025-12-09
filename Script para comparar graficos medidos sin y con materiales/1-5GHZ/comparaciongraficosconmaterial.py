import skrf as rf
import matplotlib.pyplot as plt

# --- CARGAR ARCHIVOS S2P ---
try:
    ntwk_rogers_medido = rf.Network('Rogers_1-5GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'Rogers_1GHz_puntos.s2p'.")
    exit()

try:
    ntwk_fr4_medido = rf.Network('FR4_1-5GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'FR4_1GHz_puntos.s2p'.")
    exit()

try:
    ntwk_rogers_simulado = rf.Network('1-5Ghz SodaLime Rogers Taco Madera_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'S21-Rogers-1G.s2p'.")
    exit()

try:
    ntwk_fr4_simulado = rf.Network('1-5Ghz SodaLime FR4 Taco Madera_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'S21-FR4-1G.s2p'.")
    exit()

print("Archivos cargados exitosamente.")

# --- Graficar ---
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# --- ROGERS ---
ntwk_rogers_medido.plot_s_db(m=1, n=0, ax=ax[0], label='Medido', color='tab:blue')
ntwk_rogers_simulado.plot_s_db(m=1, n=0, ax=ax[0], label='Con SodaLime encima', color='tab:orange', linestyle='--')
ax[0].set_xlim(0, 5e9)
ax[0].set_title('ROGERS 4350B $|S_{21}|$')
ax[0].set_xlabel('Frecuencia (Hz)')
ax[0].set_ylabel('Magnitud (dB)')
ax[0].grid(True, which="both", ls="--", alpha=0.6)
ax[0].legend()

# --- FR4 ---
ntwk_fr4_medido.plot_s_db(m=1, n=0, ax=ax[1], label='Medido', color='tab:blue')
ntwk_fr4_simulado.plot_s_db(m=1, n=0, ax=ax[1], label='Con SodaLime encima', color='tab:orange', linestyle='--')
ax[1].set_xlim(0, 5e9)
ax[1].set_title('FR4 $|S_{21}|$')
ax[1].set_xlabel('Frecuencia (Hz)')
ax[1].set_ylabel('Magnitud (dB)')
ax[1].grid(True, which="both", ls="--", alpha=0.6)
ax[1].legend()

# --- Título general ---
plt.suptitle('Comparación de Parámetros $|S_{21}|$ a 1GHZ (Transmisión)')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()