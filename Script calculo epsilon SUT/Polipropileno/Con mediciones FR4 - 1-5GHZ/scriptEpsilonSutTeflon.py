import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ====================================================================
# I. CONSTANTES Y PARÁMETROS (AJUSTAR AQUÍ)
# ====================================================================

EPSILON_0 = 8.854e-12
INF_SIMULADO = 100.0

E_1_MRR = 4.3
H_MRR = 0.00156
W_LINE = 0.00285

S_SUT = 0.0022
E_2_REF = 2.28

F_0_MEASURED = 1.531859e9 # Resonancia SIN SUT (descargada) en Hz
F_1_MEASURED = 1.47642e9 # Resonancia CON SUT (cargada) en Hz
#F_1_MEASURED += 0.01e9

# ====================================================================
# II. IMPLEMENTACIÓN DEL MÉTODO VARIACIONAL
# ====================================================================

def calculate_h_beta(beta, t=0):
    return 1.0

def calculate_g_beta(beta, E1, E2, E3, H, S, D):
    beta_abs = abs(beta)

    h_coth = 1.0 / np.tanh(beta_abs * H) if H < INF_SIMULADO else 1.0
    s_coth = 1.0 / np.tanh(beta_abs * S) if S < INF_SIMULADO else 1.0
    d_coth = 1.0 / np.tanh(beta_abs * D) if D < INF_SIMULADO else 1.0

    Num = E3 * d_coth + E2 * s_coth
    Den = beta_abs * (E3 * d_coth * (E1 * h_coth + E2 * s_coth) + E2 * (E2 + E1 * h_coth * s_coth))

    return Num / Den

def calculate_f_beta_sq(beta, w):
    u = beta * w / 2.0
    term1 = (8/5) * (np.sin(u) / u) if u != 0 else 8/5

    if u != 0:
        term2_u = u
    else:
        term2_u = 1.0e-12

    term2 = (12.0 / 5.0) / (term2_u**2) * (np.cos(u) - (2.0 * np.sin(u) / term2_u) + (np.sin(u/2) / (u/2))**2)

    return (term1 + term2)**2

def capacitance_integrand(beta, E1, E2, E3, H, S, D, w):
    f_sq = calculate_f_beta_sq(beta, w)
    g_val = calculate_g_beta(beta, E1, E2, E3, H, S, D)
    h_val = calculate_h_beta(beta)

    return f_sq * g_val * h_val

def calculate_capacitance(E1, E2, E3, H, S, D, w):
    """Calcula la Capacitancia C usando la integración numérica (Ec. 2), ajustada para estabilidad."""

    # -------------------------------------------------------------------
    # AJUSTES DE ESTABILIDAD NUMÉRICA:
    # 1. INTEGRATION_LIMIT: Límite superior finito en lugar de np.inf.
    # 2. MAX_ITERATIONS: Aumenta las subdivisiones de integración para mejor convergencia.
    # Si sigue fallando, prueba con 20000 y 200.
    # -------------------------------------------------------------------
    INTEGRATION_LIMIT = 20000
    MAX_ITERATIONS = 200

    # La integral es simétrica, integramos de 0 a INTEGRATION_LIMIT y multiplicamos por 2
    integral_result, error = quad(
        capacitance_integrand, 0, INTEGRATION_LIMIT,
        args=(E1, E2, E3, H, S, D, w),
        limit=MAX_ITERATIONS
    )

    # Ec. 2: 1/C = (1 / (pi * Q^2 * E0)) * Integral
    C = (np.pi * EPSILON_0) / (2.0 * integral_result)

    # Manejo de error: Si la integral es demasiado pequeña o grande, el resultado es dudoso.
    if error > 1e-4:
        print(f"ADVERTENCIA: La integral tuvo un alto error estimado: {error:.2e}")

    return C

# ====================================================================
# III. ALGORITMO (FIGURA 6)
# ====================================================================

def step1_2_calculate_c0_c1():
    C0 = calculate_capacitance(
        E1=1, E2=1, E3=1,
        H=H_MRR, S=INF_SIMULADO, D=INF_SIMULADO,
        w=W_LINE
    )

    C1 = calculate_capacitance(
        E1=E_1_MRR, E2=1, E3=1,
        H=H_MRR, S=INF_SIMULADO, D=INF_SIMULADO,
        w=W_LINE
    )

    return C0, C1

def step3_4_calculate_eff_perms(C0, C1, F0, F1):
    E_f0 = C1 / C0
    E_f1 = E_f0 * (F0 / F1)**2
    return E_f0, E_f1

def step5_map_e_f_to_e_2(C0, E_f1):
    """
    Mapea la permitividad efectiva medida (E_f1) a la permitividad del SUT (E_2).

    Ajuste: Se reduce el rango de mapeo (np.linspace) para forzar mayor precisión
    cerca del valor de referencia E_2_REF=2.3.
    """

    # Rango original: np.linspace(1.0, 10.0, 100)
    # Nuevo Rango Centrado en E_2_REF=2.3
    #E_2_REF = 2.3  # Usamos el valor de referencia de tu prueba de Teflon
    e2_vector = np.linspace(E_2_REF - 1.0, E_2_REF + 3.0, 500) # Rango más fino y centrado

    ef_vector = []

    for E2_test in e2_vector:
        C_sut = calculate_capacitance(
            E1=E_1_MRR, E2=E2_test, E3=1.0,
            H=H_MRR, S=S_SUT, D=INF_SIMULADO,
            w=W_LINE
        )
        E_f_test = C_sut / C0
        ef_vector.append(E_f_test)

    ef_vector = np.array(ef_vector)

    # Buscar E_2_calculated por interpolación inversa
    # Usamos interp para encontrar qué valor de e2_vector corresponde a E_f1
    # Aseguramos que los vectores estén ordenados para la interpolación
    idx = np.argsort(ef_vector)
    E_2_calculated = np.interp(E_f1, ef_vector[idx], e2_vector[idx])

    return E_2_calculated, e2_vector, ef_vector

# ====================================================================
# IV. FUNCIÓN PRINCIPAL Y COMPARACIÓN
# ====================================================================

def run_analysis():
    print("--- ANÁLISIS DE PERMITIVIDAD CON MRR ---")
    print(f"Parámetros MRR: E1={E_1_MRR}, H={H_MRR*1e3:.2f} mm, W={W_LINE*1e3:.2f} mm")
    print(f"Parámetros SUT: S={S_SUT*1e3:.2f} mm, E2_Ref={E_2_REF}")
    print(f"Frecuencias: F0={F_0_MEASURED/1e6:.2f} MHz, F1={F_1_MEASURED/1e6:.2f} MHz")
    print("-" * 55)

    C0, C1 = step1_2_calculate_c0_c1()

    E_f0, E_f1 = step3_4_calculate_eff_perms(C0, C1, F_0_MEASURED, F_1_MEASURED)

    print(f"Permitividad Efectiva Descargada (E_f0): {E_f0:.4f}")
    print(f"Permitividad Efectiva Cargada (E_f1): {E_f1:.4f}")
    print("-" * 55)

    E_2_calculated, e2_vector, ef_vector = step5_map_e_f_to_e_2(C0, E_f1)

    error_abs = abs(E_2_calculated - E_2_REF)
    error_percent = (error_abs / E_2_REF) * 100

    print(f"Valor E2 de Referencia (Teórico): {E_2_REF:.4f}")
    print(f"Valor E2 Calculado (por el script): {E_2_calculated:.4f}")
    print(f"Desviación Absoluta: {error_abs:.4f}")
    print(f"Error Porcentual: {error_percent:.2f}%")

    plt.figure(figsize=(8, 5))
    plt.plot(ef_vector, e2_vector, 'b-', label='Relación Funcional E_f vs E_2')

    plt.plot(E_f1, E_2_calculated, 'ro', markersize=8, label=f'Medición (E2={E_2_calculated:.2f})')

    plt.xlabel('Permitividad Efectiva E_f')
    plt.ylabel('Permitividad Relativa de la Muestra E_2')
    plt.title('Mapeo de Permitividad Efectiva a Permitividad del SUT')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    run_analysis()