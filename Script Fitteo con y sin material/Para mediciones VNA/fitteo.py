# ============================================================
# 📦 IMPORTACIÓN DE LIBRERÍAS ESENCIALES
# ============================================================
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import signal              # usado para find_peaks
from scipy.optimize import minimize   # usado en algunos ajustes
import lmfit                          # núcleo de los ajustes no lineales
from lmfit import Model               # para definir modelos de ajuste

print("✅ Librerías importadas correctamente")
print(f"📦 lmfit version: {lmfit.__version__}")
print(f"📊 numpy version: {np.__version__}")
print(f"🐼 pandas version: {pd.__version__}")

# ============================================================
# 🎨 CONFIGURACIÓN DE ESTILO, CONSTANTES Y COLORES
# ============================================================
warnings.filterwarnings('ignore', category=RuntimeWarning)

plt.style.use('default')
plt.rcParams.update({
    'figure.figsize': (12, 8),
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.formatter.useoffset': False,
})
np.set_printoptions(precision=4, suppress=True)

COLORES = {
    'resonancia':   '#E74C3C',  # rojo para marcar f0
    'ajuste':       '#8FBC8F',  # verde para curvas de modelo/fit
    'datos':        '#3498DB',  # azul para datos
    'fondo':        '#95A5A6',  # gris
    'destacado':    '#F39C12',  # naranja para resaltar
    'comparacion':  '#9B59B6',  # morado para comparativas
    'incertidumbre':'#E67E22',  # naranja oscuro para barras de error
}

ARCHIVO_S2P = '1-5Ghz SodaLime Rogers Taco Madera.s2p'  # archivo Touchstone de trabajo

print("🎨 Estilo aplicado y paleta definida")
print("📁 Archivo objetivo:", ARCHIVO_S2P)

# ============================================================
# 📂 FUNCIÓN DE CARGA DE DATOS S21 - ROGERS 1GHz
# ============================================================

def cargar_datos_rogers_s2p():
    """
    Carga datos S21 del archivo Rogers-1GHz.s2p en formato Touchstone (RI).
    Maneja correctamente comas decimales y cabeceras de NanoVNA.
    """
    if not os.path.exists(ARCHIVO_S2P):
        print(f"❌ Archivo no encontrado: {ARCHIVO_S2P}")
        return None

    try:
        data = []
        with open(ARCHIVO_S2P, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                # Saltar comentarios y encabezados
                if not line or line.startswith(('!', '#')):
                    continue
                # Cambiar coma por punto decimal
                line = line.replace(',', '.')
                values = [float(x) for x in line.split()]
                if len(values) == 9:
                    data.append(values)

        if not data:
            print("⚠️ No se encontraron datos válidos en el archivo.")
            return None

        data = np.array(data)

        # Frecuencia en GHz
        f_Hz = data[:, 0]
        f_GHz = f_Hz / 1e9

        # S21 (columnas 3 y 4 -> índice 3 y 4)
        s21_real = data[:, 3]
        s21_imag = data[:, 4]
        S21 = s21_real + 1j * s21_imag

        # Magnitud y fase
        S21_mag = np.abs(S21)
        S21_db = 20 * np.log10(S21_mag + 1e-15)
        S21_ang = np.angle(S21, deg=True)

        df = pd.DataFrame({
            "f_GHz": f_GHz,
            "S21": S21,
            "S21_mag": S21_mag,
            "S21_db": S21_db,
            "S21_ang": S21_ang
        })

        print(f"📁 Formato detectado: RI (Real–Imag). {len(df)} puntos cargados correctamente.")
        return df

    except Exception as e:
        print(f"⚠️ Error al procesar {ARCHIVO_S2P}: {e}")
        return None

# ============================================================
# 📂 CARGA Y VISUALIZACIÓN INICIAL DE DATOS S21 - ROGERS 1GHz
# ============================================================

df = cargar_datos_rogers_s2p()

if df is not None:
    # Variables maestras
    f = df["f_GHz"].values
    S21 = df["S21"].values
    S21_db = df["S21_db"].values
    S21_mag = df["S21_mag"].values

    # Métricas básicas
    n_puntos = len(f)
    f_min, f_max = f[0], f[-1]
    s21_min, s21_max = S21_db.min(), S21_db.max()
    s21_mean = np.mean(S21_db)
    datos_validos = not (np.any(np.isnan(S21_db)) or np.any(np.isinf(S21_db)))

    # Resumen en consola
    print(f"✅ {ARCHIVO_S2P} procesado correctamente")
    print(f"   📊 {n_puntos:,} puntos en el rango {f_min:.3f} – {f_max:.3f} GHz")
    print(f"   📈 |S21| en dB: min = {s21_min:.1f}, max = {s21_max:.1f}, promedio = {s21_mean:.1f}")
    print(f"   🔍 Integridad de datos: {'✓ OK' if datos_validos else '⚠️ Problemas detectados'}")

    # Visualización
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Magnitud
    ax1.plot(f, S21_db, color=COLORES["datos"], lw=1.2)
    ax1.set_ylabel("|S21| [dB]")
    ax1.set_title("Magnitud de S21 en función de la frecuencia")
    ax1.grid(True, alpha=0.3)

    # Fase
    fase_deg = np.angle(S21, deg=True)
    ax2.plot(f, fase_deg, color=COLORES["comparacion"], lw=1.2)
    ax2.set_ylabel("Fase [°]")
    ax2.set_xlabel("Frecuencia [GHz]")
    ax2.set_title("Fase de S21 en función de la frecuencia")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
else:
    print("❌ No se pudieron cargar los datos.")

# ============================================================
# 📐 MODELOS LORENTZIANOS
# ============================================================

def modelo_lorentziano_lineal(f, f0, Q, A_real, A_imag, B_real, B_imag):
    A = A_real + 1j * A_imag
    B = B_real + 1j * B_imag
    return A + B / (1 + 2j * Q * (f - f0) / f0)

def modelo_lorentziano_exacto(f, f0, Q, A_real, A_imag, B_real, B_imag):
    A = A_real + 1j * A_imag
    B = B_real + 1j * B_imag
    return A + B / (1 + 2j * Q * ((f/f0) - (f0/f)))

def ajustar_lorentziano(f, S21, f0_inicial, Q_inicial=50,
                        ventana=0.15, modo="lineal"):
    """
    Ajusta una resonancia con modelo Lorentziano (lineal o exacto).
    """
    # --- Selección de ventana ---
    mask = (f >= f0_inicial - ventana/2) & (f <= f0_inicial + ventana/2)
    f_win = f[mask]
    S21_win = S21[mask]

    if len(f_win) < 10:
        return {"exito": False, "error": "Pocos puntos en la ventana"}

    # --- Modelo según modo ---
    if modo == "lineal":
        modelo_func = modelo_lorentziano_lineal
    elif modo == "exacto":
        modelo_func = modelo_lorentziano_exacto
    else:
        raise ValueError("modo debe ser 'lineal' o 'exacto'")

    modelo = Model(modelo_func)

    # Parámetros iniciales
    params = modelo.make_params(
        f0=f0_inicial,
        Q=Q_inicial,
        A_real=0.0, A_imag=0.0,
        B_real=np.abs(S21_win).max(),
        B_imag=0.0
    )

    # Límites
    params["f0"].set(min=f0_inicial*0.98, max=f0_inicial*1.02)
    params["Q"].set(min=5, max=5000)

    # Residuales (real + imag)
    datos = np.concatenate([S21_win.real, S21_win.imag])
    def residuales(pars, freqs, datos):
        S21_model = modelo_func(freqs, **pars.valuesdict())
        return datos - np.concatenate([S21_model.real, S21_model.imag])

    # Ajuste
    resultado = lmfit.minimize(residuales, params, args=(f_win, datos))
    if not resultado.success:
        return {"exito": False, "error": "Ajuste no convergió"}

    # Modelo ajustado
    S21_fit = modelo_func(f_win, **resultado.params.valuesdict())

    # Métricas
    ss_res = np.sum((datos - np.concatenate([S21_fit.real, S21_fit.imag]))**2)
    ss_tot = np.sum((datos - np.mean(datos))**2)
    R2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0

    return {
        "exito": True,
        "f0": resultado.params["f0"].value,
        "Q": resultado.params["Q"].value,
        "params": resultado.params,
        "frecuencias": f_win,
        "S21_datos": S21_win,
        "S21_fit": S21_fit,
        "R2": R2,
        "chi2": resultado.chisqr,
        "chi2_red": resultado.redchi,
        "modo": modo
    }

def analizar_lorentziano(df_simple, f, S21, ventana=0.15, modo="lineal"):
    """
    Analiza todas las resonancias con modelo Lorentziano (lineal o exacto).
    """
    resultados = []
    print("="*70)
    print(f"🎯 ANÁLISIS LORENTZIANO ({modo.upper()})")
    print("="*70)

    for _, row in df_simple.iterrows():
        numero = int(row["numero"])
        f0_ini = row["f_res_GHz"]
        Q_ini = row["Q_est"]

        print(f"\n🔸 PROCESANDO RESONANCIA #{numero} ({f0_ini:.4f} GHz)")

        ajuste = ajustar_lorentziano(
            f, S21,
            f0_inicial=f0_ini,
            Q_inicial=Q_ini if not np.isnan(Q_ini) else 50,
            ventana=ventana,
            modo=modo
        )

        if ajuste["exito"]:
            resultados.append(ajuste)
            p = ajuste["params"]
            print("   ✅ AJUSTE EXITOSO:")
            print(f"      ├─ f₀: {ajuste['f0']:.6f} ± {p['f0'].stderr or 0:.6f} GHz")
            print(f"      ├─ Q: {ajuste['Q']:.1f} ± {p['Q'].stderr or 0:.1f}")
            print(f"      ├─ R²: {ajuste['R2']:.4f}")
            print(f"      ├─ χ²: {ajuste['chi2']:.2e}")
            print(f"      └─ Modo: {modo}")
        else:
            print(f"   ❌ AJUSTE FALLIDO: {ajuste['error']}")
    print("\n")
    return resultados

# ============================================================
# 📊 MÉTODO SIMPLE (f_res y Q por -3 dB)
# ============================================================

def metodo_simple(f, S21_db, altura_minima=-40, separacion_minima=30, prominencia_min=3, plot=True):
    """
    Método simple para determinar frecuencia de resonancia y Q por -3 dB.

    Parámetros
    ----------
    f : array
        Frecuencias en GHz
    S21_db : array
        Magnitud de S21 en dB
    altura_minima : float
        Nivel mínimo para considerar un pico (dB)
    separacion_minima : int
        Separación mínima entre resonancias (en puntos)
    prominencia_min : float
        Prominencia mínima (criterio interno de find_peaks)
    plot : bool
        Si True, genera gráfico con las resonancias detectadas

    Returns
    -------
    df_res : DataFrame con resonancias y Q:
        - numero
        - f_res_GHz
        - S21_res_dB
        - f_izq_3dB
        - f_der_3dB
        - BW_3dB
        - Q_est
    """
    # --- detectar resonancias ---
    idx, _ = signal.find_peaks(
        S21_db,
        height=altura_minima,
        prominence=prominencia_min,
        distance=separacion_minima
    )

    resultados = []
    for i, k in enumerate(idx):
        f_res = f[k]
        nivel_max = S21_db[k]
        nivel_3dB = nivel_max - 3

        # buscar cruces -3 dB
        i_left = np.where(S21_db[:k] <= nivel_3dB)[0]
        i_right = np.where(S21_db[k:] <= nivel_3dB)[0]

        f_izq = f[i_left[-1]] if len(i_left) > 0 else np.nan
        f_der = f[k + i_right[0]] if len(i_right) > 0 else np.nan

        BW = f_der - f_izq if (np.isfinite(f_izq) and np.isfinite(f_der)) else np.nan
        Q_est = f_res / BW if BW and BW > 0 else np.nan

        resultados.append([i+1, f_res, nivel_max, f_izq, f_der, BW, Q_est])

    df_res = pd.DataFrame(resultados, columns=[
        "numero", "f_res_GHz", "S21_res_dB", "f_izq_3dB", "f_der_3dB", "BW_3dB", "Q_est"
    ])

    # --- salida en consola ---
    print("📊 MÉTODO SIMPLE (−3 dB): RESULTADOS")
    print("=" * 60)
    if len(df_res) > 0:
        for _, r in df_res.iterrows():
            print(f"🔸 Resonancia #{int(r['numero'])}")
            print(f"   ├─ Frecuencia: {r['f_res_GHz']:.4f} GHz")
            print(f"   ├─ Magnitud: {r['S21_res_dB']:.2f} dB")
            print(f"   ├─ f(-3dB) izq: {r['f_izq_3dB']:.4f} GHz")
            print(f"   ├─ f(-3dB) der: {r['f_der_3dB']:.4f} GHz")
            print(f"   ├─ Δf (-3dB): {r['BW_3dB']:.6f} GHz")
            print(f"   └─ ✨ Q estimado: {r['Q_est']:.1f}\n")
    else:
        print("⚠️ No se detectaron resonancias con los criterios actuales")

    print(f"✅ Análisis completado para {len(df_res)} resonancias")

    # --- gráfico opcional ---
    if plot and len(df_res) > 0:
        plt.figure(figsize=(12, 6))
        plt.plot(f, S21_db, color=COLORES["datos"], lw=1.2, label="|S21| datos")
        plt.scatter(df_res["f_res_GHz"], df_res["S21_res_dB"],
                    color=COLORES["resonancia"], s=80, marker="o", label="Resonancias")
        for _, row in df_res.iterrows():
            plt.text(row["f_res_GHz"], row["S21_res_dB"]+1,
                     f"{row['f_res_GHz']:.3f} GHz", ha="center", va="bottom",
                     fontsize=9, color=COLORES["resonancia"])
        plt.xlabel("Frecuencia [GHz]")
        plt.ylabel("|S21| [dB]")
        plt.title("Método Simple: detección de f_res y Q")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.show()

    return df_res

# ============================================================
# EJECUCIÓN DEL MÉTODO SIMPLE Y GUARDADO DE RESULTADOS
# ============================================================
df_simple = metodo_simple(
    f, S21_db,
    altura_minima=-40,
    separacion_minima=30,
    prominencia_min=3,
    plot=True
)

# ============================================================
# Ajuste lineal
# ============================================================
resultados_lineal = analizar_lorentziano(
    df_simple,
    f,
    S21,
    ventana=0.8,
    modo="lineal"
)
# ============================================================
# Ajuste exacto
resultados_exacto = analizar_lorentziano(
    df_simple,
    f,
    S21,
    ventana=0.2,
    modo="exacto"
)
# ============================================================

# ============================================================
# 📊 VISUALIZADOR GENÉRICO - LORENTZIANO (LINEAL O EXACTO)
# ============================================================

def visualizar_lorentziano(resultados):
    """
    Visualiza los resultados de un ajuste Lorentziano (lineal o exacto).
    Muestra magnitud (dB) y fase (grados) para cada resonancia.

    Parámetros
    ----------
    resultados : list of dict
        Salida de analizar_lorentziano (contiene 'modo')
    """
    if len(resultados) == 0:
        print("⚠️ No hay resultados para visualizar")
        return

    modo = resultados[0].get("modo", "desconocido").upper()
    n_res = len(resultados)
    n_cols = 2  # magnitud y fase
    n_rows = n_res

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    fig.suptitle(f"Método Lorentziano ({modo})", fontsize=16, fontweight="bold")

    if n_res == 1:
        axes = np.array([axes])  # caso especial

    for i, res in enumerate(resultados):
        f_win = res["frecuencias"]
        S21_datos = res["S21_datos"]
        S21_fit = res["S21_fit"]

        # Magnitud en dB
        mag_datos = 20*np.log10(np.abs(S21_datos))
        mag_fit = 20*np.log10(np.abs(S21_fit))

        # Fase en grados (unwrap para evitar saltos de 360°)
        # fase_datos = np.unwrap(np.angle(S21_datos)) * 180/np.pi
        # fase_fit = np.unwrap(np.angle(S21_fit)) * 180/np.pi

        fase_datos = np.unwrap(np.angle(S21_datos)) * 180/np.pi
        fase_fit = np.unwrap(np.angle(S21_fit)) * 180/np.pi

        ax_mag = axes[i, 0] if n_res > 1 else axes[0]
        ax_fase = axes[i, 1] if n_res > 1 else axes[1]

        # === Magnitud ===
        ax_mag.plot(f_win, mag_datos, 'o', color=COLORES["datos"],
                    markersize=4, alpha=0.7, label="Datos")
        ax_mag.plot(f_win, mag_fit, '-', color=COLORES["ajuste"],
                    linewidth=2, label=f"Ajuste {modo}")
        ax_mag.axvline(x=res["f0"], color="gray", ls="--", lw=1)
        ax_mag.set_ylabel("|S21| [dB]", fontweight="bold")
        ax_mag.set_title(f"Resonancia {i+1} - Q={res['Q']:.1f}",
                         fontweight="bold", color=COLORES["ajuste"])
        ax_mag.grid(True, alpha=0.3)
        ax_mag.legend(fontsize=9)

        # === Fase ===
        ax_fase.plot(f_win, fase_datos, 'o', color=COLORES["datos"],
                     markersize=4, alpha=0.7, label="Datos")
        ax_fase.plot(f_win, fase_fit, '-', color=COLORES["ajuste"],
                     linewidth=2, label=f"Ajuste {modo}")
        ax_fase.axvline(x=res["f0"], color="gray", ls="--", lw=1)
        ax_fase.set_ylabel("Fase [°]", fontweight="bold")
        ax_fase.set_xlabel("Frecuencia [GHz]", fontweight="bold")
        ax_fase.grid(True, alpha=0.3)
        ax_fase.legend(fontsize=9)

    plt.tight_layout()
    plt.show()

    print(f"📊 Visualización Lorentziano ({modo}) completada para {n_res} resonancias")

# Para modelo lineal
#visualizar_lorentziano(resultados_lineal)

# Para modelo exacto
visualizar_lorentziano(resultados_exacto)