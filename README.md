# 1 - Requisitos para correr los codigos

* Tener instalado Python 3
* Instalar las librerias scikit-rf y lmfit:

```bash
  pip install scikit-rf
  pip install lmfit
```

El formato de los datos de los archivos .s2p que utiliza la libreria scikit-rf es con punto para los decimales:
* El CST ya exporta los archivos .s2p en este formato con lo cual se pueden importar 
* El LiteVNA por ejemplo los genera con coma, por lo cual se deben buscar y reemplazar todas las comas por puntos con un programa como el bloc de notas

# 2 - Scripts presentes en este repositorio
En todos los casos se dejan codigos de ejemplo

En los que importan archivos .s2p usando la libreria scikit-rf se encuentran los archivos .s2p que utilizamos nosotros durante el proyecto, de mediciones del VNA o simulaciones de CST. En estos casos, se debe colocar en la misma carpeta el archivo .s2p de interes y colocar el nombre de mismo en la linea de codigo que importa dicho archivo

## 2.1 - Script para comparar graficos medidos y simulados
Hay un script para cada resonador que se puede encontrar en la carpeta correspondiente.

Carga 4 archivos .s2p de los parametros S21 en Rogers y FR4 tanto simulados como medidos, los grafica y los compara.

Los nombres de los archivos .s2p de interes se deben reemplazar en las lineas que utilizan la funcion rf.Network

Por ejemplo, esta linea carga los parametros .s2p exportados desde el VNA para el resonador de 1GHZ fabricado en Rogers 4350B:

```
ntwk_rogers_medido = rf.Network('Rogers_1GHz_puntos.s2p')
```
Al correr el script genera un grafico como el siguiente:

<img width="1400" height="600" alt="Figure_2" src="https://github.com/user-attachments/assets/d5d44bbf-816d-4211-a2b3-dc97ba4e0634" />


## 2.2 - Script para comparar graficos medidos sin y con materiales

Hay un script para los resonadores de 1GHZ y 1.5GHZ que permite comparar los parametros S21 del resonador solo y al colocarle un material dielectrico encima, tanto para FR4 como para Rogers, para ver el corrimiento que se produce en los picos de resonancia

```
try:
    ntwk_rogers_medido = rf.Network('Rogers_1GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'Rogers_1GHz_puntos.s2p'.")
    exit()

try:
    ntwk_fr4_medido = rf.Network('FR4_1GHz_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'FR4_1GHz_puntos.s2p'.")
    exit()

try:
    ntwk_rogers_simulado = rf.Network('1Ghz SodaLime Rogers Taco Madera_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'S21-Rogers-1G.s2p'.")
    exit()

try:
    ntwk_fr4_simulado = rf.Network('1Ghz SodaLime FR4 Taco Madera_puntos.s2p')
except FileNotFoundError:
    print("ERROR: No se encontró 'S21-FR4-1G.s2p'.")
    exit()
```
El nombre de las variables quedo igual que en el del script 2.1, pero se deben reemplazar el nombre de los archivos que se cargan en  ntwk_rogers_simulado y en ntwk_fr4_simulado por los archivos .s2p de medicion con el material dielectrico que se haya colocado encima de los resonadores. Los otros 2 son los .s2p de los resonadores sin nada encima
El grafico generado es como el siguiente:

<img width="1400" height="600" alt="Figure_3" src="https://github.com/user-attachments/assets/377c0b90-7851-4001-bfde-1efe77a1388d" />


## 2.3 - Script Fitteo con y sin material
Este es el script que realiza el curve fitting a los parametros S21, genera varios graficos con los resultados de cada pico de resonancia antes y despues del fitteo y al final de todo imprime los resultados de las frecuencias de resonancia luego del fitting mediante 3 metodos, siendo el ultimo, el Lorentziano, el mas exacto:

```
🔸 PROCESANDO RESONANCIA #1 (1.3219 GHz)
   ✅ AJUSTE EXITOSO:
      ├─ f₀: 1.321506 ± 0.000269 GHz
      ├─ Q: 27.5 ± 0.6
      ├─ R²: 0.9697
      ├─ χ²: 3.20e-02
      └─ Modo: exacto

🔸 PROCESANDO RESONANCIA #2 (2.6438 GHz)
   ✅ AJUSTE EXITOSO:
      ├─ f₀: 2.643616 ± 0.000547 GHz
      ├─ Q: 35.5 ± 1.0
      ├─ R²: 0.9497
      ├─ χ²: 1.91e-01
      └─ Modo: exacto

🔸 PROCESANDO RESONANCIA #3 (3.9516 GHz)
   ✅ AJUSTE EXITOSO:
      ├─ f₀: 3.951010 ± 0.000774 GHz
      ├─ Q: 41.7 ± 1.4
      ├─ R²: 0.9397
      ├─ χ²: 5.72e-01
      └─ Modo: exacto
 ```

Hay un script que esta hecho para los parametros .s2p generados por el CST y otro para los exportados desde el LiteVNA debido a esta diferencia de formatos, cada uno se encuentra en la carpeta correspondiente

El archivo .s2p de interes a fittear se debe cargar en la linea 46:

```
ARCHIVO_S2P = '1-5Ghz SodaLime Rogers Taco Madera.s
```
## 2.4 - Script para calculo epsilon relativo del sustrato
Hay 2 carpetas, con un script para calcular el epsilon relativo del FR4 y Rogers, tanto para parametros s21 obtenidos de simulaciones del CST como para las mediciones post fitting

Dentro de ellas hay una carpeta con un script para cada resonador (1 1.5 y 2GHZ)

Ambos ya tienen cargados los parametros geometricos de los resonadores

En el script para simulaciones, se debe cargar el archivo .s2p correspondiente en las lineas 31 y 32:

```
ntwk_rogers = rf.Network('S21-Rogers-1G.s2p')
ntwk_fr4 = rf.Network('S21-FR4-1G.s2p')
```

En el script para mediciones en cambio, se deben colocar las frecuencias de resonancia obtenidas directamente del script de fitting, y hardcodearlas en las lineas 30 y 38:

```
frecuencias_resonancia_rogers = [0.980079, 1.958852, 2.925703]
frecuencias_resonancia_fr4 = [1.036478, 2.071893, 3.097262]
```

Al correr el script, se generara un grafico con los resultados del calculo del epsilon efectivo y luego del relativo (el de interes) para cada modo resonante:

<img width="1000" height="600" alt="Figure_4" src="https://github.com/user-attachments/assets/931a9b68-e8ad-4bdc-809d-3408758a668d" />


## 2.5 - Script calculo epsilon del material SUT
Hay 4 scripts, para los resonadores de 1GHZ y 1.5GHZ de FR4 y Rogers
En el ejemplo esta hecho para el teflon

En las lineas 17, 18, 20 y 21:

```
S_SUT = 0.00152
E_2_REF = 2.1

F_0_MEASURED = 1.531859e9   # Resonancia SIN SUT (descargada) en Hz
F_1_MEASURED = 1.412021e9 # Resonancia CON SUT (cargada) en Hz
```

Se deben cargar en el orden en que estan:

* Espesor del material dielectrico a medir (SUT)
* Epsilon relativo aproximado del material (obtenido de tablas)
* Frecuencia de resonancia del resonador sin SUT, la obtenida luego de hacer el fitting
* Frecuencia de resonancia del resonador con el SUT, la obtenida luego de hacer el fitting

Al correr el programa se imprime lo siguiente:

```
--- ANÁLISIS DE PERMITIVIDAD CON MRR ---
Parámetros MRR: E1=4.3, H=1.56 mm, W=2.60 mm
Parámetros SUT: S=21.52 mm, E2_Ref=2.1
Frecuencias: F0=1036.48 MHz, F1=986.26 MHz
-------------------------------------------------------
Permitividad Efectiva Descargada (E_f0): 3.2206
Permitividad Efectiva Cargada (E_f1): 3.5569
-------------------------------------------------------
Valor E2 de Referencia (Teórico): 2.1000
Valor E2 Calculado (por el script): 1.9810
Desviación Absoluta: 0.1190
Error Porcentual: 5.67%
```
Y se abre una ventana con un la recta generada por el metodo variacional que utiliza esta tecnica de calculo, indicando con un punto rojo el valor final obtenido del epsilon relativo del material a medir 

<img width="800" height="500" alt="Figure_1" src="https://github.com/user-attachments/assets/556f595e-1547-4f01-8bb6-0c8a504c5795" />
