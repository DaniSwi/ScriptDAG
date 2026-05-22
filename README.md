# Motor de Pre-procesamiento para NCSP: Extracción de Subsistemas Lineales

Este proyecto es una herramienta de pre-procesamiento matemático modular diseñada para analizar Problemas de Satisfacción de Restricciones Numéricas (NCSP). Toma instancias de ecuaciones complejas generadas en texto plano, construye un Grafo Acíclico Dirigido (DAG) para fusionar subexpresiones comunes y extrae los subsistemas de sumas n-arias en un formato de álgebra lineal estandarizado (`.npz`).

El sistema separa el componente lineal ($A \cdot y = b$) de la estructura no lineal ("cascarón") y cuenta con un auditor matemático incorporado que demuestra la equivalencia estricta del procesamiento.

---

## 🛠 Instalación y Configuración del Entorno (Recomendado)

Para garantizar la reproducibilidad del código y evitar conflictos con otras librerías en tu sistema, se recomienda ejecutar este proyecto dentro de un entorno virtual (`.venv`).

El proyecto está construido en **Python 3** y requiere `numpy` y `sympy`. Sigue estos pasos para configurar tu entorno:

### Paso 1: Crear el Entorno Virtual
Abre tu terminal (o la consola integrada de tu editor, como Visual Studio Code), asegúrate de estar en la carpeta raíz del proyecto y ejecuta:

```bash
python -m venv .venv
```

### Paso 2: Activar el Entorno Virtual
Debes activar el entorno para que las instalaciones se guarden aquí y no de forma global.

En Windows (PowerShell):

```PowerShell
.\.venv\Scripts\Activate.ps1
```
En Windows (Command Prompt - CMD):

```DOS
.\.venv\Scripts\activate.bat
```
En macOS y Linux:

```Bash
source .venv/bin/activate
```

### Paso 3: Instalar las Dependencias
Con el entorno activado, instala las librerías matemáticas requeridas:

```Bash
pip install numpy sympy
```

# Arquitectura del Proyecto

El código fuente está desacoplado en 5 módulos con responsabilidades únicas:

- parser_instancias.py: Interfaz de lectura. Traduce los archivos .txt ignorando comentarios y utilizando sympy para convertir las cadenas de texto en Árboles de Sintaxis Abstracta (AST).

- dag_builder.py: Núcleo de grafos. Implementa la lógica de unificación; detecta subexpresiones repetidas y las mapea al mismo nodo en memoria (solucionando el problema de localidad).

- matrix_extractor.py: Motor de extracción quirúrgica. Localiza las sumas n-arias, extrae sus coeficientes para la matriz lineal, y reemplaza la suma en la ecuación original por una variable auxiliar w_i.

- main.py: El orquestador. Conecta las capas de lectura y extracción, y empaqueta los resultados en un archivo binario .npz listo para ser consumido por solvers externos.

- validator.py: Auditor matemático. Realiza ingeniería inversa tomando el .npz resultante, reconstruye el sistema lineal y usa el motor de simplificación de SymPy para demostrar empíricamente que la diferencia entre el sistema procesado y el original es exactamente 0.

# Guía de uso

## 1. Generación de las Matrices
Coloca tu archivo de instancia (en la carpeta de instancia) y ejecuta el orquestador principal:

```Bash
python main.py
```

## 2. Auditoría y Validación (Prueba de Correctitud)
Para garantizar que el pre-procesamiento no alteró la integridad matemática del sistema original, ejecuta el validador:


```Bash
python validator.py
```

Salida esperada: Un reporte en consola confirmando, ecuación por ecuación, que el sistema reconstruido es perfectamente equivalente a la fuente original.

# Estructura del archivo .npz 
El archivo de salida está optimizado con NumPy. Si deseas cargar estos datos en tus propios modelos o algoritmos de optimización, el archivo contiene los siguientes arreglos accesibles mediante llaves:

```Python
import numpy as np
datos = np.load('subsistema_n_ario.npz')

# Acceso a los componentes:
A = datos['A']
y = datos['y']
b = datos['b']
cascarones = datos['shell_eqs']
```

- datos['A']: Matriz 2D de coeficientes numéricos (flotantes) de las sumas n-arias.
- datos['y']: Vector 1D de cadenas de texto. Actúa como leyenda, indicando qué subexpresión matemática (ej. x12) corresponde a cada columna de la Matriz A.
- datos['b']: Vector numérico 1D con las constantes de las sumas extraídas.
- datos['shell_eqs']: Vector de cadenas de texto con las ecuaciones "cascarón" originales, donde las sumas lineales extraídas han sido reemplazadas por las variables auxiliares $w_0, w_1, \dots, w_n$.
