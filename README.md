# ScriptDAG: identificación de subsistemas n-arios en NCSP

Preprocesador y contractor para problemas de satisfacción de restricciones numéricas
(NCSP) que implementa la Sección 3 de:

> Araya, I. & Reyes, V. (2019). *Enhancing interval constraint propagation by
> identifying and filtering n-ary subsystems.*

Lee instancias `.txt` en formato tipo Minibex, construye un DAG que fusiona las
subexpresiones comunes, identifica las sumas n-arias y extrae el subsistema lineal
asociado, junto con dominios sound para cada nodo (Secciones 3.1–3.2). Encima de eso,
`nsc.py` implementa el contractor **NSC** de la Sección 3.3.

**Alcance.** El proyecto implementa el NSC basado en el conditioner de Gauss–Jordan
(§3.3.2). La variante con simplex (§3.3.1) no está implementada, porque necesita un
LP solver con cotas certificadas. Tampoco hay un solver branch & bound: NSC contrae
una caja, pero no bisecta.

---

## Paso a paso

Todos los comandos se ejecutan desde la carpeta raíz del repo (la que contiene
`main.py`).

### 1. Requisitos

- **Python 3.9 o superior.** Para ver tu versión, ejecuta `python --version`. En
  Windows, si `python` no existe, prueba con `py --version` y usa `py` en lugar de
  `python` en todos los comandos.
- **git**, para clonar el repo. También puedes descargarlo como ZIP desde GitHub.

### 2. Descargar el repo

```bash
git clone https://github.com/DaniSwi/ScriptDAG.git
```

```bash
cd ScriptDAG
```

### 3. Crear el entorno virtual e instalar numpy

Este paso se hace una sola vez. Primero crea el entorno:

```bash
python -m venv .venv
```

Después actívalo; el comando depende de tu terminal.

- **Windows, PowerShell:**

  ```powershell
  .\.venv\Scripts\Activate.ps1
  ```

  Si aparece un error de *"ejecución de scripts deshabilitada"*, usa CMD o la opción
  sin activar que se explica más abajo.

- **Windows, CMD:**

  ```bat
  .venv\Scripts\activate.bat
  ```

- **macOS / Linux:**

  ```bash
  source .venv/bin/activate
  ```

Con el entorno activo, el prompt empieza con `(.venv)`. Entonces instala numpy, que es
la única dependencia:

```bash
pip install numpy
```

**Sin activar el entorno:** puedes llamar directamente al Python del entorno, por
ejemplo `.venv\Scripts\python main.py --lote` en Windows o
`.venv/bin/python main.py --lote` en macOS/Linux.

Cada vez que abras una terminal nueva, vuelve a activar el entorno; no hace falta
reinstalar.

### 4. Comprobar que todo funciona

```bash
python test_paper.py
```

```bash
python test_nsc.py
```

Los dos deben terminar en `Todo OK.` (`test_nsc.py` tarda unos segundos).

### 5. Procesar las instancias incluidas

```bash
python main.py --lote
```

Este comando procesa todos los `.txt` de `instances/` y escribe un `.npz` por instancia
en `outputs/`. Además, audita cada una. Al final debe mostrar
`RESUMEN: 3 exitosos | 0 fallidos`.

### 6. Procesar una sola instancia

Esta es la salida real con `inst001`; la línea `escrito` indica dónde quedó el
resultado:

```
python main.py instances/inst001.txt

--- Pre-procesamiento de instances/inst001.txt ---
  DAG      : {'nodos': 33, ..., 'compartidos_CSE': 19}
  sistema  : {'filas_w': 10, 'columnas_y': 12, ..., 'variables_auxiliares_nuevas': 0}
  A        : 10 filas (w) x 12 columnas (y)
  escrito  : outputs/inst001.npz
--- AUDITORIA instances/inst001.txt ---
  [PASS] V1 conteo de restricciones  archivo=8  parseadas=8
  ...
  => TODO OK
```

### 7. Comparar HC4 contra NSC (opcional)

Agrega `--nsc` a cualquiera de los comandos anteriores. Para cada instancia, imprime
cuánto contrae la caja inicial HC4 sobre el DAG (REF) y cuánto NSC:

```bash
python main.py instances/inst001.txt --nsc
```

```
  REF (HC4 en el DAG) : 4/8 variables contraidas, reduccion de volumen 10^9.42, 112 iter, 0.08s
  NSC                 : 4/8 variables contraidas, reduccion de volumen 10^9.42, 112 iter, 0.11s
```

### Todas las opciones

| Opción | Qué hace |
|---|---|
| `--lote` | Procesa toda la carpeta `instances/` (es lo que pasa si no se da ningún archivo). |
| `--nsc` | Compara además REF contra NSC. |
| `--txt` | Escribe también `outputs/<nombre>_reescrita.txt`, el sistema reescrito con variables auxiliares (estilo Ceberio). |
| `--sin-dominios` | No calcula dominios; todos quedan en `(-inf, inf)`. |
| `--sin-validar` | Omite la auditoría (solo junto con `--lote`). |

## Usar tus propias instancias

Copia tu archivo `.txt` a `instances/` y ejecuta `python main.py --lote`. También puedes
procesarlo directamente con `python main.py ruta/a/tu_instancia.txt`. El formato es el
siguiente:

```
// los comentarios empiezan con // o #
Variables
x1 in [-1, 5];
x2 in [-1e8, 1e8];
Constraints
5*x2 + x1*x2 - x2 = 3;
2*x2*sin(x1 + x2) + x1^2 <= 2;
exp(x1) - 1/x2 = 0;
end
```

- Una variable por línea: `nombre in [inferior, superior];`.
- Una restricción por línea, terminada en `;`, con `=`, `<=`, `>=`, `<` o `>`.
- Operadores `+ - * / ^`. Funciones: `sin`, `cos`, `exp`, `log`, `sqrt`, `atan`,
  `tanh`, `Abs`.
- Las líneas `minimize ...` o `maximize ...` se ignoran.

**Todavía no se soporta:** la sección `constants`, la sección `function`, arreglos de
variables (`x[3]`), desigualdades dobles (`-1 <= x^2 <= 4`) ni restricciones que
ocupen varias líneas. Si el archivo trae alguna de estas, el programa se detiene con un
`ParseError` que indica el archivo y la línea; nunca la ignora en silencio.

## El sistema extraído

Por cada nodo suma n-aria del DAG hay una fila:

$$w_i = c_i + \sum_j A_{ij}\, y_j$$

- `c` es el vector de **constantes aditivas** de cada suma.
- `w` son las **raíces** de las sumas, es decir, los nodos suma del DAG.
- `y` son las **columnas**: variables, potencias, productos, funciones y también
  otros `w_k` cuando una suma aparece dentro de otra (las columnas mixtas de la
  Fig. 3 del paper).

Las `w_i` y las `y_j` son **etiquetas de nodos del DAG, no variables nuevas**. Tal
como dice el paper (*"we do not modify the original system"*), la salida oficial es
un anexo al problema original y no lo reescribe. Productos, potencias y funciones
son columnas de `A`, y sus coeficientes numéricos se sacan hacia afuera. Por eso
`x2^2`, `-x2^2` y `5*x2^2` comparten una sola columna, y `A` conserva los
coeficientes reales en vez de quedar como una matriz de incidencia de ±1.

Los coeficientes se calculan en aritmética exacta (`Fraction`). Los dominios son
sobre-aproximaciones sound, con redondeo dirigido hacia afuera. Cuando un nodo no se
puede acotar, su dominio es `(-inf, +inf)`.

## El contractor NSC (Sección 3.3)

`nsc.contraer` repite hasta punto fijo:

1. **HC4-revise sobre el DAG**: una fase forward y una backward sobre todos los
   nodos. Los dominios de los nodos intermedios persisten entre pasadas, así que un
   nodo compartido conserva la información de proyección que reciba desde cualquiera
   de sus padres.
2. **Proyección sobre `PA·y = P·(w − c)`**: cada fila se trata como una restricción
   lineal y se proyecta sobre sus `y`. `P` es el conditioner de Gauss–Jordan con
   pivote de máxima magnitud por columna, calculado una sola vez en racionales
   exactos.

Toda la aritmética es de intervalos, con redondeo hacia afuera. Una intersección
vacía significa que la caja no tiene soluciones.

```python
from exportator import construir
from nsc import Preparado, conditioner_de, contraer
inst, sis, _ = construir('instances/inst001.txt')
dom, iteraciones = contraer(Preparado(inst, sis, conditioner_de(sis)))
```

`test_nsc.py` comprueba tres cosas:
- Con la matriz del paper, el conditioner reproduce exactamente la `P` y la `PA` de
  la §3.3.2, incluida la fila dependiente con −1/6.
- `P·A = PA` en racionales para cada instancia.
- Soundness: se toma un punto, se encierran con intervalos los valores de las raíces
  en ese punto y se usan como lados derechos, de modo que el punto es una solución
  real del sistema. Ni REF ni NSC pueden eliminarlo.

Esta prueba detecta errores plantados, como usar `a` en lugar de `1/a`, dividir de
forma ingenua por un intervalo que contiene 0 u olvidar la rama negativa de una raíz
par.

**Qué tanto contrae.** En el ejemplo (1) del paper (dominios `[-1, 5]`), REF reduce
el volumen de la caja en 10^0.15 y NSC en 10^1.36. En las tres instancias de
`instances/`, partiendo de la caja inicial, NSC y REF contraen lo mismo. La
ganancia que reporta el paper se mide durante el branch & bound, en cajas más
pequeñas, y ese solver no está en el repo.

## Arquitectura

```
parser_instancias.py   .txt -> Instancia (vars, dominios, ops, rhs, raices, dag)
  dag_builder.py       ast de Python -> DAG normalizado (hash-consing, Algoritmo 1)
matrix_extractor.py    DAG -> SistemaLineal (A exacta, c, y, w, subsistema de productos)
  intervals.py         aritmetica de intervalos con redondeo dirigido
  dominios.py          propagacion de dominios por el DAG + muestreadores
exportator.py          pipeline completo; exportar_npz / exportar_txt
validator.py           auditoria (V1-V5)
nsc.py                 contractor NSC (Sec. 3.3): conditioner, HC4-revise, punto fijo
batch_processor.py     recorre una carpeta
main.py                CLI
test_paper.py          regresion sobre el ejemplo (1) del paper
test_nsc.py            conditioner del paper, P*A = PA, soundness de NSC
```

El parseo usa el módulo `ast` de Python con un recorrido iterativo. No usa SymPy,
porque SymPy distribuye `-6*(x1+x3)` y así destruye las subexpresiones comunes que
motivan la técnica. `GUIA-refactor-ScriptDAG.md` explica paso a paso por qué cada
módulo quedó así.

## Auditoría

`validator.py` corre estas comprobaciones, que no dependen del parser del pipeline:

| | Qué comprueba |
|---|---|
| V1 | Las restricciones del `.txt` son las mismas que se parsearon (nada se descarta en silencio). |
| V2 | Cada fila de `A` reproduce exactamente, en racionales, su nodo suma. |
| V3 | Equivalencia numérica contra el texto **original**, evaluado de forma independiente en puntos aleatorios. |
| V4a | Soundness de los dominios: el valor de cada nodo en puntos muestreados cae dentro de su intervalo. |
| V4b | Isotonía por inclusión: en sub-cajas, los intervalos se achican. |
| V4c | Soundness local de cada operador, muestreando directamente en los dominios de sus argumentos. |
| V5 | No se introdujeron variables auxiliares nuevas. |

Estas comprobaciones son **evidencia empírica, no una demostración**. V2 es exacta;
V3 y V4 son muestreos. Para confirmar que el validator detecta errores reales, se
puede introducir a propósito un bug en `intervals.reciproco` y verificar que V4c
falla.

## Formato del `.npz`

```python
import numpy as np
z = np.load('outputs/inst001.npz')
A, c, y, w = z['A'], z['c'], z['y_names'], z['w_names']
```

```
A, c                       (m,n) y (m,)  float64    w = c + A·y
A_aug, rhs_aug             [A | -I] y -c            A_aug·[y;w] = rhs_aug
A_fila/A_col/A_num/A_den   A exacta, COO de racionales (strings)
y_names, y_lo, y_hi        columnas y sus dominios
w_names, w_lo, w_hi        filas y sus dominios
w_es_columna               indice de columna de cada w, o -1
var_names, var_lo, var_hi  variables ORIGINALES del problema
prod_nombre/tipo/dato/args subsistema no lineal (Proposicion 1)
restricciones              ecuaciones cascaron, 'nodo op rhs'
P_fila/P_col/P_num/P_den   conditioner P exacto (columnas = filas de A)
PA_fila/.../PA_den         P·A exacto
PA_pivote                  columna del pivote de cada fila de PA, o -1 si es dependiente
meta_*                     conteos para auditoria
```

`A_exacta` (las cuatro columnas COO) es la fuente de verdad. `A` en float64 está
solo por conveniencia.

## Limitaciones conocidas

- **El parser cubre solo parte de Minibex.** No acepta: la sección `constants`,
  desigualdades dobles (`-1 <= x^2 <= 4`), restricciones que ocupan varias líneas,
  arreglos de variables (`x[3] in [0,1]`) ni la sección `function`. Cuando una
  línea no se puede parsear, el parser levanta `ParseError` en vez de saltársela.
- `abs`, `min` y `max` se parsean, pero su dominio queda en `(-inf, inf)`: es sound,
  pero no sirve para filtrar.
- Si el coeficiente de la raíz no es 1, el `rhs` se divide en float. NSC lo
  compensa ensanchando ese valor 2 ulps. Los literales
  decimales se interpretan como el racional escrito (`0.1` → `1/10`).
- Solo hay 3 instancias de prueba, y `inst001`/`inst002` vienen del mismo generador.
