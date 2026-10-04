# ScriptDAG: identificación de subsistemas n-arios en NCSP

Preprocesador para problemas de satisfacción de restricciones numéricas (NCSP) que
implementa las Secciones 3.1–3.2 de:

> Araya, I. & Reyes, V. (2019). *Enhancing interval constraint propagation by
> identifying and filtering n-ary subsystems.*

Lee instancias `.txt` en formato tipo Minibex, construye un DAG que fusiona las
subexpresiones comunes, identifica las sumas n-arias y extrae el subsistema lineal
asociado, junto con dominios sound para cada nodo.

**Alcance.** El proyecto cubre solo el preprocesamiento. La Sección 3.3 (filtrado
NSC: conditioner por Gauss–Jordan, HC4-revise, iteración a punto fijo) **todavía no
está implementada**. El `.npz` de salida contiene todo lo que hace falta para
construirla encima.

---

## Instalación

Necesitas **Python ≥ 3.9** (por `math.nextafter`) y **numpy**. Ninguna otra
dependencia.

```bash
python -m venv .venv
```

Activa el entorno según tu shell:

```powershell
.\.venv\Scripts\Activate.ps1
```

```bash
source .venv/bin/activate
```

Y luego instala numpy:

```bash
pip install numpy
```

## Uso

Para procesar todas las instancias de `instances/` y escribir los resultados en
`outputs/` (incluye la auditoría):

```bash
python main.py --lote
```

Para procesar una sola instancia:

```bash
python main.py instances/inst001.txt
```

Si además quieres el `.txt` reescrito con variables auxiliares (opcional, ver más
abajo):

```bash
python main.py instances/inst001.txt --txt
```

Para la regresión contra el ejemplo (1) del paper:

```bash
python test_paper.py
```

Otras opciones: `--sin-dominios` y `--sin-validar`.

**Qué deberías ver:** `main.py --lote` imprime `3 exitosos | 0 fallidos`, con 21
`[PASS]` y ningún `[FAIL]`. `test_paper.py` termina en `Todo OK.`

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

## Arquitectura

```
parser_instancias.py   .txt -> Instancia (vars, dominios, ops, rhs, raices, dag)
  dag_builder.py       ast de Python -> DAG normalizado (hash-consing, Algoritmo 1)
matrix_extractor.py    DAG -> SistemaLineal (A exacta, c, y, w, subsistema de productos)
  intervals.py         aritmetica de intervalos con redondeo dirigido
  dominios.py          propagacion de dominios por el DAG + muestreadores
exportator.py          pipeline completo; exportar_npz / exportar_txt
validator.py           auditoria (V1-V5)
batch_processor.py     recorre una carpeta
main.py                CLI
test_paper.py          regresion sobre el ejemplo (1) del paper
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
- Si el coeficiente de la raíz no es 1, el `rhs` se divide en float. Los literales
  decimales se interpretan como el racional escrito (`0.1` → `1/10`).
- Solo hay 3 instancias de prueba, y `inst001`/`inst002` vienen del mismo generador.
