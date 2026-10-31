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
(§3.3.2) y un branch & bound por intervalos (`bb.py`) para medirlo como en el paper:
cuántas cajas y cuánto tiempo necesita cada estrategia. La variante con simplex
(§3.3.1) no está implementada, porque necesita un LP solver con cotas certificadas.
El solver es Python puro, así que no alcanza a resolver sistemas del tamaño de
`inst001` en un tiempo razonable (ver *Limitaciones*).

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

```bash
python test_bb.py
```

Los tres deben terminar en `Todo OK.` `test_nsc.py` tarda unos segundos y
`test_bb.py`, unos 15 segundos.

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

### 8. Resolver con branch & bound (opcional)

Agrega `--bb` para resolver el sistema completo dos veces, primero con REF y después
con NSC. Imprime cuántas cajas y cuánto tiempo usó cada estrategia, y la ganancia
`gt` (tiempo) y `gb` (cajas) definida igual que en la Tabla 1 del paper. Con el
ejemplo (1) del paper:

```
  REF :     203 cajas,    2 cajas solucion (2 certificadas),     0.35s  [completo]
  NSC :      45 cajas,    2 cajas solucion (2 certificadas),     0.17s  [completo]
  ganancia NSC vs REF: gt = 2.13  gb = 4.51
```

Cada estrategia tiene un límite de 600 s y de 100 000 cajas. Si se alcanza, la línea
dice `TRUNCADO` y no se calcula la ganancia, porque un árbol incompleto no permite
comparar. Con las instancias incluidas pasa eso: son de 8 variables con dominios
`[-1e8, 1e8]` y el árbol no termina.

### 9. Reproducir el experimento de la Sección 4.2 (opcional)

`experimento.py` genera instancias aleatorias con la misma estructura que las del
paper, las resuelve con REF y NSC, y resume la ganancia relativa promedio de NSC
(Definición 3 del paper) para cada número de sumas P:

```bash
python experimento.py --tipo mult
```

```bash
python experimento.py --tipo func
```

`mult` corresponde a las ecuaciones (6) del paper (productos de sumas, Fig. 6) y
`func` a las (7) (sumas de funciones de sumas, Fig. 7). Por defecto usa los mismos
parámetros del paper (n = 8, Q = 100, P de 8 a 16, 20 instancias por P) salvo los
dominios, que son `[-2, 2]` en lugar de `[-1e8, 1e8]`. Con todo por defecto tarda
horas; para una prueba rápida:

```bash
python experimento.py --tipo mult --P 8 10 --instancias 3 --limite 60
```

Los resultados quedan en `experimentos/resultados_<tipo>.csv`, y las instancias
generadas, en `experimentos/instancias/`.

### Todas las opciones

| Opción | Qué hace |
|---|---|
| `--lote` | Procesa toda la carpeta `instances/` (es lo que pasa si no se da ningún archivo). |
| `--nsc` | Compara además cuánto contraen REF y NSC la caja inicial. |
| `--bb` | Resuelve con branch & bound usando REF y NSC, y compara cajas y tiempo. |
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
pequeñas; ver la sección siguiente.

## El branch & bound (`bb.py`)

`bb.resolver` explora el árbol en profundidad. En cada caja:
- la contrae con REF o con NSC;
- si queda vacía, la descarta;
- si todas las variables miden menos que `eps` (1e-8 por defecto), la guarda como
  caja solución;
- si no, la bisecta por la mitad, en la variable que elige SmearSumRel (el
  bisector del paper).

Después de contraer con REF o NSC se aplica un **Newton intervalar** (`newton.py`,
Hansen–Sengupta), igual para las dos estrategias, como hace Ibex en el paper. Newton
además **certifica** las cajas solución: una caja certificada contiene exactamente
una solución. Newton solo se aplica a sistemas cuadrados de igualdades. Las
decisiones que el paper deja a Ibex están documentadas como S1–S7 en `bb.py` y
T1–T5 en `newton.py`.

**Resultados.** En el ejemplo (1) del paper, REF necesita 203 cajas y NSC 45
(`gb = 4.5`, `gt = 2.1`). Ambos encuentran las mismas 2 soluciones, y las dos quedan
certificadas. El paper reporta 83 cajas para HC4 sobre el DAG y 33 para NSC. En
`inst001` ninguna de las dos estrategias termina en 30 minutos (ver *Limitaciones*).

`test_bb.py` comprueba cuatro cosas:
- `sin` y `cos` por intervalos son sound y ajustados.
- El jacobiano por intervalos cumple el valor medio: `f(x) − f(y) ∈ J(X)·(x − y)`.
  Es la propiedad de la que depende que Newton sea sound.
- En el ejemplo (1), las dos estrategias encuentran las mismas 2 soluciones, las
  certifican y NSC usa menos cajas.
- Soundness del árbol completo: con una solución fabricada como en `test_nsc.py`,
  alguna caja solución la contiene.

La prueba del jacobiano detecta derivadas con bug (por ejemplo, derivar `sin` como
`sin`). La del árbol es más débil: una bisección que pierde el 10% de cada mitad la
hace fallar en 9 de 24 árboles, y si pierde solo el 1% no la detecta.

## Instancias aleatorias y experimento (Sección 4.2)

`generador.py` replica la estructura del generador del paper
([vareyesr/ncsp_generator](https://github.com/vareyesr/ncsp_generator), del que salieron
`inst001` e `inst002`), con los dominios como parámetro:

- **`mult`, ecuación (6):** cada restricción es `Π_k (a + Σ_i b·t)^d = R`, con
  `d ∈ {−1, 1}`.
- **`func`, ecuación (7):** cada restricción es `Σ_k u(a + Σ_i b·t) = R`, con
  `u ∈ {log, exp, x⁻¹, identidad}`.
- Los términos `t` salen de un pool con las n variables y algunas potencias al
  cuadrado o al cubo. `a` y `b` son enteros en `[-5, 5]`.
- En total hay P sumas y Q términos. `R` se calcula en un punto aleatorio de
  `[-1, 1]ⁿ`, que queda anotado como solución conocida.

`experimento.py` resuelve cada instancia con REF y NSC y calcula `gt` y `gb` con la
Definición 3. Solo promedia las instancias que **ambas** estrategias terminaron, y
comprueba en cada una que la solución conocida quedó dentro de una caja solución.
Las decisiones están documentadas como G1–G2 y E1–E4 en los docstrings. Con
`--reanudar` se retoma una corrida interrumpida sin repetir los árboles ya resueltos.

**Resultados.** 10 instancias por P y tipo, límite de 300 s por árbol, dominios
`[-2, 2]`. `n` es el número de instancias que ambas estrategias terminaron (las
únicas que entran al promedio); gt y gb > 1 significa que NSC gana:

| P | `mult` n | `mult` gt | `mult` gb | `func` n | `func` gt | `func` gb |
|---|---|---|---|---|---|---|
| 8 | 7 | **2.46** | **3.62** | 4 | **2.48** | **3.18** |
| 10 | 8 | 1.02 | 1.34 | 7 | 1.14 | 1.58 |
| 12 | 7 | 0.79 | 1.04 | 8 | 0.79 | 1.14 |
| 14 | 3 | 0.79 | 1.00 | 5 | 0.72 | 1.02 |
| 16 | 2 | 0.85 | 1.03 | 4 | 0.75 | 1.05 |

- **Se reproduce la tendencia de las Fig. 6 y 7 del paper:** NSC gana claramente con
  pocas sumas (P = 8: entre 3 y 3.6 veces menos cajas, 2.5 veces más rápido), y la
  ganancia cae a ~1 a medida que crece P. El paper lo explica igual: con más sumas,
  la caja de `b` se ensancha y la proyección filtra menos.
- **La magnitud es menor que en el paper**, que reporta ganancias de 30 a 90 con
  P = 8. La diferencia más probable son los dominios: con `[-1e8, 1e8]` el filtrado
  inicial pesa mucho más.
- **Con P ≥ 12, gt queda en ~0.8:** NSC no reduce el árbol y su costo extra se nota.
  En el paper queda en ~1; en Python, la proyección pesa más en proporción.
- **Ninguna solución perdida:** en todos los árboles terminados, la solución conocida
  quedó dentro de una caja solución.
- **Muchos árboles truncados,** sobre todo con P grande, por los polos de los
  términos `(suma)^(-1)` (ver *Limitaciones*). Los promedios de P = 14 y 16 se basan
  en pocas instancias (2 a 5).

Los datos completos están en `experimentos/resultados_mult.csv` y
`experimentos/resultados_func.csv`.

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
bb.py                  branch & bound por intervalos sobre REF o NSC
newton.py              Newton intervalar, jacobiano por intervalos, SmearSumRel
generador.py           instancias aleatorias como las de la Sec. 4.2 del paper
experimento.py         experimento de la Sec. 4.2: REF vs NSC por cada P
batch_processor.py     recorre una carpeta
main.py                CLI
test_paper.py          regresion sobre el ejemplo (1) del paper
test_nsc.py            conditioner del paper, P*A = PA, soundness de NSC
test_bb.py             sin/cos, jacobiano, ejemplo (1) con B&B, soundness, generador
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
- **El branch & bound no termina en `inst001`.** En 30 minutos, REF revisó 328 115
  cajas y NSC 237 573, y ninguno encontró una solución todavía: el árbol se pierde en
  los dominios de `[-1e8, 1e8]`. Va a unos 3–5 ms por caja, uno o dos órdenes de
  magnitud más lento que Ibex en C++.
- `sin` y `cos` no contraen hacia atrás (no hay regla de inversión); solo hacia
  adelante.
- Solo hay 3 instancias de prueba, y `inst001`/`inst002` vienen del mismo generador.
