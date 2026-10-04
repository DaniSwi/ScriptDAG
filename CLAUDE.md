# ScriptDAG — contexto para Claude Code

Implementa la Sección 3 de **Araya & Reyes (2019), *Enhancing interval constraint
propagation by identifying and filtering n-ary subsystems*** (J. Glob. Optim. 74:1–20).
Lee instancias `.txt` tipo Minibex, construye un DAG fusionando subexpresiones comunes,
identifica las sumas n-arias y extrae el subsistema lineal asociado (Secciones 3.1–3.2).
Encima de eso, `nsc.py` implementa el contractor **NSC** (Sección 3.3).

**Qué falta:** la variante NSC+Simplex (§3.3.1), porque necesita un LP solver con cotas
certificadas y el repo solo depende de numpy; y un solver branch & bound que bisecte.
Sin ese solver no se puede reproducir la Tabla 1 del paper: NSC contrae una caja, pero
no explora el árbol de búsqueda.

## Comandos

```bash
python3 main.py --lote            # procesa instances/ -> outputs/ + auditoría
python3 main.py --lote --nsc      # además compara REF (HC4 en el DAG) vs NSC
python3 main.py instances/inst001.txt
python3 main.py instances/inst001.txt --txt   # además el .txt reescrito (opcional)
python3 test_paper.py             # regresión contra el ejemplo (1) del paper
python3 test_nsc.py               # conditioner del paper, P·A = PA, soundness de NSC
```

Sin dependencias salvo **numpy**. Requiere **Python ≥ 3.9** (`math.nextafter`).
No hay `requirements.txt`; `pip install numpy` basta. Hay un `.venv` local (Python 3.14,
no versionado).

**Verificación esperada:**
- `main.py --lote` imprime `3 exitosos | 0 fallidos` con 21 `[PASS]` y cero `[FAIL]`.
- `test_paper.py` y `test_nsc.py` terminan en `Todo OK.` (`test_nsc.py` tarda unos 6 s).

Si tocas algo del núcleo, córrelos todos antes de dar nada por bueno.

El PDF del paper está fuera del repo, en
`C:\Users\PC RST GALAX\Downloads\2019_-_Enhancing_interval_constraint_propagation_by_identifying.pdf`.

## Arquitectura

```
parser_instancias.py   .txt -> Instancia (vars, dominios, ops, rhs, raices, dag)
  dag_builder.py       ast de Python -> DAG normalizado (hash-consing)   <-- el núcleo
matrix_extractor.py    DAG -> SistemaLineal (A exacta, c, y, w, subsistema de productos)
  intervals.py         aritmética de intervalos sound con redondeo dirigido
  dominios.py          propagación de dominios por el DAG + muestreadores
nsc.py                 contractor NSC: conditioner, HC4-revise, punto fijo (Sec. 3.3)
exportator.py          construir() = pipeline completo; exportar_npz / exportar_txt
validator.py           7 comprobaciones (V1, V2, V3, V4a, V4b, V4c, V5)
batch_processor.py     recorre una carpeta
main.py                CLI
test_paper.py          regresión sobre el ejemplo (1)
test_nsc.py            pruebas de NSC
```

El sistema extraído es **`w_i = c_i + Σ_j A[i][j]·y_j`**, con una fila por cada nodo
suma n-aria del DAG y columnas `y_j` que son nodos del DAG: variables, potencias,
productos, funciones, **y también otros `w_k`** (columnas mixtas, el `y5 = b5` de la
Fig. 3 del paper). Ojo con la notación: `c` es el vector de **constantes aditivas**;
el vector de raíces de las sumas es `w`. El `b` del paper corresponde a `w − c`.

### NSC (`nsc.py`)

`contraer(Preparado(inst, sis, cond))` repite hasta punto fijo (ratio 1e-3):

1. `hc4_revise`: forward + backward sobre todo el DAG, con los dominios de los nodos
   intermedios persistentes entre pasadas.
2. `proyectar_lineal`: proyecta cada fila de `PA·y = P·(w − c)` sobre sus `y`.

Con `cond=None` y `usar_lineal=False` es REF (HC4 sobre el DAG), la referencia del paper.
Las decisiones que el paper no fija están documentadas como N1–N7 en el docstring.
La más importante es N1: el pivote de Gauss–Jordan se elige por columna, en orden, y es
el de mayor magnitud; en caso de empate se elige la fila con menos no-ceros, y los
pivotes no se normalizan. Con esa regla se reproduce exactamente la `P` de la §3.3.2.

## Invariantes — no romper

Estos puntos son el resultado de una revisión que encontró que la versión anterior
**no implementaba la técnica del paper**. Si un cambio los viola, es una regresión,
por muy limpio que se vea el diff.

1. **Cero variables auxiliares nuevas.** Productos, potencias y funciones son *nodos
   del DAG*, o sea columnas de `A`. Nunca variables `v_i` del problema. La versión
   anterior creaba una `v_i` por cada producto repetido y eso absorbía los coeficientes:
   `4*x2` se volvía `v_0` con coeficiente `1`, dejando una matriz de incidencia de ±1
   sobre la cual Gauss–Jordan no tiene nada que hacer.
2. **Un nodo `mul` o `pow` jamás lleva coeficiente numérico.** El coeficiente se saca
   siempre hacia afuera al construir (`(c·n)^k = c^k · n^k` para `k` entero). Por eso
   `x2^2`, `-x2^2` y `5*x2^2` comparten **una sola** columna.
3. **No reintroducir SymPy en el pipeline.** Distribuye `-6*(x1+x3)` automáticamente y
   con eso destruye la CSE que motiva el paper. Además `sympify(..., evaluate=False)`
   revienta con `RecursionError` en la línea de 33 953 caracteres de `instprofe.txt`.
   El parseo es con el módulo `ast` y una normalización explícita que es el Algoritmo 1.
4. **El recorrido del AST tiene que ser iterativo.** `instprofe.txt` tiene cadenas
   `BinOp` de profundidad 1492; cualquier visitante recursivo revienta.
5. **Coeficientes en `Fraction`, no en float.** `A_exacta` es la fuente de verdad;
   `sis.A` (float64) es una conveniencia. Las comparaciones de equivalencia se hacen
   en racionales. `P` y `PA` también se calculan en `Fraction`.
6. **Los dominios solo pueden sobre-aproximar.** Si no se sabe acotar, el valor es
   `(-inf, +inf)`, nunca `[-1e8, 1e8]`. Recortar un dominio pierde soluciones.
   Todo redondeo va hacia afuera (`intervals.dn` / `intervals.up`). Lo mismo vale
   para NSC: una proyección que no se sabe invertir **no contrae**, y una división por
   un intervalo que contiene 0 no contrae.
7. **El parser falla ruidosamente.** Una línea de restricción que no se puede parsear
   levanta `ParseError` con archivo y línea. Nunca se descarta en silencio: esa era la
   razón por la que el validator viejo podía decir "ÉXITO ABSOLUTO" con restricciones
   perdidas.

### La prueba de que el validator sirve

No basta con que las 7 comprobaciones pasen. Si tocas `intervals.py` o `dominios.py`,
verifica que el validator **muerde** reintroduciendo un bug a propósito:

```python
import intervals as iv
from exportator import construir
from validator import validar
def buggy(I):
    a, b = I
    return (-0.5, 0.5) if a <= 0 <= b else (min(1/b, 1/a), max(1/b, 1/a))
iv.reciproco = buggy
validar(*construir('instances/inst001.txt'))   # debe salir [FAIL] en V4c
```

V4c (soundness local por operador) es la única que lo detecta: muestreando desde las `x`
el bug es invisible, porque para que `1/w` se salga de `[-0.5, 0.5]` hace falta `|w| < 2`
y ningún punto aleatorio del espacio de las `x` llega ahí.

### La prueba de que `test_nsc.py` sirve

Si tocas `nsc.py`, comprueba que `test_nsc.soundness(inst, sis, 40)` devuelve
violaciones cuando reemplazas por monkeypatch alguna de estas funciones con un bug:
- `_proyectar_suma`: usar `a` en lugar de `1/a`.
- `_division_proyeccion`: hacer `num * [1/hi, 1/lo]` aunque el divisor contenga 0.
- `_proyectar_raiz_entera`: olvidar la rama negativa de una raíz par.

Las tres se verificaron. Los dos primeros bugs se detectan en el ejemplo (1) y en
inst001; el tercero, solo en inst001.

## Formato del `.npz`

```
A, c                     (m,n) y (m,)  float64   w = c + A·y
A_aug, rhs_aug           [A | -I] y -c           A_aug·[y;w] = rhs_aug
A_fila/A_col/A_num/A_den A exacta, COO de racionales (strings)
y_names, y_lo, y_hi      columnas y sus dominios
w_names, w_lo, w_hi      filas y sus dominios
w_es_columna             índice de columna de cada w, o -1
var_names, var_lo, var_hi    variables ORIGINALES del problema
prod_nombre/tipo/dato/args   subsistema no lineal (Proposición 1)
restricciones            ecuaciones cascarón, 'nodo op rhs'
P_fila/P_col/P_num/P_den conditioner exacto (columnas de P = filas de A)
PA_fila/.../PA_den       P·A exacto
PA_pivote                columna del pivote de cada fila de PA, o -1 si es dependiente
meta_*                   conteos para auditoría
```

Las filas de `P` y `PA` van primero las de pivote, en orden de pivote, y después las
dependientes.

Decisión de diseño: **la salida oficial no modifica el sistema original.** El paper dice
*"we do not modify the original system"* y *"we do not need concretization"*; las `w_i`
y `y_j` son etiquetas de nodos, no variables nuevas. El `.txt` reescrito (estilo
I-CSE/Ceberio) existe como opción con `--txt`, marcando las auxiliares como no
bisectables.

## Estado conocido y límites

**NSC no gana nada en las 3 instancias de `instances/`.** Partiendo de la caja inicial
([-1e8, 1e8]), NSC contrae exactamente lo mismo que REF. En el ejemplo (1) del paper,
con dominios `[-1, 5]`, sí gana: reduce el volumen en 10^1.36 contra 10^0.15 de REF.
El paper mide la ganancia dentro del branch & bound, en cajas pequeñas.

**Cobertura del parser: 5 de 10 variantes realistas de Minibex.** Medido. Rechaza:
sección `constants`, desigualdad doble (`-1 <= x^2 <= 4`), restricción multilínea,
arreglos de variables (`x[3] in [0,1]`), sección `function`. Acepta `min`/`max`/`abs`,
exponente fraccionario, notación científica, variables a ambos lados, `exp`/`log`.
**Este es el principal bloqueador para que el repo sea usable por terceros.**

**`abs`, `min` y `max` parsean pero dan dominio `(-inf, inf)`** — sound pero inútil
para filtrar. Faltan sus reglas en `intervals.funcion` (solo existe `Abs`, con
mayúscula). `exp` y `log` sí están bien. En NSC, `sin`, `cos` y las demás funciones
sin regla de inversión no contraen hacia atrás.

**Aproximación.** La extracción del subsistema lineal es exacta en racionales; los
dominios son sobre-aproximaciones. Quedan dos grietas menores:
- Cuando el coeficiente de la raíz no es 1, el parser divide el `rhs` en float. NSC lo
  compensa ensanchando ese valor 2 ulps (N6 en `nsc.py`); la auditoría no lo compensa.
- Los literales decimales se interpretan como el racional que la persona escribió
  (`0.1` → `1/10`), no como el double.

**El ejemplo (1) de `test_paper.py` sigue reconstruido**, no es el literal del paper.
`test_nsc.py` sí usa las tres ecuaciones literales (página 2 del PDF). Dos decisiones de
diseño quedaron sin árbitro:
- Si `c·(suma)` conserva la subsuma como nodo compartido o la distribuye. Hoy es un
  flag, `distribuir_sumas_escaladas`, con default conservador.
- Si una suma escalada dentro de una función genera su propia fila. Hoy sí, y con eso
  se reproduce la fila `6|6` y la fila dependiente −1/6 del paper.

**Solo hay 3 instancias** en `instances/`, y `inst001`/`inst002` vienen del mismo
generador. No hay un benchmark real contra el cual medir. Los benchmarks del paper son
los de COPRIN y los de https://github.com/vareyesr/ncsp_generator.

### Pendientes

- No hay licencia. Sin ella, formalmente nadie puede reutilizar el código.
- NSC+Simplex (§3.3.1) y un branch & bound para medir contra la Tabla 1.

## Convenciones

- Comentarios y docstrings en español, **sin acentos** en los `.py` (evita líos de
  encoding al mover archivos entre Windows y Linux). Los `.md` sí llevan acentos.
- Finales de línea **LF**. Ojo en Windows: si editas con un script de Python en modo
  texto, escribe con `newline=''` o el archivo queda en CRLF y el diff completo cambia.
- Cada módulo abre con un docstring que dice **qué bug o desviación arregla**, con el
  código (D1…D6 para desviaciones respecto del paper, B1…B7 para bugs, N1…N7 para las
  decisiones de NSC). Mantén esa trazabilidad al modificar.
- Nombres de funciones y variables en español, igual que el resto del repo.
- `GUIA-refactor-ScriptDAG.md` documenta paso a paso por qué cada archivo quedó así.
  Léelo antes de cambiar el núcleo (`dag_builder.py`, `matrix_extractor.py`).

## Referencia

Araya, I. & Reyes, V. (2019). *Enhancing interval constraint propagation by identifying
and filtering n-ary subsystems.* J. Glob. Optim. 74:1–20. doi:10.1007/s10898-019-00738-5.
Secciones 3.1, 3.2 y 3.3.2 implementadas; 3.3.1 (simplex) pendiente.
