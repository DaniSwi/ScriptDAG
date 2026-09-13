# ScriptDAG → NSC: paso a paso, archivo por archivo

Objetivo de esta pasada (lo que acordamos):

- **Alcance:** fiel al paper completo — fix D1/D2, bugs unsound B1–B7, y parseo
  sin SymPy para no perder la CSE de `-6*(x1+x3)` (D3).
- **Salida oficial:** `.npz` **anexo**, sin tocar el sistema original (D4 resuelto
  a favor del estilo paper). El `.txt` reescrito queda como opción (`--txt`).
- **Fuera de alcance:** la Sección 3.3 (conditioner `P`, simplex, punto fijo).
  Queda anotada al final como el paso siguiente.

Resultado verificado sobre el repo actual:

| | antes | después |
|---|---|---|
| `inst001` | `A` 10×**36**, 33 columnas `v_i`, entradas ±1 | `A` 10×**12**, coeficientes −5…5 |
| `inst002` | igual | `A` 10×**12** |
| `instprofe` | `A` 3×**4343** + **4339 variables nuevas**, ~33 s | `A` 3×**1513**, **0 variables nuevas**, 96.7 % densa, **0.40 s** |
| ejemplo (1) | coeficientes 4 y 2 perdidos | `A` 5×8, filas `4\|1\|1` · `1\|1` · `1\|1\|2` · `6\|6` · `4\|1\|1\|-6\|1` |
| auditoría | 1 check (equivalencia algebraica) | 7 checks, y **muerde** si reintroduces B1 |

Árbol final (10 archivos, 3 nuevos):

```
intervals.py        NUEVO   aritmética de intervalos sound       (B1, B2, B6)
dag_builder.py      REESCRITO  DAG + Algoritmo 1                 (D1, D2, D3, D6, B7)
parser_instancias.py REESCRITO parser + dominios                 (B3, B5, D3)
matrix_extractor.py REESCRITO  Sección 3.2                       (D1, D2, D5)
dominios.py         NUEVO   propagación de dominios por el DAG   (B1, B2)
exportator.py       REESCRITO  .npz completo, .txt opcional      (D4, D5)
validator.py        REESCRITO  auditoría de verdad               (B5)
batch_processor.py  RETOCADO   validación reactivada
main.py             RETOCADO   CLI, se borra asd()               (B4)
test_paper.py       NUEVO   regresión contra el ejemplo (1)
```

---

## Paso 0 — preparación

```bash
git checkout -b nsc-fiel-al-paper
python3 -c "import numpy, sympy"   # sympy ya NO es necesario, pero no estorba
```

El pipeline nuevo **no importa SymPy en ninguna parte**. Esa es la decisión de
fondo: SymPy era a la vez lo bueno del repo (te regalaba el hash-consing) y lo
que rompía la técnica (distribuye `-6*(x1+x3)` antes de que puedas verlo). Se
reemplaza por el módulo `ast` de Python + una normalización explícita, que es
literalmente `generate_sum`/`generate_product` del Algoritmo 1.

Orden de trabajo obligatorio: **1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10**.
Cada paso deja el repo ejecutable, así que puedes commitear en cada uno.

---

## Paso 1 — `intervals.py` (archivo NUEVO)

**Qué arregla:** B1, B2, B6.

**Qué borras:** las funciones `calcular_dominios_w` y `calcular_dominios_v` de
`exportator.py` — entera la lógica de intervalos sale de ahí (se hace en el
Paso 5). No dejes ninguna aritmética de intervalos dentro del exportador.

**Qué colocas:** un módulo con redondeo dirigido y una regla sound por
operación. Los tres puntos críticos:

```python
def dn(x): return math.nextafter(x, -math.inf)   # redondeo hacia afuera (B6)
def up(x): return math.nextafter(x,  math.inf)

def reciproco(I):
    a, b = I
    if a > 0.0 or b < 0.0:      return (dn(1.0/b), up(1.0/a))
    if a == 0.0 and b > 0.0:    return (dn(1.0/b), INF)      # [1/b, +inf)
    if b == 0.0 and a < 0.0:    return (NEG_INF, up(1.0/a))  # (-inf, 1/a]
    return ENTERO                                            # 0 en el interior

def potencia_entera(I, n):
    if n == 0: return (1.0, 1.0)
    if n < 0:  return reciproco(potencia_entera(I, -n))       # <- B1
    ...
```

La versión vieja hacía `if n % 2 != 0` con `n = -1`, y como en Python
`-1 % 2 == 1` tomaba la rama "impar" y devolvía `[lo**n, hi**n]`. Para
`w ∈ [-2,2]` eso daba `[-0.5, 0.5]`: un intervalo que **no contiene ningún
valor real** de `1/w`. Con `calcular_dominios=True` eso eliminaba todas las
soluciones.

Regla de oro del módulo: **si no sabes acotar, devuelve `(-inf, +inf)`**. Nunca
`[-1e8, 1e8]`.

**Verificación:**

```bash
python3 -c "import intervals as iv; print(iv.reciproco((-2,2)), iv.potencia_entera((-2,2),-2))"
# (-inf, inf) (0.25, inf)      <- antes: (-0.5,0.5) y (0,0.25)
```

---

## Paso 2 — `dag_builder.py` (reescritura completa)

**Qué arregla:** D1, D2, D3, D6, B7, y de paso elimina las variables globales.

**Qué borras:**

- `dag_nodes_cache` y `node_counter` globales, y `clear_cache()`. El DAG pasa a
  vivir en una instancia de `DagBuilder`, así que ya no hay estado que limpiar
  entre instancias (y ya no puedes olvidarte de limpiarlo).
- `DagNode.ref_count` y `DagNode.children` tal como estaban (D6: `children` se
  llenaba y nunca se usaba; el DAG "real" era el de SymPy).

**Qué colocas:** un DAG con hash-consing propio y seis tipos de nodo:

```
var   data = nombre                 args = ()
const data = Fraction               args = ()
add   data = (const, coefs)         args = (nodos sumandos)     -> una FILA de A
mul   data = None                   args = (factores)           -> una COLUMNA
pow   data = Fraction (exponente)   args = (base,)              -> una COLUMNA
func  data = nombre                 args = (argumentos)         -> una COLUMNA
```

Las dos reglas que cambian todo:

1. **Un nodo `mul` nunca lleva coeficiente numérico.** El coeficiente se saca
   siempre hacia afuera durante la construcción. Por eso `x2**2`, `-x2**2` y
   `5*x2**2` comparten una sola columna, y desaparece D2 (antes `v_0=x2**2`,
   `v_1=-v_0` y `v_19=5*v_0` eran tres columnas distintas y sus relaciones eran
   invisibles para el sistema lineal).

2. **Para potencias con exponente entero, el coeficiente también sale:**
   `(c·n)^k = c^k · n^k`. Así `(2x)^2` y `x^2` comparten nodo.

El normalizador devuelve una "forma" en vez de un nodo:

```python
('n', coef, nodo)        # coef * nodo
('s', const, terminos)   # suma PENDIENTE, todavía sin materializar
```

Las sumas se mantienen **pendientes** mientras el padre siga siendo `+` o `-`.
Sin eso, `a+b+c` (que Python parsea como `(a+b)+c`) generaría una fila espuria
por cada suma binaria intermedia. Se materializa con `_materializar()` en
cuanto el padre es un producto, una potencia o una función.

`generate_sum` queda así:

```python
def _suma(self, izq, der):
    const, terminos = Fraction(0), []
    for forma in (izq, der):
        if forma[0] == 's':                      # suma pendiente: se fusiona
            const += forma[1]; terminos.extend(forma[2]); continue
        c, nodo = forma[1], forma[2]
        if nodo is None:
            const += c
        elif nodo.kind == 'add' and (c == 1 or self.distribuir_sumas_escaladas):
            sub_const, sub_coefs = nodo.data    # aplanado n-ario
            const += c * sub_const
            terminos.extend((c*sc, sa) for sc, sa in zip(sub_coefs, nodo.args))
        else:
            terminos.append((c, nodo))          # <- aquí vive el -6*(x1+x3)
    return ('s', const, terminos)
```

Fíjate en la condición `c == 1`: una suma hija de una suma **se aplana**
(es lo que dice el paper), pero con coeficiente ≠ 1 **no**, porque aplanarla
obligaría a distribuir y eso es exactamente lo que destruía la CSE (D3). El
flag `distribuir_sumas_escaladas=True` recupera el comportamiento distribuido
si algún día quieres comparar.

`generate_product` es el mismo patrón con aplanamiento n-ario de factores, y
las divisiones entran como `Pow(-1)`:

```python
def _inverso(self, par):
    c, nodo = par
    if nodo is None: return Fraction(1,1)/c, None
    return Fraction(1,1)/c, self.pow(nodo, -1)
```

**Un detalle que te va a morder si no lo haces:** el recorrido tiene que ser
**iterativo**. `instprofe.txt` tiene una línea de 33 953 caracteres con una
cadena `BinOp` de profundidad 1492; cualquier visitante recursivo revienta con
`RecursionError`. Y `sympy.sympify(..., evaluate=False)` revienta ahí también
(su `EvaluateFalseTransformer` es recursivo) — esa es la razón práctica, además
de la conceptual, para dejar SymPy fuera.

```python
def _linealizar(self, raiz, variables):
    res, pila = {}, [(raiz, False)]
    while pila:
        nodo, procesado = pila.pop()
        if not procesado:
            pila.append((nodo, True))
            for h in self._hijos_ast(nodo): pila.append((h, False))
            continue
        res[id(nodo)] = self._combinar(nodo, res, variables)
    return res[id(raiz)]
```

Sobre D6: `ref_count` pasa a llamarse `in_degree` y se incrementa en el punto
correcto (cada vez que un padre referencia al hijo, también cuando el nodo ya
estaba en caché). `registrar_raiz()` suma la arista de la restricción.

`escalar(coef, nodo)` tiene un caso especial: si `nodo` es una suma, el
coeficiente se **reparte dentro de la suma**, de modo que `sin(6*(x1+x3))`
produce la fila `6x1+6x3` (la fila `6|6` de la referencia) en vez de un
producto por una constante.

**Verificación:**

```bash
python3 -c "
from dag_builder import DagBuilder
b=DagBuilder(); V={'x1':1,'x2':1,'x3':1}
c,n=b.desde_string('4*x2 + sin(x1+x3) + 2*x2*sin(x1+x3)', V)
print(n.to_str()); print(b.resumen())"
# sin(x3 + x1) + 4*(x2) + 2*((sin(x3 + x1))*(x2))
# {'nodos': 8, ..., 'compartidos_CSE': 5}
```

Lo que tienes que ver: los coeficientes `4` y `2` **fuera** de los nodos, y
`x1+x3` como un único nodo compartido.

---

## Paso 3 — `parser_instancias.py` (reescritura)

**Qué arregla:** B3, B5, y conecta con el DAG nuevo.

**Qué borras:**

- `import sympy` y las dos llamadas a `sp.sympify`.
- `rhs_val = float(rhs_str)` (revienta si el lado derecho no es constante).
- El `if len(parts) == 3:` **sin `else`**, que descartaba restricciones en
  silencio. Ese era el detonante de B5: el validator usaba el mismo parser, así
  que una restricción perdida era invisible y la auditoría igual decía
  "ÉXITO ABSOLUTO".

**Qué colocas:**

```python
_NUM = r'[+-]?\s*(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?|inf|oo|INF)'
_RE_DOMINIO = re.compile(r'\[\s*(' + _NUM + r')\s*,\s*(' + _NUM + r')\s*\]')
```

El `\s*` después del signo es el arreglo de B3: `'[- 2, 2]'` (que está en
`instprofe.txt`) se leía como `[2.0, 2.0]`, o sea `y` fijada en 2. Unsound.

Conteo duro de restricciones:

```python
inst.lineas_restriccion += 1          # se cuenta TODA línea de la sección
...
if len(partes) != 3:
    raise ParseError(f"{ruta}:{nro}: restriccion no parseable: ...")
...
if inst.n_restricciones() != inst.lineas_restriccion:
    raise ParseError(...)
```

Y el lado derecho ya no tiene que ser un `float`: si trae variables se pasa
todo a la izquierda (`lhs - rhs`) y la restricción queda contra `0`. Si el
coeficiente de la raíz no es 1 y el operador es una desigualdad, se divide y se
**invierte el sentido cuando el coeficiente es negativo** — es un detalle chico
que también pierde soluciones si se ignora.

La función devuelve un objeto `Instancia` (no una tupla) con `var_names`,
`var_domains`, `ops`, `rhs`, `raices`, `coef_raiz`, `dag`,
`lineas_restriccion`. Los dominios ya salen como `(lo, hi)` de floats, no como
strings: nadie vuelve a parsear `'[-1e8,1e8]'` río abajo.

**Verificación:**

```bash
python3 -c "
from parser_instancias import parse_dominio, parse_instance
print(parse_dominio('[- 2, 2]'))            # (-2.0, 2.0)
i = parse_instance('instances/instprofe.txt')
print(i.n_restricciones(), i.dag.resumen())"
# 3 {'nodos': 1518, 'por_tipo': {'var': 6, 'pow': 21, 'mul': 1488, 'add': 3}, ...}
```

Ojo al tiempo: esto tarda 0.3 s. La versión con SymPy tardaba ~33 s.

---

## Paso 4 — `matrix_extractor.py` (reescritura)

**Qué arregla:** D1, D2, D5. Es el archivo donde estaba el problema central.

**Qué borras — esto es lo más importante de toda la guía:**

```python
# TODO ESTO SE VA:
v_cache = {}
v_equations = []
v_counter = 0
...
if isinstance(expr, sp.Mul):
    vars_in_mul = [arg for arg in new_expr.args if arg.free_symbols]
    is_n_ary_mult = len(vars_in_mul) > 1
    node = dag_nodes_cache.get(expr)
    is_repeated = node and node.ref_count > 1
    if is_n_ary_mult or is_repeated:
        v_sym = sp.Symbol(f'v_{v_counter}')     # <-- D1: aquí se comían los coeficientes
        ...
elif isinstance(expr, (sp.Pow, sp.Function)):
    ...                                          # <-- lo mismo
```

Las variables auxiliares `v_i` se creaban **antes** de que la rama `Add` hiciera
la separación coeficiente/nodo con `as_coeff_Mul()`. Como `4*x2` es un `Mul`
repetido, se volvía `v_0 = 4*x2` y entraba a la matriz con coeficiente `1.0`.
Por eso `A` salía como matriz de incidencia de ceros y unos: sobre eso el
conditioner de Gauss–Jordan (que elige como pivote el de mayor magnitud) no
tiene nada que hacer, y el simplex no tiene de dónde sacar cotas.

En el diseño nuevo **la separación coeficiente/nodo ya la hizo el DagBuilder**,
así que el extractor es casi trivial: recorre el DAG, una fila por cada nodo
`add`, y las columnas son los nodos hijos.

```python
for nodo in dag.nodes:                  # una fila por suma n-aria
    if nodo.kind == 'add':
        nodo.label = f"w_{len(sis.w_nodes)}"
        ...

for nodo in sis.w_nodes:
    const, coefs = nodo.data
    fila = {}
    for coef, hijo in zip(coefs, nodo.args):
        j = columna(hijo)
        fila[j] = fila.get(j, Fraction(0)) + coef    # "coefficients are summed up"
    sis.A_exacta.append(fila)
    sis.c_exacta.append(const)
```

Productos, potencias y funciones siguen siendo **nodos del DAG**, es decir
columnas de `A`; su filtrado le toca a HC4, no al subsistema lineal. Y las
columnas mixtas (un `w_k` que además es columna, el `y5 = b5` de la Fig. 3)
salen solas porque un nodo `add` puede ser hijo de otro `add`.

**Qué colocas además (D5):** se guarda `A_exacta` en `Fraction` (no solo el
float), y se emite el **segundo subsistema** (productos / potencias /
funciones) que la Proposición 1 pide y que antes no se exportaba:

```python
for nodo in dag.nodes:
    if nodo.kind in ('mul', 'pow', 'func'):
        sis.subsistema_productos.append((_nombre(nodo), nodo.kind, dato,
                                         [_nombre(a) for a in nodo.args]))
```

Y un método `aumentada()` que devuelve `[A | -I]` con `[A|-I]·[y;w] = -c`,
que es la forma en que un consumidor externo puede armar el sistema sin
re-parsear strings.

Fíjate en el cambio de semántica de `b`: en el código viejo `b` era el vector
de **constantes aditivas**; en el paper `b` es el vector de **auxiliares de las
raíces**. Por eso aquí el vector de constantes se llama `c` y las raíces se
llaman `w`, explícitamente. El sistema es `w = c + A·y`.

**Verificación:**

```bash
python3 -c "
from exportator import construir
i,s,d = construir('instances/inst001.txt')
print(s.shape()); print(s.resumen()); print(list(s.y_names))"
# (10, 12)
# {'filas_w': 10, 'columnas_y': 12, 'densidad': 0.8333, 'coef_min': -5.0,
#  'coef_max': 5.0, 'variables_auxiliares_nuevas': 0, 'nodos_producto': 15}
# ['x7', '(x7)^(2)', 'x4', '(x6)^(3)', '(x3)^(3)', 'x2', '(x0)^(2)', ...]
```

Si ves `10×12` con coeficientes −5…5, el bug D1 está muerto. Si ves `10×36`
con ±1, todavía queda una rama de productos viva.

---

## Paso 5 — `dominios.py` (archivo NUEVO)

**Qué arregla:** B1 y B2 a nivel de sistema.

**Qué colocas:** propagación de intervalos por el DAG en orden topológico.
Como los hijos siempre se crean antes que los padres (hash-consing), recorrer
`dag.nodes` en orden de `id` ya es un orden topológico válido — no hace falta
ordenar nada.

```python
for nodo in dag.nodes:
    if   nodo.kind == 'var':  dom[nodo.id] = var_domains.get(nodo.data, iv.ENTERO)
    elif nodo.kind == 'add':
        const, coefs = nodo.data
        acum = iv.from_fraction(const)
        for c, hijo in zip(coefs, nodo.args):
            acum = iv.suma(acum, iv.mult(iv.from_fraction(c), dom[hijo.id]))
        dom[nodo.id] = acum
    ...
```

Nota B6: el coeficiente exacto se convierte a intervalo con
`iv.from_fraction(c)`, no con `float(c)`. Para coeficientes enteros da lo
mismo; para cualquier racional no representable, `float()` introduce un recorte
inválido.

El default para una variable sin dominio es `(-inf, inf)`, **no** `[-1e8,1e8]`.
Ese era B2: `w_1` de `inst001` tiene rango real `[-2.6e9, 2.6e9]` y se
declaraba `[-1e8, 1e8]`; ese recorte puede perder soluciones.

Aquí van también los tres muestreadores que usa el validator
(`muestrear`, `sub_caja`, `muestra_en_intervalo`) y `evaluar_en_punto` /
`valor_local`. El muestreo **no es uniforme a propósito**: mezcla uniforme,
log-uniforme y extremos/ceros, porque el muestreo uniforme casi nunca produce
denominadores chicos, que es justo donde vive B1.

---

## Paso 6 — `exportator.py` (reescritura)

**Qué arregla:** D4 (decisión de salida) y D5 (contenido del `.npz`).

**Qué borras:** `parse_domain_string`, `calcular_dominios_w`,
`calcular_dominios_v` (se fueron a `intervals.py` + `dominios.py`), y el bloque
de `exportar_a_txt` que escribía las `v_i` como variables reales.

**Qué colocas:** `construir()` (parse → DAG → sistema → dominios),
`exportar_npz()` y `exportar_txt()`. El `.npz` ahora trae:

```
A, c, A_aug=[A|-I], rhs_aug
A_fila/A_col/A_num/A_den        A exacta, COO de racionales           (B6)
y_names, y_lo, y_hi             columnas y sus dominios sound
w_names, w_lo, w_hi             filas y sus dominios sound            (D5)
w_es_columna                    índice de columna de cada w, o -1     (columnas mixtas)
var_names, var_lo, var_hi       variables ORIGINALES
prod_nombre/tipo/dato/args      subsistema no lineal (Proposición 1)  (D5)
restricciones                   ecuaciones cascarón
meta_*                          conteos para auditoría
```

Sobre D4: la salida por defecto **no reescribe el problema**. El paper es
explícito — *"we do not modify the original system"*, *"we do not need
concretization"* — y las `w_i`/`y_j` son etiquetas de nodos del DAG, no
variables nuevas. En `instprofe` la versión vieja emitía **4339 variables
nuevas** para un problema de 6 variables; para un B&B por intervalos eso
multiplica el espacio de búsqueda, y si el solver bisecta sobre las auxiliares
el resultado sale peor que el `REF` de la Tabla 1.

El `.txt` reescrito sigue disponible con `exportar_texto=True` (o `--txt`), pero
ahora marca las auxiliares como no bisectables y sus dominios son sound
(`-inf/inf` cuando no se pueden acotar), nunca el `[-1e8,1e8]` de antes.

---

## Paso 7 — `validator.py` (reescritura)

**Qué arregla:** B5. Este es el archivo que más cambia de *propósito*.

El validator viejo daba confianza falsa por tres razones: usaba el mismo parser
que el pipeline, solo comprobaba equivalencia algebraica de la reescritura (o
sea, pasaba **igual de bien con y sin los bugs D1/D2**), y `sp.simplify` se
colgaba en `instprofe` — por eso `batch_processor.py:32` ya lo tenía comentado.

**Qué borras:** `sp.simplify`, `sp.cancel`, el bucle de `subs` hasta punto fijo,
y el `import parse_instance` como fuente de verdad.

**Qué colocas:** siete checks, ninguno dependiente del parser del pipeline.

| check | qué comprueba |
|---|---|
| **V1** | líneas de restricción del `.txt` == restricciones parseadas (relee el archivo por su cuenta) |
| **V2** | cada fila de `A` reproduce **exactamente** (en `Fraction`) el nodo suma asociado |
| **V3** | equivalencia numérica contra el texto **original**, evaluado con `eval()` de Python en puntos aleatorios — camino totalmente independiente |
| **V4a** | los dominios contienen el valor real de cada nodo (Monte Carlo sobre las `x`) |
| **V4b** | isotonía por inclusión: `X' ⊂ X ⟹ F(X') ⊂ F(X)` sobre sub-cajas aleatorias |
| **V4c** | soundness **local** de cada operador: se muestrea directo en los dominios de los hijos |
| **V5** | cero variables auxiliares nuevas |

V4c es el que de verdad importa y vale la pena que entiendas por qué. Con
muestreo desde las `x` (V4a) el bug B1 **no se detecta**: para que `1/w` se
salga de `[-0.5, 0.5]` hace falta `|w| < 2`, y ningún punto aleatorio del
espacio de las `x` produce eso. V4c muestrea `w` **directamente** en su propio
dominio declarado, sesgado hacia magnitudes chicas, y lo pilla al primer
intento.

La prueba de que el validator sirve es la prueba negativa: reintroduce B1 a
propósito y mira si muerde.

```python
import intervals as iv
from exportator import construir
from validator import validar
def buggy(I):
    a, b = I
    return (-0.5, 0.5) if a <= 0 <= b else (min(1/b,1/a), max(1/b,1/a))
iv.reciproco = buggy
validar(*construir('instances/inst001.txt'))
```

Debe salir:

```
[FAIL] V4c soundness local de cada operador  violacion: ('pow', '(-5 + -5*(x7) + ...',
       [2.789651503627017e-06], 358467.71494569536, (-0.5, 0.5))
```

Si no falla, tu validator sigue siendo decorativo.

---

## Paso 8 — `batch_processor.py` (retoque)

**Qué borras:** la línea comentada

```python
#validar_sistema(ruta_instancia, ruta_npz, verbose=False)
```

**Qué colocas:** la validación activada por defecto (ya no se cuelga), más el
resumen del DAG y del sistema por instancia — que es lo que te deja ver de un
vistazo si la matriz es la del paper o la de incidencia:

```python
print(f"    DAG     : {inst.dag.resumen()}")
print(f"    sistema : {sis.resumen()}")
```

Y `raise` si la auditoría falla, para que el resumen final no diga "3 exitosos"
cuando uno está mal.

---

## Paso 9 — `main.py` (retoque)

**Qué borras:** la función `asd()` completa. Era código muerto **y** estaba rota
(B4): hacía `vars_dict, constr_list = parse_instance(...)` cuando la función
devolvía 3 valores, y `A, y, b, shell = extract_linear_system(...)` cuando
devolvía 5. Ejecutarla lanzaba `ValueError: too many values to unpack`.

**Qué colocas:** un CLI mínimo.

```bash
python main.py instances/inst001.txt        # una instancia
python main.py instances/inst001.txt --txt  # además el .txt reescrito
python main.py --lote                       # toda la carpeta instances/
```

Deja el docstring con el contenido real del `.npz`: el que estaba decía
`datos['b'] = vector de resultados`, y ahora `b` no existe (es `c`, las
constantes) y el vector de raíces es `w`.

---

## Paso 10 — `test_paper.py` (archivo NUEVO)

**Qué colocas:** cinco aserciones de regresión sobre el ejemplo (1). Están
escritas para no depender de los nombres, solo de las propiedades que el paper
exige, así que puedes pegar las ecuaciones literales del paper en `EJEMPLO_1`
sin tocar nada más:

```python
def test_coeficientes_en_la_matriz():   # D1: 4, 2 y -6 están DENTRO de A
def test_sin_variables_auxiliares():    # D4: cero v_i
def test_una_columna_por_nodo():        # D2: sin columnas duplicadas
def test_cse_de_la_suma_compartida():   # D3: (x1+x3) sigue siendo un nodo
def test_modo_distribuido_...():        # el flag reproduce la forma -6|-6
```

Salida esperada:

```
A = 5 filas x 8 columnas
  w_0 = 0 + {'(x3)^(2)': '1', 'x2': '4', '(x2)*(x1)': '1'}
  w_1 = 0 + {'x3': '1', 'x1': '1'}
  w_2 = 0 + {'(x3)^(2)': '1', '(x2)*(x1)': '1', '(x2)*(sin(w_1))': '2'}
  w_3 = 0 + {'x3': '6', 'x1': '6'}
  w_4 = 0 + {'(x3)^(2)': '1', 'x2': '4', '(x2)*(x1)': '1', 'w_1': '-6', 'sin(w_3)': '1'}
```

Cinco filas, ocho columnas, los coeficientes `4`, `2`, `6` y `-6` dentro de la
matriz, y `w_1 = x1+x3` apareciendo **a la vez como fila y como columna**: la
columna mixta de la Fig. 3.

---

## Verificación final

```bash
python3 test_paper.py       # las 5 aserciones del ejemplo (1)
python3 main.py --lote      # las 3 instancias + auditoría completa
```

Esperado:

```
>>> inst001     A 10x12, coef -5..5, 0 auxiliares, 7/7 PASS
>>> inst002     A 10x12, coef -5..5, 0 auxiliares, 7/7 PASS
>>> instprofe   A 3x1513, densidad 0.967, 0 auxiliares, 0.40 s, 7/7 PASS
 RESUMEN: 3 exitosos | 0 fallidos
```

---

## Mapa: hallazgo → archivo → prueba que lo cubre

| | dónde se arregla | prueba |
|---|---|---|
| D1 coeficientes absorbidos por `v_i` | `dag_builder` + `matrix_extractor` | `test_paper.py::test_coeficientes_en_la_matriz`, V2 |
| D2 mismo nodo en varias columnas | `dag_builder.escalar` / `_potencia` | `test_paper.py::test_una_columna_por_nodo` |
| D3 SymPy destruye la CSE | `dag_builder` (sin SymPy) | `test_paper.py::test_cse_de_la_suma_compartida` |
| D4 se modificaba el sistema original | `exportator` (.npz anexo) | `test_paper.py::test_sin_variables_auxiliares`, V5 |
| D5 el `.npz` no era consumible | `exportator.exportar_npz` | inspección de `np.load` |
| D6 `ref_count` ≠ in-degree | `dag_builder.in_degree` | `resumen()['compartidos_CSE']` |
| B1 potencias/divisiones negativas | `intervals.reciproco`, `potencia_entera` | V4c (+ prueba negativa) |
| B2 dominios truncados a ±1e8 | `dominios` (default `(-inf,inf)`) | V4a, V4b |
| B3 regex pierde el signo | `parser_instancias._RE_DOMINIO` | `parse_dominio('[- 2, 2]')` |
| B4 `asd()` rota | `main.py` (eliminada) | — |
| B5 validator decorativo | `validator` (7 checks) | prueba negativa de B1 |
| B6 redondeo | `intervals.from_fraction`, `A_exacta` | V2 (exacto en `Fraction`) |
| B7 rendimiento | `dag_builder` (hash-consing + iterativo) | 33 s → 0.40 s |

---

## Lo que queda fuera (y hay que decirlo en el informe)

Nada de la **Sección 3.3** está implementado: ni HC4-revise, ni el conditioner
`P` por Gauss–Jordan con pivote de máxima magnitud, ni la proyección
`P·A·y = P·b`, ni el simplex, ni la iteración a punto fijo. Hoy el repo es
**solo el preprocesamiento** — pero ahora es el preprocesamiento correcto, y el
`.npz` sí contiene todo lo que hace falta para escribir NSC encima.

El paso siguiente natural, si decides hacerlo:

1. `nsc.py` con `conditioner(A)` → Gauss–Jordan eligiendo como pivote el de
   mayor magnitud, sobre `A_exacta` (racionales) para no arrastrar error.
2. `filtrar(P·A, P·c, dominios)` → HC4-revise fila a fila.
3. Bucle a punto fijo sobre los dos subsistemas (lineal + productos).
4. Comparar contra la Tabla 1 del paper.

Solo ahí tiene sentido medir. Con la matriz de incidencia de antes no había
nada que medir: sobre ceros y unos, Gauss–Jordan no produce ningún conditioner
útil.
