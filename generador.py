"""
Capa 11 (NUEVA): generador de instancias aleatorias como las de la Sec. 4.2.

Replica la estructura de github.com/vareyesr/ncsp_generator (el generador del
paper; inst001/inst002 salieron de ahi) con un cambio: los dominios de las
variables son un parametro.  Con los [-1e8, 1e8] del paper el branch & bound
en Python no termina (ver CLAUDE.md), asi que el experimento usa dominios mas
chicos.

  tipo 'mult'  ec. (6):  prod_k ( a_jk + sum_i b_ijk t_ijk )^d_jk  = R_j
  tipo 'func'  ec. (7):  sum_k  u( a_jk + sum_i b_ijk t_ijk )       = R_j
               con u en {log, exp, potencia -1, identidad}

  - n variables y n restricciones (sistema cuadrado, como en el paper)
  - t_ijk sale de un pool: las n variables + (pool - n) potencias al
    cuadrado o al cubo de variables al azar
  - a_jk en [-5, 5], b_ijk en [-5, 5] sin 0, d_jk en {-1, 1}, enteros
  - sum_j p_j = P (sumas en el sistema), sum q_jk = Q (terminos)
  - R = valor en un punto aleatorio de [-1, 1]^n, que queda anotado como
    solucion conocida (mismo formato que inst001)

Decision G1: el paper no dice como reparte P y Q; aqui cada restriccion
recibe al menos una suma y cada suma al menos un termino, y el resto se
reparte uniformemente al azar.  Decision G2: si en el punto la expresion no
esta definida (log de algo <= 0, division por 0, overflow), se sortea de nuevo
la restriccion.
"""

import math
import os
import random


def _repartir(total, partes, rng):
    """Composicion aleatoria de total en partes, cada una >= 1 (G1)."""
    cuenta = [1] * partes
    for _ in range(total - partes):
        cuenta[rng.randrange(partes)] += 1
    return cuenta


def _pool(n, tam, rng):
    terminos = [(f"x{i}", 1) for i in range(n)]
    usados = set()
    while len(terminos) < tam:
        t = (f"x{rng.randrange(n)}", rng.choice((2, 3)))
        if t not in usados:
            usados.add(t)
            terminos.append(t)
    return terminos


def _texto_termino(var, e):
    return var if e == 1 else f"{var}^{e}"


def _suma(q, pool, rng):
    """(texto, funcion de evaluacion) de  a + sum b*t."""
    a = rng.randint(-5, 5)
    partes = [(rng.choice([b for b in range(-5, 6) if b != 0]), rng.choice(pool))
              for _ in range(q)]
    txt = str(a) + "".join(f"{'+' if b > 0 else '-'}{abs(b)}*{_texto_termino(v, e)}"
                          for b, (v, e) in partes)

    def valor(x):
        return a + sum(b * x[v] ** e for b, (v, e) in partes)
    return f"({txt})", valor


def _restriccion(tipo, qs, pool, rng, x):
    sumas = [_suma(q, pool, rng) for q in qs]
    if tipo == 'mult':
        piezas, val = [], 1.0
        for txt, f in sumas:
            d = rng.choice((-1, 1))
            piezas.append(txt if d == 1 else f"{txt}^(-1)")
            val *= f(x) ** d
        return "*".join(piezas), val
    piezas, val = [], 0.0
    for txt, f in sumas:
        u = rng.choice(('log', 'exp', 'inv', 'id'))
        s = f(x)
        if u == 'log':
            piezas.append(f"log{txt}")
            val += math.log(s)
        elif u == 'exp':
            piezas.append(f"exp{txt}")
            val += math.exp(s)
        elif u == 'inv':
            piezas.append(f"{txt}^(-1)")
            val += 1.0 / s
        else:
            piezas.append(txt)
            val += s
    return "+".join(piezas), val


def generar(tipo, n, P, Q, dominio, semilla, tam_pool=None):
    """Devuelve el texto de una instancia."""
    rng = random.Random(semilla)
    pool = _pool(n, tam_pool or n + n // 2, rng)
    x = {f"x{i}": rng.uniform(-1.0, 1.0) for i in range(n)}
    p = _repartir(P, n, rng)
    q = _repartir(Q, P, rng)
    lineas, k = [], 0
    for j in range(n):
        qs = q[k:k + p[j]]
        k += p[j]
        while True:                                    # G2
            try:
                expr, R = _restriccion(tipo, qs, pool, rng, x)
            except (ValueError, ZeroDivisionError, OverflowError):
                continue
            if math.isfinite(R) and abs(R) < 1e12:
                break
        lineas.append(f"{expr} = {R!r};")
    lo, hi = dominio
    return ("//One known solution for this problem is:\n//"
            + ",".join(repr(x[f"x{i}"]) for i in range(n)) + ",\n"
            + "Variables\n\n"
            + "".join(f"x{i} in [{lo},{hi}];\n" for i in range(n))
            + "\nConstraints\n\n" + "\n".join(lineas) + "\nend\n")


def escribir(carpeta, tipo, n, P, Q, dominio, semillas):
    """Escribe una instancia por semilla y devuelve las rutas."""
    os.makedirs(carpeta, exist_ok=True)
    rutas = []
    for s in semillas:
        ruta = os.path.join(carpeta, f"{tipo}_n{n}_P{P}_Q{Q}_s{s}.txt")
        with open(ruta, 'w', newline='\n') as f:
            f.write(generar(tipo, n, P, Q, dominio, s))
        rutas.append(ruta)
    return rutas
