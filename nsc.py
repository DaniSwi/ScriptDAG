"""
Capa 8 (NUEVA): el contractor NSC (Seccion 3.3 del paper).

Hasta aqui el repo solo hacia el preprocesamiento (Secciones 3.1-3.2).  Este
modulo implementa el filtrado, tal como lo describe la Seccion 3 del paper:

  1. HC4-revise sobre el DAG: una fase forward y una backward sobre TODOS los
     nodos, con dominios persistentes en los nodos intermedios (por eso los
     nodos compartidos conservan la informacion de proyeccion).
  2. Contraccion del subsistema lineal  A*y = w - c  por proyeccion sobre
     PA*y = P*(w - c), con P el conditioner de Gauss-Jordan (Seccion 3.3.2).
     Los dominios mejorados son los de los nodos del DAG, asi que "mapearlos
     de vuelta al DAG" es automatico.
  3. Repetir 1-2 hasta punto fijo.

Decisiones (documentadas porque el paper no las fija):
  N1  Pivote de Gauss-Jordan: por columna, en orden, el de mayor magnitud
      (lo que dice el paper).  Empates: la fila con menos no-ceros (menos
      relleno).  Con esa regla se reproduce EXACTAMENTE la P de la Sec. 3.3.2
      (ver test_nsc.py).  Los pivotes NO se normalizan a 1: PA conserva el 2 y
      el 6 que muestra el paper.
  N2  P y PA se calculan en Fraction sobre A_exacta (invariante 5), y solo se
      pasan a intervalos al proyectar, con redondeo hacia afuera.
  N3  Las filas de PA que quedan en cero son relaciones lineales entre las w
      (la fila 5 del ejemplo: b5 - b4/6 = 0).  El paper solo proyecta sobre y,
      asi que se guardan en Conditioner.dependientes pero no se usan.
  N4  Se proyecta solo sobre y (el paper: "projecting the domains of b over
      y").  Una w que tambien es columna (columna mixta) se contrae igual,
      porque es el mismo nodo del DAG.
  N5  Contraccion SOUND: toda operacion es aritmetica de intervalos con
      redondeo hacia afuera; una division por un intervalo que contiene 0 no
      contrae.  Una interseccion vacia devuelve None (la caja no tiene
      soluciones).
  N6  Cierra la grieta del rhs dividido en float: si la raiz de una
      restriccion tenia coeficiente != 1, el parser guardo rhs/coef
      redondeado al mas cercano; aqui ese valor se ensancha 2 ulps.
  N7  No se implementa la variante con simplex (Sec. 3.3.1): requiere un LP
      solver con cotas certificadas y el repo solo depende de numpy.
"""

import math
import time
from fractions import Fraction

import intervals as iv
from dominios import evaluar_nodo, propagar_dominios


# ===================================================================
#  Conditioner por Gauss-Jordan (Seccion 3.3.2)
# ===================================================================

class Conditioner:
    """PA = P * A en aritmetica exacta.

    filas        indices de las filas de PA que se usan para proyectar
                 (las de pivote, en orden de pivote)
    PA[r]        {columna j: Fraction}
    P[r]         {fila original i de A: Fraction}
    pivotes      [(r, j)]  fila de PA y columna donde quedo el pivote
    dependientes filas r con PA[r] == 0 (relaciones entre las w, N3)
    """

    def __init__(self, P, PA, pivotes, dependientes):
        self.P = P
        self.PA = PA
        self.pivotes = pivotes
        self.dependientes = dependientes
        self.filas = [r for r, _ in pivotes]

    def rango(self):
        return len(self.pivotes)


def _restar(dst, src, f):
    """dst <- dst - f*src, sobre filas dispersas {indice: Fraction}."""
    for k, v in src.items():
        nv = dst.get(k, 0) - f * v
        if nv:
            dst[k] = nv
        else:
            dst.pop(k, None)


def conditioner(A_exacta, n_columnas):
    """Gauss-Jordan con pivote de maxima magnitud por columna (N1, N2)."""
    m = len(A_exacta)
    PA = [dict(f) for f in A_exacta]
    P = [{i: Fraction(1)} for i in range(m)]
    libres = set(range(m))
    pivotes = []
    for j in range(n_columnas):
        if not libres:
            break
        candidatas = [i for i in libres if PA[i].get(j)]
        if not candidatas:
            continue
        p = max(candidatas, key=lambda i: (abs(PA[i][j]), -len(PA[i]), -i))
        libres.discard(p)
        pivotes.append((p, j))
        piv = PA[p][j]
        for i in range(m):
            if i != p and PA[i].get(j):
                f = PA[i][j] / piv
                _restar(PA[i], PA[p], f)
                _restar(P[i], P[p], f)
    dependientes = sorted(libres)
    assert all(not PA[i] for i in dependientes), "fila libre no nula tras Gauss-Jordan"
    return Conditioner(P, PA, pivotes, dependientes)


def conditioner_de(sis):
    return conditioner(sis.A_exacta, len(sis.y_nodes))


# ===================================================================
#  Proyecciones (fase backward de HC4-revise)
# ===================================================================

def _division_proyeccion(num, den):
    """Envoltura sound de {x : x*d = n, n en num, d en den}."""
    if den[0] <= 0.0 <= den[1]:
        if num[0] <= 0.0 <= num[1]:
            return iv.ENTERO                  # 0*x = 0 para todo x
        if den[0] == 0.0 and den[1] == 0.0:
            return None                       # x*0 = n != 0: imposible
    return iv.mult(num, iv.reciproco(den))


def _dn2(x):
    return iv.dn(iv.dn(x))


def _up2(x):
    return iv.up(iv.up(x))


def _raiz(x, n, hacia_arriba):
    """x^(1/n) para x >= 0, redondeado hacia el lado pedido.

    math.sqrt es correctamente redondeada; x**(1/n) no lo es (1/n tampoco es
    exacto), asi que se ensancha con un margen relativo holgado."""
    if x == iv.INF:
        return iv.INF
    if x <= 0.0:
        return 0.0
    if n == 2:
        r = math.sqrt(x)
        return iv.up(r) if hacia_arriba else iv.dn(r)
    r = x ** (1.0 / n)
    return iv.up(r * (1 + 1e-12)) if hacia_arriba else iv.dn(r * (1 - 1e-12))


def _proyectar_raiz_entera(base, w, n):
    """base inter {b : b^n en w}, n entero > 0."""
    if n % 2 == 1:
        lo = -_raiz(-w[0], n, True) if w[0] < 0 else _raiz(w[0], n, False)
        hi = -_raiz(-w[1], n, False) if w[1] < 0 else _raiz(w[1], n, True)
        return iv.interseccion(base, (lo, hi))
    t = iv.interseccion(w, (0.0, iv.INF))
    if t is None:
        return None
    r = (_raiz(t[0], n, False), _raiz(t[1], n, True))
    pos = iv.interseccion(base, r)
    neg = iv.interseccion(base, (-r[1], -r[0]))
    if pos is None:
        return neg
    if neg is None:
        return pos
    return iv.hull(pos, neg)


def _proyectar_pow(nodo, dom):
    e = Fraction(nodo.data)
    h = nodo.args[0].id
    w = dom[nodo.id]
    if e.denominator != 1:
        if e != Fraction(1, 2):
            return True                       # no se sabe invertir: no contrae
        t = iv.interseccion(w, (0.0, iv.INF))
        if t is None:
            return False
        nuevo = iv.interseccion(dom[h], (iv.dn(t[0] * t[0]), iv.up(t[1] * t[1])))
    else:
        n = int(e)
        if n == 0:
            return True
        if n < 0:
            w = iv.reciproco(w)               # b^|n| = 1/w
            n = -n
        nuevo = _proyectar_raiz_entera(dom[h], w, n)
    if nuevo is None:
        return False
    dom[h] = nuevo
    return True


def _proyectar_func(nodo, dom):
    if len(nodo.args) != 1:
        return True
    h = nodo.args[0].id
    w, a = dom[nodo.id], dom[h]
    nombre = nodo.data
    if nombre == 'exp':
        if w[1] <= 0.0:
            return False
        lo = _dn2(math.log(w[0])) if w[0] > 0.0 else iv.NEG_INF
        hi = _up2(math.log(w[1])) if math.isfinite(w[1]) else iv.INF
        nuevo = iv.interseccion(a, (lo, hi))
    elif nombre == 'log':
        lo = _dn2(math.exp(w[0])) if w[0] > -745 else 0.0
        hi = _up2(math.exp(w[1])) if w[1] < 709 else iv.INF
        nuevo = iv.interseccion(a, (max(lo, 0.0), hi))
    elif nombre == 'sqrt':
        t = iv.interseccion(w, (0.0, iv.INF))
        if t is None:
            return False
        nuevo = iv.interseccion(a, (iv.dn(t[0] * t[0]), iv.up(t[1] * t[1])))
    elif nombre == 'Abs':
        if w[1] < 0.0:
            return False
        nuevo = iv.interseccion(a, (-w[1], w[1]))
    else:
        return True                           # sin, cos, ...: no contrae
    if nuevo is None:
        return False
    dom[h] = nuevo
    return True


def _prefijos(terminos, operar, neutro):
    pre = [neutro]
    for t in terminos:
        pre.append(operar(pre[-1], t))
    suf = [neutro]
    for t in reversed(terminos):
        suf.append(operar(t, suf[-1]))
    suf.reverse()
    return pre, suf


def _proyectar_suma(rhs, items, dom):
    """Proyecta  sum_k a_k * y_k  in rhs  sobre cada y_k.

    items = [(node_id, a_k intervalo, 1/a_k intervalo)].  Se usan sumas
    prefijo/sufijo en vez de 'total - termino' para que un unico termino no
    acotado no impida contraer a los demas."""
    terminos = [iv.mult(a, dom[h]) for h, a, _ in items]
    pre, suf = _prefijos(terminos, iv.suma, (0.0, 0.0))
    for k, (h, _, inv) in enumerate(items):
        resto = iv.suma(pre[k], suf[k + 1])
        nuevo = iv.interseccion(dom[h], iv.mult(iv.resta(rhs, resto), inv))
        if nuevo is None:
            return False
        dom[h] = nuevo
    return True


# ===================================================================
#  HC4-revise sobre el DAG (paso 1 de NSC)
# ===================================================================

class Preparado:
    """Todo lo que no cambia entre iteraciones, convertido una sola vez a
    intervalos: coeficientes de las sumas, sus inversos, P, PA y c."""

    def __init__(self, inst, sis, cond):
        self.inst, self.sis, self.cond = inst, sis, cond
        self.sumas = {}
        for nodo in inst.dag.nodes:
            if nodo.kind == 'add':
                const, coefs = nodo.data
                self.sumas[nodo.id] = (
                    iv.from_fraction(const),
                    [(h.id, iv.from_fraction(c), iv.from_fraction(1 / c))
                     for c, h in zip(coefs, nodo.args)])
        self.restricciones = restricciones_de(inst)
        self.lineal = []
        if cond is not None:
            c_iv = [iv.from_fraction(q) for q in sis.c_exacta]
            for r in cond.filas:
                P_r = [(sis.w_nodes[i].id, iv.from_fraction(p), c_iv[i])
                       for i, p in sorted(cond.P[r].items())]
                PA_r = [(sis.y_nodes[j].id, iv.from_fraction(q), iv.from_fraction(1 / q))
                        for j, q in sorted(cond.PA[r].items())]
                self.lineal.append((P_r, PA_r))


def restricciones_de(inst):
    """[(node_id, intervalo)] que cada raiz debe cumplir (N6)."""
    out = []
    for raiz, coef, op, rhs in zip(inst.raices, inst.coef_raiz, inst.ops, inst.rhs):
        r = (rhs, rhs) if coef == 1 else (_dn2(rhs), _up2(rhs))
        if op in ('<=', '<'):
            r = (iv.NEG_INF, r[1])
        elif op in ('>=', '>'):
            r = (r[0], iv.INF)
        out.append((raiz.id, r))
    return out


def hc4_revise(prep, dom):
    """Una pasada forward + backward sobre todo el DAG.  False si vacia."""
    nodos = prep.inst.dag.nodes
    for nodo in nodos:
        if nodo.kind in ('var', 'const'):
            continue
        I = iv.interseccion(dom[nodo.id], evaluar_nodo(nodo, dom))
        if I is None:
            return False
        dom[nodo.id] = I
    for nid, r in prep.restricciones:
        I = iv.interseccion(dom[nid], r)
        if I is None:
            return False
        dom[nid] = I
    for nodo in reversed(nodos):
        k = nodo.kind
        if k == 'add':
            const, items = prep.sumas[nodo.id]
            ok = _proyectar_suma(iv.resta(dom[nodo.id], const), items, dom)
        elif k == 'mul':
            hijos = [a.id for a in nodo.args]
            pre, suf = _prefijos([dom[h] for h in hijos], iv.mult, (1.0, 1.0))
            ok = True
            for k2, h in enumerate(hijos):
                otros = iv.mult(pre[k2], suf[k2 + 1])
                proy = _division_proyeccion(dom[nodo.id], otros)
                nuevo = None if proy is None else iv.interseccion(dom[h], proy)
                if nuevo is None:
                    ok = False
                    break
                dom[h] = nuevo
        elif k == 'pow':
            ok = _proyectar_pow(nodo, dom)
        elif k == 'func':
            ok = _proyectar_func(nodo, dom)
        else:
            ok = True
        if not ok:
            return False
    return True


# ===================================================================
#  Proyeccion sobre PA*y = P*(w - c)  (paso 2 de NSC)
# ===================================================================

def proyectar_lineal(prep, dom):
    for P_r, PA_r in prep.lineal:
        rhs = (0.0, 0.0)
        for wid, p, c in P_r:
            rhs = iv.suma(rhs, iv.mult(p, iv.resta(dom[wid], c)))
        if not _proyectar_suma(rhs, PA_r, dom):
            return False
    return True


# ===================================================================
#  Punto fijo (paso 3)
# ===================================================================

def _ancho(I):
    return I[1] - I[0]


def _cambio_significativo(antes, dom, ratio):
    for nid, I in dom.items():
        a, b = _ancho(antes[nid]), _ancho(I)
        if b < a and (math.isinf(a) or a - b > ratio * a):
            return True
    return False


def contraer(prep, var_domains=None, usar_lineal=True, ratio=1e-3, max_iter=200):
    """NSC (usar_lineal=True) o HC4 sobre el DAG (usar_lineal=False) hasta
    punto fijo, partiendo de la caja var_domains.

    Devuelve (dom, iteraciones); dom es None si la caja no tiene soluciones."""
    inst = prep.inst
    dom = propagar_dominios(inst.dag, var_domains or inst.var_domains)
    for it in range(1, max_iter + 1):
        antes = dict(dom)
        if not hc4_revise(prep, dom):
            return None, it
        if usar_lineal and not proyectar_lineal(prep, dom):
            return None, it
        if not _cambio_significativo(antes, dom, ratio):
            return dom, it
    return dom, max_iter


# ===================================================================
#  Reporte REF (HC4 en el DAG) vs NSC
# ===================================================================

def _caja_vars(dom, inst):
    ids = {n.data: n.id for n in inst.dag.nodes if n.kind == 'var'}
    return {v: dom[ids[v]] for v in inst.var_names if v in ids}


def comparar(inst, sis, var_domains=None, verbose=True):
    """Corre REF y NSC sobre la caja inicial e imprime cuanto contrae cada uno."""
    cond = conditioner_de(sis)
    resultados = {}
    for nombre, lineal in (('REF (HC4 en el DAG)', False), ('NSC', True)):
        prep = Preparado(inst, sis, cond if lineal else None)
        t0 = time.time()
        dom, it = contraer(prep, var_domains, usar_lineal=lineal)
        resultados[nombre] = (dom, it, time.time() - t0)
    if verbose:
        ini = var_domains or inst.var_domains
        print(f"  conditioner: rango {cond.rango()}/{len(sis.w_nodes)} filas, "
              f"{len(cond.dependientes)} dependientes")
        for nombre, (dom, it, t) in resultados.items():
            if dom is None:
                print(f"  {nombre:20s}: caja VACIA (sin soluciones) en {it} iter, {t:.2f}s")
                continue
            caja = _caja_vars(dom, inst)
            contraidas = sum(1 for v, I in caja.items()
                             if _ancho(I) < _ancho(ini[v]))
            ganancia = sum(math.log10(_ancho(ini[v]) / max(_ancho(I), 1e-300))
                           for v, I in caja.items()
                           if math.isfinite(_ancho(ini[v])) and _ancho(ini[v]) > 0)
            print(f"  {nombre:20s}: {contraidas}/{len(caja)} variables contraidas, "
                  f"reduccion de volumen 10^{ganancia:.2f}, {it} iter, {t:.2f}s")
    return cond, resultados
