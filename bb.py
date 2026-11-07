"""
Capa 9 (NUEVA): branch & bound por intervalos sobre el contractor NSC.

El paper mide NSC DENTRO de un solver (Ibex): cuenta cajas del arbol de
busqueda y tiempo (Tabla 1).  Contraer una sola caja no basta para ver la
ganancia, asi que este modulo arma el solver minimo:

  pila <- [caja inicial]
  mientras haya cajas:
      contraer la caja con REF o NSC (nsc.contraer)
      vacia                      -> descartar
      ancho de toda var < eps    -> caja solucion
      si no                      -> bisectar y apilar las dos mitades

Decisiones (el paper usa Ibex; aqui hay que fijarlas explicitamente):
  S1  Busqueda en profundidad (pila).  No cambia el numero de cajas cuando se
      exploran todas las soluciones, y usa poca memoria.
  S2  Bisector: SmearSumRel (el del paper, Newton.smear_sum_rel) cuando hay
      jacobiano, si no la variable de mayor ancho; siempre por el punto
      medio y siempre igual para REF y NSC.  bisector='ancho' fuerza la
      variable mas ancha.  Un dominio no acotado se corta en 0 o en un punto
      finito, para no producir mitades (-inf, -inf).
  S3  Una caja es solucion cuando todas sus variables miden menos que eps
      (absoluto).  Por si sola, la caja solo garantiza que no se perdio
      ninguna solucion; la existencia la certifica Newton cuando aplica (S6).
  S4  Cada caja se contrae desde los dominios de sus VARIABLES: los nodos
      intermedios se recalculan por la fase forward, como en Ibex.
  S5  Limites de cajas y de tiempo: si se alcanzan, el resultado lo dice
      (completo=False) en vez de reportar numeros de un arbol truncado como
      si fueran finales.
  S6  Newton intervalar (newton.py) despues de cada contraccion, para REF y
      NSC por igual (Ibex lo usa en todas las estrategias del paper).  Se
      alterna contraer -> Newton hasta que Newton deje de mejorar.  Cada caja
      solucion se intenta certificar con epsilon-inflacion
      (Newton.certificar): certificada = una caja apenas mas grande contiene
      exactamente una solucion.
  S8  Tres estrategias, las de la Sec. 3 del paper: ARB = HC4 sobre los
      arboles de expresion (sin preprocesamiento: inst construido con
      compartir=False), REF = HC4 sobre el DAG, NSC = DAG + subsistema lineal.
      Las tres usan el mismo Newton, bisector y ratio.
  S7  ratio del punto fijo REF/NSC dentro del arbol = 0.1 (en nsc.contraer el
      default es 1e-3).  Medido en el ejemplo (1): mismas cajas (+-10%) y
      casi la mitad del tiempo por caja.  Es el mismo para REF y NSC.
"""

import math
import time

from newton import Newton
from nsc import Preparado, conditioner_de, contraer


class ResultadoBB:
    def __init__(self, estrategia):
        self.estrategia = estrategia
        self.cajas = 0             # nodos del arbol (cada contraccion cuenta 1)
        self.soluciones = []       # cajas solucion: {var: (lo, hi)}
        self.certificadas = []     # bool por caja solucion (Newton, S6)
        self.descartadas = 0
        self.tiempo = 0.0
        self.completo = True
        self.motivo = ''

    def resumen(self):
        estado = 'completo' if self.completo else f'TRUNCADO ({self.motivo})'
        return (f"{self.estrategia:4s}: {self.cajas:7d} cajas, "
                f"{len(self.soluciones):4d} cajas solucion "
                f"({sum(self.certificadas)} certificadas), "
                f"{self.tiempo:8.2f}s  [{estado}]")


def _bisectar(caja, var_names, v=None):
    """S2: corta la variable v (por defecto la mas ancha) por el punto medio
    (o por un punto finito si el dominio no es acotado)."""
    if v is None or not caja[v][1] > caja[v][0]:
        v = max(var_names, key=lambda n: caja[n][1] - caja[n][0])
    lo, hi = caja[v]
    if math.isinf(lo) and math.isinf(hi):
        m = 0.0
    elif math.isinf(lo):
        m = min(hi - 1.0, -1.0) if hi > 0 else hi * 2.0 - 1.0
    elif math.isinf(hi):
        m = max(lo + 1.0, 1.0) if lo < 0 else lo * 2.0 + 1.0
    else:
        m = lo / 2.0 + hi / 2.0            # evita overflow con dominios enormes
    if not lo < m < hi:                     # caja de 1-2 floats: no se puede cortar
        return None
    izq, der = dict(caja), dict(caja)
    izq[v] = (lo, m)
    der[v] = (m, hi)
    return izq, der


def _mejora(antes, despues, ratio):
    return any((a[1] - a[0]) - (d[1] - d[0]) > ratio * (a[1] - a[0])
               for a, d in zip(antes.values(), despues.values()))


def resolver(inst, sis, usar_nsc=True, eps=1e-8, ratio=0.1, usar_newton=True,
             bisector='smear', max_cajas=100000, max_segundos=600.0,
             restricciones=None, nombre=None):
    """Branch & bound completo.  restricciones permite reemplazar los lados
    derechos (lo usa la prueba de soundness con soluciones fabricadas).
    nombre etiqueta el resultado (por ejemplo 'ARB' para un inst construido
    sin CSE); por defecto es 'NSC' o 'REF'."""
    res = ResultadoBB(nombre or ('NSC' if usar_nsc else 'REF'))
    prep = Preparado(inst, sis, conditioner_de(sis) if usar_nsc else None)
    if restricciones is not None:
        prep.restricciones = restricciones
    ids = {n.data: n.id for n in inst.dag.nodes if n.kind == 'var'}
    var_names = [v for v in inst.var_names if v in ids]
    newton = Newton(inst, prep.restricciones)
    usar_newton = usar_newton and newton.aplicable()

    t0 = time.time()
    pila = [dict(inst.var_domains)]
    while pila:
        if res.cajas >= max_cajas:
            res.completo, res.motivo = False, f'limite de {max_cajas} cajas'
            break
        if time.time() - t0 > max_segundos:
            res.completo, res.motivo = False, f'limite de {max_segundos:.0f}s'
            break
        caja = pila.pop()
        res.cajas += 1
        dom, _ = contraer(prep, caja, usar_lineal=usar_nsc, ratio=ratio)
        if dom is None:
            res.descartadas += 1
            continue
        caja = {v: dom[ids[v]] for v in var_names}
        while usar_newton:
            nueva, _ = newton.paso(caja)
            if nueva is None or not _mejora(caja, nueva, 0.01):
                break
            dom, _ = contraer(prep, nueva, usar_lineal=usar_nsc, ratio=ratio)
            if dom is None:
                nueva = None
                break
            caja = {v: dom[ids[v]] for v in var_names}
        if usar_newton and nueva is None:
            res.descartadas += 1
            continue
        if all(caja[v][1] - caja[v][0] < eps for v in var_names):
            res.soluciones.append(caja)
            res.certificadas.append(usar_newton and newton.certificar(caja))
            continue
        v = newton.smear_sum_rel(caja) if usar_newton and bisector == 'smear' else None
        mitades = _bisectar(caja, var_names, v)
        if mitades is None:
            res.soluciones.append(caja)
            res.certificadas.append(usar_newton and newton.certificar(caja))
            continue
        pila.extend(reversed(mitades))
    res.tiempo = time.time() - t0
    return res


def solucion_conocida(ruta):
    """Lee el punto del comentario '//One known solution...' (ncsp_generator).

    Devuelve una lista de floats en el orden de declaracion, o None."""
    with open(ruta) as f:
        lineas = [l.strip() for l in f]
    for k, l in enumerate(lineas[:-1]):
        if 'known solution' in l.lower():
            try:
                return [float(t) for t in lineas[k + 1].lstrip('/').split(',') if t.strip()]
            except ValueError:
                return None
    return None


def distancia_a_soluciones(punto, soluciones):
    """Distancia (norma infinito) del punto a la caja solucion mas cercana."""
    mejor = math.inf
    for caja in soluciones:
        d = max(max(lo - x, 0.0, x - hi) for x, (lo, hi) in zip(punto, caja.values()))
        mejor = min(mejor, d)
    return mejor


def _ganancia(base, otra):
    if base.completo and otra.completo and otra.cajas and otra.tiempo > 0:
        return f"gt = {base.tiempo / otra.tiempo:.2f}  gb = {base.cajas / otra.cajas:.2f}"
    return "no comparable (algun arbol truncado)"


def comparar_bb(inst, sis, ruta=None, **kw):
    """Resuelve con ARB (si hay ruta), REF y NSC e imprime la comparacion
    (gt y gb del paper).  ARB necesita la ruta para reconstruir sin CSE."""
    resultados = []
    if ruta is not None:
        from exportator import construir
        inst_arb, sis_arb, _ = construir(ruta, compartir=False)
        resultados.append(resolver(inst_arb, sis_arb, usar_nsc=False, nombre='ARB', **kw))
    resultados += [resolver(inst, sis, usar_nsc=False, **kw),
                   resolver(inst, sis, usar_nsc=True, **kw)]
    for r in resultados:
        print(f"  {r.resumen()}")
    por = {r.estrategia: r for r in resultados}
    print(f"  ganancia NSC vs REF: {_ganancia(por['REF'], por['NSC'])}")
    if 'ARB' in por:
        print(f"  ganancia REF vs ARB: {_ganancia(por['ARB'], por['REF'])}")
        print(f"  ganancia NSC vs ARB: {_ganancia(por['ARB'], por['NSC'])}")
    punto = solucion_conocida(ruta) if ruta else None
    if punto is not None and len(punto) == len(inst.var_names):
        for r in resultados:
            if not r.completo:
                print(f"  {r.estrategia}: arbol truncado, no se puede comprobar la solucion conocida")
                continue
            d = distancia_a_soluciones(punto, r.soluciones)
            print(f"  {r.estrategia}: solucion conocida a {d:.1e} de la caja solucion mas cercana")
    return resultados
