"""
Capa 6: auditoria del sistema extraido.

El validator anterior daba una confianza falsa (B5):
  - usaba EL MISMO parser que el pipeline, asi que una restriccion descartada
    en silencio era invisible;
  - solo comprobaba equivalencia algebraica de la reescritura, asi que pasaba
    igual de bien CON y SIN los bugs D1/D2;
  - no verificaba nada de los dominios;
  - sp.simplify() se colgaba en instprofe.

Esta version comprueba cinco cosas, y ninguna depende del parser del pipeline:

  V1  conteo de restricciones: lineas del .txt == restricciones parseadas.
  V2  cada fila de A reproduce EXACTAMENTE (en Fraction) el nodo suma asociado.
  V3  equivalencia numerica contra el texto ORIGINAL, evaluado con eval() de
      Python en puntos aleatorios: un camino independiente del pipeline.
  V4  soundness de los dominios: Monte Carlo sobre las variables, se comprueba
      que el valor real de CADA nodo cae dentro del intervalo declarado.
  V5  no se introdujeron variables auxiliares nuevas (invariante del paper).
"""

import math
import random
import re
from fractions import Fraction

import intervals as iv
from dominios import (evaluar_en_punto, muestrear, propagar_dominios, sub_caja,
                      valor_local, muestra_en_intervalo)

_RE_OP = re.compile(r'(>=|<=|==|=|>|<)')

_NS_MATE = {k: getattr(math, k) for k in
            ('sin', 'cos', 'tan', 'exp', 'log', 'sqrt', 'atan', 'asin', 'acos',
             'sinh', 'cosh', 'tanh')}
_NS_MATE['abs'] = abs
_NS_MATE['Abs'] = abs
_NS_MATE['pi'] = math.pi


class Resultado:
    def __init__(self):
        self.checks = []
        self.ok = True

    def add(self, nombre, ok, detalle=""):
        self.checks.append((nombre, ok, detalle))
        self.ok = self.ok and ok

    def imprimir(self, titulo=""):
        print(f"--- AUDITORIA {titulo} ---")
        for nombre, ok, det in self.checks:
            print(f"  [{'PASS' if ok else 'FAIL'}] {nombre}" + (f"  {det}" if det else ""))
        print(f"  => {'TODO OK' if self.ok else 'HAY FALLOS'}")


def _lineas_restriccion_crudas(ruta):
    """Relee el archivo SIN usar el parser del pipeline (independencia, B5)."""
    fuera = []
    modo = None
    for linea in open(ruta):
        t = linea.strip()
        if not t or t.startswith('//') or t.startswith('#'):
            continue
        b = t.lower()
        if b.startswith('minimize') or b.startswith('maximize'):
            continue
        if b in ('variables', 'variables:'):
            modo = 'v'; continue
        if b in ('constraints', 'constraints:'):
            modo = 'c'; continue
        if b == 'end':
            break
        if modo == 'c':
            fuera.append(t.rstrip(';').strip())
    return fuera


def validar(inst, sis, dom, n_muestras=25, semilla=0, verbose=True):
    r = Resultado()
    rng = random.Random(semilla)
    crudas = _lineas_restriccion_crudas(inst.ruta)

    # ---------------------------------------------------------------- V1
    r.add("V1 conteo de restricciones",
          len(crudas) == inst.n_restricciones() == inst.lineas_restriccion,
          f"archivo={len(crudas)}  parseadas={inst.n_restricciones()}")

    # ---------------------------------------------------------------- V2
    malas = []
    for i, nodo in enumerate(sis.w_nodes):
        const, coefs = nodo.data
        esperado = {}
        for c, hijo in zip(coefs, nodo.args):
            j = sis.col_de_nodo[hijo.id]
            esperado[j] = esperado.get(j, Fraction(0)) + c
        esperado = {j: v for j, v in esperado.items() if v != 0}
        real = {j: v for j, v in sis.A_exacta[i].items() if v != 0}
        if esperado != real or Fraction(sis.c_exacta[i]) != Fraction(const):
            malas.append(i)
    r.add("V2 cada fila de A reproduce su nodo suma (exacto)", not malas,
          f"filas malas: {malas[:5]}" if malas else f"{len(sis.w_nodes)} filas")

    # ---------------------------------------------------------------- V3
    fallas_v3, saltadas = [], 0
    for k, cruda in enumerate(crudas):
        partes = _RE_OP.split(cruda)
        lhs_txt = partes[0].replace('^', '**')
        raiz, coef = inst.raices[k], inst.coef_raiz[k]
        for _ in range(max(3, n_muestras // 3)):
            punto = muestrear(inst.var_domains, rng, cota=10.0)
            ns = dict(_NS_MATE); ns.update(punto)
            try:
                ref = eval(lhs_txt, {"__builtins__": {}}, ns)
            except Exception:
                saltadas += 1
                continue
            val = evaluar_en_punto(inst.dag, punto).get(raiz.id)
            if val is None:
                saltadas += 1
                continue
            rec = float(coef) * val
            escala = max(1.0, abs(ref), abs(rec))
            if abs(ref - rec) > 1e-7 * escala:
                fallas_v3.append((k, ref, rec))
                break
    r.add("V3 equivalencia numerica con el texto original", not fallas_v3,
          f"fallas: {fallas_v3[:2]}" if fallas_v3 else
          f"{len(crudas)} restricciones, {saltadas} evaluaciones indefinidas omitidas")

    # ---------------------------------------------------------------- V4
    if dom is None:
        r.add("V4 soundness de dominios", True, "omitido (dominios no calculados)")
    else:
        violaciones = []
        for _ in range(n_muestras):
            punto = muestrear(inst.var_domains, rng, cota=1e3)
            vals = evaluar_en_punto(inst.dag, punto)
            for nodo in inst.dag.nodes:
                v = vals.get(nodo.id)
                if v is None or not math.isfinite(v):
                    continue
                if not iv.contiene(dom[nodo.id], v, tol=1e-6):
                    violaciones.append((nodo.id, nodo.kind, v, dom[nodo.id]))
            if violaciones:
                break
        r.add("V4a soundness de dominios (puntos)", not violaciones,
              f"violacion: {violaciones[0]}" if violaciones else
              f"{n_muestras} puntos x {len(inst.dag.nodes)} nodos")

        # V4b: isotonia por inclusion.  Es MUCHO mas sensible que V4a: detecta
        # B1 (1/[a,b] con 0 dentro) aunque ningun punto aleatorio lo toque.
        fallas = []
        for _ in range(max(5, n_muestras // 2)):
            sub = sub_caja(inst.var_domains, rng, cota=1e3)
            dsub = propagar_dominios(inst.dag, sub)
            for nodo in inst.dag.nodes:
                (slo, shi), (glo, ghi) = dsub[nodo.id], dom[nodo.id]
                margen = 1e-6 * max(1.0, abs(glo) if math.isfinite(glo) else 0.0,
                                    abs(ghi) if math.isfinite(ghi) else 0.0)
                if slo < glo - margen or shi > ghi + margen:
                    fallas.append((nodo.id, nodo.kind, nodo.to_str()[:40],
                                   (slo, shi), (glo, ghi)))
                    break
            if fallas:
                break
        r.add("V4b isotonia por inclusion (sub-cajas)", not fallas,
              f"violacion: {fallas[0]}" if fallas else
              f"{max(5, n_muestras // 2)} sub-cajas")

        # V4c: soundness LOCAL de cada operador.  Se muestrea directamente en
        # los dominios de los hijos, asi que si la regla de intervalos de una
        # operacion esta mal (B1: 1/[a,b] con 0 dentro), salta aqui aunque
        # ningun punto del espacio de las x llegue a esa region.
        locales = []
        for nodo in inst.dag.nodes:
            if not nodo.args:
                continue
            for _ in range(40):
                vals = [muestra_en_intervalo(dom[h.id], rng) for h in nodo.args]
                if any(v is None for v in vals):
                    continue
                real = valor_local(nodo, vals)
                if real is None or not math.isfinite(real):
                    continue
                if not iv.contiene(dom[nodo.id], real, tol=1e-6):
                    locales.append((nodo.kind, nodo.to_str()[:45], vals[:2],
                                    real, dom[nodo.id]))
                    break
            if locales:
                break
        r.add("V4c soundness local de cada operador", not locales,
              f"violacion: {locales[0]}" if locales else
              f"{len(inst.dag.nodes)} nodos x 40 muestras")

    # ---------------------------------------------------------------- V5
    nombres_permitidos = set(inst.var_names)
    sospechosas = [n for n in sis.y_names
                   if re.fullmatch(r'v_\d+', n) and n not in nombres_permitidos]
    r.add("V5 cero variables auxiliares nuevas", not sospechosas,
          f"{len(sospechosas)} sospechosas" if sospechosas else
          f"{len(sis.y_names)} columnas, todas nodos del DAG")

    if verbose:
        r.imprimir(inst.ruta)
    return r
