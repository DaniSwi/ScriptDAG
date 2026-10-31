"""
Pruebas del branch & bound (bb.py), del Newton intervalar (newton.py) y de
sin/cos por intervalos (B8).

  1. sin/cos por intervalos son sound (Monte Carlo denso, con intervalos
     chicos y centrados en maximos/minimos).  El validator NO cubre esto: en
     las instancias el argumento de sin mide mas de 2*pi y la rama de los
     extremos nunca se ejecuta.
  2. Jacobiano por intervalos: cumple el valor medio,
     f(x) - f(y) en J(X) (x - y) para x, y en X.  Es la propiedad de la que
     depende que Newton sea sound.
  3. Ejemplo (1) del paper: REF y NSC terminan, encuentran las MISMAS dos
     soluciones (residuo < 1e-6), Newton certifica ambas y NSC usa menos
     cajas que REF.
  4. Soundness del arbol completo: con una solucion fabricada (ver
     test_nsc._caso_con_solucion) alguna caja solucion la contiene.
  5. Generador (generador.py): las instancias pasan la auditoria y el arbol
     encuentra la solucion anotada.

    python test_bb.py
"""

import math
import random

import intervals as iv
import test_nsc as tn
from bb import resolver
from dominios import evaluar_en_punto, muestrear, propagar_dominios
from exportator import construir
from newton import jacobiano


def test_trig_sound():
    rng = random.Random(1)
    malos = total = 0
    extremos = [k * math.pi / 2 for k in range(-40, 40)]
    for _ in range(5000):
        c = rng.uniform(-100, 100) if rng.random() < 0.7 else rng.choice(extremos)
        w = 10 ** rng.uniform(-12, 1)
        a, b = c - w * rng.random(), c + w * rng.random()
        for nombre, f in (('sin', math.sin), ('cos', math.cos)):
            lo, hi = iv.funcion(nombre, [(a, b)])
            for x in [a, b] + [rng.uniform(a, b) for _ in range(10)]:
                total += 1
                malos += not lo <= f(x) <= hi
    assert malos == 0, f"sin/cos por intervalos no son sound: {malos}/{total}"
    ancho = iv.funcion('sin', [(2.0, 2.0 + 1e-8)])
    assert ancho[1] - ancho[0] < 1e-7, f"sin no es ajustado: {ancho}"
    print(f"[PASS] sin/cos sound en {total} muestras y ajustados en cajas chicas")


def valor_medio(inst, n, cota, semilla=0):
    """Cuenta violaciones de f(x) - f(y) en J(X)(x - y) en cajas aleatorias."""
    rng = random.Random(semilla)
    col = {v: j for j, v in enumerate(inst.var_names)}
    malos = total = 0
    for _ in range(n):
        centro = muestrear(inst.var_domains, rng, cota=cota)
        X = {}
        for v, x in centro.items():
            r, (lo, hi) = 10 ** rng.uniform(-4, 0), inst.var_domains[v]
            X[v] = (max(lo, x - r), min(hi, x + r))
        G = jacobiano(inst.dag, propagar_dominios(inst.dag, X), col)
        for _ in range(5):
            p = {v: rng.uniform(*X[v]) for v in X}
            q = {v: rng.uniform(*X[v]) for v in X}
            fp, fq = evaluar_en_punto(inst.dag, p), evaluar_en_punto(inst.dag, q)
            for raiz in inst.raices:
                if fp[raiz.id] is None or fq[raiz.id] is None:
                    continue
                acc = (0.0, 0.0)
                for v, j in col.items():
                    d = p[v] - q[v]
                    acc = iv.suma(acc, iv.mult(G[raiz.id].get(j, (0.0, 0.0)), (d, d)))
                total += 1
                malos += not iv.contiene(acc, fp[raiz.id] - fq[raiz.id], tol=1e-9)
    return malos, total


def test_jacobiano_valor_medio():
    for (inst, _, _), nombre, cota in ((tn._ejemplo_1(), 'ejemplo (1)', 5.0),
                                       (construir('instances/inst001.txt'), 'inst001', 3.0)):
        malos, total = valor_medio(inst, 200, cota)
        assert malos == 0, f"{nombre}: el jacobiano no encierra f(x)-f(y): {malos}/{total}"
        print(f"[PASS] jacobiano cumple el valor medio en {nombre} ({total} pares)")


def _residuos(p):
    x1, x2, x3 = p['x1'], p['x2'], p['x3']
    return [5 * x2 + x1 * x3 + x2 * x3 - x2 - 3,
            2 * x2 * math.sin(x1 + x3) + x1 * x3 + x2 * x3 - 2,
            math.sin(6 * x1 + 6 * x3) + 4 * x2 - 6 * (x1 + x3) + x1 * x3 + x2 * x3 + 10]


def _centro(caja):
    return {v: (lo + hi) / 2 for v, (lo, hi) in caja.items()}


def test_ejemplo_1():
    inst, sis, _ = tn._ejemplo_1()
    ref = resolver(inst, sis, usar_nsc=False)
    nsc = resolver(inst, sis, usar_nsc=True)
    assert ref.completo and nsc.completo
    for r in (ref, nsc):
        for s in r.soluciones:
            assert max(abs(x) for x in _residuos(_centro(s))) < 1e-6, \
                f"{r.estrategia}: caja solucion que no es solucion: {s}"
    # mismas soluciones: cada caja de REF esta cerca de una de NSC y viceversa
    def cerca(a, b):
        return max(abs(a[v] - b[v]) for v in a) < 1e-6
    cn = [_centro(s) for s in nsc.soluciones]
    cr = [_centro(s) for s in ref.soluciones]
    assert all(any(cerca(p, q) for q in cn) for p in cr)
    assert all(any(cerca(p, q) for q in cr) for p in cn)
    grupos = []
    for p in cn:
        if not any(cerca(p, q) for q in grupos):
            grupos.append(p)
    assert len(grupos) == 2, f"esperaba 2 soluciones, hay {len(grupos)}"
    assert nsc.cajas < ref.cajas, "NSC no redujo el arbol de busqueda"
    assert all(ref.certificadas) and all(nsc.certificadas), "Newton no certifico las soluciones"
    print(f"[PASS] ejemplo (1): 2 soluciones certificadas; REF {ref.cajas} cajas, "
          f"NSC {nsc.cajas} cajas (gb = {ref.cajas / nsc.cajas:.1f})")


def soundness_bb(inst, sis, n, semilla=0):
    """Devuelve (fallas, arboles completos revisados).

    Algunos sistemas fabricados son dificiles (mas de 20000 cajas): esos
    arboles quedan truncados y no se revisan, porque una caja no explorada no
    prueba nada.  Medido (con Newton): una biseccion que pierde el 10% de cada
    mitad falla en 9 de 24 arboles; si pierde solo el 1%, en 0 de 24.  Es
    decir, una perdida muy chica puede pasar desapercibida (ver CLAUDE.md)."""
    rng = random.Random(semilla)
    fallas, revisados = [], 0
    for _ in range(n):
        _, restricciones, p = tn._caso_con_solucion(inst, rng)
        for usar_nsc in (False, True):
            r = resolver(inst, sis, usar_nsc=usar_nsc, restricciones=restricciones,
                         max_cajas=20000)
            if not r.completo:
                continue
            revisados += 1
            if not any(all(c[v][0] <= x <= c[v][1] for v, x in p.items())
                       for c in r.soluciones):
                fallas.append((r.estrategia, p))
    return fallas, revisados


def test_soundness_bb():
    inst, sis, _ = tn._ejemplo_1()
    fallas, revisados = soundness_bb(inst, sis, 8)
    assert not fallas, f"el arbol perdio una solucion conocida: {fallas[:2]}"
    assert revisados >= 10, f"muy pocos arboles completos ({revisados}): la prueba no dice nada"
    print(f"[PASS] soundness del arbol: {revisados} arboles completos con solucion "
          "fabricada, ninguno la pierde (REF y NSC)")


def test_generador():
    """Instancias chicas de los dos tipos: pasan la auditoria y el arbol
    encuentra la solucion anotada en el comentario."""
    import os
    import tempfile
    import generador
    from bb import distancia_a_soluciones, solucion_conocida
    from validator import validar
    carpeta = tempfile.mkdtemp()
    for tipo in ('mult', 'func'):
        for ruta in generador.escribir(carpeta, tipo, 4, 6, 16, (-10, 10), [1, 2, 3]):
            inst, sis, dom = construir(ruta)
            assert validar(inst, sis, dom, verbose=False).ok, f"auditoria fallo: {ruta}"
            p = solucion_conocida(ruta)
            for usar_nsc in (False, True):
                r = resolver(inst, sis, usar_nsc=usar_nsc, max_segundos=60)
                assert r.completo, f"no termino: {ruta}"
                assert distancia_a_soluciones(p, r.soluciones) < 1e-6, \
                    f"{r.estrategia} no encontro la solucion conocida de {ruta}"
            os.remove(ruta)
    print("[PASS] generador: 6 instancias auditadas, solucion conocida encontrada (REF y NSC)")


if __name__ == "__main__":
    test_trig_sound()
    test_jacobiano_valor_medio()
    test_ejemplo_1()
    test_soundness_bb()
    test_generador()
    print("\nTodo OK.")
