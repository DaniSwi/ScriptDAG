"""
Pruebas del contractor NSC (Seccion 3.3 del paper, modulo nsc.py).

  1. El conditioner reproduce EXACTAMENTE la P y la PA de la Seccion 3.3.2.
  2. P*A == PA en racionales, para cada instancia.
  3. Soundness: en cajas aleatorias que contienen una solucion CONOCIDA, ni
     REF ni NSC la eliminan.  La solucion se fabrica asi: se toma un punto p
     y el lado derecho de cada restriccion se reemplaza por la envoltura por
     intervalos del valor de su raiz en p.  Asi p es solucion real del
     sistema, sin depender de evaluaciones en float.
  4. Sobre el ejemplo (1) literal del paper, NSC contrae mas que HC4 en el
     DAG (REF) y nunca menos.

    python test_nsc.py
"""

import glob
import os
import random
import tempfile
from fractions import Fraction as F

import intervals as iv
from dominios import muestrear, propagar_dominios
from exportator import construir
from nsc import Preparado, conditioner, conditioner_de, contraer

EJEMPLO_1_TXT = """Variables
x1 in [-1,5];
x2 in [-1,5];
x3 in [-1,5];
Constraints
5*x2 + x1*x3 + x2*x3 - x2 = 3;
2*x2*sin(x1 + x3) + x1*x3 + x2*x3 = 2;
sin(6*x1 + 6*x3) + 4*x2 - 6*(x1 + x3) + x1*x3 + x2*x3 = -10;
end
"""


def _ejemplo_1():
    fd, ruta = tempfile.mkstemp(suffix='.txt')
    with os.fdopen(fd, 'w') as f:
        f.write(EJEMPLO_1_TXT)
    try:
        return construir(ruta)
    finally:
        os.remove(ruta)


def _igual_salvo_signo(a, b):
    return a == b or a == {k: -v for k, v in b.items()}


def test_conditioner_reproduce_paper():
    """La matriz A y la P de la Sec. 3.3.2, copiadas del paper."""
    A = [{0: F(2), 1: F(1), 2: F(1)},
         {1: F(1), 2: F(1), 3: F(1), 4: F(-6), 5: F(4)},
         {1: F(1), 2: F(1), 5: F(4)},
         {6: F(6), 7: F(6)},
         {6: F(1), 7: F(1)}]
    P_paper = [{0: F(1), 2: F(-1)}, {2: F(1)}, {1: F(-1), 2: F(1)}, {3: F(1)}]
    PA_paper = [{0: F(2), 5: F(-4)}, {1: F(1), 2: F(1), 5: F(4)},
                {3: F(-1), 4: F(6)}, {6: F(6), 7: F(6)}]
    c = conditioner(A, 8)
    for k, r in enumerate(c.filas):
        assert _igual_salvo_signo(c.P[r], P_paper[k]), (k, c.P[r])
        assert _igual_salvo_signo(c.PA[r], PA_paper[k]), (k, c.PA[r])
    assert len(c.dependientes) == 1
    assert c.P[c.dependientes[0]] == {4: F(1), 3: F(-1, 6)}, "falta la fila -0.17"
    print("[PASS] el conditioner reproduce P y PA de la Seccion 3.3.2")


def _instancias():
    return sorted(glob.glob(os.path.join('instances', '*.txt')))


def test_P_por_A_exacto():
    for ruta in _instancias():
        inst, sis, _ = construir(ruta, calcular_dominios=False)
        c = conditioner_de(sis)
        for r in range(len(sis.A_exacta)):
            prod = {}
            for i, p in c.P[r].items():
                for j, q in sis.A_exacta[i].items():
                    prod[j] = prod.get(j, 0) + p * q
            prod = {j: v for j, v in prod.items() if v != 0}
            assert prod == c.PA[r], f"{ruta}: fila {r} de P*A != PA"
        print(f"[PASS] P*A == PA exacto en {os.path.basename(ruta)} "
              f"(rango {c.rango()}/{len(sis.A_exacta)})")


def _caso_con_solucion(inst, rng):
    """(caja inicial, restricciones, punto) con el punto solucion real."""
    p = muestrear(inst.var_domains, rng, cota=10.0)
    env = propagar_dominios(inst.dag, {v: (x, x) for v, x in p.items()})
    restricciones = [(raiz.id, (iv.dn(env[raiz.id][0]), iv.up(env[raiz.id][1])))
                     for raiz in inst.raices]
    caja = {}
    for v, x in p.items():
        lo, hi = inst.var_domains[v]
        r = 10 ** rng.uniform(-3, 1)
        caja[v] = (max(lo, x - r * rng.random()), min(hi, x + r * rng.random()))
    return caja, restricciones, p


def soundness(inst, sis, n, semilla=0):
    """Devuelve la lista de violaciones (vacia si todo es sound)."""
    rng = random.Random(semilla)
    cond = conditioner_de(sis)
    ids = {nd.data: nd.id for nd in inst.dag.nodes if nd.kind == 'var'}
    violaciones = []
    for _ in range(n):
        caja, restricciones, p = _caso_con_solucion(inst, rng)
        for lineal in (False, True):
            prep = Preparado(inst, sis, cond if lineal else None)
            prep.restricciones = restricciones
            dom, _ = contraer(prep, caja, usar_lineal=lineal)
            if dom is None:
                violaciones.append(('caja vacia', lineal, p))
                continue
            for v, x in p.items():
                lo, hi = dom[ids[v]]
                if not lo <= x <= hi:
                    violaciones.append((v, lineal, x, (lo, hi)))
    return violaciones


def test_soundness():
    casos = [(_ejemplo_1(), 'ejemplo (1)', 40)]
    casos += [(construir(r), os.path.basename(r), 3 if 'profe' in r else 20)
              for r in _instancias()]
    for (inst, sis, _), nombre, n in casos:
        malas = soundness(inst, sis, n)
        assert not malas, f"{nombre}: NSC/REF eliminaron una solucion: {malas[:2]}"
        print(f"[PASS] soundness en {nombre}: {n} cajas con solucion conocida, "
              f"ninguna se pierde (REF y NSC)")


def test_nsc_contrae_mas_que_ref():
    inst, sis, _ = _ejemplo_1()
    cond = conditioner_de(sis)
    ref, _ = contraer(Preparado(inst, sis, None), usar_lineal=False)
    nsc, _ = contraer(Preparado(inst, sis, cond), usar_lineal=True)
    mejor = False
    for nd in inst.dag.nodes:
        if nd.kind != 'var':
            continue
        (a, b), (c, d) = ref[nd.id], nsc[nd.id]
        assert a - 1e-9 <= c and d <= b + 1e-9, f"NSC contrajo menos que REF en {nd.data}"
        mejor = mejor or (d - c) < 0.99 * (b - a)
    assert mejor, "NSC no contrajo mas que REF en el ejemplo (1)"
    print("[PASS] ejemplo (1): NSC contrae estrictamente mas que REF")


if __name__ == "__main__":
    test_conditioner_reproduce_paper()
    test_P_por_A_exacto()
    test_soundness()
    test_nsc_contrae_mas_que_ref()
    print("\nTodo OK.")
