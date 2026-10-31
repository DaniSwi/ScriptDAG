"""
Capa 12 (NUEVA): experimento de la Sec. 4.2 (Fig. 6 y Fig. 7 del paper).

Para cada P genera instancias con generador.py, las resuelve con REF y NSC
(bb.py) y resume la ganancia relativa promedio de NSC respecto de REF, en
tiempo (gt) y en cajas (gb), con la Definicion 3 del paper:

    t_r(a, b, pi) = T(a, pi) / (T(a, pi) + T(b, pi))
    ganancia(a vs b) = sum_pi t_r(b, a, pi) / sum_pi t_r(a, b, pi)

Uso (los valores por defecto son los del paper salvo el dominio):

    python experimento.py --tipo mult
    python experimento.py --tipo func --P 8 10 12 --instancias 5 --limite 120

Decisiones:
  E1  Dominio [-2, 2] por defecto (el paper usa [-1e8, 1e8]; en Python el
      branch & bound no termina con eso).  n = 8, pool de 12 terminos y
      Q = 100, como en el paper.
  E2  Solo entran al promedio las instancias que AMBAS estrategias
      terminaron dentro del limite: comparar un arbol truncado no tiene
      sentido.  Las truncadas se cuentan y se reportan.
  E3  Las tareas (instancia, estrategia) corren en paralelo en varios
      procesos.  Los tiempos se miden con la maquina cargada; afecta a REF y
      NSC por igual, pero gb (cajas) es la metrica que no depende de eso.
  E4  Cada resultado registra si la solucion conocida quedo a menos de 1e-6
      de una caja solucion (R se escribio redondeado, asi que el punto es
      solucion solo hasta ese redondeo): es una comprobacion de soundness
      gratis sobre cada instancia del experimento.
"""

import argparse
import csv
import multiprocessing as mp
import os
import time

import bb
import generador
from exportator import construir

CAMPOS = ['tipo', 'P', 'semilla', 'ruta', 'estrategia', 'cajas', 'tiempo',
          'completo', 'soluciones', 'certificadas', 'sol_conocida_dentro']


def _tarea(args):
    tipo, P, semilla, ruta, usar_nsc, limite = args
    inst, sis, _ = construir(ruta)
    r = bb.resolver(inst, sis, usar_nsc=usar_nsc, max_segundos=limite,
                    max_cajas=10 ** 7)
    dentro = ''
    if r.completo:
        p = bb.solucion_conocida(ruta)
        dentro = bb.distancia_a_soluciones(p, r.soluciones) < 1e-6
    return {'tipo': tipo, 'P': P, 'semilla': semilla, 'ruta': ruta,
            'estrategia': r.estrategia, 'cajas': r.cajas,
            'tiempo': round(r.tiempo, 4), 'completo': r.completo,
            'soluciones': len(r.soluciones), 'certificadas': sum(r.certificadas),
            'sol_conocida_dentro': dentro}


def ganancia(pares):
    """Definicion 3: pares = [(valor REF, valor NSC)] -> ganancia de NSC."""
    a = sum(ref / (ref + nsc) for ref, nsc in pares if ref + nsc > 0)
    b = sum(nsc / (ref + nsc) for ref, nsc in pares if ref + nsc > 0)
    return a / b if b > 0 else float('nan')


def resumir(filas):
    """{(tipo, P): dict con n, truncadas, gt, gb, ...}"""
    por = {}
    for f in filas:
        por.setdefault((f['tipo'], int(f['P']), int(f['semilla'])), {})[f['estrategia']] = f
    resumen = {}
    for (tipo, P, _), d in sorted(por.items()):
        r = resumen.setdefault((tipo, P), {'total': 0, 'truncadas': 0, 'tiempos': [],
                                           'cajas': [], 'sol_perdida': 0})
        r['total'] += 1
        if 'REF' not in d or 'NSC' not in d:
            continue
        ref, nsc = d['REF'], d['NSC']
        for x in (ref, nsc):
            if str(x['sol_conocida_dentro']) == 'False':
                r['sol_perdida'] += 1
        if str(ref['completo']) != 'True' or str(nsc['completo']) != 'True':
            r['truncadas'] += 1
            continue
        r['tiempos'].append((float(ref['tiempo']), float(nsc['tiempo'])))
        r['cajas'].append((int(ref['cajas']), int(nsc['cajas'])))
    for r in resumen.values():
        r['n'] = len(r['cajas'])
        r['gt'] = ganancia(r['tiempos'])
        r['gb'] = ganancia(r['cajas'])
    return resumen


def imprimir(resumen):
    print(f"{'tipo':5s} {'P':>3s} {'n':>3s} {'trunc':>5s} {'gt':>7s} {'gb':>7s} "
          f"{'cajas REF':>10s} {'cajas NSC':>10s}  sol. perdidas")
    for (tipo, P), r in sorted(resumen.items()):
        cr = sum(c[0] for c in r['cajas']) / max(r['n'], 1)
        cn = sum(c[1] for c in r['cajas']) / max(r['n'], 1)
        print(f"{tipo:5s} {P:3d} {r['n']:3d} {r['truncadas']:5d} {r['gt']:7.2f} "
              f"{r['gb']:7.2f} {cr:10.0f} {cn:10.0f}  {r['sol_perdida']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--tipo', choices=('mult', 'func'), default='mult')
    ap.add_argument('--P', type=int, nargs='+', default=[8, 10, 12, 14, 16])
    ap.add_argument('--instancias', type=int, default=20)
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--Q', type=int, default=100)
    ap.add_argument('--dominio', type=float, default=2.0)
    ap.add_argument('--limite', type=float, default=300.0, help='segundos por arbol')
    ap.add_argument('--procesos', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--carpeta', default='experimentos')
    ap.add_argument('--reanudar', action='store_true',
                    help='agrega al CSV existente y salta los arboles ya resueltos; '
                         'usar con los MISMOS parametros que la corrida original')
    a = ap.parse_args()

    carpeta_inst = os.path.join(a.carpeta, 'instancias')
    tareas = []
    for P in a.P:
        semillas = range(1, a.instancias + 1)
        rutas = generador.escribir(carpeta_inst, a.tipo, a.n, P, a.Q,
                                   (-a.dominio, a.dominio), semillas)
        for s, ruta in zip(semillas, rutas):
            for usar_nsc in (False, True):
                tareas.append((a.tipo, P, s, ruta, usar_nsc, a.limite))

    salida = os.path.join(a.carpeta, f"resultados_{a.tipo}.csv")
    filas = []
    if a.reanudar and os.path.exists(salida):
        with open(salida, newline='') as f:
            filas = list(csv.DictReader(f))
        hechas = {(int(r['P']), int(r['semilla']), r['estrategia']) for r in filas}
        tareas = [t for t in tareas
                  if (t[1], t[2], 'NSC' if t[4] else 'REF') not in hechas]
        print(f"reanudando: {len(filas)} arboles ya estaban en {salida}")
    print(f"{len(tareas)} arboles en {a.procesos} procesos -> {salida}")
    t0 = time.time()
    with open(salida, 'a' if filas else 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=CAMPOS)
        if not filas:
            w.writeheader()
        with mp.Pool(a.procesos) as pool:
            for k, fila in enumerate(pool.imap_unordered(_tarea, tareas), 1):
                w.writerow(fila)
                f.flush()
                filas.append(fila)
                if k % 10 == 0 or k == len(tareas):
                    print(f"  {k}/{len(tareas)} arboles, {time.time() - t0:.0f}s", flush=True)
    print()
    imprimir(resumir(filas))


if __name__ == '__main__':
    main()
