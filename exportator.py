"""
Capa 5: exportacion.

Decision de diseno (D4 de la revision): la salida OFICIAL es un .npz ANEXO al
sistema original.  NO se reescribe el problema con variables nuevas, porque el
paper dice explicitamente "we do not modify the original system" y "we do not
need concretization": las w_i y las y_j son ETIQUETAS de nodos del DAG, no
variables nuevas del problema.  Reescribir el .txt (estilo I-CSE / Ceberio)
sigue disponible con exportar_txt=True, pero ya no es el camino por defecto y
marca las auxiliares como no bisectables.

Arregla D5: el .npz de antes no traia ni las w, ni sus dominios, ni los
dominios de y, ni el subsistema de productos; era imposible armar A*y = b
sin re-parsear strings.

Contenido del .npz
------------------
  A            (m,n) float64   matriz de coeficientes
  c            (m,)  float64   termino constante de cada fila:  w_i = c_i + A[i]*y
  A_aug        (m,n+m)         [A | -I]  tal que  A_aug * [y;w] = -c
  A_fila/col/num/den            A en aritmetica EXACTA (COO de racionales)
  y_names, y_lo, y_hi           columnas y sus dominios (sound)
  w_names, w_lo, w_hi           filas y sus dominios (sound)
  w_es_columna                  indice de columna de cada w, o -1 (columnas mixtas)
  var_names, var_lo, var_hi     variables ORIGINALES del problema
  prod_nombre/tipo/dato/args    subsistema no lineal (Proposicion 1)
  restricciones                 'nodo op rhs' (ecuaciones cascaron)
  meta_*                        conteos para auditoria
"""

import numpy as np

import dominios
import intervals as iv
from matrix_extractor import extraer_sistema_lineal, registrar_restricciones
from parser_instancias import parse_instance


def construir(ruta_instancia, calcular_dominios=True, distribuir_sumas_escaladas=False):
    """parse -> DAG -> sistema lineal -> dominios.  Devuelve (inst, sis, dom)."""
    inst = parse_instance(ruta_instancia, distribuir_sumas_escaladas)
    sis = extraer_sistema_lineal(inst.dag)
    registrar_restricciones(sis, inst.raices, inst.ops, inst.rhs)
    dom = dominios.propagar_dominios(inst.dag, inst.var_domains) if calcular_dominios else None
    return inst, sis, dom


def _cajas(nodos, dom):
    if dom is None:
        return (np.full(len(nodos), -np.inf), np.full(len(nodos), np.inf))
    lo = np.array([dom[n.id][0] for n in nodos], dtype=float)
    hi = np.array([dom[n.id][1] for n in nodos], dtype=float)
    return lo, hi


def exportar_npz(ruta_npz, inst, sis, dom):
    filas, cols, nums, dens = [], [], [], []
    for i, fila in enumerate(sis.A_exacta):
        for j, q in fila.items():
            filas.append(i); cols.append(j)
            nums.append(str(q.numerator)); dens.append(str(q.denominator))

    y_lo, y_hi = _cajas(sis.y_nodes, dom)
    w_lo, w_hi = _cajas(sis.w_nodes, dom)
    var_lo = np.array([inst.var_domains[v][0] for v in inst.var_names], dtype=float)
    var_hi = np.array([inst.var_domains[v][1] for v in inst.var_names], dtype=float)
    A_aug, rhs_aug = sis.aumentada()

    np.savez(
        ruta_npz,
        A=sis.A, c=sis.c, A_aug=A_aug, rhs_aug=rhs_aug,
        A_fila=np.array(filas, dtype=np.int64), A_col=np.array(cols, dtype=np.int64),
        A_num=np.array(nums), A_den=np.array(dens),
        y_names=np.array(sis.y_names), y_lo=y_lo, y_hi=y_hi,
        w_names=np.array(sis.w_names), w_lo=w_lo, w_hi=w_hi,
        w_es_columna=np.array([sis.col_de_nodo.get(n.id, -1) for n in sis.w_nodes],
                              dtype=np.int64),
        var_names=np.array(inst.var_names), var_lo=var_lo, var_hi=var_hi,
        prod_nombre=np.array([p[0] for p in sis.subsistema_productos]),
        prod_tipo=np.array([p[1] for p in sis.subsistema_productos]),
        prod_dato=np.array([p[2] for p in sis.subsistema_productos]),
        prod_args=np.array([" | ".join(p[3]) for p in sis.subsistema_productos]),
        restricciones=np.array([f"{a} {b} {c!r}" for a, b, c in sis.restricciones]),
        meta_n_restricciones=np.int64(inst.n_restricciones()),
        meta_lineas_restriccion=np.int64(inst.lineas_restriccion),
        meta_n_nodos=np.int64(len(inst.dag.nodes)),
        meta_aux_nuevas=np.int64(0),
    )


def exportar_txt(ruta_txt, inst, sis, dom):
    """Salida OPCIONAL estilo Ceberio: sistema reescrito con w_i explicitas.

    Solo para alimentar un solver que no se puede tocar por dentro.  Las
    auxiliares se marcan como no bisectables y sus dominios son sound
    (-inf/inf cuando no se pueden acotar), nunca el [-1e8,1e8] de antes (B2).
    """
    def dstr(par):
        lo, hi = par
        f = lambda v: ('-inf' if v == -np.inf else ('inf' if v == np.inf else repr(float(v))))
        return f"[{f(lo)}, {f(hi)}]"

    y_lo, y_hi = _cajas(sis.y_nodes, dom)
    w_lo, w_hi = _cajas(sis.w_nodes, dom)

    with open(ruta_txt, 'w') as f:
        f.write("// Reescritura estilo I-CSE/Ceberio del sistema original.\n")
        f.write("// Las variables w_i son AUXILIARES: no deben bisectarse.\n")
        f.write("variables\n")
        for v in inst.var_names:
            f.write(f"  {v} in {dstr(inst.var_domains[v])};\n")
        for k, nombre in enumerate(sis.w_names):
            f.write(f"  {nombre} in {dstr((w_lo[k], w_hi[k]))};   // no-bisectable\n")
        f.write("\nconstraints\n")
        for i, nombre in enumerate(sis.w_names):
            piezas = []
            const = sis.c_exacta[i]
            if const != 0:
                piezas.append(str(float(const)))
            for j, q in sorted(sis.A_exacta[i].items()):
                col = sis.y_names[j].replace('^', '^')
                piezas.append(f"({float(q)})*({col})" if q != 1 else f"({col})")
            f.write(f"  {nombre} = {' + '.join(piezas) if piezas else '0'};\n")
        for nombre, op, rhs in sis.restricciones:
            f.write(f"  {nombre} {op} {rhs!r};\n")
        f.write("end\n")


def procesar_instancia(ruta_instancia, ruta_npz, ruta_txt=None,
                       calcular_dominios=True, exportar_texto=False,
                       distribuir_sumas_escaladas=False):
    inst, sis, dom = construir(ruta_instancia, calcular_dominios,
                               distribuir_sumas_escaladas)
    exportar_npz(ruta_npz, inst, sis, dom)
    if exportar_texto and ruta_txt:
        exportar_txt(ruta_txt, inst, sis, dom)
    return inst, sis, dom
