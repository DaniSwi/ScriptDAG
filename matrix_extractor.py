"""
Capa 3: extraccion del subsistema lineal n-ario (Seccion 3.2 del paper).

Cambios respecto de la version anterior:
  D1  Se elimina por completo el mecanismo de variables auxiliares v_i.
      Los productos, potencias y funciones son NODOS del DAG, es decir
      COLUMNAS de A, tal como dice el paper.  Los coeficientes constantes
      quedan dentro de A y no escondidos dentro de un v_i.
  D2  Un nodo aparece en UNA sola columna: la separacion coeficiente/nodo la
      hizo el DagBuilder, asi que x2^2, -x2^2 y 5*x2^2 comparten columna.
  D5  Se exporta el sistema completo:  A*y + c = w , con los nombres y la
      correspondencia fila <-> w, y ademas el subsistema de productos
      (Proposicion 1) que antes no se emitia.
  B6  A se conserva tambien en aritmetica exacta (Fraction) ademas del float.

Forma del sistema (Fig. 3 del paper):

        w_i  =  c_i + sum_j  A[i][j] * y_j          i = 0 .. m-1

  - cada fila i corresponde a un nodo 'add' del DAG (una suma n-aria);
  - las columnas y_j son nodos del DAG: variables, potencias, productos,
    funciones, y tambien otros w_k (columnas mixtas, el y5=b5 de la Fig. 3);
  - la matriz aumentada [A | -I] da el sistema homogeneo  [A|-I]*[y;w] = -c.
"""

from fractions import Fraction

import numpy as np

from dag_builder import DagBuilder


class SistemaLineal:
    def __init__(self):
        self.A = None              # np.ndarray float (m x n), redondeo al mas cercano
        self.A_exacta = []         # [ {j: Fraction} ]  filas exactas
        self.c = None              # np.ndarray float (m,)   termino constante
        self.c_exacta = []
        self.y_nodes = []
        self.y_names = []
        self.w_nodes = []
        self.w_names = []
        self.col_de_nodo = {}      # node.id -> j
        self.fila_de_nodo = {}     # node.id -> i
        self.subsistema_productos = []   # [(nombre, tipo, dato, [nombres operandos])]
        self.restricciones = []          # [(nombre_raiz, op, rhs)]
        self.var_names = []

    # ---------------------------------------------------------------- info
    def shape(self):
        return self.A.shape

    def aumentada(self):
        """[A | -I] tal que  [A|-I] * [y ; w] = -c ."""
        m, n = self.A.shape
        M = np.zeros((m, n + m))
        M[:, :n] = self.A
        w_col = {nd.id: k for k, nd in enumerate(self.w_nodes)}
        for k in range(m):
            M[k, n + k] = -1.0
        # si un w tambien es columna de y, la columna duplicada es intencional:
        # el paper permite columnas mixtas.
        return M, np.array([-v for v in self.c])

    def resumen(self):
        m, n = self.A.shape
        densidad = float(np.count_nonzero(self.A)) / (m * n) if m * n else 0.0
        coefs = self.A[self.A != 0]
        return {
            'filas_w': m, 'columnas_y': n,
            'densidad': round(densidad, 4),
            'coef_min': float(coefs.min()) if coefs.size else 0.0,
            'coef_max': float(coefs.max()) if coefs.size else 0.0,
            'variables_auxiliares_nuevas': 0,
            'nodos_producto': len(self.subsistema_productos),
        }


def _nombre(nodo):
    return nodo.to_str(usar_labels=True)


def extraer_sistema_lineal(dag: DagBuilder):
    """Construye A, c, y, w a partir del DAG ya normalizado."""
    sis = SistemaLineal()

    # 1. una fila por cada nodo 'add' (suma n-aria).  Etiqueta w_i.
    for nodo in dag.nodes:
        if nodo.kind == 'add':
            nodo.label = f"w_{len(sis.w_nodes)}"
            sis.fila_de_nodo[nodo.id] = len(sis.w_nodes)
            sis.w_nodes.append(nodo)
            sis.w_names.append(nodo.label)

    # 2. columnas: todo nodo que aparece como sumando de alguna fila
    def columna(nodo):
        j = sis.col_de_nodo.get(nodo.id)
        if j is None:
            j = len(sis.y_nodes)
            sis.col_de_nodo[nodo.id] = j
            sis.y_nodes.append(nodo)
            sis.y_names.append(_nombre(nodo))
        return j

    for nodo in sis.w_nodes:
        const, coefs = nodo.data
        fila = {}
        for coef, hijo in zip(coefs, nodo.args):
            j = columna(hijo)
            fila[j] = fila.get(j, Fraction(0)) + coef   # "coefficients are summed up"
        sis.A_exacta.append(fila)
        sis.c_exacta.append(const)

    # 3. version numerica
    m, n = len(sis.w_nodes), len(sis.y_nodes)
    A = np.zeros((m, n))
    for i, fila in enumerate(sis.A_exacta):
        for j, q in fila.items():
            A[i, j] = float(q)
    sis.A = A
    sis.c = np.array([float(q) for q in sis.c_exacta])

    # 4. subsistema de productos / potencias / funciones (Proposicion 1)
    for nodo in dag.nodes:
        if nodo.kind in ('mul', 'pow', 'func'):
            operandos = [_nombre(a) for a in nodo.args]
            dato = (str(nodo.data) if nodo.kind == 'func'
                    else (str(nodo.data) if nodo.kind == 'pow' else ''))
            sis.subsistema_productos.append((_nombre(nodo), nodo.kind, dato, operandos))

    sis.var_names = [nd.data for nd in dag.nodes if nd.kind == 'var']
    return sis


def registrar_restricciones(sis, raices, ops, rhs_vals):
    """Las 'ecuaciones cascaron': cada restriccion original queda como
    <nombre del nodo raiz> <op> <constante>, sin tocar el sistema original."""
    for nodo, op, rhs in zip(raices, ops, rhs_vals):
        sis.restricciones.append((_nombre(nodo), op, float(rhs)))
    return sis
