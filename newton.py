"""
Capa 10 (NUEVA): Newton intervalar (Hansen-Sengupta) sobre el DAG.

Sin este contractor el branch & bound de bb.py no termina en sistemas del
tamano de inst001: cerca de una solucion HC4/NSC contraen linealmente y el
arbol se llena de cajas.  Ibex usa Newton en TODAS las estrategias del paper
(REF, NSC, NSC+Simplex), asi que agregarlo no favorece a ninguna.

Paso de Newton sobre la caja X, para f(x) = raiz(x) - R  (R = lado derecho):
    x~ = punto medio de X
    J  = jacobiano por intervalos sobre TODO X
    C  = inversa de mid(J)                       (precondicionador, float)
    Gauss-Seidel sobre  (C J) (x - x~) = -C f(x~):
        X_i <- X_i inter ( x~_i - (b_i + sum_{j!=i} A_ij (X_j - x~_j)) / A_ii )

Decisiones:
  T1  El jacobiano se calcula por diferenciacion automatica forward sobre el
      DAG.  Los VALORES de los nodos se evaluan sobre toda X (forward puro),
      no con los dominios contraidos por HC4/NSC: esos solo encierran los
      valores en las soluciones, y el teorema del valor medio necesita
      encerrarlos en todo X.
  T2  Si un A_ii contiene 0, esa variable no se contrae (no se usa division
      extendida).  Si mid(J) es singular o J tiene entradas no acotadas, el
      paso no hace nada.  Siempre sound.
  T3  Certificacion: si cada componente de la imagen de Gauss-Seidel cae en el
      INTERIOR de X_i (y ningun A_ii contiene 0), el teorema de
      Hansen-Sengupta garantiza que X contiene una unica solucion.
  T4  Solo se aplica a sistemas cuadrados de igualdades (n ecuaciones '=' y n
      variables).  En otro caso aplicable() es False y el solver sigue sin
      Newton.
  T5  Derivadas sin regla conocida (funciones raras, potencias fraccionarias
      con base que toca 0, log/sqrt con argumento que toca 0) dan (-inf, inf):
      la columna queda inutil y Newton simplemente no contrae.
"""

import math

import numpy as np

import intervals as iv
from dominios import propagar_dominios


# ===================================================================
#  Jacobiano por intervalos (diferenciacion automatica forward)
# ===================================================================

def _acumular(destino, gradiente, factor):
    for j, g in gradiente.items():
        t = iv.mult(factor, g)
        destino[j] = iv.suma(destino[j], t) if j in destino else t


def _pow_real(I, e):
    """I^e con I > 0 y e no entero (monotona), redondeo hacia afuera holgado."""
    a, b = I
    if not a > 0.0:
        return iv.ENTERO
    try:
        pa = a ** e
        pb = b ** e if math.isfinite(b) else (iv.INF if e > 0 else 0.0)
    except OverflowError:
        return iv.ENTERO
    lo, hi = min(pa, pb), max(pa, pb)
    return (iv.dn(lo * (1 - 1e-12)), iv.up(hi * (1 + 1e-12)))


def _derivada_func(nombre, a):
    """f'(a) por intervalos, para funciones de un argumento (T5)."""
    if nombre == 'sin':
        return iv.funcion('cos', [a])
    if nombre == 'cos':
        s = iv.funcion('sin', [a])
        return (-s[1], -s[0])
    if nombre == 'exp':
        return iv.funcion('exp', [a])
    if nombre == 'log':
        return iv.reciproco(a) if a[0] > 0.0 else iv.ENTERO
    if nombre == 'sqrt':
        if not a[0] > 0.0:
            return iv.ENTERO
        return iv.mult((0.5, 0.5), iv.reciproco(iv.funcion('sqrt', [a])))
    if nombre == 'atan':
        return iv.reciproco(iv.suma((1.0, 1.0), iv.potencia_entera(a, 2)))
    if nombre == 'tanh':
        t2 = iv.potencia_entera(iv.funcion('tanh', [a]), 2)
        return iv.interseccion(iv.resta((1.0, 1.0), t2), (0.0, 1.0)) or (0.0, 1.0)
    if nombre == 'Abs':
        if a[0] >= 0.0:
            return (1.0, 1.0)
        if a[1] <= 0.0:
            return (-1.0, -1.0)
        return (-1.0, 1.0)
    return iv.ENTERO


def jacobiano(dag, dom, columna_de_var):
    """{node.id: {j: intervalo}}: derivada de cada nodo respecto de cada
    variable original j, con los valores de los nodos en dom (T1)."""
    grad = {}
    for nodo in dag.nodes:
        k = nodo.kind
        g = {}
        if k == 'var':
            j = columna_de_var.get(nodo.data)
            if j is not None:
                g[j] = (1.0, 1.0)
        elif k == 'add':
            _, coefs = nodo.data
            for c, h in zip(coefs, nodo.args):
                _acumular(g, grad[h.id], iv.from_fraction(c))
        elif k == 'mul':
            vals = [dom[h.id] for h in nodo.args]
            pre = [(1.0, 1.0)]
            for v in vals:
                pre.append(iv.mult(pre[-1], v))
            suf = [(1.0, 1.0)]
            for v in reversed(vals):
                suf.append(iv.mult(v, suf[-1]))
            suf.reverse()
            for idx, h in enumerate(nodo.args):
                _acumular(g, grad[h.id], iv.mult(pre[idx], suf[idx + 1]))
        elif k == 'pow':
            e = nodo.data
            b = dom[nodo.args[0].id]
            if e.denominator == 1:
                d = iv.mult(iv.from_fraction(e), iv.potencia_entera(b, int(e) - 1))
            else:
                d = iv.mult(iv.from_fraction(e), _pow_real(b, float(e - 1)))
            _acumular(g, grad[nodo.args[0].id], d)
        elif k == 'func':
            if len(nodo.args) == 1:
                d = _derivada_func(nodo.data, dom[nodo.args[0].id])
                _acumular(g, grad[nodo.args[0].id], d)
            else:
                for h in nodo.args:
                    for j in grad[h.id]:
                        g[j] = iv.ENTERO
        grad[nodo.id] = g
    return grad


# ===================================================================
#  Paso de Hansen-Sengupta
# ===================================================================

class Newton:
    """Datos fijos del sistema: que nodos son las ecuaciones y sus lados
    derechos (los mismos intervalos que usa nsc.Preparado)."""

    def __init__(self, inst, restricciones):
        self.dag = inst.dag
        presentes = {n.data for n in inst.dag.nodes if n.kind == 'var'}
        self.var_names = [v for v in inst.var_names if v in presentes]
        self.col = {v: j for j, v in enumerate(self.var_names)}
        self.ecuaciones = restricciones
        self.es_cuadrado = (len(restricciones) == len(self.var_names)
                            and all(op in ('=', '==') for op in inst.ops))

    def aplicable(self):
        return self.es_cuadrado and len(self.var_names) > 0

    def paso(self, caja):
        """Devuelve (caja_nueva o None si vacia, certificada)."""
        n = len(self.var_names)
        X = [caja[v] for v in self.var_names]
        if not all(math.isfinite(a) and math.isfinite(b) for a, b in X):
            return caja, False

        dom = propagar_dominios(self.dag, caja)
        grad = jacobiano(self.dag, dom, self.col)
        J = [[grad[nid].get(j, (0.0, 0.0)) for j in range(n)] for nid, _ in self.ecuaciones]
        if not all(math.isfinite(a) and math.isfinite(b) for fila in J for a, b in fila):
            return caja, False

        xm = [a / 2.0 + b / 2.0 for a, b in X]
        val = propagar_dominios(self.dag, {v: (x, x) for v, x in zip(self.var_names, xm)})
        F = [iv.resta(val[nid], R) for nid, R in self.ecuaciones]

        Jm = np.array([[a / 2.0 + b / 2.0 for a, b in fila] for fila in J])
        try:
            C = np.linalg.inv(Jm)
        except np.linalg.LinAlgError:
            return caja, False
        if not np.all(np.isfinite(C)):
            return caja, False

        A = [[(0.0, 0.0)] * n for _ in range(n)]
        b = [(0.0, 0.0)] * n
        for i in range(n):
            for k in range(n):
                c = float(C[i, k])
                if c == 0.0:
                    continue
                ci = (c, c)
                b[i] = iv.suma(b[i], iv.mult(ci, F[k]))
                for j in range(n):
                    A[i][j] = iv.suma(A[i][j], iv.mult(ci, J[k][j]))

        certificada = True
        X = list(X)
        for i in range(n):
            a_ii = A[i][i]
            if a_ii[0] <= 0.0 <= a_ii[1]:
                certificada = False
                continue
            s = b[i]
            for j in range(n):
                if j != i:
                    s = iv.suma(s, iv.mult(A[i][j], iv.resta(X[j], (xm[j], xm[j]))))
            y = iv.resta((xm[i], xm[i]), iv.mult(s, iv.reciproco(a_ii)))
            if not (X[i][0] < y[0] and y[1] < X[i][1]):
                certificada = False
            nuevo = iv.interseccion(X[i], y)
            if nuevo is None:
                return None, False
            X[i] = nuevo
        return dict(zip(self.var_names, X)), certificada

    def smear_sum_rel(self, caja):
        """Variable a bisectar segun SmearSumRel (el bisector del paper).

        Para cada ecuacion i, el impacto de x_j es |J_ij| * ancho(X_j),
        normalizado por la suma de la fila; se suma sobre las ecuaciones y se
        elige el maximo.  Devuelve None si no se puede calcular (y bb.py usa
        la variable mas ancha)."""
        dom = propagar_dominios(self.dag, caja)
        grad = jacobiano(self.dag, dom, self.col)
        puntaje = [0.0] * len(self.var_names)
        for nid, _ in self.ecuaciones:
            imp = []
            for j, v in enumerate(self.var_names):
                a, b = grad[nid].get(j, (0.0, 0.0))
                ancho = caja[v][1] - caja[v][0]
                m = max(abs(a), abs(b)) * ancho if ancho > 0 else 0.0
                imp.append(min(m, 1e300) if m == m else 1e300)
            total = sum(imp)
            if total > 0.0 and math.isfinite(total):
                for j, x in enumerate(imp):
                    puntaje[j] += x / total
        j = max(range(len(puntaje)), key=puntaje.__getitem__)
        if puntaje[j] <= 0.0:
            return None
        return self.var_names[j]

    def certificar(self, caja, intentos=6):
        """T3 con epsilon-inflacion: en una caja ya diminuta el redondeo impide
        que la imagen caiga ESTRICTAMENTE dentro, asi que se prueba sobre una
        caja algo inflada alrededor de ella (lo que hace Ibex).  True significa
        que la caja inflada (radio ~1.5x el original + 1e-10 relativo) contiene
        exactamente una solucion."""
        X = dict(caja)
        for _ in range(intentos):
            inflada = {}
            for v, (a, b) in X.items():
                m, r = a / 2.0 + b / 2.0, (b - a) / 2.0
                d = 1.5 * r + 1e-10 * max(1.0, abs(m))
                inflada[v] = (iv.dn(m - d), iv.up(m + d))
            nueva, cert = self.paso(inflada)
            if cert:
                return True
            if nueva is None:
                return False
            X = nueva
        return False
