"""
Capa 1: parseo de instancias .txt (formato tipo Minibex/Ibex).

Cambios respecto de la version anterior:
  B3  La regex de dominios aceptaba solo '-2' pegado.  '[- 2, 2]' (que esta en
      instprofe.txt) se leia como [2, 2]: unsound.  Ahora se acepta espacio
      entre el signo y el numero, notacion cientifica, 'inf' y '+oo'/'-oo'.
  B5  Antes, una linea de restriccion que no calzara con 'len(parts)==3' se
      descartaba EN SILENCIO y el validator no podia notarlo.  Ahora se cuentan
      las lineas y cualquier linea no parseable levanta ParseError.
  D3  Ya no se usa sympify(): se construye directamente el DAG normalizado.
      SymPy distribuye -6*(x1+x3) y con eso destruye la CSE del paper.
  --  El lado derecho ya no tiene que ser un float: si trae variables se pasa
      al lado izquierdo (lhs - rhs) y la restriccion queda contra 0.
"""

import re

from dag_builder import DagBuilder, ParseError

# numero con signo opcional separado por espacios, notacion cientifica o infinito
_NUM = r'[+-]?\s*(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?|inf|oo|INF)'
_RE_DOMINIO = re.compile(r'\[\s*(' + _NUM + r')\s*,\s*(' + _NUM + r')\s*\]')
_RE_VAR = re.compile(r'^([A-Za-z_]\w*)\s+in\s+(\[[^\]]*\])\s*;?\s*$')
_RE_OP = re.compile(r'(>=|<=|==|=|>|<)')


def parse_numero(texto):
    """'- 2' -> -2.0 ,  '1e8' -> 1e8 ,  '-inf' -> -inf   (arregla B3)."""
    t = texto.replace(' ', '')
    signo = 1.0
    while t and t[0] in '+-':
        if t[0] == '-':
            signo = -signo
        t = t[1:]
    if t.lower() in ('inf', 'oo'):
        return signo * float('inf')
    return signo * float(t)


def parse_dominio(cadena):
    """'[- 2, 2]' -> (-2.0, 2.0).  Si no calza, se levanta ParseError."""
    m = _RE_DOMINIO.search(cadena)
    if not m:
        raise ParseError(f"dominio no reconocido: {cadena!r}")
    lo, hi = parse_numero(m.group(1)), parse_numero(m.group(2))
    if lo > hi:
        raise ParseError(f"dominio invalido (lo > hi): {cadena!r}")
    return (lo, hi)


class Instancia:
    """Todo lo que el resto del pipeline necesita saber del archivo de entrada."""

    def __init__(self, ruta):
        self.ruta = ruta
        self.var_names = []              # orden de declaracion
        self.var_domains = {}            # nombre -> (lo, hi)
        self.ops = []                    # '=', '<=', ...
        self.rhs = []                    # float
        self.raices = []                 # DagNode raiz de cada restriccion
        self.coef_raiz = []              # Fraction: la raiz vale coef * nodo
        self.dag = None
        self.lineas_restriccion = 0      # cuantas lineas de restriccion tenia el .txt
        self.objetivo_omitido = False

    def n_restricciones(self):
        return len(self.ops)


def parse_instance(ruta, distribuir_sumas_escaladas=False):
    """Lee el archivo y devuelve una Instancia con el DAG ya construido."""
    inst = Instancia(ruta)
    dag = DagBuilder(distribuir_sumas_escaladas=distribuir_sumas_escaladas)
    inst.dag = dag

    with open(ruta, 'r') as f:
        lineas = f.readlines()

    modo = None
    for nro, cruda in enumerate(lineas, 1):
        linea = cruda.strip()
        if not linea or linea.startswith('//') or linea.startswith('#'):
            continue
        baja = linea.lower()

        if baja.startswith('minimize') or baja.startswith('maximize'):
            inst.objetivo_omitido = True
            continue
        if baja in ('variables', 'variables:'):
            modo = 'vars'
            continue
        if baja in ('constraints', 'constraints:'):
            modo = 'constraints'
            continue
        if baja == 'end':
            break

        if modo == 'vars':
            m = _RE_VAR.match(linea)
            if not m:
                raise ParseError(f"{ruta}:{nro}: declaracion de variable invalida: {linea!r}")
            nombre, dom = m.group(1), m.group(2)
            if nombre in inst.var_domains:
                raise ParseError(f"{ruta}:{nro}: variable duplicada: {nombre}")
            inst.var_names.append(nombre)
            inst.var_domains[nombre] = parse_dominio(dom)

        elif modo == 'constraints':
            inst.lineas_restriccion += 1          # B5: se cuenta TODA linea
            cuerpo = linea.rstrip(';').strip()
            partes = _RE_OP.split(cuerpo)
            if len(partes) != 3:
                raise ParseError(
                    f"{ruta}:{nro}: restriccion no parseable "
                    f"({len(partes)} piezas tras separar por el operador): {cuerpo[:80]!r}")
            izq_txt, op, der_txt = partes[0], partes[1], partes[2]
            coef_i, nodo_i = dag.desde_string(izq_txt, inst.var_domains)
            coef_d, nodo_d = dag.desde_string(der_txt, inst.var_domains)

            if nodo_d is None:
                # caso normal:  f(x) op constante
                if nodo_i is None:
                    raise ParseError(f"{ruta}:{nro}: restriccion sin variables")
                raiz, coef, rhs = nodo_i, coef_i, float(coef_d)
            else:
                # lado derecho con variables: se pasa todo a la izquierda
                forma = dag._suma(('n', coef_i, nodo_i), ('n', -coef_d, nodo_d))
                coef, raiz = dag._materializar(forma)
                rhs = 0.0
                if raiz is None:
                    raise ParseError(f"{ruta}:{nro}: restriccion trivial")
            dag.registrar_raiz(raiz)
            inst.raices.append(raiz)
            inst.coef_raiz.append(coef)
            inst.ops.append('=' if op == '==' else op)
            inst.rhs.append(rhs / float(coef) if coef != 1 and op in ('=', '==') else rhs)
            if coef != 1 and op not in ('=', '=='):
                # con desigualdad, dividir por un coeficiente negativo invierte el
                # sentido; se corrige explicitamente para no perder soluciones.
                inst.rhs[-1] = rhs / float(coef)
                if coef < 0:
                    inst.ops[-1] = {'<': '>', '<=': '>=', '>': '<', '>=': '<='}[inst.ops[-1]]

    if not inst.var_names:
        raise ParseError(f"{ruta}: no se declararon variables")
    if inst.n_restricciones() != inst.lineas_restriccion:
        raise ParseError(f"{ruta}: se leyeron {inst.lineas_restriccion} lineas de "
                         f"restriccion pero solo se parsearon {inst.n_restricciones()}")
    return inst
