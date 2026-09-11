"""
Capa 0 (NUEVA): aritmetica de intervalos *sound* con redondeo dirigido.

Arregla los bugs B1, B2 y B6 de la revision:
  B1  potencias negativas y divisiones cuyo denominador contiene 0
      devolvian intervalos que NO contienen la expresion (eliminaban soluciones).
  B2  los dominios de w_i/v_i se truncaban a [-1e8, 1e8].
  B6  toda la aritmetica usaba float sin redondeo dirigido.

Convencion: un intervalo es una tupla (lo, hi) de floats con lo <= hi.
Se permiten -inf / +inf.  Nunca se devuelve un intervalo que no contenga
al conjunto real de valores de la expresion: si no sabemos acotar, se
devuelve (-inf, +inf), que es sound aunque no filtre.
"""

import math
from fractions import Fraction

INF = math.inf
NEG_INF = -math.inf
ENTERO = (-INF, INF)


# ---------------------------------------------------------------- redondeo
def dn(x):
    """Redondeo hacia abajo (siguiente float hacia -inf)."""
    if x != x:                      # NaN
        return NEG_INF
    if x == NEG_INF or x == INF:
        return x
    return math.nextafter(x, NEG_INF)


def up(x):
    """Redondeo hacia arriba (siguiente float hacia +inf)."""
    if x != x:                      # NaN
        return INF
    if x == NEG_INF or x == INF:
        return x
    return math.nextafter(x, INF)


def from_fraction(q):
    """Encierra un racional exacto en un intervalo de floats (B6)."""
    q = Fraction(q)
    f = float(q)
    if Fraction(f) == q:
        return (f, f)
    return (dn(f), up(f))


def hull(*intervalos):
    lo = min(i[0] for i in intervalos)
    hi = max(i[1] for i in intervalos)
    return (lo, hi)


def contiene(intervalo, valor, tol=1e-9):
    lo, hi = intervalo
    margen = tol * max(1.0, abs(lo) if math.isfinite(lo) else 0.0,
                       abs(hi) if math.isfinite(hi) else 0.0)
    return (lo - margen) <= valor <= (hi + margen)


def es_finito(intervalo):
    return math.isfinite(intervalo[0]) and math.isfinite(intervalo[1])


# ------------------------------------------------------------- operaciones
def suma(I, J):
    return (dn(I[0] + J[0]), up(I[1] + J[1]))


def suma_escalar(I, c):
    return (dn(I[0] + c), up(I[1] + c))


def _prod_seguro(a, b):
    """a*b tratando 0*inf como 0 (convencion de IA para coeficientes exactos)."""
    if a == 0.0 or b == 0.0:
        return 0.0
    return a * b


def mult(I, J):
    ps = [_prod_seguro(I[0], J[0]), _prod_seguro(I[0], J[1]),
          _prod_seguro(I[1], J[0]), _prod_seguro(I[1], J[1])]
    if any(p != p for p in ps):          # NaN -> no sabemos acotar
        return ENTERO
    return (dn(min(ps)), up(max(ps)))


def mult_escalar(I, c):
    """c * I con c float exacto (o racional ya convertido)."""
    if c == 0.0:
        return (0.0, 0.0)
    a, b = _prod_seguro(I[0], c), _prod_seguro(I[1], c)
    if a != a or b != b:
        return ENTERO
    return (dn(min(a, b)), up(max(a, b)))


def reciproco(I):
    """1/I.  Si 0 pertenece al interior de I el resultado real es no acotado (B1).

    Casos:
        0 no pertenece a I          -> [1/hi, 1/lo]
        I = [0, b], b > 0           -> [1/b, +inf)
        I = [a, 0], a < 0           -> (-inf, 1/a]
        0 en el interior de I       -> (-inf, +inf)   (envoltura de la union)
    """
    a, b = I
    if a > 0.0 or b < 0.0:
        return (dn(1.0 / b), up(1.0 / a))
    if a == 0.0 and b > 0.0:
        return (dn(1.0 / b), INF)
    if b == 0.0 and a < 0.0:
        return (NEG_INF, up(1.0 / a))
    return ENTERO


def division(I, J):
    return mult(I, reciproco(J))


def potencia_entera(I, n):
    """I**n con n entero, incluyendo n negativo (B1)."""
    a, b = I
    if n == 0:
        return (1.0, 1.0)
    if n < 0:
        return reciproco(potencia_entera(I, -n))
    try:
        an, bn = a ** n, b ** n
    except (OverflowError, ValueError):
        return ENTERO
    if n % 2 == 1:                        # impar: monotona creciente
        return (dn(an), up(bn))
    # par
    if a <= 0.0 <= b:
        return (0.0, up(max(an, bn)))
    return (dn(min(an, bn)), up(max(an, bn)))


def potencia(I, J):
    """I**J en el caso general (exponente no constante): no se acota."""
    return ENTERO


_SIEMPRE_ACOTADAS = {'sin': (-1.0, 1.0), 'cos': (-1.0, 1.0)}


def funcion(nombre, args):
    """Evaluacion por intervalos de funciones elementales.  Sound por defecto."""
    if nombre in _SIEMPRE_ACOTADAS:
        return _SIEMPRE_ACOTADAS[nombre]
    I = args[0] if args else ENTERO
    a, b = I
    try:
        if nombre == 'exp':
            return (dn(math.exp(a)) if a > -745 else 0.0,
                    up(math.exp(b)) if b < 709 else INF)
        if nombre == 'Abs':
            if a <= 0.0 <= b:
                return (0.0, up(max(abs(a), abs(b))))
            return (dn(min(abs(a), abs(b))), up(max(abs(a), abs(b))))
        if nombre == 'sqrt':
            if b < 0.0:
                return ENTERO                 # sin parte real
            lo = math.sqrt(max(a, 0.0))
            return (dn(lo), up(math.sqrt(b)) if math.isfinite(b) else INF)
        if nombre == 'log':
            if b <= 0.0:
                return ENTERO
            lo = NEG_INF if a <= 0.0 else dn(math.log(a))
            return (lo, up(math.log(b)) if math.isfinite(b) else INF)
        if nombre == 'atan':
            return (dn(math.atan(a)) if math.isfinite(a) else -math.pi / 2,
                    up(math.atan(b)) if math.isfinite(b) else math.pi / 2)
        if nombre == 'tanh':
            return (dn(math.tanh(a)), up(math.tanh(b)))
    except (OverflowError, ValueError):
        return ENTERO
    return ENTERO
