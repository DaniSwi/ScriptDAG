"""
Capa 4 (NUEVA): propagacion de dominios por el DAG, con intervalos sound.

Arregla B1 y B2: antes los dominios de w_i/v_i se declaraban [-1e8, 1e8]
aunque el rango real fuera [-2.6e9, 2.6e9] (recorte => perdida de soluciones),
y las potencias negativas devolvian intervalos que no contenian a la expresion.

Aqui NO se recorta nada: si no se sabe acotar, el dominio es (-inf, +inf).
El unico compromiso permitido es sobre-aproximar, nunca sub-aproximar.
"""

import math
import random
from fractions import Fraction

import intervals as iv


def propagar_dominios(dag, var_domains):
    """Devuelve {node.id: (lo, hi)} para TODOS los nodos del DAG.

    Los hijos siempre tienen id menor que el padre (hash-consing), asi que
    recorrer dag.nodes en orden es un orden topologico valido.
    """
    dom = {}
    for nodo in dag.nodes:
        if nodo.kind == 'var':
            dom[nodo.id] = tuple(var_domains.get(nodo.data, iv.ENTERO))
        else:
            dom[nodo.id] = evaluar_nodo(nodo, dom)
    return dom


def evaluar_nodo(nodo, dom):
    """Fase forward de UN nodo: su intervalo a partir de los dominios de sus
    hijos.  La usan propagar_dominios y el HC4-revise de nsc.py."""
    k = nodo.kind
    if k == 'var':
        return dom[nodo.id]
    if k == 'const':
        return iv.from_fraction(nodo.data)
    if k == 'add':
        const, coefs = nodo.data
        acum = iv.from_fraction(const)
        for c, hijo in zip(coefs, nodo.args):
            acum = iv.suma(acum, iv.mult(iv.from_fraction(c), dom[hijo.id]))
        return acum
    if k == 'mul':
        acum = (1.0, 1.0)
        for hijo in nodo.args:
            acum = iv.mult(acum, dom[hijo.id])
        return acum
    if k == 'pow':
        e = Fraction(nodo.data)
        base = dom[nodo.args[0].id]
        if e.denominator == 1:
            return iv.potencia_entera(base, int(e))
        if base[0] >= 0:
            return iv.funcion('sqrt', [base]) if e == Fraction(1, 2) else iv.ENTERO
        return iv.ENTERO
    if k == 'func':
        return iv.funcion(nodo.data, [dom[a.id] for a in nodo.args])
    return iv.ENTERO


def evaluar_en_punto(dag, punto):
    """Evalua todos los nodos en un punto concreto {nombre_var: float}.

    Sirve para la validacion Monte Carlo: comprobar que el valor real de cada
    nodo cae dentro del dominio declarado.  Devuelve {node.id: float} (o None
    si la evaluacion no esta definida, p.ej. 1/0).
    """
    val = {}
    for nodo in dag.nodes:
        k = nodo.kind
        try:
            if k == 'var':
                val[nodo.id] = float(punto[nodo.data])
            elif k == 'const':
                val[nodo.id] = float(nodo.data)
            elif k == 'add':
                const, coefs = nodo.data
                s = float(const)
                for c, h in zip(coefs, nodo.args):
                    s += float(c) * val[h.id]
                val[nodo.id] = s
            elif k == 'mul':
                p = 1.0
                for h in nodo.args:
                    p *= val[h.id]
                val[nodo.id] = p
            elif k == 'pow':
                val[nodo.id] = float(val[nodo.args[0].id]) ** float(nodo.data)
            elif k == 'func':
                f = getattr(math, nodo.data, None)
                if nodo.data == 'Abs':
                    f = abs
                if f is None:
                    val[nodo.id] = None
                else:
                    val[nodo.id] = f(*[val[a.id] for a in nodo.args])
            else:
                val[nodo.id] = None
        except (ZeroDivisionError, ValueError, OverflowError, TypeError):
            val[nodo.id] = None
    return val


def muestrear(var_domains, rng, cota=1e3):
    """Punto aleatorio dentro de las cajas.

    NO es uniforme a proposito: el muestreo uniforme casi nunca produce
    denominadores chicos, que es justo donde vive el bug B1.  Se mezclan tres
    regimenes: uniforme, log-uniforme (magnitudes chicas) y extremos/ceros.
    """
    punto = {}
    for nombre, (lo, hi) in var_domains.items():
        lo, hi = max(lo, -cota), min(hi, cota)
        if lo > hi:
            lo = hi = 0.0
        modo = rng.random()
        if modo < 0.34:
            punto[nombre] = rng.uniform(lo, hi)
        elif modo < 0.67:
            mag = 10 ** rng.uniform(-9, math.log10(max(abs(lo), abs(hi), 1e-9)))
            v = mag * rng.choice((-1.0, 1.0))
            punto[nombre] = min(max(v, lo), hi)
        else:
            punto[nombre] = rng.choice((lo, hi, 0.0, (lo + hi) / 2))
            punto[nombre] = min(max(punto[nombre], lo), hi)
    return punto


def sub_caja(var_domains, rng, cota=1e3):
    """Sub-caja aleatoria (posiblemente muy estrecha) dentro de las cajas.

    Sirve para comprobar la ISOTONIA POR INCLUSION, que es la propiedad que
    define un filtrado sound:  X' subset X  =>  F(X') subset F(X).
    Una sola evaluacion puntual casi nunca detecta B1; una sub-caja centrada
    cerca de una raiz del denominador lo detecta de inmediato.
    """
    caja = {}
    for nombre, (lo, hi) in var_domains.items():
        lo, hi = max(lo, -cota), min(hi, cota)
        if lo > hi:
            lo = hi = 0.0
        c = muestrear({nombre: (lo, hi)}, rng, cota)[nombre]
        radio = (hi - lo) * 10 ** rng.uniform(-12, -1)
        caja[nombre] = (max(lo, c - radio), min(hi, c + radio))
    return caja


def valor_local(nodo, vals_hijos):
    """Evalua UN nodo a partir de valores concretos de sus hijos.

    Se usa en la auditoria local (V4c): permite muestrear directamente en los
    dominios de los hijos, sin tener que llegar ahi propagando desde las x.
    Es lo que hace detectable el bug B1: para el nodo (w)^(-1) se muestrean
    valores de w cercanos a 0, que un muestreo desde las x nunca produce.
    """
    k = nodo.kind
    try:
        if k == 'add':
            const, coefs = nodo.data
            return float(const) + sum(float(c) * v for c, v in zip(coefs, vals_hijos))
        if k == 'mul':
            p = 1.0
            for v in vals_hijos:
                p *= v
            return p
        if k == 'pow':
            return float(vals_hijos[0]) ** float(nodo.data)
        if k == 'func':
            f = abs if nodo.data == 'Abs' else getattr(math, nodo.data, None)
            return None if f is None else f(*vals_hijos)
    except (ZeroDivisionError, ValueError, OverflowError, TypeError):
        return None
    return None


def muestra_en_intervalo(intervalo, rng, cota=1e12):
    """Valor aleatorio dentro de un intervalo, sesgado hacia magnitudes chicas
    y hacia los extremos (donde viven los bugs de aritmetica de intervalos)."""
    lo, hi = intervalo
    lo, hi = max(lo, -cota), min(hi, cota)
    if lo > hi:
        return None
    modo = rng.random()
    if modo < 0.3:
        return rng.uniform(lo, hi)
    if modo < 0.7:
        mag = 10 ** rng.uniform(-12, math.log10(max(abs(lo), abs(hi), 1e-12)))
        return min(max(mag * rng.choice((-1.0, 1.0)), lo), hi)
    return rng.choice((lo, hi, 0.0, (lo + hi) / 2.0))
