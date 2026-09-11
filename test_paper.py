"""
Pruebas de regresion contra la tecnica del paper (Araya & Reyes 2019, Sec. 3.2).

NOTA: el sistema de abajo es el ejemplo (1) del paper RECONSTRUIDO a partir de
la forma de la matriz reportada en la revision (5 filas, la fila de eq3 con
coeficientes 4|1|1|-6|...).  Si tienes el paper a mano, reemplaza EJEMPLO_1 por
las tres ecuaciones literales: las aserciones no dependen de los nombres, solo
de las propiedades que el paper exige.

    python test_paper.py
"""

from fractions import Fraction

from dag_builder import DagBuilder
from matrix_extractor import extraer_sistema_lineal

EJEMPLO_1 = [
    "4*x2 + x1*x2 + x3^2",                         # eq1
    "x1*x2 + x3^2 + 2*x2*sin(x1+x3)",              # eq2
    "4*x2 + x1*x2 + x3^2 - 6*(x1+x3) + sin(6*(x1+x3))",   # eq3
]
VARS = {v: None for v in ('x1', 'x2', 'x3')}


def construir(exprs, distribuir=False):
    dag = DagBuilder(distribuir_sumas_escaladas=distribuir)
    for e in exprs:
        c, n = dag.desde_string(e, VARS)
        dag.registrar_raiz(n)
    return dag, extraer_sistema_lineal(dag)


def mostrar(sis):
    print("  columnas y:", sis.y_names)
    for i, nombre in enumerate(sis.w_names):
        fila = {sis.y_names[j]: str(q) for j, q in sorted(sis.A_exacta[i].items())}
        print(f"  {nombre} = {sis.c_exacta[i]} + {fila}")


def test_coeficientes_en_la_matriz():
    """D1: los coeficientes 4 y 2 deben estar DENTRO de A, no escondidos."""
    dag, sis = construir(EJEMPLO_1)
    coefs = {q for fila in sis.A_exacta for q in fila.values()}
    assert Fraction(4) in coefs, "el coeficiente 4 de 4*x2 se perdio (bug D1)"
    assert Fraction(2) in coefs, "el coeficiente 2 de 2*x2*sin(..) se perdio (bug D1)"
    assert Fraction(-6) in coefs, "el coeficiente -6 se perdio (bug D1)"
    print("[PASS] D1: los coeficientes 4, 2 y -6 estan en A")


def test_sin_variables_auxiliares():
    """D4: el paper no modifica el sistema original; cero variables nuevas."""
    dag, sis = construir(EJEMPLO_1)
    assert not [n for n in sis.y_names if n.startswith('v_')], \
        "se crearon variables auxiliares v_i (bug D1/D4)"
    print(f"[PASS] D4: 0 variables auxiliares nuevas ({len(sis.y_names)} columnas)")


def test_una_columna_por_nodo():
    """D2: el mismo nodo no puede ocupar dos columnas."""
    dag, sis = construir(EJEMPLO_1)
    assert len(sis.y_names) == len(set(sis.y_names)), "hay columnas duplicadas (bug D2)"
    assert len(sis.col_de_nodo) == len(sis.y_nodes)
    print("[PASS] D2: una columna por nodo, sin duplicados")


def test_cse_de_la_suma_compartida():
    """D3: (x1+x3) debe seguir siendo UN nodo compartido por sin(x1+x3) y
    por -6*(x1+x3); SymPy con evaluate=True lo distribuye y lo destruye."""
    dag, sis = construir(EJEMPLO_1)
    compartidas = [n for n in dag.nodes if n.kind == 'add' and n.in_degree > 1]
    assert compartidas, "ninguna suma quedo compartida: se perdio la CSE (bug D3)"
    n = compartidas[0]
    print(f"[PASS] D3: la suma '{n.to_str()}' es un nodo compartido "
          f"(in_degree={n.in_degree}) y aparece como columna mixta")


def test_modo_distribuido_reproduce_la_referencia():
    """Con distribuir_sumas_escaladas=True, -6*(x1+x3) se reparte en -6,-6,
    que es la forma citada en 5.1 de la revision."""
    dag, sis = construir(EJEMPLO_1, distribuir=True)
    ultima = sis.A_exacta[-1]
    vals = sorted(str(q) for q in ultima.values())
    assert vals.count('-6') == 2, f"esperaba dos -6, obtuve {vals}"
    print(f"[PASS] modo distribuido: la fila de eq3 tiene {vals}")


if __name__ == "__main__":
    dag, sis = construir(EJEMPLO_1)
    print(f"A = {sis.A.shape[0]} filas x {sis.A.shape[1]} columnas")
    mostrar(sis)
    print()
    test_coeficientes_en_la_matriz()
    test_sin_variables_auxiliares()
    test_una_columna_por_nodo()
    test_cse_de_la_suma_compartida()
    test_modo_distribuido_reproduce_la_referencia()
    print("\nTodo OK.")
