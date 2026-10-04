"""
Punto de entrada.  Uso:

    python main.py instances/inst001.txt            # una instancia
    python main.py instances/inst001.txt --txt      # ademas el .txt reescrito
    python main.py --lote                           # toda la carpeta instances/
    python main.py --lote --nsc                     # ademas REF vs NSC (Sec. 3.3)

Se elimino la funcion asd(), que estaba rota (B4: desempaquetaba 2 y 4 valores
de funciones que devuelven 3 y 5) y era codigo muerto.

Que contiene el .npz (D5):  A, c, A_aug=[A|-I], A exacta en COO racional,
y_names/y_lo/y_hi, w_names/w_lo/w_hi, var_names/var_lo/var_hi, el subsistema
de productos y las ecuaciones cascaron.  Con eso se arma  A*y + c = w  sin
tener que re-parsear ningun string:

    import numpy as np
    z = np.load('outputs/inst001.npz')
    A, c, y, w = z['A'], z['c'], z['y_names'], z['w_names']
"""

import os
import sys

from batch_processor import procesar_carpeta
from nsc import comparar
from exportator import construir, exportar_npz, exportar_txt
from validator import validar


def procesar_uno(ruta, carpeta_salida='outputs', exportar_texto=False,
                 calcular_dominios=True, correr_nsc=False):
    os.makedirs(carpeta_salida, exist_ok=True)
    base = os.path.splitext(os.path.basename(ruta))[0]
    print(f"--- Pre-procesamiento de {ruta} ---")
    inst, sis, dom = construir(ruta, calcular_dominios)
    print(f"  DAG      : {inst.dag.resumen()}")
    print(f"  sistema  : {sis.resumen()}")
    print(f"  A        : {sis.A.shape[0]} filas (w) x {sis.A.shape[1]} columnas (y)")
    npz = os.path.join(carpeta_salida, f"{base}.npz")
    exportar_npz(npz, inst, sis, dom)
    print(f"  escrito  : {npz}")
    if exportar_texto:
        txt = os.path.join(carpeta_salida, f"{base}_reescrita.txt")
        exportar_txt(txt, inst, sis, dom)
        print(f"  escrito  : {txt}")
    validar(inst, sis, dom)
    if correr_nsc:
        comparar(inst, sis)
    return inst, sis, dom


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = {a for a in sys.argv[1:] if a.startswith('--')}
    if '--lote' in flags or not args:
        procesar_carpeta('instances', 'outputs',
                         calcular_dominios='--sin-dominios' not in flags,
                         exportar_texto='--txt' in flags,
                         validar_todo='--sin-validar' not in flags,
                         correr_nsc='--nsc' in flags)
    else:
        for ruta in args:
            procesar_uno(ruta, exportar_texto='--txt' in flags,
                         calcular_dominios='--sin-dominios' not in flags,
                         correr_nsc='--nsc' in flags)
