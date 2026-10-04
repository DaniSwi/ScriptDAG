"""
Capa 7: procesamiento masivo de una carpeta de instancias.

Cambios: la validacion vuelve a estar ACTIVADA (antes estaba comentada porque
sp.simplify se colgaba en instprofe; el validator nuevo no usa simplify).
Ademas se imprime el resumen del sistema de cada instancia, que es lo que
permite ver de un vistazo si la matriz es la del paper o la de incidencia.
"""

import glob
import os
import time
import traceback

from exportator import construir, exportar_npz, exportar_txt
from nsc import comparar
from validator import validar


def procesar_carpeta(carpeta_entrada, carpeta_salida, calcular_dominios=True,
                     exportar_texto=False, validar_todo=True,
                     distribuir_sumas_escaladas=False, correr_nsc=False):
    os.makedirs(carpeta_salida, exist_ok=True)
    archivos = sorted(glob.glob(os.path.join(carpeta_entrada, '*.txt')))
    if not archivos:
        print(f"[!] No hay .txt en '{carpeta_entrada}'.")
        return

    print("=" * 72)
    print(f" PROCESAMIENTO MASIVO: {len(archivos)} instancias")
    print(f" dominios={'calculados (sound)' if calcular_dominios else 'omitidos (-inf,inf)'}"
          f" | txt reescrito={'si' if exportar_texto else 'no'}"
          f" | validacion={'si' if validar_todo else 'no'}")
    print("=" * 72)

    ok = fallidos = 0
    for ruta in archivos:
        base = os.path.splitext(os.path.basename(ruta))[0]
        print(f"\n>>> {base}")
        t0 = time.time()
        try:
            inst, sis, dom = construir(ruta, calcular_dominios,
                                       distribuir_sumas_escaladas)
            exportar_npz(os.path.join(carpeta_salida, f"{base}.npz"), inst, sis, dom)
            if exportar_texto:
                exportar_txt(os.path.join(carpeta_salida, f"{base}_reescrita.txt"),
                             inst, sis, dom)
            print(f"    DAG     : {inst.dag.resumen()}")
            print(f"    sistema : {sis.resumen()}")
            print(f"    tiempo  : {time.time() - t0:.2f}s")
            if validar_todo:
                r = validar(inst, sis, dom, verbose=False)
                for nombre, bien, det in r.checks:
                    print(f"    [{'PASS' if bien else 'FAIL'}] {nombre}  {det}")
                if not r.ok:
                    raise RuntimeError("la auditoria fallo")
            if correr_nsc:
                comparar(inst, sis)
            ok += 1
        except Exception as e:
            fallidos += 1
            print(f"    [X] ERROR: {e}")
            traceback.print_exc()

    print("\n" + "=" * 72)
    print(f" RESUMEN: {ok} exitosos | {fallidos} fallidos")
    print("=" * 72)


if __name__ == "__main__":
    procesar_carpeta('instances', 'outputs',
                     calcular_dominios=True,
                     exportar_texto=False,
                     validar_todo=True)
