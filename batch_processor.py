import os
import glob
from exportator import procesar_instancia
from validator import validar_sistema

def procesar_carpeta(carpeta_entrada, carpeta_salida):
    """
    Ejecuta el pipeline completo y valida cada instancia automáticamente.
    """
    if not os.path.exists(carpeta_salida):
        os.makedirs(carpeta_salida)

    archivos_txt = glob.glob(os.path.join(carpeta_entrada, '*.txt'))
    
    if not archivos_txt:
        print(f"[!] No se encontraron archivos .txt en la carpeta '{carpeta_entrada}'.")
        return

    print(f"===========================================================")
    print(f" INICIANDO PROCESAMIENTO MASIVO: {len(archivos_txt)} instancias detectadas")
    print(f"===========================================================\n")

    exitosos = 0
    fallidos = 0

    for ruta_instancia in archivos_txt:
        nombre_base = os.path.splitext(os.path.basename(ruta_instancia))[0]
        
        ruta_npz = os.path.join(carpeta_salida, f"{nombre_base}.npz")
        ruta_txt_out = os.path.join(carpeta_salida, f"{nombre_base}_procesada.txt")
        
        print(f">>> Analizando instancia: {nombre_base}")
        
        try:
            procesar_instancia(ruta_instancia, ruta_npz, ruta_txt_out)
            validar_sistema(ruta_instancia, ruta_npz, verbose=False) 
            exitosos += 1
        except Exception as e:
            print(f"  [X] ERROR CRÍTICO procesando {nombre_base}: {str(e)}")
            fallidos += 1
            
        print("-" * 60)
        
    print("\n===========================================================")
    print(f" RESUMEN DEL LOTE COMPLETADO")
    print(f"   -> Procesados con éxito : {exitosos}")
    print(f"   -> Fallidos / Errores   : {fallidos}")
    print(f"   -> Resultados en        : '{carpeta_salida}'")
    print(f"===========================================================\n")

if __name__ == "__main__":
    CARPETA_ENTRADA = 'instances' 
    CARPETA_SALIDA = 'outputs'
    
    procesar_carpeta(CARPETA_ENTRADA, CARPETA_SALIDA)