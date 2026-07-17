import os
import glob

# Importamos los módulos de nuestra arquitectura de pre-procesamiento
# Corregido: Importamos desde 'exportator' que contiene la función actualizada de 3 argumentos
from exportator import procesar_instancia
from validator import validar_sistema

def procesar_carpeta(carpeta_entrada, carpeta_salida):
    """
    Lee todos los .txt en la carpeta de entrada, ejecuta el pipeline completo 
    (pre-procesamiento y validación) y guarda los .npz y .txt en la salida.
    """
    if not os.path.exists(carpeta_salida):
        os.makedirs(carpeta_salida)

    # Buscamos todos los archivos .txt dentro de la carpeta 'instances/'
    archivos_txt = glob.glob(os.path.join(carpeta_entrada, '*.txt'))
    
    if not archivos_txt:
        print(f"[!] No se encontraron archivos .txt en la carpeta '{carpeta_entrada}'.")
        return

    print(f"===========================================================")
    print(f" INICIANDO PROCESAMIENTO MASIVO: {len(archivos_txt)} instancias detectadas")
    print(f"===========================================================\n")

    exitosos = 0
    fallidos = 0

    # Iteramos sobre cada instancia encontrada
    for ruta_instancia in archivos_txt:
        # Extraemos el nombre base sin la extensión ni la carpeta (ej. 'instprofe')
        nombre_base = os.path.splitext(os.path.basename(ruta_instancia))[0]
        
        # Generamos las rutas de salida automáticas con el mismo nombre base
        ruta_npz = os.path.join(carpeta_salida, f"{nombre_base}.npz")
        ruta_txt_out = os.path.join(carpeta_salida, f"{nombre_base}_procesada.txt")
        
        print(f">>> Analizando instancia: {nombre_base}")
        
        try:
            # Paso 1: Extracción, construcción del DAG y creación de NPZ y TXT legible
            procesar_instancia(ruta_instancia, ruta_npz, ruta_txt_out)
            
            # Paso 2: Validación Matemática
            # Importante: enviamos verbose=False para que no sature la consola con cada ecuación
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
    # Usamos tu estructura actual de carpetas
    CARPETA_ENTRADA = 'instances' 
    CARPETA_SALIDA = 'outputs'
    
    procesar_carpeta(CARPETA_ENTRADA, CARPETA_SALIDA)