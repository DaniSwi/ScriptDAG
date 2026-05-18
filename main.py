import numpy as np
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def procesar_instancia(file_path, output_name):
    print(f"--- Procesando: {file_path} ---")
    
    # 1. Limpiar memoria del DAG (útil si corres esto en un bucle para muchos archivos)
    dag_builder.clear_cache()
    
    # 2. Parseo
    print("1. Parseando el archivo...")
    vars_dict, constr_list = parse_instance(file_path)
    
    # 3. Construcción del DAG y Extracción
    print("2. Construyendo DAG y extrayendo matrices...")
    A, y_vars, b = extract_linear_system(constr_list)
    
    print(f"   -> Dimensiones de la Matriz A: {A.shape}")
    print(f"   -> Número de variables auxiliares 'y': {len(y_vars)}")
    
    # 4. Exportación en formato .npz (el formato optimizado para estos casos)
    print(f"3. Guardando resultados en {output_name}...")
    np.savez(output_name, A=A, y=y_vars, b=b)
    print("¡Proceso completado con éxito!\n")

if __name__ == "__main__":
    #aca se coloca la ruta del archivo 
    ruta_instancia = 'instances/inst001.txt'
    ruta_salida = 'outputs/subsistema_n_ario.npz' #cambiar nombres por acomodo
    
    procesar_instancia(ruta_instancia, ruta_salida)

    #la instancia que entrega? un .npz que contiene:
    """
    datos['A'] = matriz de coeficientes
    datos['y'] = vector de variables ['x22', 'x32'] etc.
    datos['b'] = vector de resultados 
    """

    """
    Entonces las personas que quieran usar el .npz para desempaquetar
    los datos, importan numpy y usan el metodo np.load(elnombredelarchivo) y así lo tienen
    """