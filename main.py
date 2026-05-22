import numpy as np
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def procesar_instancia(file_path, output_npz_name):
    print(f"--- Iniciando Pre-procesamiento de: {file_path} ---")
    
    # 1. Limpiar memoria del caché del DAG
    # Esto es crucial para garantizar consistencia si procesas múltiples archivos en un bucle
    dag_builder.clear_cache()
    
    # 2. Capa de Parseo
    print("1. Cargando archivo y construyendo árboles matemáticos en SymPy...")
    vars_dict, constr_list = parse_instance(file_path)
    
    # 3. Capa de Construcción del DAG y Extracción
    print("2. Ejecutando análisis del DAG y extrayendo componentes lineales y no lineales...")
    # Ahora recibimos los 4 elementos, incluyendo las ecuaciones cascarón
    A, y_vars, b, shell_equations = extract_linear_system(constr_list)
    
    print(f"   [INFO] Dimensiones de la Matriz A (Coeficientes lineales): {A.shape}")
    print(f"   [INFO] Cantidad de términos variables en vector 'y': {len(y_vars)}")
    print(f"   [INFO] Cantidad de ecuaciones cascarón estructurales 'w': {len(shell_equations)}")
    
    # Convertimos las ecuaciones simbólicas de SymPy a strings estándar.
    # Esto nos permite serializarlas de forma masiva y segura dentro de la arquitectura de NumPy.
    shell_eqs_strs = np.array([f"{str(lhs)} = {str(rhs)}" for lhs, rhs in shell_equations])
    
    # 4. Capa de Exportación y Almacenamiento
    print(f"3. Empaquetando y guardando datos consolidados en: {output_npz_name}...")
    # Agregamos 'shell_eqs' como un arreglo de texto dentro del mismo contenedor comprimido .npz
    np.savez(output_npz_name, A=A, y=y_vars, b=b, shell_eqs=shell_eqs_strs)
    
    # Vista previa informativa en la consola de comandos
    print("\n======================================================================")
    print("   VISTA PREVIA DE LAS ECUACIONES CASCARÓN (ESTRUCTURA NO LINEAL)")
    print("======================================================================")
    for idx, eq in enumerate(shell_eqs_strs):
        print(f" Restricción {idx}: {eq}")
    print("======================================================================")
    print("¡Pre-procesamiento completado! El archivo está listo para ser auditado.\n")

if __name__ == "__main__":
    #aca se coloca la ruta del archivo 
    ruta_instancia = 'instances/inst001.txt'
    ruta_salida = 'outputs/subsistema_n_ario_test.npz' #cambiar nombres por acomodo
    
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