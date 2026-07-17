import numpy as np
import sympy as sp
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def exportar_a_txt(A, y_vars, b, shell_eqs_strs, vars_dict, domains_dict, output_txt_name):
    """
    Genera un archivo .txt legible por humanos con la instancia procesada.
    Incluye las variables originales, las auxiliares w_i, y los subsistemas extraídos.
    """
    y_syms = [sp.sympify(y) for y in y_vars]
    num_w = A.shape[0]
    
    with open(output_txt_name, 'w') as f:
        f.write("variables\n")
        
        # 1. Escribimos las variables originales
        for var_name in vars_dict.keys():
            dominio_original = domains_dict.get(var_name, "[-1e8, 1e8]")
            f.write(f"  {var_name} in {dominio_original};\n")
            
        # 2. Escribimos las variables auxiliares w
        for i in range(num_w):
            f.write(f"  w_{i} in [-1e8, 1e8];\n")
            
        f.write("\nconstraints\n")
        
        # 3. Reconstruimos y escribimos las ecuaciones lineales (w = A*y + b)
        for i in range(num_w):
            expr_lineal = float(b[i])
            for j in range(len(y_vars)):
                coeff = A[i, j]
                if coeff != 0:
                    expr_lineal += coeff * y_syms[j]
                    
            # SymPy formatea la matemática perfecto, solo cambiamos ** por ^
            expr_str = str(expr_lineal).replace('**', '^').replace(' ', '')
            f.write(f"  w_{i}={expr_str};\n")
            
        # 4. Escribimos los cascarones (Estructura no lineal)
        for eq_str in shell_eqs_strs:
            eq_format = str(eq_str).replace('**', '^')
            f.write(f"  {eq_format};\n")
            
        f.write("end\n")

def procesar_instancia(file_path, output_npz_name, output_txt_name):
    print(f"--- Iniciando Pre-procesamiento de: {file_path} ---")
    
    # 1. Limpiar memoria del caché del DAG
    dag_builder.clear_cache()
    
    # 2. Capa de Parseo
    print("1. Cargando archivo y construyendo árboles en SymPy...")
    vars_dict, domains_dict, constr_list = parse_instance(file_path)
    
    # 3. Capa de Construcción del DAG y Extracción
    print("2. Extrayendo componentes lineales y no lineales...")
    # *Nota: Asumimos que matrix_extractor ya devuelve (A, y_names, b, shell_equations)
    A, y_vars, b, shell_equations = extract_linear_system(constr_list)
    
    # Convertimos los cascarones a texto (incluyendo los operadores >=, <=, etc.)
    shell_eqs_strs = np.array([f"{str(lhs)} {op} {str(rhs)}" for lhs, op, rhs in shell_equations])
    
    # 4. Capa de Exportación Computacional (.npz)
    print(f"3. Guardando datos numéricos en: {output_npz_name}...")
    np.savez(output_npz_name, A=A, y=y_vars, b=b, shell_eqs=shell_eqs_strs)
    
    # 5. Capa de Exportación Legible (.txt)
    print(f"4. Generando reporte legible en: {output_txt_name}...")
    exportar_a_txt(A, y_vars, b, shell_eqs_strs, vars_dict, domains_dict, output_txt_name)
    
    print("\n¡Pre-procesamiento completado con doble exportación!\n")

if __name__ == "__main__":
    #aca se coloca la ruta del archivo 
    ruta_instancia = 'instances/instprofe.txt'
    ruta_salida_npz = 'outputs/subsistema_n_ario_profe.npz' #cambiar nombres por acomodo
    ruta_salida_txt = 'outputs/subsistema_n_ario_profe.txt' #ruta para el archivo legible
    
    procesar_instancia(ruta_instancia, ruta_salida_npz, ruta_salida_txt)

    #la instancia que entrega? un .npz que contiene:
    """
    datos['A'] = matriz de coeficientes
    datos['y'] = vector de variables ['x22', 'x32'] etc.
    datos['b'] = vector de resultados 
    datos['shell_eqs'] = las ecuaciones originales reemplazadas
    """