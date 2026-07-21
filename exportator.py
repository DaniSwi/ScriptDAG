import numpy as np
import sympy as sp
import re
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def calcular_dominios_w(A, y_vars, b, domains_dict):
    """
    Calcula los límites exactos [min, max] para cada variable auxiliar w_i
    basado en la matriz A y los dominios de las variables originales.
    """
    w_domains = []
    
    # 1. Limpiamos y convertimos los dominios originales de texto a números [min, max]
    limites_y = []
    for var in y_vars:
        dom_str = domains_dict.get(var, "[-1e8, 1e8]")
        # Extraemos los números dentro de los corchetes
        numeros = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', dom_str)
        if len(numeros) == 2:
            limites_y.append([float(numeros[0]), float(numeros[1])])
        else:
            limites_y.append([-1e8, 1e8]) # Fallback seguro
            
    # 2. Aplicamos Aritmética de Intervalos fila por fila
    for i in range(A.shape[0]):
        min_w = float(b[i])
        max_w = float(b[i])
        
        for j in range(len(y_vars)):
            coef = A[i, j]
            min_y, max_y = limites_y[j]
            
            if coef > 0:
                min_w += coef * min_y
                max_w += coef * max_y
            elif coef < 0:
                min_w += coef * max_y # Si el coeficiente es negativo, el max_y aporta al mínimo
                max_w += coef * min_y # y el min_y aporta al máximo
                
        w_domains.append(f"[{min_w}, {max_w}]")
        
    return w_domains

def exportar_a_txt(A, y_vars, b, shell_eqs_strs, vars_dict, domains_dict, w_domains, output_txt_name):
    """Genera un archivo .txt con el formato original estricto."""
    y_syms = [sp.sympify(y) for y in y_vars]
    num_w = A.shape[0]
    
    with open(output_txt_name, 'w') as f:
        f.write("variables\n")
        
        for var_name in vars_dict.keys():
            dominio_original = domains_dict.get(var_name, "[-1e8, 1e8]")
            f.write(f"  {var_name} in {dominio_original};\n")
            
        # NUEVO: Imprimimos los dominios calculados en lugar de [-1e8, 1e8]
        for i in range(num_w):
            f.write(f"  w_{i} in {w_domains[i]};\n")
            
        f.write("\nconstraints\n")
        
        for i in range(num_w):
            expr_lineal = float(b[i])
            for j in range(len(y_vars)):
                coeff = A[i, j]
                if coeff != 0:
                    expr_lineal += coeff * y_syms[j]
                    
            expr_str = str(expr_lineal).replace('**', '^').replace(' ', '')
            f.write(f"  w_{i}={expr_str};\n")
            
        for eq_str in shell_eqs_strs:
            eq_format = str(eq_str).replace('**', '^')
            f.write(f"  {eq_format};\n")
            
        f.write("end\n")

def procesar_instancia(file_path, output_npz_name, output_txt_name):
    print(f"--- Iniciando Pre-procesamiento de: {file_path} ---")
    dag_builder.clear_cache()
    
    print("1. Cargando archivo y construyendo árboles en SymPy...")
    vars_dict, domains_dict, constr_list = parse_instance(file_path)
    
    print("2. Extrayendo componentes lineales y no lineales...")
    A, y_vars, b, shell_equations = extract_linear_system(constr_list)
    shell_eqs_strs = np.array([f"{str(lhs)} {op} {str(rhs)}" for lhs, op, rhs in shell_equations])
    
    # NUEVA CAPA: Calculamos los dominios reales de las 'w'
    print("3. Pre-calculando dominios reales con Aritmética de Intervalos...")
    w_domains = calcular_dominios_w(A, y_vars, b, domains_dict)
    
    print(f"4. Guardando datos numéricos en: {output_npz_name}...")
    np.savez(output_npz_name, A=A, y=y_vars, b=b, shell_eqs=shell_eqs_strs) # Aquí también podríamos guardar los dominios de W si la IA/Solver los necesita
    
    print(f"5. Generando reporte legible en: {output_txt_name}...")
    # Pasamos la nueva lista de dominios a la función de exportación
    exportar_a_txt(A, y_vars, b, shell_eqs_strs, vars_dict, domains_dict, w_domains, output_txt_name)
    
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