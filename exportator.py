import numpy as np
import sympy as sp
import re
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def parse_domain_string(dom_str):
    """Extrae los valores numéricos [min, max] de un string como '[-15, 15]'"""
    numeros = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', dom_str)
    if len(numeros) == 2:
        return [float(numeros[0]), float(numeros[1])]
    return [-1e8, 1e8]

def calcular_dominios_w(A, y_vars, b, domains_dict):
    """Calcula límites de las sumas (w_i) usando aritmética lineal"""
    w_domains = {}
    limites_y = [parse_domain_string(domains_dict.get(var, "[-1e8, 1e8]")) for var in y_vars]
            
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
                min_w += coef * max_y 
                max_w += coef * min_y 
                
        w_domains[f"w_{i}"] = [min_w, max_w]
        
    return w_domains

def calcular_dominios_v(v_equations, known_domains):
    """
    Calcula límites de variables no-lineales (v_i) usando Aritmética de Intervalos.
    known_domains ya incluye los dominios de x_i y w_i.
    """
    v_domains = {}
    
    for v_sym, expr in v_equations:
        def eval_interval(ex):
            # 1. Caso Base: Es un número
            if isinstance(ex, sp.Number):
                val = float(ex)
                return [val, val]
                
            # 2. Caso Base: Es una variable ya conocida (x, y, z, w_i, o un v_i anterior)
            if isinstance(ex, sp.Symbol):
                return known_domains.get(str(ex), [-1e8, 1e8])
                
            # 3. Multiplicación (Evaluamos las 4 esquinas del rectángulo)
            if isinstance(ex, sp.Mul):
                args_intervals = [eval_interval(arg) for arg in ex.args]
                res = args_intervals[0]
                for i in range(1, len(args_intervals)):
                    I1, I2 = res, args_intervals[i]
                    p = [I1[0]*I2[0], I1[0]*I2[1], I1[1]*I2[0], I1[1]*I2[1]]
                    res = [min(p), max(p)]
                return res
                
            # 4. Potencia (Manejo especial para exponentes pares/impares)
            if isinstance(ex, sp.Pow):
                base_i = eval_interval(ex.base)
                exp_i = eval_interval(ex.exp)
                # Si el exponente es un número entero exacto:
                if exp_i[0] == exp_i[1] and int(exp_i[0]) == exp_i[0]:
                    n = int(exp_i[0])
                    if n % 2 != 0: # Impar
                        return [base_i[0]**n, base_i[1]**n]
                    else: # Par (El mínimo podría ser 0 si cruza el eje)
                        if base_i[0] <= 0 <= base_i[1]:
                            return [0, max(base_i[0]**n, base_i[1]**n)]
                        else:
                            return [min(base_i[0]**n, base_i[1]**n), max(base_i[0]**n, base_i[1]**n)]
                            
            # 5. Funciones Trigonométricas (Siempre oscilan entre -1 y 1)
            if isinstance(ex, (sp.sin, sp.cos)):
                return [-1.0, 1.0]
                
            # Fallback de seguridad
            return [-1e8, 1e8]
            
        # Calculamos el dominio de la ecuación v_i actual
        dom = eval_interval(expr)
        v_domains[str(v_sym)] = dom
        
        # Lo inyectamos en known_domains porque un futuro v_k podría usar este v_i
        known_domains[str(v_sym)] = dom 
        
    return v_domains

def exportar_a_txt(A, y_vars, b, shell_eqs_strs, v_equations, vars_dict, domains_dict, w_domains, v_domains, output_txt_name):
    y_syms = [sp.sympify(y) for y in y_vars]
    num_w = A.shape[0]
    
    with open(output_txt_name, 'w') as f:
        f.write("variables\n")
        
        # 1. Variables originales
        for var_name in vars_dict.keys():
            f.write(f"  {var_name} in {domains_dict.get(var_name, '[-1e8, 1e8]')};\n")
            
        # 2. Variables W con dominios matemáticos reales
        for i in range(num_w):
            min_w, max_w = w_domains[f"w_{i}"]
            f.write(f"  w_{i} in [{min_w}, {max_w}];\n")
            
        # 3. Variables V con dominios no-lineales calculados
        for v_sym, _ in v_equations:
            min_v, max_v = v_domains[str(v_sym)]
            f.write(f"  {v_sym} in [{min_v}, {max_v}];\n")
            
        f.write("\nconstraints\n")
        
        # 4. Subsistema W Lineal
        for i in range(num_w):
            expr_lineal = float(b[i])
            for j in range(len(y_vars)):
                coeff = A[i, j]
                if coeff != 0:
                    expr_lineal += coeff * y_syms[j]
            expr_str = str(expr_lineal).replace('**', '^').replace(' ', '')
            f.write(f"  w_{i}={expr_str};\n")
            
        # 5. Subsistema V No-Lineal
        for v_sym, expr in v_equations:
            expr_str = str(expr).replace('**', '^').replace(' ', '')
            f.write(f"  {v_sym}={expr_str};\n")
            
        # 6. Cascarones
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
    A, y_vars, b, shell_equations, v_equations = extract_linear_system(constr_list)
    shell_eqs_strs = np.array([f"{str(lhs)} {op} {str(rhs)}" for lhs, op, rhs in shell_equations])
    v_eqs_strs = np.array([f"{str(sym)}={str(expr)}" for sym, expr in v_equations])
    
    print("3. Pre-calculando dominios reales con Aritmética de Intervalos...")
    # Consolidamos todos los dominios conocidos en un solo diccionario maestro
    known_domains = {var: parse_domain_string(dom) for var, dom in domains_dict.items()}
    
    # Dominios de sumas W
    w_domains = calcular_dominios_w(A, y_vars, b, domains_dict)
    known_domains.update(w_domains)
    
    # Dominios de multiplicaciones V
    v_domains = calcular_dominios_v(v_equations, known_domains)
    
    print(f"4. Guardando datos numéricos en: {output_npz_name}...")
    np.savez(output_npz_name, A=A, y=y_vars, b=b, shell_eqs=shell_eqs_strs, v_eqs=v_eqs_strs)
    
    print(f"5. Generando reporte legible en: {output_txt_name}...")
    exportar_a_txt(A, y_vars, b, shell_eqs_strs, v_equations, vars_dict, domains_dict, w_domains, v_domains, output_txt_name)
    
    print("\n¡Pre-procesamiento completado con doble exportación!\n")