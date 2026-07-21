import numpy as np
import sympy as sp
from dag_builder import build_dag_node

def extract_linear_system(constraints):
    """
    Extrae el sistema lineal A * y = b de las sumas n-arias,
    y preserva las ecuaciones originales usando variables auxiliares w_i.
    Incluye Unificación Real de Variables (CSE).
    """
    A_rows = []
    b_vector = []
    y_variables = {}
    
    shell_equations = []
    w_counter = 0

    # NUEVO: Diccionario caché para la Unificación (CSE)
    # Llave = Expresión matemática exacta (ej. x + y + z)
    # Valor = Variable auxiliar asignada (ej. w_0)
    sum_cache = {}

    def replace_and_extract_sums(expr):
        nonlocal w_counter
        
        if isinstance(expr, sp.Add):
            
            # ========================================================
            # 1. LA MAGIA DE LA UNIFICACIÓN: 
            # Si SymPy detecta que esta suma ya se procesó antes, 
            # simplemente devuelve el w_i antiguo y corta la ejecución aquí.
            # ========================================================
            if expr in sum_cache:
                return sum_cache[expr]
            
            # ========================================================
            # 2. Si es una suma totalmente nueva, creamos un nuevo w_i
            # ========================================================
            w_sym = sp.Symbol(f'w_{w_counter}')
            w_counter += 1
            
            row_dict = {}
            constant_term = 0.0
            
            for term in expr.args:
                coeff, var_part = term.as_coeff_Mul()
                
                if var_part == 1:
                    constant_term += float(coeff)
                else:
                    if var_part not in y_variables:
                        y_variables[var_part] = len(y_variables)
                    
                    col_index = y_variables[var_part]
                    row_dict[col_index] = row_dict.get(col_index, 0.0) + float(coeff)
            
            A_rows.append(row_dict)
            b_vector.append(constant_term)
            
            # 3. Guardamos la nueva suma en nuestro caché para el futuro
            sum_cache[expr] = w_sym
            
            return w_sym
            
        elif isinstance(expr, sp.Mul):
            return sp.Mul(*[replace_and_extract_sums(arg) for arg in expr.args])
            
        elif isinstance(expr, sp.Pow):
            return sp.Pow(replace_and_extract_sums(expr.base), replace_and_extract_sums(expr.exp))
            
        elif isinstance(expr, sp.Function):
            return expr.func(*[replace_and_extract_sums(arg) for arg in expr.args])
            
        else:
            return expr

    # LÓGICA PRINCIPAL DEL EXTRACTOR
    for lhs_expr, op, rhs_val in constraints:
        build_dag_node(lhs_expr)
        shell_lhs = replace_and_extract_sums(lhs_expr)
        shell_equations.append((shell_lhs, op, rhs_val))

    num_rows = len(A_rows)
    num_cols = len(y_variables)
    A_matrix = np.zeros((num_rows, num_cols))
    
    for i, row_dict in enumerate(A_rows):
        for j, coeff in row_dict.items():
            A_matrix[i, j] = coeff
            
    y_names = np.array([str(expr) for expr, idx in sorted(y_variables.items(), key=lambda item: item[1])])
    b_array = np.array(b_vector)
    
    return A_matrix, y_names, b_array, shell_equations