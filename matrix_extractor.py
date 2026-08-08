import numpy as np
import sympy as sp
from dag_builder import dag_nodes_cache, build_dag_node

def extract_linear_system(constraints):
    A_rows = []
    b_vector = []
    y_variables = {}
    shell_equations = []
    
    sum_cache = {}
    w_counter = 0
    
    v_cache = {}
    v_equations = [] 
    v_counter = 0

    def process_expr(expr):
        nonlocal w_counter, v_counter
        
        if isinstance(expr, (sp.Symbol, sp.Number)):
            return expr
            
        processed_args = [process_expr(arg) for arg in expr.args]
        new_expr = expr.func(*processed_args) if processed_args else expr
            
        if isinstance(expr, sp.Mul):
            # Modalidad Multiplicación N-Aria pedida por el profesor (ej. A * B * V_11)
            vars_in_mul = [arg for arg in new_expr.args if arg.free_symbols]
            is_n_ary_mult = len(vars_in_mul) > 1
            
            node = dag_nodes_cache.get(expr)
            is_repeated = node and node.ref_count > 1
            
            # Englobamos si es un producto n-ario O si es una subexpresión repetida
            if is_n_ary_mult or is_repeated:
                if expr not in v_cache:
                    v_sym = sp.Symbol(f'v_{v_counter}')
                    v_counter += 1
                    v_cache[expr] = v_sym
                    v_equations.append((v_sym, new_expr))
                return v_cache[expr]
            return new_expr
                
        elif isinstance(expr, (sp.Pow, sp.Function)):
            node = dag_nodes_cache.get(expr)
            if node and node.ref_count > 1:
                if expr not in v_cache:
                    v_sym = sp.Symbol(f'v_{v_counter}')
                    v_counter += 1
                    v_cache[expr] = v_sym
                    v_equations.append((v_sym, new_expr))
                return v_cache[expr]
            return new_expr
                
        elif isinstance(expr, sp.Add):
            # Modalidad Sumas N-Arias -> Matriz Lineal Mixta
            if expr not in sum_cache:
                w_sym = sp.Symbol(f'w_{w_counter}')
                w_counter += 1
                row_dict = {}
                constant_term = 0.0
                
                # Aquí v_1 y v_2 entran naturalmente a la matriz lineal
                for term in new_expr.args:
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
                sum_cache[expr] = w_sym
            return sum_cache[expr]
        else:
            return new_expr

    for lhs_expr, op, rhs_val in constraints:
        build_dag_node(lhs_expr)
        
    for lhs_expr, op, rhs_val in constraints:
        shell_lhs = process_expr(lhs_expr)
        shell_equations.append((shell_lhs, op, rhs_val))

    num_rows = len(A_rows)
    num_cols = len(y_variables)
    A_matrix = np.zeros((num_rows, num_cols))
    
    for i, row_dict in enumerate(A_rows):
        for j, coeff in row_dict.items():
            A_matrix[i, j] = coeff
            
    y_names = np.array([str(expr) for expr, idx in sorted(y_variables.items(), key=lambda item: item[1])])
    b_array = np.array(b_vector)
    
    return A_matrix, y_names, b_array, shell_equations, v_equations