import numpy as np
from dag_builder import build_dag_node

def extract_linear_system(constraints):
    """
    Extrae el sistema lineal A * y = b de las sumas n-arias del DAG.
    """
    A_rows = []
    b_vector = []
    y_variables = {}
    
    for lhs_expr, rhs_val in constraints:
        root_node = build_dag_node(lhs_expr)
        
        if root_node.is_n_ary_sum:
            row_dict = {}
            constant_term = 0.0
            
            for child_node in root_node.children:
                term = child_node.expr
                coeff, var_part = term.as_coeff_Mul()
                
                if var_part == 1:
                    constant_term += float(coeff)
                else:
                    if var_part not in y_variables:
                        y_variables[var_part] = len(y_variables)
                    
                    col_index = y_variables[var_part]
                    row_dict[col_index] = row_dict.get(col_index, 0.0) + float(coeff)
            
            A_rows.append(row_dict)
            b_vector.append(rhs_val - constant_term)

    num_rows = len(A_rows)
    num_cols = len(y_variables)
    A_matrix = np.zeros((num_rows, num_cols))
    
    for i, row_dict in enumerate(A_rows):
        for j, coeff in row_dict.items():
            A_matrix[i, j] = coeff
            
    y_names = np.array([str(expr) for expr, idx in sorted(y_variables.items(), key=lambda item: item[1])])
    b_array = np.array(b_vector)
    
    return A_matrix, y_names, b_array