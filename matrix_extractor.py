import numpy as np
import sympy as sp
from dag_builder import build_dag_node

def extract_linear_system(constraints):
    """
    Extrae el sistema lineal A * y = b de las sumas n-arias,
    y preserva las ecuaciones originales usando variables auxiliares w_i.
    """
    A_rows = []
    b_vector = []
    y_variables = {}
    
    # Aquí guardaremos las ecuaciones "cascarón"
    shell_equations = []
    
    # Contador para generar nombres únicos de variables w (w_0, w_1, ...)
    w_counter = 0

    def replace_and_extract_sums(expr):
        """
        Función recursiva (el 'bisturí'). Recorre la expresión matemática.
        Si encuentra una suma, la extrae para la matriz A y devuelve una variable w_i en su lugar.
        """
        nonlocal w_counter
        
        # 1. Si encontramos una suma n-aria (El objetivo)
        if isinstance(expr, sp.Add):
            # Creamos la nueva variable auxiliar
            w_sym = sp.Symbol(f'w_{w_counter}')
            w_counter += 1
            
            row_dict = {}
            constant_term = 0.0
            
            # Analizamos los hijos de la suma
            for term in expr.args:
                coeff, var_part = term.as_coeff_Mul()
                
                if var_part == 1:
                    constant_term += float(coeff)
                else:
                    if var_part not in y_variables:
                        y_variables[var_part] = len(y_variables)
                    
                    col_index = y_variables[var_part]
                    row_dict[col_index] = row_dict.get(col_index, 0.0) + float(coeff)
            
            # Guardamos la fila en nuestra matriz A
            A_rows.append(row_dict)
            
            # Guardamos la constante. 
            # Como la ecuación es w_i = A*y + b, el término 'b' es la constante de la suma.
            b_vector.append(constant_term)
            
            # Devolvemos la variable w_i para que SymPy la ponga en lugar de la suma gigante
            return w_sym
            
        # 2. Si es una multiplicación, revisamos sus partes internamente
        elif isinstance(expr, sp.Mul):
            return sp.Mul(*[replace_and_extract_sums(arg) for arg in expr.args])
            
        # 3. Si es una potencia (ej. Suma^-1), revisamos la base y el exponente
        elif isinstance(expr, sp.Pow):
            return sp.Pow(replace_and_extract_sums(expr.base), replace_and_extract_sums(expr.exp))
            
        # 4. Si es una función (sin, cos, exp), revisamos su argumento
        elif isinstance(expr, sp.Function):
            return expr.func(*[replace_and_extract_sums(arg) for arg in expr.args])
            
        # 5. Si es un término simple (variable o número), lo dejamos intacto
        else:
            return expr

    # ==========================================
    # LÓGICA PRINCIPAL DEL EXTRACTOR
    # ==========================================
    for lhs_expr, rhs_val in constraints:
        # Construimos el DAG (Mantenemos la lógica de la sección 3.1)
        build_dag_node(lhs_expr)
        
        # Aplicamos nuestro "bisturí" a la expresión del lado izquierdo
        shell_lhs = replace_and_extract_sums(lhs_expr)
        
        # Guardamos la ecuación cascarón final (ej. w_0**-1 = -0.12389)
        shell_equations.append((shell_lhs, rhs_val))

    # Construimos la matriz numpy A final
    num_rows = len(A_rows)
    num_cols = len(y_variables)
    A_matrix = np.zeros((num_rows, num_cols))
    
    for i, row_dict in enumerate(A_rows):
        for j, coeff in row_dict.items():
            A_matrix[i, j] = coeff
            
    y_names = np.array([str(expr) for expr, idx in sorted(y_variables.items(), key=lambda item: item[1])])
    b_array = np.array(b_vector)
    
    # Ahora devolvemos también los cascarones
    return A_matrix, y_names, b_array, shell_equations