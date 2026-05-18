import re
import sympy as sp

def parse_instance(file_path):
    """
    Lee el archivo .txt de la instancia y extrae las variables y restricciones.
    Retorna un diccionario de variables y una lista de ecuaciones SymPy.
    """
    variables = {}
    constraints = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    mode = None
    for line in lines:
        line = line.strip()
        if not line or line.startswith('//'):
            continue
            
        if line == 'Variables':
            mode = 'vars'
            continue
        elif line == 'Constraints':
            mode = 'constraints'
            continue
        elif line == 'end':
            break
            
        if mode == 'vars':
            match = re.match(r'([a-zA-Z0-9_]+)\s+in\s+\[.*\];', line)
            if match:
                var_name = match.group(1)
                variables[var_name] = sp.Symbol(var_name)
                
        elif mode == 'constraints':
            line = line.rstrip(';')
            if '=' in line:
                lhs_str, rhs_str = line.split('=')
                lhs_str = lhs_str.replace('^', '**')
                
                lhs_expr = sp.sympify(lhs_str, locals=variables)
                rhs_val = float(rhs_str)
                constraints.append((lhs_expr, rhs_val))
                
    return variables, constraints