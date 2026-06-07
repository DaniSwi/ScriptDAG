import re
import sympy as sp

def parse_instance(file_path):
    """
    Lee el archivo de texto y extrae variables y restricciones.
    Soporta mayúsculas/minúsculas y operadores de inecuación (>=, <=, =).
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
            
        line_lower = line.lower()
        
        # 1. Ignoramos la función objetivo (nuestro sistema procesa el DAG de restricciones)
        if line_lower.startswith('minimize') or line_lower.startswith('maximize'):
            print(" [!] Nota: Función objetivo detectada y omitida de la matriz.")
            continue
            
        # 2. Detección flexible de secciones
        if line_lower == 'variables':
            mode = 'vars'
            continue
        elif line_lower == 'constraints':
            mode = 'constraints'
            continue
        elif line_lower == 'end':
            break
            
        # 3. Extracción de datos
        if mode == 'vars':
            match = re.match(r'([a-zA-Z0-9_]+)\s+in\s+\[.*\];', line)
            if match:
                var_name = match.group(1)
                variables[var_name] = sp.Symbol(var_name)
                
        elif mode == 'constraints':
            line = line.rstrip(';')
            
            # Rompemos la línea usando CUALQUIER operador relacional
            parts = re.split(r'(>=|<=|==|=|>|<)', line)
            
            if len(parts) == 3:
                lhs_str, op, rhs_str = parts
                lhs_str = lhs_str.replace('^', '**')
                
                lhs_expr = sp.sympify(lhs_str, locals=variables)
                rhs_val = float(rhs_str)
                
                # AHORA GUARDAMOS 3 COSAS: (Lado_Izquierdo, Operador, Lado_Derecho)
                constraints.append((lhs_expr, op, rhs_val))
                
    return variables, constraints