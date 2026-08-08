import re
import sympy as sp

def parse_instance(file_path):
    """
    Lee el archivo de texto y extrae variables, dominios y restricciones.
    """
    variables = {}
    domains = {}
    constraints = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    mode = None
    for line in lines:
        line = line.strip()
        if not line or line.startswith('//'):
            continue
            
        line_lower = line.lower()
        if line_lower.startswith('minimize') or line_lower.startswith('maximize'):
            print(" [!] Nota: Función objetivo detectada y omitida.")
            continue
            
        if line_lower == 'variables':
            mode = 'vars'
            continue
        elif line_lower == 'constraints':
            mode = 'constraints'
            continue
        elif line_lower == 'end':
            break
            
        if mode == 'vars':
            match = re.match(r'([a-zA-Z0-9_]+)\s+in\s+(\[.*?\]);', line)
            if match:
                var_name = match.group(1)
                var_domain = match.group(2)
                variables[var_name] = sp.Symbol(var_name)
                domains[var_name] = var_domain
                
        elif mode == 'constraints':
            line = line.rstrip(';')
            parts = re.split(r'(>=|<=|==|=|>|<)', line)
            if len(parts) == 3:
                lhs_str, op, rhs_str = parts
                lhs_str = lhs_str.replace('^', '**')
                lhs_expr = sp.sympify(lhs_str, locals=variables)
                rhs_val = float(rhs_str)
                constraints.append((lhs_expr, op, rhs_val))
                
    return variables, domains, constraints