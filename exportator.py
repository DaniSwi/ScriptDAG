import numpy as np
import sympy as sp
import re
from parser_instancias import parse_instance
from matrix_extractor import extract_linear_system
import dag_builder

def parse_domain_string(dom_str):
    numeros = re.findall(r'-?\d+\.?\d*(?:e-?\d+)?', dom_str)
    if len(numeros) == 2: return [float(numeros[0]), float(numeros[1])]
    return [-1e8, 1e8]

def calcular_dominios_w(A, y_vars, b, domains_dict):
    w_domains = {}
    limites_y = [parse_domain_string(domains_dict.get(var, "[-1e8, 1e8]")) for var in y_vars]
    for i in range(A.shape[0]):
        min_w = float(b[i]); max_w = float(b[i])
        for j in range(len(y_vars)):
            coef = A[i, j]; min_y, max_y = limites_y[j]
            if coef > 0: min_w += coef * min_y; max_w += coef * max_y
            elif coef < 0: min_w += coef * max_y; max_w += coef * min_y 
        w_domains[f"w_{i}"] = [min_w, max_w]
    return w_domains

def calcular_dominios_v(v_equations, known_domains):
    v_domains = {}
    for v_sym, expr in v_equations:
        def eval_interval(ex):
            if isinstance(ex, sp.Number): val = float(ex); return [val, val]
            if isinstance(ex, sp.Symbol): return known_domains.get(str(ex), [-1e8, 1e8])
            if isinstance(ex, sp.Mul):
                args_intervals = [eval_interval(arg) for arg in ex.args]
                res = args_intervals[0]
                for i in range(1, len(args_intervals)):
                    I1, I2 = res, args_intervals[i]
                    p = [I1[0]*I2[0], I1[0]*I2[1], I1[1]*I2[0], I1[1]*I2[1]]
                    res = [min(p), max(p)]
                return res
            if isinstance(ex, sp.Pow):
                base_i = eval_interval(ex.base); exp_i = eval_interval(ex.exp)
                if exp_i[0] == exp_i[1] and int(exp_i[0]) == exp_i[0]:
                    n = int(exp_i[0])
                    if n % 2 != 0: return [base_i[0]**n, base_i[1]**n]
                    else:
                        if base_i[0] <= 0 <= base_i[1]: return [0, max(base_i[0]**n, base_i[1]**n)]
                        else: return [min(base_i[0]**n, base_i[1]**n), max(base_i[0]**n, base_i[1]**n)]
            if isinstance(ex, (sp.sin, sp.cos)): return [-1.0, 1.0]
            return [-1e8, 1e8]
        dom = eval_interval(expr)
        v_domains[str(v_sym)] = dom
        known_domains[str(v_sym)] = dom 
    return v_domains

def exportar_a_txt(A, y_vars, b, shell_eqs_strs, v_equations, vars_dict, domains_dict, w_domains, v_domains, output_txt_name, usar_dominios):
    y_syms = [sp.sympify(y) for y in y_vars]
    num_w = A.shape[0]
    
    with open(output_txt_name, 'w') as f:
        f.write("variables\n")
        # Variables originales
        for var_name in vars_dict.keys():
            f.write(f"  {var_name} in {domains_dict.get(var_name, '[-1e8, 1e8]')};\n")
            
        # Variables W
        for i in range(num_w):
            dom_str = f"[{w_domains[f'w_{i}'][0]}, {w_domains[f'w_{i}'][1]}]" if usar_dominios else "[-1e8, 1e8]"
            f.write(f"  w_{i} in {dom_str};\n")
            
        # Variables V
        for v_sym, _ in v_equations:
            dom_str = f"[{v_domains[str(v_sym)][0]}, {v_domains[str(v_sym)][1]}]" if usar_dominios else "[-1e8, 1e8]"
            f.write(f"  {v_sym} in {dom_str};\n")
            
        f.write("\nconstraints\n")
        
        # Ecuaciones W (Matriz)
        for i in range(num_w):
            expr_lineal = float(b[i])
            for j in range(len(y_vars)):
                coeff = A[i, j]
                if coeff != 0: expr_lineal += coeff * y_syms[j]
            expr_str = str(expr_lineal).replace('**', '^').replace(' ', '')
            f.write(f"  w_{i}={expr_str};\n")
            
        # Ecuaciones V (No lineales / Multiplicaciones)
        for v_sym, expr in v_equations:
            expr_str = str(expr).replace('**', '^').replace(' ', '')
            f.write(f"  {v_sym}={expr_str};\n")
            
        # Cascarones Topológicos
        for eq_str in shell_eqs_strs:
            eq_format = str(eq_str).replace('**', '^')
            f.write(f"  {eq_format};\n")
        f.write("end\n")

def procesar_instancia(file_path, output_npz_name, output_txt_name, calcular_dominios=False):
    dag_builder.clear_cache()
    vars_dict, domains_dict, constr_list = parse_instance(file_path)
    A, y_vars, b, shell_equations, v_equations = extract_linear_system(constr_list)
    shell_eqs_strs = np.array([f"{str(lhs)} {op} {str(rhs)}" for lhs, op, rhs in shell_equations])
    v_eqs_strs = np.array([f"{str(sym)}={str(expr)}" for sym, expr in v_equations])
    
    w_domains = {}
    v_domains = {}
    if calcular_dominios:
        known_domains = {var: parse_domain_string(dom) for var, dom in domains_dict.items()}
        w_domains = calcular_dominios_w(A, y_vars, b, domains_dict)
        known_domains.update(w_domains)
        v_domains = calcular_dominios_v(v_equations, known_domains)
    
    np.savez(output_npz_name, A=A, y=y_vars, b=b, shell_eqs=shell_eqs_strs, v_eqs=v_eqs_strs)
    exportar_a_txt(A, y_vars, b, shell_eqs_strs, v_equations, vars_dict, domains_dict, w_domains, v_domains, output_txt_name, calcular_dominios)