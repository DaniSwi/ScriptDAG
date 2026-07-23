import numpy as np
import sympy as sp
import re
from parser_instancias import parse_instance

def validar_sistema(ruta_instancia, ruta_npz, verbose=True):
    print("======================================================")
    print(" INICIANDO AUDITORÍA MATEMÁTICA (VALIDACIÓN CON SYMPY)")
    print("======================================================\n")

    # Extraemos restricciones originales
    _, _, original_constraints = parse_instance(ruta_instancia)
    
    datos = np.load(ruta_npz)
    A_matrix = datos['A']
    y_names = datos['y']
    b_vector = datos['b']
    shell_eqs_strs = datos['shell_eqs']
    
    # Soporte para versiones anteriores del npz
    v_eqs_strs = datos.get('v_eqs', [])

    y_syms = [sp.sympify(y_str) for y_str in y_names]
    w_dict = {}
    
    if verbose: print("--- FASE 1: RECONSTRUCCIÓN DE VARIABLES w_i y v_i ---")
        
    for i in range(len(A_matrix)):
        expr_lineal = float(b_vector[i])
        for j in range(len(y_names)):
            coeff = A_matrix[i, j]
            if coeff != 0.0:
                expr_lineal += coeff * y_syms[j]
        w_simbolo = sp.Symbol(f'w_{i}')
        w_dict[w_simbolo] = expr_lineal
        
    v_list = []
    for eq_str in v_eqs_strs:
        lhs, rhs = str(eq_str).split('=')
        v_list.append((sp.sympify(lhs.strip()), sp.sympify(rhs.strip())))

    print("\n--- FASE 2: SUSTITUCIÓN Y COMPROBACIÓN ---\n")
    validacion_exitosa = True
    
    for i in range(len(original_constraints)):
        orig_lhs, orig_op, orig_rhs = original_constraints[i]
        
        shell_str = str(shell_eqs_strs[i])
        parts = re.split(r'(>=|<=|==|=|>|<)', shell_str)
        shell_lhs_sym = sp.sympify(parts[0].strip())
        shell_rhs_val = float(parts[2].strip())
        
        # Primero reemplazamos W, luego desenrollamos V en reversa
        rec_lhs = shell_lhs_sym.subs(w_dict)
        for v_sym, v_expr in reversed(v_list):
            rec_lhs = rec_lhs.subs(v_sym, v_expr)
            
        if not np.isclose(orig_rhs, shell_rhs_val):
            print(f"  [!] ERROR: Los resultados (RHS) no coinciden.\n")
            validacion_exitosa = False
            continue

        diferencia = sp.simplify(orig_lhs - rec_lhs)
        
        if diferencia == 0:
            if verbose: print(f"  [PASS] Equivalencia perfecta en ecuación {i}.")
        else:
            print(f"  [FAIL] Diferencia detectada en ecuación {i}: {diferencia}\n")
            validacion_exitosa = False

    print("======================================================")
    if validacion_exitosa:
        print(" RESULTADO FINAL: ÉXITO ABSOLUTO.")
    else:
        print(" RESULTADO FINAL: FALLO.")
    print("======================================================")