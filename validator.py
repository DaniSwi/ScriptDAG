import numpy as np
import sympy as sp
import re
from parser_instancias import parse_instance

def validar_sistema(ruta_instancia, ruta_npz, verbose=True):
    print("======================================================")
    print(" INICIANDO AUDITORÍA MATEMÁTICA (VALIDACIÓN CON SYMPY)")
    print("======================================================\n")

    # Extraemos las restricciones e ignoramos las variables y los dominios
    _, _, original_constraints = parse_instance(ruta_instancia)
    
    datos = np.load(ruta_npz)
    A_matrix = datos['A']
    y_names = datos['y']
    b_vector = datos['b']
    shell_eqs_strs = datos['shell_eqs']

    y_syms = [sp.sympify(y_str) for y_str in y_names]
    w_dict = {}
    
    if verbose:
        print("--- FASE 1: RECONSTRUCCIÓN DE VARIABLES w_i ---")
        
    for i in range(len(A_matrix)):
        expr_lineal = float(b_vector[i])
        
        if verbose:
            print(f"\n[w_{i}] Construyendo desde la Fila {i} de la Matriz A:")
            print(f"      -> Constante inicial (b): {expr_lineal}")
            
        for j in range(len(y_names)):
            coeff = A_matrix[i, j]
            if coeff != 0.0:
                expr_lineal += coeff * y_syms[j]
                if verbose:
                    print(f"      -> Agregando término: ({coeff}) * {y_syms[j]}")
                
        w_simbolo = sp.Symbol(f'w_{i}')
        w_dict[w_simbolo] = expr_lineal
        
        if verbose:
            print(f"      => RESULTADO LINEAL: w_{i} = {expr_lineal}")

    print("\n--- FASE 2: SUSTITUCIÓN Y COMPROBACIÓN ---\n")
    
    validacion_exitosa = True
    
    for i in range(len(original_constraints)):
        orig_lhs, orig_op, orig_rhs = original_constraints[i]
        
        shell_str = str(shell_eqs_strs[i])
        parts = re.split(r'(>=|<=|==|=|>|<)', shell_str)
        shell_lhs_sym = sp.sympify(parts[0].strip())
        shell_rhs_val = float(parts[2].strip())
        
        rec_lhs = shell_lhs_sym.subs(w_dict)
        
        if verbose:
            print(f"ECUACIÓN {i}:")
            print(f"  1. Original parseada : {orig_lhs} = {orig_rhs}")
            print(f"  2. Cascarón guardado : {shell_lhs_sym} = {shell_rhs_val}")
            print(f"  3. Sustituyendo w_i  : {rec_lhs} = {shell_rhs_val}")
        
        if not np.isclose(orig_rhs, shell_rhs_val):
            print(f"  [!] ERROR: Los resultados (RHS) no coinciden.\n")
            validacion_exitosa = False
            continue

        diferencia = sp.simplify(orig_lhs - rec_lhs)
        
        if verbose:
            print(f"  4. Resta Matemática  : ({orig_lhs}) - ({rec_lhs})")
            print(f"  5. Simplificación    : {diferencia}")
            
        if diferencia == 0:
            print(f"  [PASS] Equivalencia perfecta confirmada.\n")
        else:
            print(f"  [FAIL] Diferencia detectada: {diferencia}\n")
            validacion_exitosa = False

    print("======================================================")
    if validacion_exitosa:
        print(" RESULTADO FINAL: ÉXITO ABSOLUTO.")
    else:
        print(" RESULTADO FINAL: FALLO.")
    print("======================================================")

if __name__ == "__main__":
    ruta_instancia = 'instances/instprofe.txt'
    ruta_npz = 'outputs/subsistema_n_ario_profe.npz'
    # Cambia a False cuando quieras correr cientos de instancias rápido
    validar_sistema(ruta_instancia, ruta_npz, verbose=True)