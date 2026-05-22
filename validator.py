import numpy as np
import sympy as sp
from parser_instancias import parse_instance

def validar_sistema(ruta_instancia, ruta_npz):
    print("======================================================")
    print(" INICIANDO AUDITORÍA MATEMÁTICA (VALIDACIÓN CON SYMPY)")
    print("======================================================\n")

    # 1. Cargar la fuente de verdad (Las ecuaciones originales)
    print("[1/4] Parseando el archivo original para obtener la 'verdad absoluta'...")
    _, original_constraints = parse_instance(ruta_instancia)

    # 2. Cargar nuestro pre-procesamiento
    print(f"[2/4] Cargando el archivo procesado: {ruta_npz}...")
    datos = np.load(ruta_npz)
    A_matrix = datos['A']
    y_names = datos['y']
    b_vector = datos['b']
    shell_eqs_strs = datos['shell_eqs']

    # 3. Reconstruir las variables lineales w_i (La operación w = A*y + b)
    print("[3/4] Reconstruyendo los sistemas lineales (Matriz A * Vector y + Vector b)...")
    
    # Convertimos los nombres de las columnas 'y' a símbolos matemáticos
    y_syms = [sp.sympify(y_str) for y_str in y_names]
    
    w_dict = {}  # Diccionario para almacenar {w_0: ecuacion_lineal_0, w_1: ecuacion_lineal_1...}
    
    for i in range(len(A_matrix)):
        # Iniciamos con la constante b_i
        expr_lineal = float(b_vector[i])
        
        # Hacemos el producto punto de la fila i con el vector y
        for j in range(len(y_names)):
            coeff = A_matrix[i, j]
            if coeff != 0.0:
                expr_lineal += coeff * y_syms[j]
                
        # Registramos la variable w_i con su ecuación matemática
        w_simbolo = sp.Symbol(f'w_{i}')
        w_dict[w_simbolo] = expr_lineal

    # 4. Sustitución y Validación cruzada
    print("[4/4] Ejecutando sustitución en cascarones y comprobando equivalencia matemática...\n")
    
    validacion_exitosa = True
    
    for i in range(len(original_constraints)):
        # Datos Originales
        orig_lhs, orig_rhs = original_constraints[i]
        
        # Datos Reconstruidos (Parseamos el cascarón guardado)
        shell_str = str(shell_eqs_strs[i])
        shell_lhs_str, shell_rhs_str = shell_str.split(' = ')
        shell_lhs_sym = sp.sympify(shell_lhs_str)
        shell_rhs_val = float(shell_rhs_str)
        
        # El Bisturí Inverso: Sustituimos las 'w_i' por sus expresiones lineales reconstruidas
        rec_lhs = shell_lhs_sym.subs(w_dict)
        
        # Verificamos la constante del lado derecho
        if not np.isclose(orig_rhs, shell_rhs_val):
            print(f"[!] ERROR en la Ecuación {i}: Los resultados (RHS) no coinciden.")
            validacion_exitosa = False
            continue

        # La Prueba Suprema: Restamos la ecuación original menos la reconstruida
        # Matemáticamente, si son idénticas, (A - B) debe ser igual a 0.
        diferencia = sp.simplify(orig_lhs - rec_lhs)
        
        if diferencia == 0:
            print(f" [PASS] Ecuación {i}: Perfectamente equivalente.")
        else:
            print(f" [FAIL] Ecuación {i}: Diferencia detectada -> {diferencia}")
            validacion_exitosa = False

    print("\n======================================================")
    if validacion_exitosa:
        print(" RESULTADO FINAL: ÉXITO ABSOLUTO.")
        print(" El sistema pre-procesado en el archivo .npz")
        print(" preserva la integridad matemática del problema.")
    else:
        print(" RESULTADO FINAL: FALLO.")
        print(" Se detectaron inconsistencias matemáticas.")
    print("======================================================")

if __name__ == "__main__":
    ruta_instancia = 'instances/inst001.txt'
    ruta_npz = 'outputs/subsistema_n_ario_test.npz'
    
    validar_sistema(ruta_instancia, ruta_npz)