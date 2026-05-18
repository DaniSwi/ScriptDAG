import sympy as sp

# Caché global para este módulo
dag_nodes_cache = {}
node_counter = 0

class DagNode:
    """Clase que representa un nodo en el DAG."""
    def __init__(self, expr, is_n_ary_sum=False):
        global node_counter
        self.id = node_counter
        node_counter += 1
        self.expr = expr
        self.is_n_ary_sum = is_n_ary_sum
        self.children = []

def build_dag_node(expr):
    """
    Construye el DAG fusionando subexpresiones comunes.
    """
    if expr in dag_nodes_cache:
        return dag_nodes_cache[expr]
    
    is_sum = isinstance(expr, sp.Add)
    node = DagNode(expr, is_n_ary_sum=is_sum)
    dag_nodes_cache[expr] = node
    
    if is_sum:
        for arg in expr.args:
            child_node = build_dag_node(arg)
            node.children.append(child_node)
            
    return node

def clear_cache():
    """Limpia el caché por si procesamos múltiples instancias en una sola ejecución."""
    global dag_nodes_cache, node_counter
    dag_nodes_cache.clear()
    node_counter = 0