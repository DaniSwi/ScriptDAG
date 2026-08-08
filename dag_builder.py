import sympy as sp

dag_nodes_cache = {}
node_counter = 0

class DagNode:
    """Clase que representa un nodo en el DAG."""
    def __init__(self, expr):
        global node_counter
        self.id = node_counter
        node_counter += 1
        self.expr = expr
        self.ref_count = 1  # Contador para detectar subexpresiones comunes (CSE)
        self.children = []

def build_dag_node(expr):
    """Construye el DAG fusionando subexpresiones comunes."""
    if expr in dag_nodes_cache:
        dag_nodes_cache[expr].ref_count += 1
        return dag_nodes_cache[expr]
    
    node = DagNode(expr)
    dag_nodes_cache[expr] = node
    
    if hasattr(expr, 'args'):
        for arg in expr.args:
            child_node = build_dag_node(arg)
            node.children.append(child_node)
            
    return node

def clear_cache():
    """Limpia el caché por si procesamos múltiples instancias."""
    global dag_nodes_cache, node_counter
    dag_nodes_cache.clear()
    node_counter = 0