"""
Capa 2: construccion del DAG (Algoritmo 1 del paper, pasos de normalizacion).

Cambios respecto de la version anterior:
  D2/D1  La granularidad de un nodo es el *nodo variable*, nunca 'coef * nodo'.
         Todo coeficiente numerico se saca hacia afuera durante la construccion,
         de modo que x2**2, -x2**2 y 5*x2**2 comparten una sola columna.
  D3     Ya no se usa sympify(): SymPy distribuye -6*(x1+x3) y destruye la CSE
         que motiva el paper.  Se parsea con el modulo `ast` de Python y se
         normaliza a mano (aplanado n-ario de sumas y productos, resta -> coef -1,
         division -> Pow(-1)), que es exactamente generate_sum/generate_product.
  D6     ref_count se reemplaza por in_degree real (aristas del DAG).
  B6     los coeficientes se guardan como Fraction exactos, no como float.
  B7     hash-consing con memoizacion por subarbol -> cada nodo se procesa 1 vez.
  --     se eliminan las variables globales: el DAG vive en una instancia de
         DagBuilder, asi que ya no hace falta clear_cache() entre instancias.

Representacion canonica de un nodo:
    var   data = nombre                      args = ()
    const data = Fraction                    args = ()
    add   data = (const, coefs)              args = (nodos sumandos)
    mul   data = None                        args = (factores, con multiplicidad)
    pow   data = Fraction (exponente)        args = (base,)
    func  data = nombre                      args = (argumentos)

Un nodo 'add' representa   const + sum(coefs[i] * args[i])   -> una fila de A.
Un nodo 'mul' NUNCA lleva coeficiente numerico (va afuera) -> una columna de A.
"""

import ast
from fractions import Fraction


class DagNode:
    __slots__ = ('id', 'kind', 'data', 'args', 'in_degree', 'label')

    def __init__(self, nid, kind, data, args):
        self.id = nid
        self.kind = kind
        self.data = data
        self.args = args
        self.in_degree = 0      # numero de aristas entrantes (padres del DAG)
        self.label = None       # 'w_i' si es fila, 'y_j'/nombre si es columna

    def __repr__(self):
        return f"<{self.kind}#{self.id} {self.to_str()}>"

    # ------------------------------------------------------------- impresion
    def to_str(self, usar_labels=False):
        if usar_labels and self.label is not None and self.kind == 'add':
            return self.label
        k = self.kind
        if k == 'var':
            return self.data
        if k == 'const':
            return _num_str(self.data)
        if k == 'pow':
            return f"({self.args[0].to_str(usar_labels)})^({_num_str(self.data)})"
        if k == 'func':
            return f"{self.data}(" + ", ".join(a.to_str(usar_labels) for a in self.args) + ")"
        if k == 'mul':
            return "*".join(f"({a.to_str(usar_labels)})" for a in self.args)
        if k == 'add':
            const, coefs = self.data
            piezas = []
            if const != 0:
                piezas.append(_num_str(const))
            for c, a in zip(coefs, self.args):
                s = a.to_str(usar_labels)
                piezas.append(s if c == 1 else f"{_num_str(c)}*({s})")
            return " + ".join(piezas) if piezas else "0"
        return "?"


def _num_str(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


class ParseError(Exception):
    pass


class DagBuilder:
    """DAG con hash-consing.  Una instancia por sistema; sin estado global."""

    def __init__(self, distribuir_sumas_escaladas=False):
        self._cache = {}        # clave canonica -> DagNode
        self.nodes = []         # indexado por id; los hijos siempre tienen id menor
        self.raices = []        # nodos raiz (uno por restriccion)
        self.distribuir_sumas_escaladas = distribuir_sumas_escaladas

    # ------------------------------------------------------- fabrica de nodos
    def _get(self, clave, kind, data, args):
        n = self._cache.get(clave)
        if n is not None:
            for _ in args:
                pass
            return n
        n = DagNode(len(self.nodes), kind, data, args)
        self._cache[clave] = n
        self.nodes.append(n)
        for a in args:
            a.in_degree += 1
        return n

    def _touch(self, nodo, args):
        """Suma las aristas cuando se reutiliza un nodo ya existente."""
        for a in args:
            a.in_degree += 1

    def var(self, nombre):
        return self._get(('var', nombre), 'var', nombre, ())

    def const(self, valor):
        q = Fraction(valor)
        return self._get(('const', q), 'const', q, ())

    def pow(self, base, exp):
        q = Fraction(exp)
        clave = ('pow', base.id, q)
        existe = clave in self._cache
        n = self._get(clave, 'pow', q, (base,))
        if existe:
            self._touch(n, (base,))
        return n

    def func(self, nombre, args):
        clave = ('func', nombre, tuple(a.id for a in args))
        existe = clave in self._cache
        n = self._get(clave, 'func', nombre, tuple(args))
        if existe:
            self._touch(n, args)
        return n

    def mul(self, factores):
        """Producto n-ario SIN coeficiente numerico (el coeficiente va afuera)."""
        factores = tuple(sorted(factores, key=lambda n: n.id))
        if len(factores) == 1:
            return factores[0]
        clave = ('mul', tuple(f.id for f in factores))
        existe = clave in self._cache
        n = self._get(clave, 'mul', None, factores)
        if existe:
            self._touch(n, factores)
        return n

    def add(self, const, terminos):
        """Suma n-aria.  terminos = [(Fraction coef, DagNode)] con nodo != const."""
        agrupado = {}
        orden = []
        for c, nodo in terminos:
            if nodo.id not in agrupado:
                agrupado[nodo.id] = [Fraction(0), nodo]
                orden.append(nodo.id)
            agrupado[nodo.id][0] += Fraction(c)      # "coefficients are summed up"
        items = [(agrupado[i][0], agrupado[i][1]) for i in orden if agrupado[i][0] != 0]
        items.sort(key=lambda t: t[1].id)
        args = tuple(t[1] for t in items)
        coefs = tuple(t[0] for t in items)
        clave = ('add', Fraction(const), tuple((c, a.id) for c, a in zip(coefs, args)))
        existe = clave in self._cache
        n = self._get(clave, 'add', (Fraction(const), coefs), args)
        if existe:
            self._touch(n, args)
        return n

    # ------------------------------------------------------------- envoltorio
    def escalar(self, coef, nodo):
        """Devuelve UN nodo que representa coef*nodo.

        Solo se usa cuando el resultado tiene que ser un nodo si o si (argumento
        de una funcion, base de una potencia fraccionaria).  Si el nodo es una
        suma, el coeficiente se reparte dentro de la suma: 6*(x1+x3) se vuelve
        el nodo 6*x1+6*x3, que es una fila mas de la matriz (la fila '6|6' de la
        referencia), en vez de un producto por una constante.
        """
        coef = Fraction(coef)
        if nodo is None:
            return self.const(coef)
        if coef == 1:
            return nodo
        if nodo.kind == 'add':
            const, coefs = nodo.data
            return self.add(coef * const,
                            [(coef * c, a) for c, a in zip(coefs, nodo.args)])
        return self.mul([self.const(coef), nodo])

    # ------------------------------------------------------- parseo + normaliz.
    def desde_string(self, expresion, variables):
        """Parsea 'expresion' (sintaxis Ibex, con ^) y devuelve (coef, nodo)."""
        codigo = expresion.replace('^', '**').strip()
        if not codigo:
            raise ParseError('expresion vacia')
        try:
            arbol = ast.parse(codigo, mode='eval').body
        except SyntaxError as e:
            raise ParseError(f"no se pudo parsear la expresion: {expresion[:80]}...") from e
        return self._materializar(self._linealizar(arbol, variables))

    def _linealizar(self, raiz, variables):
        """Recorrido post-orden ITERATIVO (las instancias tienen miles de terminos,
        una version recursiva revienta con RecursionError)."""
        res = {}
        pila = [(raiz, False)]
        while pila:
            nodo, procesado = pila.pop()
            if not procesado:
                pila.append((nodo, True))
                for h in self._hijos_ast(nodo):
                    pila.append((h, False))
                continue
            res[id(nodo)] = self._combinar(nodo, res, variables)
        return res[id(raiz)]

    @staticmethod
    def _hijos_ast(n):
        if isinstance(n, ast.BinOp):
            return (n.left, n.right)
        if isinstance(n, ast.UnaryOp):
            return (n.operand,)
        if isinstance(n, ast.Call):
            return tuple(n.args)
        return ()

    # Una "forma" es el resultado parcial de linealizar un subarbol:
    #   ('n', coef, nodo)        ->  coef * nodo   (nodo puede ser None)
    #   ('s', const, terminos)   ->  suma PENDIENTE, aun sin materializar
    # Las sumas se mantienen pendientes mientras el padre siga siendo + o -,
    # para que 'a+b+c' (que Python parsea como (a+b)+c) genere UNA sola fila
    # n-aria y no una fila espuria por cada suma binaria intermedia.

    def _materializar(self, forma):
        """Convierte una forma en el par (coef, nodo) definitivo."""
        if forma[0] == 'n':
            return forma[1], forma[2]
        _, const, terminos = forma
        if not terminos:
            return const, None
        if const == 0 and len(terminos) == 1:
            return terminos[0]
        return Fraction(1), self.add(const, terminos)

    def _combinar(self, n, res, variables):
        # --- hojas ---------------------------------------------------------
        if isinstance(n, ast.Constant):
            if isinstance(n.value, bool) or not isinstance(n.value, (int, float)):
                raise ParseError(f"constante no numerica: {n.value!r}")
            q = (Fraction(n.value).limit_denominator(10 ** 15)
                 if isinstance(n.value, float) else Fraction(n.value))
            return ('n', q, None)
        if isinstance(n, ast.Name):
            if variables is not None and n.id not in variables:
                raise ParseError(f"variable no declarada: {n.id}")
            return ('n', Fraction(1), self.var(n.id))

        # --- unarios -------------------------------------------------------
        if isinstance(n, ast.UnaryOp):
            f = res[id(n.operand)]
            if isinstance(n.op, ast.UAdd):
                return f
            if isinstance(n.op, ast.USub):
                return self._negar(f)
            raise ParseError(f"operador unario no soportado: {type(n.op).__name__}")

        # --- llamadas a funcion --------------------------------------------
        if isinstance(n, ast.Call):
            if not isinstance(n.func, ast.Name):
                raise ParseError("llamada a funcion no soportada")
            args = [self.escalar(*self._materializar(res[id(a)])) for a in n.args]
            return ('n', Fraction(1), self.func(n.func.id, args))

        # --- binarios ------------------------------------------------------
        if isinstance(n, ast.BinOp):
            izq, der, op = res[id(n.left)], res[id(n.right)], n.op
            if isinstance(op, ast.Add):
                return self._suma(izq, der)
            if isinstance(op, ast.Sub):
                return self._suma(izq, self._negar(der))
            if isinstance(op, ast.Mult):
                return self._producto(self._materializar(izq), self._materializar(der))
            if isinstance(op, ast.Div):
                return self._producto(self._materializar(izq),
                                      self._inverso(self._materializar(der)))
            if isinstance(op, ast.Pow):
                return self._potencia(self._materializar(izq), self._materializar(der))
            raise ParseError(f"operador no soportado: {type(op).__name__}")

        raise ParseError(f"nodo AST no soportado: {type(n).__name__}")

    def _negar(self, forma):
        if forma[0] == 'n':
            return ('n', -forma[1], forma[2])
        return ('s', -forma[1], [(-c, nd) for c, nd in forma[2]])

    # ------------------------------------------------- generate_sum (paso 2-3)
    def _suma(self, izq, der):
        const = Fraction(0)
        terminos = []
        for forma in (izq, der):
            if forma[0] == 's':
                const += forma[1]
                terminos.extend(forma[2])
                continue
            c, nodo = forma[1], forma[2]
            if nodo is None:
                const += c
            elif nodo.kind == 'add' and (c == 1 or self.distribuir_sumas_escaladas):
                # aplanado n-ario: una suma hija de una suma se fusiona en el padre.
                # con coeficiente != 1 NO se aplana, para no destruir la CSE de
                # -6*(x1+x3)  (bug D3).
                sub_const, sub_coefs = nodo.data
                const += c * sub_const
                terminos.extend((c * sc, sa) for sc, sa in zip(sub_coefs, nodo.args))
            else:
                terminos.append((c, nodo))
        return ('s', const, terminos)

    # --------------------------------------------- generate_product (paso 2-3)
    def _producto(self, izq, der):
        coef = izq[0] * der[0]
        factores = []
        for _, nodo in (izq, der):
            if nodo is None:
                continue
            if nodo.kind == 'mul':
                factores.extend(nodo.args)        # aplanado n-ario de productos
            else:
                factores.append(nodo)
        if not factores:
            return ('n', coef, None)
        return ('n', coef, self.mul(factores))

    def _inverso(self, par):
        c, nodo = par
        if c == 0:
            raise ParseError("division por cero literal")
        if nodo is None:
            return Fraction(1, 1) / c, None
        return Fraction(1, 1) / c, self.pow(nodo, -1)

    def _potencia(self, base, exp):
        cb, nb = base
        ce, ne = exp
        if ne is not None:
            # exponente variable: no es un Pow del paper, es una funcion binaria
            return ('n', Fraction(1),
                    self.func('pow', [self.escalar(cb, nb), self.escalar(ce, ne)]))
        if nb is None:
            try:
                return ('n', cb ** ce, None)
            except (ValueError, ZeroDivisionError) as e:
                raise ParseError(f"potencia constante invalida: {cb}^{ce}") from e
        if ce.denominator == 1:
            # (cb*nb)^k = cb^k * nb^k  -> el coeficiente sale afuera (arregla D2)
            return ('n', cb ** ce, self.pow(nb, ce))
        # exponente fraccionario: sacar el coeficiente no es sound para bases
        # negativas, asi que la base queda envuelta.
        return ('n', Fraction(1), self.pow(self.escalar(cb, nb), ce))

    # ------------------------------------------------------------- utilidades
    def registrar_raiz(self, nodo):
        if nodo is not None:
            nodo.in_degree += 1
            self.raices.append(nodo)

    def orden_topologico(self):
        """Los hijos siempre se crean antes que los padres -> id ascendente."""
        return list(self.nodes)

    def resumen(self):
        from collections import Counter
        c = Counter(n.kind for n in self.nodes)
        compartidos = sum(1 for n in self.nodes if n.in_degree > 1)
        return {'nodos': len(self.nodes), 'por_tipo': dict(c),
                'compartidos_CSE': compartidos}
