
using SymPy

const sympy_parsing_mathematica = SymPy.PyCall.pyimport("sympy.parsing.mathematica")

# define a string which is converted into a Mathematica expression
s = "1/Sqrt[(x+y)^3]"
ex = sympy_parsing_mathematica.parse_mathematica(s)
f_expr = SymPy.SymPyCore.walk_expression(Sym(ex))

# generates a function
@generated ftest(x,y) = f_expr
ftest(1,0)