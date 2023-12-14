using SymPy
const sympy_parsing_mathematica = SymPy.PyCall.pyimport("sympy.parsing.mathematica")

s = "(Sqrt[Pi])/(2*EllipticE[m])"
ex = sympy_parsing_mathematica.mathematica(s, Dict("EllipticE[x]"=>"elliptic_e(x)"))


SymPy.convert_expr(ex, use_julia_code=true) # :(sqrt(pi) ./ (2 * elliptic_e(m)))