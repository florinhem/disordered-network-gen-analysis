
using SymPy
const sympy_parsing_mathematica = SymPy.PyCall.pyimport("sympy.parsing.mathematica")
mathematica2julia(s::AbstractString, substitutions::Pair{<:AbstractString,<:AbstractString}...) =
           SymPy.walk_expression(sympy_parsing_mathematica."mathematica"(s, Dict(substitutions...)))

mathematica2julia("((x3 (-y1 + y2) + x2 (y1 - y3) + 
x1 (-y2 + y3))/((x1 - x2) (x1 - x3) (x2 - x3)))")