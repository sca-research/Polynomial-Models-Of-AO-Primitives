def generate_field_equations(variables):
    """
    Generates the field equations for a given variable set.

    INPUT:
    "variables" -- variables of a polynomial ring over a finite field.

    OUTPUT:
    "field_equations" -- field equations for given variable set.
    """
    q = variables[0].parent().base_ring().order()
    field_equations = [var**q - var for var in variables]
    return field_equations

def highest_degree_component(poly):
    """
    Extracts the highest degree component of a polynomial.

    INPUT:
    "poly" -- a polynomial.

    OUTPUT:
    "poly_top" -- highest degree component of a polynomial.
    """
    poly_top = 0
    d = poly.degree()
    for (coeff, mon) in zip(poly.coefficients(), poly.monomials()):
        if mon.degree() == d:
            poly_top += coeff * mon
    return poly_top

def macaulay_bound(polys, n_vars=None):
    """
    Computes the Macaulay bound of a polynomial system.

    INPUT:
    "polys" -- a polynomial system

    "n_vars" -- number of variables to be considered for
    the Macaulay bound. If not specified the number of 
    variables of the base ring will be used.

    OUTPUT:
    Macaulay bound
    """
    if n_vars is None:
        n_vars = len(polys[0].parent().gens())
    degs = sorted([poly.degree() for poly in polys], reverse=True)
    l = min([n_vars + 1, len(polys)])
    return sum(degs[0:l]) - l + 1
