def cm_regularity(res):
    """
    Castelnuovo-Mumford regularity of a minimal graded free resolution.

    INPUT:
    "res" -- A (minimal) graded free resolution.

    OUTPUT:
    The regularity.
    """
    reg = 0
    for i in range(0, len(res)):
        b = max(res.betti(i)) - i + 1
        if b > reg:
            reg = b
    return reg

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
