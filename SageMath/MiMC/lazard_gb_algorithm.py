# SageMath functions to generate Macaulay matrix and perform Gaussian elimination
# The code is forked from https://asdm.gmbh/2021/03/15/d_reg/

from copy import copy
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.misc.flatten import flatten
from sage.misc.timing import walltime
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

def all_monoms_up_to_deg(ring, d):
    """
    Generates all monomials in a polynomial ring up to degree d.

    INPUT:
    "ring" -- A polynomial ring.
    "d" -- An integer.

    OUTPUT:
    List of monomials sorted via degree.
    """
    all_monoms = [ring(1)]
    for i in range(1, d + 1):
        all_monoms = all_monoms + [ring({tuple(a):1}) for a in WeightedIntegerVectors(i, len(ring.gens()) * [1])]
    return sorted(all_monoms, reverse=True)

def polynomial_division(f, divisors):
    """
    Multivariate polynomial division.

    INPUT:
    "f" -- A polynomial.
    "divisors" -- A list of polynomials.

    OUTPUT:
    Quotients and remainder.
    """
    f_original = f
    quotients = [0]*len(divisors)
    rem = 0
    while f != 0:
        i = 0
        division_occured = False
        while i < len(divisors) and not division_occured:
            divisable = False
            try:
                divisable = divisors[i].lt().divides(f.lt())
            except NotImplementedError as e:
                pass # _beautiful_ solution
            if divisable:
                q, _ = f.lt().quo_rem(divisors[i].lt())
                quotients[i] += q
                f = f - q * divisors[i]
                division_occured = True
            else:
                i += 1
        if not division_occured:
            r = f.lt()
            rem += r
            f -= r
    assert f_original == sum([q*d for q, d in zip(quotients, divisors)]) + rem
    return quotients, rem

def reduced_gb(gb):
    """
    Computes the reduced Groebner basis of
    an input Groebner basis.

    INPUT:
    "gb" -- A Groebner basis.

    OUTPUT:
    The reduced basisi.
    """
    ring = gb[0].parent()
    for i in range(0, len(gb)):
        gb[i] = gb[i] / gb[i].lc()
    lms = [poly.lm() for poly in gb]
    lms_gb = ring.ideal(lms).groebner_basis()
    gb_red = [poly for poly in gb if poly.lm() in lms_gb]
    return sorted(gb_red, reverse=True)

def s_poly(f, g):
    """
    The S-polynomial of f and g

    INPUT:
    "f" -- A polynomial.
    "g" -- A polynomial.

    OUTPUT:
    S_> (f, g).
    """
    l = f.lm().lcm(g.lm())
    factor_f = l // f.lt()
    factor_g = l // g.lt()
    return factor_f * f - factor_g * g

def buchberger_criterion(gb):
    """
    Buchberger's criterion for Groebner basis.

    INPUT:
    "gb" -- A list of polynomials.

    OUTPUT:
    True or False.
    """
    for j in range(len(gb)):
        for i in range(j):
            s = s_poly(gb[i], gb[j])
            _, rem = polynomial_division(s, gb)
            if rem:
                return False
    return True

def lazard_gb_algorithm(polys, 
                        solving_deg=None, 
                        print_intermediate_basis=False):
    """
    Generates the Macaulay matrix for increasing d,
    performs Gaussian elimination, and verifies whether
    a Groebner basis has been found via Buchberger's
    criterion.

    INPUT:
    "polys" -- A list of polynomials.
    "solving_deg" -- Generate the Macaulay matrix in a
                     given degree.
                     Default value is None.
    "print_intermediate_basis" -- Boolean to print intermediate
                                  bases.
                                  Default value is False.

    OUTPUT:
    A Groebner basis of the input system.
    """
    ring = polys[0].parent()
    print(f"Ring: {ring}")
    print(f"Input polynomials:\n{polys}")
    polys.append(ring(0))
    t = walltime()
    if solving_deg is None:
        d = 0
        is_gb = False
        while not is_gb:
            print(f"\n--- Degree {d} ---")
            print("Computing all monomials up to degree:", d)
            s = walltime()
            monoms = vector(all_monoms_up_to_deg(ring, d))
            print("Time needed:", walltime(s))
            print("Computing Macaulay matrix.")
            polys_m = flatten([[poly * mon for poly in polys if (poly * mon).degree() <= d] for mon in monoms])
            if polys_m == []:
                polys_m = [ring(0)]
            M, v = PolynomialSequence(polys_m).coefficients_monomials()
            s = walltime()
            print("Time needed:", walltime(s))
            print("Performing Gaussian Elimination.")
            s = walltime()
            M = M.echelon_form()
            print("Time needed:", walltime(s))
            gb = [poly for poly in M * v if poly != ring(0)]
            if len(gb) > 0:
                gb_red = reduced_gb(gb)
            else:
                gb_red = gb
            quos_rems = [polynomial_division(p, gb_red) for p in polys]
            rems = [r for _, r in quos_rems]
            is_gb = buchberger_criterion(gb_red) and not any(rems)
            if print_intermediate_basis:
                print(f"Without redundancies:\n{gb_red}")
            print(f"Is Groebner Basis: {is_gb}")
            if not is_gb:
                d += 1
    else:
        d = solving_deg
        print(f"\n--- Degree {d} ---")
        print("Computing all monomials up to degree:", d)
        s = walltime()
        monoms = vector(all_monoms_up_to_deg(ring, d))
        print("Time needed:", walltime(s))
        print("Computing Macaulay matrix.")
        s = walltime()
        polys_m = flatten([[poly * mon for poly in polys if (poly * mon).degree() <= d] for mon in monoms])
        if polys_m == []:
            polys_m = [ring(0)]
        M, v = PolynomialSequence(polys_m).coefficients_monomials()
        print("Time needed:", walltime(s))
        print("Performing Gaussian Elimination.")
        s = walltime()
        M = M.echelon_form()
        print("Time needed:", walltime(s))
        gb = [poly for poly in M * v if poly != ring(0)]
        gb_red = reduced_gb(gb)
        if print_intermediate_basis:
            print(f"Without redundancies:\n{gb_red}")
    print("Total time needed:", walltime(t))
    return gb_red
