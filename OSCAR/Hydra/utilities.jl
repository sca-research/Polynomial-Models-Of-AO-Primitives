using Oscar

function generate_field_equations(variables::Union{Array{FqMPolyRingElem, 1}, AbstractAlgebra.Generic.MatSpaceElem{FqMPolyRingElem}})
    """
    Generates the field equations for a given variable set.

    INPUT:
    "variables" -- variables of a polynomial ring over a finite field.

    OUTPUT:
    "field_equations" -- field equations for given variable set.
    """
    q = order(base_ring(variables[1]))
    field_equations = map(var -> var^q - var, variables)
    return field_equations
end

function highest_degree_component(poly::FqMPolyRingElem)
    """
    Extracts the highest degree component of a polynomial.

    INPUT:
    "poly" -- a polynomial.

    OUTPUT:
    "poly_top" -- highest degree component of a polynomial.
    """
    poly_top = 0
    d = total_degree(poly)
    for term in terms(poly)
        if total_degree(term) == d
            poly_top += term
        end
    end
    return poly_top
end

function circulant_matrix(first_row)
    """
    Generates a right circulant matrix.

    INPUT:
    "first_row" -- first row of the matrix.

    OUTPUT:
    "circ_mat" -- a right circulant matrix.
    """
    first_row = Vector(first_row)
    R = parent(first_row[1])
    n = length(first_row)
    circ_mat = zero_matrix(R, n, n)
    row = circshift(first_row, -1)
    for i in 1:n
        row = circshift(row, 1)
        for j in 1:n
            circ_mat[i, j] = row[j]
        end
    end
    return circ_mat
end
