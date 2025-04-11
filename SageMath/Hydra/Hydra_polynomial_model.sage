load("Hydra.sage")
load("./utilities.sage")

def generate_Hydra_variables_m_samples(rounds_head, m):
    """
    Generates variables for a Hydra instance as strings.

    INPUT:
    "rounds_head" -- Number of rounds for heads.
    "m" -- Number of samples.

    OUTPUT:
    List of variables as strings.
    """
    variables = []

    for i in range(0, 4):
        variables.append("y_b" + str(i + 1))

    for i in range(0, 4):
        variables.append("z_b" + str(i + 1))
    
    for i in range(0, rounds_head - 1):
        for j in range(0, 8):
            variables.append("x_s" + str(1) + "_b" + str(j + 1) + "_r" + str(i + 1))
    
    for i in range(1, m):
        for j in range(0, rounds_head):
            for k in range(0, 8):
                variables.append("x_s" + str(i + 1) + "_b" + str(k + 1) + "_r" + str(j))
    
    for i in range(0, 4):
        variables.append("k_b" + str(i + 1))

    return variables

def generate_Hydra_polynomials_m_samples(hydra,
                                         m=2,
                                         nonce=None,
                                         samples=None,
                                         termorder="degrevlex",
                                         field_equations=False,
                                         info_level=0):
    """
    Generates a polynomial model for m Hydra heads.

    INPUT:
    "hydra" -- A Hydra instance. If no argument is supplied,
               then a random instance is generated.
    "m" -- Number of samples.
    "nonce" -- An integer/field element."
               If no nonce is provided a random one is generated.
    "samples" -- Samples of Hydra heads.
                 Expects a list of 8 * m integers/field elements.
                 If no samples are provided, new once are generated
                 with the nonca and a random master key.
    "termorder" -- Termorder for the polynomial ring as string.
                   E.g., "lex", "degrevlex".
    "field_equations" -- Boolean value, if set to true field equations
                         are added to the Ciminion polynomials.
                         Default value "False".
    
    OUTPUT:
    Hydra polynomial system.
    """
    print_key = False
    if nonce is None:
        nonce = hydra.field.random_element()
    
    if samples is None:
        print_key = True
        key = [hydra.field.random_element() for i in range(0, 4)]
        samples = hydra.key_stream(nonce, key, samples=m)
    else:
        m = int(len(samples) / 8)
    
    if info_level > 0:
        print("Nonce:", nonce)
        if print_key:
            print("Key:", key)
        print("Number of samples:", m)
        print("Samples:", samples)
        print("Term order:", termorder)
    
    variables = generate_Hydra_variables_m_samples(hydra.rounds_head, m)
    P = PolynomialRing(hydra.field, variables, order=termorder)
    variables = list(P.gens())

    polynomials = []
    N = m * 8 * hydra.rounds_head
    key_variables = vector(P, variables[N:N + 4])
    y_z = vector(P, variables[0:8])
    variables = variables[8:N]
    
    large_key = vector(P, list(key_variables) + list(hydra.matrix_body_E * key_variables))

    # First head
    current_state = y_z
    next_state = vector(P, variables[0:8])
    v_in = current_state
    polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[0]) - next_state
    polynomials = polynomials + list(polys)
    for i in range(1, hydra.rounds_head - 1):
        current_state = next_state
        next_state = vector(P, variables[8 * i:8 * (i + 1)])
        polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[i]) - next_state
        polynomials = polynomials + list(polys)
    current_state = next_state
    next_state = vector(hydra.field, samples[0:8])
    polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[hydra.rounds_head - 1]) + v_in - next_state
    polynomials = polynomials + list(polys)

    # Other heads
    N = 8 * (hydra.rounds_head - 1)
    rolling_state = y_z
    for j in range(1, m):
        # Rolling function
        current_state = rolling_state
        next_state = vector(P, variables[N:N + 8])
        polys = hydra.rolling_function(current_state, hydra.matrix_rolling_function) - next_state
        polynomials = polynomials + list(polys)
        rolling_state = next_state
        N += 8

        # Head
        current_state = next_state
        next_state = vector(P, variables[N:N + 8])
        v_in = current_state
        polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[0]) - next_state
        polynomials = polynomials + list(polys)
        for i in range(1, hydra.rounds_head - 1):
            current_state = next_state
            next_state = vector(P, variables[N + 8 * i:N + 8 * (i + 1)])
            polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[i]) - next_state
            polynomials = polynomials + list(polys)
        current_state = next_state
        next_state = vector(hydra.field, samples[8 * j:8 * (j + 1)])
        polys = hydra.head_round_function(current_state, large_key, hydra.matrix_head, hydra.constants_head[hydra.rounds_head - 1]) + v_in - next_state
        polynomials = polynomials + list(polys)
        N += 8 * (hydra.rounds_head - 1)
    
    if field_equations:
        polynomials = polynomials + generate_field_equations(variables)

    return polynomials

def lm_transformation_matrix(field, n):
    """
    Generates a matrix to eliminate leading monomials.

    INPUT:
    "field" -- A field.
    "n" -- Dimension of the matrix.

    OUTPUT:
    A nxn matrix.
    """
    mat = identity_matrix(field, n)
    for i in range(0, n - 1):
        mat[i, n - 1] -= 1
    return mat

def transform_Hydra_polynomial_system(hydra, hydra_polys, m):
    """
    Eliminates variables in Hydra polynomial system via affine
    equations.

    INPUT:
    "hydra" -- A Hydra instance.
    "hydra_polys" -- A Hydra polynomial systems.
                     Expects them to be in the shape of the
                     "generate_Hydra_polynomials_m_samples"
                     output.
    "m" -- Number of Hydra samples.

    OUTPUT:
    Hydra polynomial system, where variables in non-linear 
    polynomials have been eliminated with affine equations.
    """
    A = lm_transformation_matrix(hydra.field, 8)
    B = lm_transformation_matrix(hydra.field, 4)
    B = block_diagonal_matrix([B, B])

    polys_transformed = []

    mat_head_inv = hydra.matrix_head.inverse()
    mat_rol_inv = hydra.matrix_rolling_function.inverse()
    for i in range(0, m * hydra.rounds_head + m - 1):
        tmp_polys = A * mat_head_inv * vector(hydra_polys[8 * i:8 * (i + 1)])
        degs = [poly.degree() for poly in tmp_polys]
        if max(degs[0:7]) != 1 or min(degs[0:7]) != 1:
            tmp_polys = B * mat_rol_inv * vector(hydra_polys[8 * i:8 * (i + 1)])
            degs = [poly.degree() for poly in tmp_polys]
        if max(degs[0:3]) != 1 or min(degs[4:7]) != 1:
            raise Exception("Transformation of Hydra polynomials failed.")
        polys_transformed = polys_transformed + list(tmp_polys)
    
    return polys_transformed
    

def is_in_generic_coordinates_Hydra_polynomials_m_samples(hydra, hydra_polys=None, m=2):
    """
    Generic coordinates verification for Hydra polynomial system.

    INPUT:
    "hydra" -- A Hydra instance.
    "hydra_polys" -- A Hydra polynomial systems.
                     Expects them to be in the shape of the
                     "generate_Hydra_polynomials_m_samples"
                     output.
    "m" -- Number of Hydra samples.

    OUTPUT:
    True or False, whether polynomial system is in generic coordinates or test fails.
    """
    if hydra_polys is None:
        hydra_polys = generate_Hydra_polynomials_m_samples(hydra=hydra, m=m)
    P = hydra_polys[0].parent()
    variables = P.gens()
    
    B = lm_transformation_matrix(hydra.field, 4)
    B = block_diagonal_matrix([B, B])

    polys_transformed = transform_Hydra_polynomial_system(hydra, hydra_polys, m)
    affine_polys = [highest_degree_component(poly) for poly in polys_transformed if poly.degree() == 1]
    
    for i in range(0, m * hydra.rounds_head):
        tmp_poly = P(0)
        for k in range(0, 8):
            tmp_poly += (-1)**floor(k / 4) * variables[8 * i:8 * (i + 1)][k]
        affine_polys.append(tmp_poly)

    M, v = Sequence(affine_polys).coefficients_monomials()

    if M.rank() == len(variables):
        return True
    
    return False

def number_of_non_linear_variables_Hydra_polynomial_system(hydra, hydra_polys, m, transformed=False):
    """
    Generic coordinates verification for Hydra polynomial system.

    INPUT:
    "hydra" -- A Hydra instance.
    "hydra_polys" -- A Hydra polynomial systems.
                     Expects them to be in the shape of the
                     "generate_Hydra_polynomials_m_samples"
                     output.
    "m" -- Number of Hydra samples.
    "transformed" -- Boolean value, whether Hydra polynomial system has already been 
                     transformed or not.
                     Default value set to "False".

    OUTPUT:
    Number of variables in Hydra polynomial system that cannot be eliminated
    via affine equations
    """
    n_vars = len(hydra_polys[0].parent().gens())
    if not transformed:
        polys_transformed = transform_Hydra_polynomial_system(hydra, hydra_polys, m)
    else:
        polys_transformed = polys
    affine_polys = [poly for poly in polys_transformed if poly.degree() == 1]
    M, v = Sequence(affine_polys).coefficients_monomials()
    return n_vars - M.rank()

def non_linear_variable_substitution_Hydra_polynomial_system(hydra, hydra_polys, m, transformed=False, info_level=1):
    """
    Variable substitution for downsized Hydra polynomial system.

    INPUT:
    "hydra" -- A Hydra instance.
    "hydra_polys" -- A Hydra polynomial systems.
                     Expects them to be in the shape of the
                     "generate_Hydra_polynomials_m_samples"
                     output.
    "m" -- Number of Hydra samples.
    "transformed" -- Boolean value, whether Hydra polynomial system has already been 
                     transformed or not.
                     Default value set to "False".
    "info_level" -- Integer, if greater than zero, then Ciminion parameters are
                        printed in console.

    OUTPUT:
    "affine_polys" -- Affine polynomials coming from the transformed Hydra polynomial system.
    "polys_subs" -- Affine polynomials used for the variable substitution.
    "polys_downsized_subs" -- Substituted non-linear polynomials in the transformed Hydra polynomial system.
    """
    if m > 2:
        raise Exception("Number of samples m > 2 not implemented.")

    P = hydra_polys[0].parent()
    variables = list(P.gens())
    n_vars = len(variables)

    n_non_lin = number_of_non_linear_variables_Hydra_polynomial_system(hydra, hydra_polys, m, transformed=transformed)
    variables_subs = generate_Hydra_variables_m_samples(hydra.rounds_head, m)
    variables_subs = variables_subs + ["x_subs_" + str(i + 1) for i in range(0, n_non_lin)]
    Q_subs = PolynomialRing(hydra.field, variables_subs, order=polys[0].parent().term_order())
    variables_subs = vector(Q_subs.gens()[n_vars:])

    if not transformed:
        polys_transformed = transform_Hydra_polynomial_system(hydra, hydra_polys, m)
    else:
        polys_transformed = hydra_polys

    affine_polys = Sequence([poly for poly in polys_transformed if poly.degree() == 1])
    M, v = affine_polys.coefficients_monomials()
    M = M.echelon_form()
    affine_polys = ideal((M * v).list())

    polys_downsized = [poly.reduce(affine_polys) for poly in polys_transformed if poly.degree() > 1]

    polys_subs = []
    for i in range(0, m * hydra.rounds_head):
        tmp_poly = P(0)
        for k in range(0, 8):
            tmp_poly += (-1)**floor(k / 4) * variables[8 * i:8 * (i + 1)][k]
        polys_subs.append(tmp_poly.reduce(affine_polys))

    M, v = Sequence(polys_subs).coefficients_monomials()

    if M[2:].rank() != n_non_lin:
       raise Exception("Variable substitution failed.")
    
    polys_subs = vector(Q_subs, vector(M[2:] * v)) - variables_subs
    M, v = Sequence(polys_subs).coefficients_monomials()
    M = M.echelon_form()
    polys_subs = ideal((M * v).list())
    
    polys_downsized_subs = [Q_subs(poly).reduce(polys_subs) for poly in polys_downsized]
    polys_downsized_subs = [poly / poly.lc() for poly in polys_downsized_subs if poly != Q_subs(0)]

    if info_level > 0:
        lms = [poly.lm() for poly in polys_downsized_subs]
        vars_squ = [var**2 for var in variables_subs]
        zero_dimensional = False
        if set(vars_squ) == set(lms):
            zero_dimensional = True
        substituion_succcess = set.union(*[set(poly.variables()) for poly in polys_downsized_subs]) <= set(variables_subs)
        print("Number of non-linear variables:", n_non_lin)
        print("Number of polynomials in substituted downsized Hydra polynomial system:", len(polys_downsized_subs))
        print("(x_subs_1^2, ..., x_subs_n^2) contained in leading terms of substituted polynomials:", zero_dimensional)
        print("All terms of substituted polynomial system contained in (x_subs_1, ..., x_subs_n):", substituion_succcess)

    affine_polys = [Q_subs(poly) for poly in affine_polys.gens()]
    polys_subs = list(polys_subs.gens())

    return affine_polys, polys_subs, polys_downsized_subs
