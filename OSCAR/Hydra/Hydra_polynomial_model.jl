using Oscar
include("Hydra.jl")
include("utilities.jl")

function generate_Hydra_variables_m_samples(rounds_head::Int64, samples::Int64)
    """
    Generates variables for a Ciminion instance as strings.

    INPUT:
    "rounds_head" -- Number of rounds for heads.
    "m" -- Number of samples.

    OUTPUT:
    Vector of variables as strings.
    """
    variables = String[]

    for i in 1:4
        push!(variables, "y_b" * string(i))
    end

    for i in 1:4
        push!(variables, "z_b" * string(i))
    end
    
    for i in 1:rounds_head - 1
        for j in 1:8
            push!(variables, "x_s" * string(1) * "_b" * string(j) * "_r" * string(i))
        end
    end

    for i in 2:samples
        for j in 1:rounds_head
            for k in 1:8
                push!(variables, "x_s" * string(i) * "_b" * string(k) * "_r" * string(j - 1))
            end
        end
    end

    for i in 1:4
        push!(variables, "k_b" * string(i))
    end

    return variables
end

function generate_Hydra_polynomials_m_samples(;hydra::Hydra=Hydra_constructor(),
                                               m=2,
                                               nonce=nothing,
                                               samples=nothing,
                                               termorder="degrevlex",
                                               field_equations=false,
                                               info_level=1)
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
    "info_level" -- Integer, if greater than zero, then Ciminion parameters are
                    printed in console.
    
    OUTPUT:
    Hydra polynomial system.
    """
    print_key = false
    if isnothing(nonce)
        nonce = rand(hydra.field)
    end

    if isnothing(samples)
        print_key = true
        key = matrix(map(x -> rand(hydra.field), 1:4))
        samples = key_stream(nonce, key, hydra, m=m)
    else
        m = Int64(length(samples) / 8)
    end

    if info_level > 0
        println("Nonce: ", nonce)
        if print_key
            println("Key: ", key)
        end
        println("Number of samples: ", m)
        println("Samples:\n", samples)
        println("Term order: ", termorder)
    end

    variables = generate_Hydra_variables_m_samples(hydra.rounds_head, m)
    
    if termorder == "degrevlex"
        P, variables = polynomial_ring(hydra.field, variables, internal_ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = polynomial_ring(hydra.field, variables, internal_ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end

    polynomials = Vector{typeof(variables[1])}()

    N = m * 8 * hydra.rounds_head
    key_variables = matrix(variables[N + 1:N + 4])
    y_z = matrix(variables[1:8])
    variables = variables[8 + 1:N]
    
    large_key = zero_matrix(P, 8, 1)
    tmp = hydra.matrix_body_E * key_variables
    for i in 1:4
        large_key[i, 1] = key_variables[i, 1]
        large_key[i + 4, 1] = tmp[i, 1]
    end

    # First head
    current_state = y_z
    next_state = matrix(variables[1:8])
    v_in = current_state
    polys = head_round_function(current_state, 
                                large_key, 
                                hydra.matrix_head, 
                                matrix(hydra.constants_head[:, 1])) - next_state
    polynomials = [polynomials; vec(polys[:, 1])]
    for i in 2:hydra.rounds_head - 1
        current_state = next_state
        next_state = matrix(variables[8 * (i - 1) + 1:8 * i])
        polys = head_round_function(current_state, 
                                    large_key, 
                                    hydra.matrix_head, 
                                    matrix(hydra.constants_head[:, i])) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
    end
    current_state = next_state
    next_state = matrix(samples[1:8])
    polys = head_round_function(current_state, 
                                large_key, 
                                hydra.matrix_head, 
                                matrix(hydra.constants_head[:, hydra.rounds_head])) + v_in - next_state
    polynomials = [polynomials; vec(polys[:, 1])]
    
    # Other heads
    N = 8 * (hydra.rounds_head - 1)
    rolling_state = y_z
    for j in 2:m
        # Rolling function
        current_state = rolling_state
        next_state = matrix(variables[N + 1:N + 8])
        polys = rolling_function(current_state, hydra.matrix_rolling_function) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        rolling_state = next_state
        N += 8

        # Head
        current_state = next_state
        next_state = matrix(variables[N + 1:N + 8])
        v_in = current_state
        polys = head_round_function(current_state, 
                                    large_key, 
                                    hydra.matrix_head, 
                                    matrix(hydra.constants_head[:, 1])) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        for i in 2:hydra.rounds_head - 1
            current_state = next_state
            next_state = matrix(variables[N + 8 * (i - 1) + 1:N + 8 * i])
            polys = head_round_function(current_state, 
                                        large_key, 
                                        hydra.matrix_head, 
                                        matrix(hydra.constants_head[:, i])) - next_state
            polynomials = [polynomials; vec(polys[:, 1])]
        end
        current_state = next_state
        next_state = matrix(samples[8 * (j - 1) + 1:8 * j])
        polys = head_round_function(current_state, 
                                    large_key, 
                                    hydra.matrix_head, 
                                    matrix(hydra.constants_head[:, hydra.rounds_head])) + v_in - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        N += 8 * (hydra.rounds_head - 1)
    end

    if field_equations
        polynomials = [polynomials; vec(generate_field_equations(gens(P)))]
    end

    return polynomials
end

function lm_transformation_matrix(field, n)
    """
    Generates a matrix to eliminate leading monomials.

    INPUT:
    "field" -- A field.
    "n" -- Dimension of the matrix.

    OUTPUT:
    A nxn matrix.
    """
    mat = zero_matrix(field, n, n)
    for i in 1:n - 1
        mat[i, i] += one(field)
        mat[i, n] -= one(field)
    end
    mat[n, n] += one(field)
    return mat
end

function transform_Hydra_polynomial_system(hydra::Hydra, hydra_polys, m)
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
    tmp = zero_matrix(hydra.field, 4, 4)
    B = [B tmp; tmp B]

    polys_transformed = Vector{typeof(hydra_polys[1])}()

    mat_head_inv = inv(hydra.matrix_head)
    mat_rol_inv = inv(hydra.matrix_rolling_function)
    for i in 1:m * hydra.rounds_head + m - 1
        tmp_polys = A * mat_head_inv * matrix(hydra_polys[8 * (i - 1) + 1:i * 8])
        degs = vec(map(total_degree, tmp_polys[:, 1]))
        if maximum(degs[1:7]) != 1 || minimum(degs[1:7]) != 1
            tmp_polys = B * mat_rol_inv * matrix(hydra_polys[8 * (i - 1) + 1:i * 8])
            degs = vec(map(total_degree, tmp_polys[:, 1]))
        end
        if maximum(degs[1:3]) != 1 || minimum(degs[5:7]) != 1
            error("Transformation of Hydra polynomials failed.")
        end
        polys_transformed = [polys_transformed; vec(tmp_polys[:, 1])]
    end
        
    return polys_transformed
end

function number_of_non_linear_variables_Hydra_polynomial_system(hydra, 
                                                                hydra_polys, 
                                                                m; 
                                                                transformed=false)
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
    P = parent(hydra_polys[1])
    n_vars = length(gens(P))
    if !transformed
        polys_transformed = transform_Hydra_polynomial_system(hydra, hydra_polys, m)
    else
        polys_transformed = hydra_polys
    end
    affine_polys = filter(poly -> total_degree(poly) == 1, polys_transformed)
    gb = groebner_basis_f4(ideal(affine_polys))
    return n_vars - length(gb)
end

function non_linear_variable_substitution_Hydra_polynomial_system(hydra::Hydra, 
                                                                  hydra_polys,
                                                                  m; 
                                                                  transformed=false, 
                                                                  info_level=1)
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
    if m > 2
        println("Number of samples m > 2 not implemented.")
        return
    end

    P = parent(hydra_polys[1])
    variables = gens(P)
    n_vars = length(variables)

    n_non_lin = number_of_non_linear_variables_Hydra_polynomial_system(hydra::Hydra,
                                                                       hydra_polys,
                                                                       m;
                                                                       transformed=transformed)
    variables_subs = generate_Hydra_variables_m_samples(hydra.rounds_head, m)
    variables_subs = [variables_subs; map(i -> "x_subs_i" * string(i), 1:n_non_lin)]
    Q_subs, variables_subs = polynomial_ring(base_ring(P), variables_subs, internal_ordering=:degrevlex)
    induce(variables_subs, degrevlex(variables_subs))
    phi = hom(P, Q_subs, variables_subs[1:n_vars])
    variables_subs = variables_subs[n_vars + 1:length(variables_subs)]

    if !transformed
        polys_transformed = transform_Hydra_polynomial_system(hydra, hydra_polys, m)
    else
        polys_transformed = hydra_polys
    end

    affine_polys = filter(poly -> total_degree(poly) == 1, polys_transformed)
    affine_polys = gens(groebner_basis_f4(ideal(affine_polys)))

    polys_downsized = filter(poly -> total_degree(poly) > 1, polys_transformed)
    polys_downsized = map(poly -> divrem(poly, affine_polys)[2], polys_downsized)

    polys_subs = Vector{typeof(variables_subs[1])}()
    for i in 1:m * hydra.rounds_head
        tmp_poly = zero(P)
        for k in 1:8
            tmp_poly += (-1)^Int64(floor((k - 1) / 4)) * variables[8 * (i - 1) + 1:8 * i][k]
        end
        push!(polys_subs, tmp_poly)
    end

    gb_subs = groebner_basis_f4(ideal(polys_subs))
    if length(gb_subs) < n_non_lin
        println("Variable substitution failed.")
    end

    affine_polys = map(phi, affine_polys)

    polys_subs = matrix(map(phi, polys_subs[3:length(polys_subs)])) - matrix(variables_subs)
    polys_subs = vec(polys_subs[:, 1])
    polys_subs = map(poly -> divrem(poly, affine_polys)[2], polys_subs)
    polys_subs = gens(groebner_basis_f4(ideal(polys_subs)))

    polys_downsized_subs = map(poly -> divrem(phi(poly), polys_subs)[2], polys_downsized)

    if info_level > 0
        lms = map(leading_monomial, polys_downsized_subs)
        vars_squ = map(var -> var^2, variables_subs)
        zero_dimensional = false
        if ideal(lms) == ideal(vars_squ)
            zero_dimensional = true
        end
        substitution_success = true
        I_subs = map(var -> var - one(K), variables_subs)
        for poly in polys_downsized_subs
            for mon in monomials(poly)
                tmp = divrem(mon, I_subs)[2]
                if tmp != one(Q_subs)
                    substitution_success = false
                end
            end
            if !substitution_success
                break
            end
        end
        println("Number of non-linear variables: ", n_non_lin)
        println("Number of polynomials in substituted downsized Hydra polynomial system: ", length(polys_downsized_subs))
        println("(x_subs_1^2, ..., x_subs_n^2) contained in leading terms of substituted polynomials: ", zero_dimensional)
        println("All terms of donwiszed polynomial system contained in (x_subs_1, ..., x_subs_n): ", substitution_success)
    end

    return affine_polys, polys_subs, polys_downsized_subs
end
