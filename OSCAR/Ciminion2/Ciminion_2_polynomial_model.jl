using Oscar
include("Ciminion_2.jl")

function generate_Ciminion_2_variables(rounds_C::Int64, rounds_E::Int64)
    """
    Generates variables for a Ciminion instance as strings.

    INPUT:
    "rounds_C" -- Rounds of C permutation.
    "rounds_E" -- Rounds of E permutation.

    OUTPUT:
    Vector of variables as strings.
    """
    variables = String[]
    for i in 1:rounds_C - 1
        for j in 1:3
            push!(variables, "x_C_" * string(j) * "__" * string(i))
        end
    end

    for i in 1:rounds_E
        for j in 1:3
            push!(variables, "x_E_" * string(j) * "__" * string(i))
        end
    end
    
    push!(variables, "y_1")
    push!(variables, "y_2")
    push!(variables, "x")
    
    return variables
end

function generate_Ciminion_2_polynomials(;ciminion_2::Ciminion_2=Ciminion_2_constructor(),
                                        plain=nothing,
                                        cipher=nothing,
                                        nonce=nothing,
                                        termorder="degrevlex",
                                        info_level=1)
    """
    Generates a polynomial model for the first key pair of Ciminion_2.

    INPUT:
    "ciminion_2" -- A Ciminion_2 instance.
    "plain" -- 2x1 matrix of integers/field elements.
               If no plaintext is provided a random one is generated.
    "cipher" -- 2x1 matrix of integers/field elements.
                If no ciphertext is provided a random one is generated.
    "nonce" -- An integer/field element.
               If no nonce is provided a random one is generated.
    "termorder" -- Termorder for the polynomial ring as string.
                   E.g., "lex", "degrevlex".
    "field_equations" -- Boolean value, if set to true field equations
                         are added to the Ciminion_2 polynomials.
                         Default value "false".
    "info_level" -- Integer, if greater than zero, then parameters are
                    printed in console.
    
    OUTPUT:
    Ciminion_2 polynomial system.
    """
    print_key = false
    if isnothing(nonce)
        nonce = rand(ciminion_2.field)
    end
    if isnothing(plain)
        plain = zero_matrix(ciminion_2.field, 2, 1)
        plain[1, 1] = rand(ciminion_2.field)
        plain[2, 1] = rand(ciminion_2.field)
    end
    if isnothing(cipher)
        print_key = true
        key = zero_matrix(ciminion_2.field, 2, 1)
        key[1, 1] = rand(ciminion_2.field)
        key[2, 1] = rand(ciminion_2.field)
        cipher = encrypt(plain, key, nonce, ciminion_2)
    end

    if info_level > 1
        println("Plaintext: ", plain)
        if print_key
            println("Key: ", key)
        end
        println("Ciphertext: ", cipher)
        println("Nonce: ", nonce)
        println("Term order: ", termorder)
    end

    variables = generate_Ciminion_2_variables(ciminion_2.rounds_C, ciminion_2.rounds_E)
    if termorder == "degrevlex"
        P, variables = polynomial_ring(ciminion_2.field, variables, internal_ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = polynomial_ring(ciminion_2.field, variables, internal_ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end
    
    polynomials = Vector{typeof(variables[1])}()
    N = 3 * (ciminion_2.rounds_C + ciminion_2.rounds_E - 1)
    key_variables = variables[N + 1: N + 2]
    auxiliary_variable = variables[N + 3]
    variables_C = variables[1:3 * (ciminion_2.rounds_C - 1)]
    variables_E = variables[3 * (ciminion_2.rounds_C - 1) + 1:N]

    current_state = zero_matrix(P, 3, 1)
    current_state[1, 1] = nonce
    current_state[2, 1] = key_variables[1]
    current_state[3, 1] = key_variables[2]
    next_state = matrix(variables_C[1:3])
    polys = round_function(current_state, 
                           matrix(ciminion_2.constants_C[:, 1])) - next_state
    polynomials = [polynomials; vec(polys[:, 1])]
    for i in 2:ciminion_2.rounds_C - 1
        current_state = next_state
        next_state = matrix(variables_C[3 * (i - 1) + 1: 3 * i])
        polys = round_function(current_state, 
                               matrix(ciminion_2.constants_C[:, i]),
                               key=key_variables) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
    end
    current_state = next_state
    next_state = matrix(variables_E[1:3])
    polys = round_function(current_state, 
                           matrix(ciminion_2.constants_C[:, ciminion_2.rounds_C]),
                           key=key_variables) - next_state
    polynomials = [polynomials; vec(polys[:, 1])]
    current_state = next_state
    next_state = matrix(variables_E[4:6])
    polys = round_function(current_state, 
                           matrix(ciminion_2.constants_E[:, 1]),
                           key=key_variables) - next_state
    polynomials = [polynomials; vec(polys[:, 1])]
    for i in 2:ciminion_2.rounds_E - 1
        current_state = next_state
        next_state = matrix(variables_E[3 * i + 1:3 * (i + 1)])
        polys = round_function(current_state, 
                               matrix(ciminion_2.constants_E[:, i]),
                               key=key_variables) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
    end
    current_state = next_state
    next_state = zero_matrix(P, 3, 1)
    next_state[1, 1] = -plain[1, 1] + cipher[1, 1]
    next_state[2, 1] = -plain[2, 1] + cipher[2, 1]
    next_state[3, 1] = auxiliary_variable
    polys = round_function(current_state, 
                           matrix(ciminion_2.constants_E[:, ciminion_2.rounds_E]),
                           key=key_variables) - next_state
                           polynomials = [polynomials; vec(polys[:, 1])]

    return polynomials
end
