using Oscar

struct Poseidon
    """
    Poseidon parameters.

    "field" -- A prime field.
    "d" -- Exponent of power permutation.
    "n" -- Number of blocks.
    "r_f" -- Number of full rounds in beginning and end of
             Poseidon.
    "r_p" -- Number of partial rounds.
    "matrix" -- Invertible matrix over prime field.
    "constants" -- Round constants of Poseidon.
    """
    field::FqField
    d::Int64
    n::Int64
    r_f::Int64
    r_p::Int64
    matrix::FqMatrix
    constants::FqMatrix
end

function Poseidon_constructor(;field=GF(2^31 - 1),
                               n=16,
                               r_f=4,
                               r_p=14,
                               mat=nothing,
                               constants=nothing,
                               info_level=1)
    """
    Initializes a Poseidon instance.
    Computes the least integer d which generates a power
    permutation over given prime field.
    The Poseidon matrix is applied in full as well as partial rounds.

    INPUT:
    "field" -- A prime field. Default prime p = 2^31 - 1.
    "n" -- Number of blocks. Default value n = 16.
    "r_f" -- Number of full rounds. Default value r_f = 4.
    "r_p" -- Number of partial rounds. Default value rp = 14.
    "mat" -- A n x n matrix over the prime field.
             If no matrix is provided a random one is generates.
    "constants" -- A n x (2 * r_f + r_p) matrix over the prime field.
                   If no matrix is provided a random one is generates.
    "info_level" -- Parameter to print information if > 0.

    OUTPUT: 
    Poseidon structure.
    """
    q = order(field)
    d = 2
    while gcd(d, q - 1) != 1
        d += 1 
    end

    if isnothing(mat)
        M = matrix_space(field, n, n)
        mat = rand(M)
    end

    if isnothing(constants)
        R = 2 * r_f + r_p
        constants = zero_matrix(field, n, R)
        for i in 1:R
            for j in 1:n
                constants[j, i] = rand(field)
            end
        end
    end

    if info_level > 0
        println("Poseidon parameters")
        println("Field: ", field)
        println("d: ", d)
        println("n: ", n)
        println("r_f: ", r_f)
        println("r_p: ", r_p)
        println("Matrix:\n", mat)
        println("Constants:\n", constants)
    end

    return Poseidon(field, d, n, r_f, r_p, mat, constants)
end

function full_round(poseidon::Poseidon, x_in, constants)
    """
    Full Poseidon round.

    INPUT:
    "poseidon" -- A Poseidon structure.
    "x_in" -- A n x 1 matrix.
    "constants" -- A n x 1 matrix.

    OUTPUT:
    A n x 1 matrix.
    """
    x_out = matrix(map(i -> x_in[i, 1]^poseidon.d, 1:poseidon.n))
    x_out = poseidon.matrix * x_out
    x_out += constants
    return x_out
end

function partial_round(poseidon::Poseidon, x_in, constants)
    """
    Partial Poseidon round.

    INPUT:
    "poseidon" -- A Poseidon structure.
    "x_in" -- A n x 1 matrix.
    "constants" -- A n x 1 matrix.

    OUTPUT:
    A n x 1 matrix.
    """
    x_out = deepcopy(x_in)
    x_out[1, 1] = x_out[1, 1]^poseidon.d
    x_out = poseidon.matrix * x_out
    x_out += constants
    return x_out
end

function permutation(poseidon::Poseidon, x_in)
    """
    Poseidon permutation.

    INPUT:
    "poseidon" -- A Poseidon structure.
    "x_in" -- A n x 1 matrix.

    OUTPUT:
    A n x 1 matrix.
    """
    x_out = poseidon.matrix * x_in
    counter = 1
    for _i in 1:poseidon.r_f
        x_out = full_round(poseidon, 
                           x_out,
                           matrix(poseidon.constants[:, counter]))
        counter += 1
    end
    for _i in 1:poseidon.r_p
        x_out = partial_round(poseidon, 
                              x_out,
                              matrix(poseidon.constants[:, counter]))
        counter += 1
    end
    for _i in 1:poseidon.r_f
        x_out = full_round(poseidon, 
                           x_out,
                           matrix(poseidon.constants[:, counter]))
        counter += 1
    end
    return x_out
end

function compression(poseidon::Poseidon, x_in)
    """
    Poseidon permutation truncation mode.
    Expects that 2 | n.

    INPUT:
    "poseidon" -- A Poseidon structure.
    "x_in" -- A n x 1 matrix.

    OUTPUT:
    A vector of n / 2 field elements.
    """
    x_out = deepcopy(x_in)
    x_out = permutation(poseidon, x_out)
    x_out += x_in
    return vec(x_out[:, 1])[1:Int64(poseidon.n / 2)]
end

function generate_variables(poseidon::Poseidon)
    """
    Generates variables for Poseidon permutation 
    truncation mode.

    INPUT:
    "poseidon" -- A Poseidon structure.

    OUTPUT:
    A vector of vectors of strings.
    """
    variables = [String[]]
    push!(variables, map(i -> "x_in_" * string(i), 1:poseidon.n))
    for i in 1:(2 * poseidon.r_f + poseidon.r_p - 1)
        push!(variables, map(j -> "x_" * string(i) * "_" * string(j), 1:poseidon.n))
    end
    push!(variables, map(i -> "x_out_" * string(i), 1:Int64(poseidon.n / 2)))
    return variables
end

function generate_compression_polynomials(poseidon::Poseidon;
                                          compr_val=nothing,
                                          substitution=true,
                                          info_level=1)
    """
    Generate iterated polynomial model for 
    Poseidon permutation truncation mode.

    INPUT:
    "poseidon" -- A Poseidon structure.
    "compr_val" -- A vector of n / 2 field elements.
                  If no vector is provided a random plaintext
                  is generated and the compression value is computed.
    "substitution" -- Boolean parameter to perform the variables substitution
                      z_in = M x_in before first round. Default is true.
    "info_level" -- Parameter to print information if > 0.
    
    OUTPUT:
    A vector of polynomials.
    """
    print_plain = false
    if isnothing(compr_val)
        print_plain = true
        plain = matrix(map(i -> rand(poseidon.field), 1:poseidon.n))
        compr_val = compression(poseidon, plain)
    end

    if info_level > 0
        if print_plain
            println("Preimage: ", plain)
        end
        println("Compression value: ", compr_val)
        println("substitution: ", substitution)
    end

    variables = generate_variables(poseidon)
    N = length(variables)
    R = 2 * poseidon.r_f + poseidon.r_p
    n_half = Int64(poseidon.n / 2)
    if substitution
        variables_2 = deepcopy(variables)
        variables_2 = [variables_2[2:N]; variables_2[1:1]]
        P, variables_2 = polynomial_ring(poseidon.field, collect(Iterators.Base.Flatten(variables_2)), internal_ordering=:degrevlex)
        intermediate_state_variables = variables_2[1:(R - 1) * poseidon.n]
        output_variables = variables_2[(R - 1) * poseidon.n + 1:(R - 1) * poseidon.n + n_half]
        input_variables = variables_2[(R - 1) * poseidon.n + n_half + 1:length(variables_2)]
    else
        P, variables = polynomial_ring(poseidon.field, collect(Iterators.Base.Flatten(variables)), internal_ordering=:degrevlex)
        input_variables = variables[1:poseidon.n]
        intermediate_state_variables = variables[poseidon.n + 1:R * poseidon.n]
        output_variables = variables[R * poseidon.n + 1:length(variables)]
    end
    compr_val = map(P, compr_val)
    polynomials = Vector{typeof(variables[1])}()

    input_state = matrix(input_variables)
    output_state = matrix([compr_val; output_variables])

    if substitution
        output_state -= inv(poseidon.matrix) * input_state
    else
        output_state -= input_state
        input_state = poseidon.matrix * input_state
    end

    counter = 1
    next_state = input_state
    for i in 1:poseidon.r_f
        current_state = next_state
        next_state = matrix(intermediate_state_variables[(counter - 1) * poseidon.n + 1:counter * poseidon.n])
        polys = full_round(poseidon,
                           current_state,
                           matrix(poseidon.constants[:, counter])) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        counter += 1
    end
    for i in 1:poseidon.r_p
        current_state = next_state
        next_state = matrix(intermediate_state_variables[(counter - 1) * poseidon.n + 1:counter * poseidon.n])
        polys = partial_round(poseidon,
                              current_state,
                              matrix(poseidon.constants[:, counter])) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        counter += 1
    end
    for i in 1:poseidon.r_f
        current_state = next_state
        if i < poseidon.r_f
            next_state = matrix(intermediate_state_variables[(counter - 1) * poseidon.n + 1:counter * poseidon.n])
        else
            next_state = output_state
        end
        polys = full_round(poseidon,
                           current_state,
                           matrix(poseidon.constants[:, counter])) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        counter += 1
    end

    return polynomials
end
