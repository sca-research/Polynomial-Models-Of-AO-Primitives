using Oscar
include("utilities.jl")

struct Hydra
    field::FqField
    d::Int64
    rounds_body_E_1::Int64
    rounds_body_E_2::Int64
    rounds_body_I::Int64
    rounds_head::Int64
    matrix_body_E::FqMatrix
    matrix_body_I::FqMatrix
    matrix_rolling_function::FqMatrix
    matrix_head::FqMatrix
    constants_body::FqMatrix
    constants_head::FqMatrix
    inital_value::FqMatrix
end

function Hydra_constructor(;field=GF(2),
                            d=nothing,
                            rounds_body_E_1=2,
                            rounds_body_E_2=4,
                            rounds_body_I=42,
                            rounds_head=39,
                            constants_body=nothing,
                            constants_head=nothing,
                            initial_value=nothing,
                            info_level=1)
    """
    Initialization of Hydra instance.

    INPUT:
    "field" -- A prime field.
    "rounds_body_E_1" -- Number of external rounds at beginning of body.
    "rounds_body_E_2" -- Number of external rounds at end of body.
    "rounds_body_I" -- Number of internal rounds of body.
    "rounds_head" -- Number of rounds in the heads.
    "d" -- Exponent of power permutation.
           If no exponent is provided, then smallest expoenent that induces
           a permutation over the field is chosen.
    "constants_body" -- Constants for the body, expects a matrix over the
                        field with dimension
                        (rounds_body_E_1 + rounds_body_I + rounds_body_E_2) x 4.
                         If no constants are provide, random constants are generated.
    "constants_head" -- Constants for the heads, expects a matrix over the
                        field of dimension rounds_head x 4.
                        If no constants are provide, random constants are generated.
    "initial_value" -- Initial value of the body, expects a 4 element list/vector
                           of field elements.
                           If no initial value is provided, then [0, 0, 0, 0] is used.
    "info_level" -- Integer, if greater than zero, then Ciminion parameters are
                    printed in console.
    """
    q = order(field)
    if isnothing(d)
        d = 2
        while gcd(d, q - 1) != 1
            d += 1 
        end
    else
        if gcd(d, q - 1) != 1
            println(d, " does not induce permutation over ", field)
            return
        end
    end

    if isnothing(constants_body)
        R = rounds_body_E_1 + rounds_body_I + rounds_body_E_2
        constants_body = zero_matrix(field, 4, R)
        for i in 1:R
            for j in 1:4
                constants_body[j, i] = rand(field)
            end
        end
    end

    if isnothing(constants_head)
        constants_head = zero_matrix(field, 8, rounds_head)
        for i in 1:rounds_head
            for j in 1:8
                constants_head[j, i] = rand(field)
            end
        end
    end

    if isnothing(initial_value)
        initial_value = zero_matrix(field, 4, 1)
    end

    tmp = [field(3), field(2), field(1), field(1)]
    matrix_body_E = circulant_matrix(tmp)

    matrix_body_I = [1 1 1 1; 
                     1 4 1 1; 
                     3 1 3 1; 
                     4 1 1 2]
    matrix_body_I = matrix(field, matrix_body_I)

    tmp = zero_matrix(field, 4, 4)
    matrix_rolling_function = [matrix_body_I tmp; tmp matrix_body_I]

    matrix_head = [3 1 1 1 1 1 1 1;
                   7 3 1 1 1 1 1 1;
                   4 1 4 1 1 1 1 1;
                   3 1 1 8 1 1 1 1;
                   7 1 1 1 7 1 1 1;
                   8 1 1 1 1 5 1 1;
                   5 1 1 1 1 1 2 1;
                   4 1 1 1 1 1 1 6]
    matrix_head = matrix(field, matrix_head)

    if info_level > 0
        println("Hydra parameters")
        println("Field: ", field)
        println("Rounds body E_1: ", rounds_body_E_1)
        println("Rounds body E_2: ", rounds_body_E_2)
        println("Rounds body I: ", rounds_body_I)
        println("Rounds head: ", rounds_head)
        println("d: ", d)
        println("Matrix body E:\n", matrix_body_E)
        println("Matrix body I:\n", matrix_body_I)
        println("Matrix head:\n", matrix_head)
        println("Constants body:\n", constants_body)
        println("Constants head:\n", constants_head)
        println("Initial value: ", initial_value)
    end

    return Hydra(field, 
                 d, 
                 rounds_body_E_1, 
                 rounds_body_E_2, 
                 rounds_body_I,
                 rounds_head, 
                 matrix_body_E, 
                 matrix_body_I,
                 matrix_rolling_function,
                 matrix_head,
                 constants_body,
                 constants_head,
                 initial_value)
end

function body_round_function_external(v_in::FqMatrix, 
                                      d::Int64, 
                                      mat::FqMatrix, 
                                      constants::FqMatrix)
    """
    External round function of the body.

    INPUT:
    "v_in" -- 4x1 matrix over field.
    "d" -- Exponent of power permutation.
    "mat" -- 4x4 matrix over field.
    "constants" -- 4x1 matrix over field.

    OUTPUT:
    4x1 matrix over field.
    """
    v_out = deepcopy(v_in)
    for i in 1:4
        v_out[i, 1] = v_out[i, 1]^d
    end

    v_out = mat * v_out + constants

    return v_out
end

function body_round_function_internal(v_in::FqMatrix, 
                                      mat::FqMatrix,
                                      constants::FqMatrix)
    """
    Internal round function of the body.

    INPUT:
    "v_in" -- 4x1 matrix over field.
    "mat" -- 4x4 matrix over field.
    "constants" -- 4x1 matrix over field.

    OUTPUT:
    4x1 matrix over field.
    """
    v_out = deepcopy(v_in)
    K = base_ring(v_in)
    s1 = zero(K)
    s2 = zero(K)

    for i in 1:4
        s1 += (-1)^(i - 1) * v_in[i, 1]
        s2 += (-1)^(Int64(floor((i - 1) / 2))) * v_in[i, 1]
    end
    s = (s1^2 + s2)^2
    
    for i in 1:4
        v_out[i, 1] += s
    end

    v_out = mat * v_out + constants

    return v_out
end

function body(nonce::FqFieldElem, 
              key::FqMatrix, 
              hydra::Hydra)
    """
    Body function of Hydra.

    INPUT:
    "nonce" -- Integer/field element.
    "key" -- 4x1 matrix over field.

    OUTPUT:
    8x1 matrix over field.
    """
    v_in = zero_matrix(parent(nonce), 4, 1)
    v_in[1, 1] += nonce
    for i in 2:4
        v_in[i, 1] += hydra.inital_value[i, 1]
    end
    v_in = hydra.matrix_body_E * (v_in + key)

    y = deepcopy(v_in)
    z = zero_matrix(hydra.field, 4, 1)
    
    for i in 1:hydra.rounds_body_E_1
        y = body_round_function_external(y, 
                                         hydra.d, 
                                         hydra.matrix_body_E, 
                                         matrix(hydra.constants_body[:, i]))
        z += y
    end

    for i in 1:hydra.rounds_body_I
        y = body_round_function_internal(y, 
                                         hydra.matrix_body_I, 
                                         matrix(hydra.constants_body[:, hydra.rounds_body_E_1 + i]))
        z += y
    end

    for i in 1:hydra.rounds_body_E_2 - 1
        y = body_round_function_external(y, 
                                         hydra.d, 
                                         hydra.matrix_body_E, 
                                         matrix(hydra.constants_body[:, hydra.rounds_body_E_1 + hydra.rounds_body_I + i]))
        z += y
    end
    y = body_round_function_external(y, 
                                     hydra.d, 
                                     hydra.matrix_body_E, 
                                     matrix(hydra.constants_body[:, hydra.rounds_body_E_1 + hydra.rounds_body_I + hydra.rounds_body_E_2]))
    y += key

    y_z = zero_matrix(hydra.field, 8, 1)
    for i in 1:4
        y_z[i, 1] = y[i, 1]
        y_z[i + 4, 1] = z[i, 1]
    end

    return y_z
end

function rolling_function(v_in::Union{FqMatrix, AbstractAlgebra.Generic.MatSpaceElem{FqMPolyRingElem}},
                          mat::FqMatrix)
    """
    Rolling function of Hydra.

    INPUT:
    "v_in" -- 8x1 matrix over field.
    "mat" -- 4x4 matrix over field.

    OUTPUT:
    8x1 matrix over field.
    """
    K = base_ring(v_in)
    v_out = zero_matrix(K, 8, 1)
    for i in 1:8
        v_out[i, 1] += v_in[i, 1]
    end
    s1 = zero(K)
    s2 = zero(K)
    t1 = zero(K)
    t2 = zero(K)
    for i in 1:4
        s1 += (-1)^(i - 1) * v_in[i, 1]
        s2 += (-1)^(Int64(floor((i - 1) / 2))) * v_in[i + 4, 1]
        t1 += (-1)^(i - 1) * v_in[i + 4, 1]
        t2 += (-1)^(Int64(floor((i - 1) / 2))) * v_in[i, 1]
    end
    v = s1 * s2
    w = t1 * t2

    for i in 1:4
        v_out[i, 1] += v
        v_out[i + 4, 1] += w
    end
    
    v_out = mat * v_out

    return v_out
end

function head_round_function(v_in::Union{FqMatrix, AbstractAlgebra.Generic.MatSpaceElem{FqMPolyRingElem}},
                             key::Union{FqMatrix, AbstractAlgebra.Generic.MatSpaceElem{FqMPolyRingElem}},
                             mat::FqMatrix,
                             constants::FqMatrix)
    """
    Round function of Hydra heads.

    INPUT:
    "v_in" -- 8x1 matrix over field.
    "key" -- 8x1 matrix over field.
    "mat" -- 8x8 matrix over field.
    "constants" -- 8x1 matrix over field.

    OUTPUT:
    8x1 matrix over field.
    """
    K = base_ring(v_in)
    v_out = zero_matrix(K, 8, 1)
    for i in 1:8
        v_out[i, 1] += v_in[i, 1]
    end
    s = zero(K)

    for i in 1:8
        s += (-1)^(Int64(floor((i - 1) / 4))) * v_in[i, 1]
    end
    s = s^2

    for i in 1:8
        v_out[i, 1] += s
    end

    v_out = mat * v_out + key + constants

    return v_out
end

function head(y_z::FqMatrix,
              key::FqMatrix,
              hydra::Hydra;
              m::Int64=1)
    """
    Head function of Hydra.

    INPUT:
    "y_z" -- 8x1 matrix over field.
    "key" -- 4x1 matrix over field.
    "m" -- Number of Hydra samples.
           Default is set to 1.
        
    OUTPUT:
    Vector of field elements of length 8 * samples.
    """
    large_key = zero_matrix(hydra.field, 8, 1)
    tmp = hydra.matrix_body_E * key
    for i in 1:4
        large_key[i, 1] += key[i, 1]
        large_key[i + 4, 1] += tmp[i, 1]
    end

    v_in = deepcopy(y_z)
    out = Vector{typeof(y_z[1, 1])}()

    for i in 1:m
        v_out = deepcopy(v_in)
        for j in 1:hydra.rounds_head
            v_out = head_round_function(v_out, 
                                        large_key, 
                                        hydra.matrix_head, 
                                        matrix(hydra.constants_head[:, j]))
        end
        v_out += v_in
        out = [out; vec(v_out[:,1])]
        v_in = rolling_function(v_in, hydra.matrix_rolling_function) 
    end

    return out
end

function key_stream(nonce::FqFieldElem, 
                    key::FqMatrix, 
                    hydra::Hydra;
                    m::Int64=1)
    """
    Generates a Hydra key stream.

    INPUT:
    "nonce" -- An integer/field element.
    "key" -- 4x1 matrix over field.
    "m" -- Number of Hydra samples.
           Default is set to 1.
        
    OUTPUT:
    Vector of field elements of length 8 * samples.
    """
    y_z = body(nonce, key, hydra)
    out = head(y_z, key, hydra, m=m)
    return out
end
