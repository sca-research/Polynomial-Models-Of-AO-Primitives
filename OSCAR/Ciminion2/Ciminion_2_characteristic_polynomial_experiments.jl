using Oscar
using Statistics
include("Ciminion_2.jl")
include("Ciminion_2_polynomial_model.jl")

# ----------------------------------------------------------------------
# Ciminion2 parameters
K = GF(2^31 - 1)
N_trials = 10
rounds_C = 2
max_rounds_E = 15

sep_1 = repeat("+", 70)
sep_2 = repeat("-", 70)

println("Ciminion2 Characteristic Polynomial")
println("Number of trials: ", N_trials)
println("Field: ", K)
println("Maximal number of rounds: ", rounds_C + max_rounds_E)
println(sep_2)

# ----------------------------------------------------------------------
# Function Declarations
function recursive_basis(variables; d=2)
    """
    Computes the K-vector space basis of
        K [x_1, ..., x_n] / (x_1^d, ..., x_n^d)
    in a recursive fashion.

    INPUT:
    "variables" -- Variables of the polynomial ring.
    "d" -- Degree, default value is 2.

    OUTPUT:
    Vector space basis of the quotient ring.
    """
    n = length(variables)
    if n > 1
        B_rec = recursive_basis(variables[1:n - 1], d=d)
        B = B_rec
        for i in 1:(d - 1)
            B = [B;
                 vec(map(mon -> variables[n]^i * mon, B_rec))]
        end
    else
        return vec(map(i -> variables[1]^i, 0:(d - 1)))
    end
    return B
end

function bench_char_poly(ciminion_2; print_time=true)
    """
    Generates a Ciminion2 polynomial system for a 
    random nonce, computes the DRL Gröbner basis,
    and benchmarks the time for the computation of
    the characteristic polynomial for the last variable.
    
    In addition, the characteristic polynomial matrix is
    initialized with random matrices, and the running time
    for the determinant is benchmarked.

    INPUT:
    ""ciminion_2" -- A Ciminion2 instance.
    "print_time" -- Print the time to console.
                    Default value is true.
    
    OUTPUT:
    OUTPUT:
    "t_mat" -- Time for multiplication matrix generation.
    "t_det" -- Time for characteristic polynomial determinant
    "t_rand" -- Time for the random polynomial determinant.
    """
    polys = generate_Ciminion_2_polynomials(ciminion_2=ciminion_2, info_level=0);

    P = parent(polys[1])
    variables = gens(P)

    variables_Q = String[]
    for var in variables
        push!(variables_Q, string(var))
    end
    variables_Q = [
                   variables_Q[3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 1:3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 2] 
                   variables_Q[1:3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E)];
                   variables_Q[3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 3:3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 3]
                  ]
    Q, variables_Q = polynomial_ring(base_ring(P), variables_Q, internal_ordering=:degrevlex)

    variables_tmp = [
                     variables_Q[2 + 1:2 + 3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E)];
                     variables_Q[1:2];
                     variables_Q[3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 3:3 * (ciminion_2.rounds_C - 1 + ciminion_2.rounds_E) + 3]
                    ]
    h = hom(P, Q, variables_tmp)

    polys = map(h, polys)
    gb = gens(groebner_basis_f4(ideal(polys)))
    polys = gens(groebner_basis_f4(ideal(polys)))

    variables_red = map(i -> variables_Q[2 + 3 * i], 2:(ciminion_2.rounds_C + ciminion_2.rounds_E - 1))
    push!(variables_red, last(variables_Q))
    variables_red = map(string, variables_red)
    P_red, variables_red = polynomial_ring(base_ring(Q), variables_red, internal_ordering=:degrevlex)

    variables_tmp = [zero(P_red), zero(P_red), zero(P_red), zero(P_red), zero(P_red)]
    for i in 1:(ciminion_2.rounds_C + ciminion_2.rounds_E - 2)
        variables_tmp = [variables_tmp; [zero(P_red), zero(P_red)]]
        push!(variables_tmp, variables_red[i])
    end
    push!(variables_tmp, last(variables_red))
    h = hom(Q, P_red, variables_tmp)

    polys = map(h, filter(poly -> total_degree(poly) == 2, polys))

    D_1 = 2^(ciminion_2.rounds_C + ciminion_2.rounds_E - 2)
    D_2 = 2 * D_1

    B = recursive_basis(variables_red);
    B_small = B[D_1 + 1:D_2];
    t_1 = time()
    mat_polys = map(mon -> divrem(last(variables_red) * mon, polys)[2], B_small)
    mat = matrix(map(poly -> map(mon -> coeff(poly, mon), B), mat_polys))
    t_2 = time()
    t_mat = t_2 - t_1
    if print_time
        println("Time needed for matrix generation: ", t_mat, "s")
    end
    
    Q, x = polynomial_ring(ciminion_2.field)

    A_0 = matrix(Q, mat[1:D_1,1:D_1]);
    A_1 = matrix(Q, mat[1:D_1,D_1 + 1:D_2]);
    
    t_1 = time()
    char_poly = x^2 * identity_matrix(Q, D_1) - x * A_1 - A_0;
    char_poly = det(char_poly)
    t_2 = time()
    t_det = t_2 - t_1
    if print_time
        println("Time needed for determinant: ", t_det, "s")
    end

    M = matrix_space(K, D_1, D_1)
    A_0 = matrix(Q, rand(M))
    A_1 = matrix(Q, rand(M))
    t_1 = time()
    char_poly_rand = x^2 * identity_matrix(Q, D_1) - x * A_1 - A_0;
    char_poly_rand = det(char_poly_rand)
    t_2 = time()
    t_rand = t_2 - t_1
    if print_time
        println("Time needed for random determinant: ", t_rand, "s")
    end

    return t_mat, t_det, t_rand
end

# ----------------------------------------------------------------------
# Trial run for precompilation
ciminion_2_trial = Ciminion_2_constructor(field=K, 
                                          rounds_C=rounds_C, 
                                          rounds_E=2,
                                          info_level=0)
_t_mat_trial, _t_det_trial, _t_rand_trial = bench_char_poly(ciminion_2_trial; print_time=false)

# ----------------------------------------------------------------------
# Execution of small scale experiments
for rounds_E in 2:max_rounds_E
    println("Rounds: ", rounds_C + rounds_E)
    times_mat = []
    times_det = []
    times_rand = []
    for _i in 1:N_trials
        ciminion_2 = Ciminion_2_constructor(field=K, 
                                            rounds_C=rounds_C, 
                                            rounds_E=rounds_E,
                                            info_level=0)
        t_mat, t_det, t_rand = bench_char_poly(ciminion_2)
        push!(times_mat, t_mat)
        push!(times_det, t_det)
        push!(times_rand, t_rand)
        println(sep_1)
    end
    println("Times Multiplication Matrix: ", times_mat)
    println("Mean Multiplication Matrix: ", mean(times_mat))
    println("Standard deviation Multiplication Matrix: ", std(times_mat))

    println("Times Determinant: ", times_det)
    println("Mean Determinant: ", mean(times_det))
    println("Standard deviation Determinant: ", std(times_det))

    println("Times Random: ", times_rand)
    println("Mean Random: ", mean(times_rand))
    println("Standard deviation Random: ", std(times_rand))
    println(sep_2)
end
