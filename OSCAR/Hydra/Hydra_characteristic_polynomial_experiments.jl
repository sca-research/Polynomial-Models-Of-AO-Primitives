using Oscar
using Statistics
include("Hydra.jl")
include("Hydra_polynomial_model.jl")

function recursive_basis(variables; d=2)
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

function bench_char_poly(hydra; print_time=true)
    m = 2
    polys = generate_Hydra_polynomials_m_samples(hydra=hydra, m=m, info_level=0);
    polys = transform_Hydra_polynomial_system(hydra, polys, m);

    affine_polys, polys_subs, polys_downsized_subs = non_linear_variable_substitution_Hydra_polynomial_system(hydra, 
                                                                                                              polys, 
                                                                                                              m; 
                                                                                                              transformed=true,
                                                                                                              info_level=0);

    P = parent(polys_downsized_subs[1])
    variables_subs = map(i -> "x_subs_i" * string(i), 1:2 * hydra.rounds_head - 2)
    P_subs, variables_subs = polynomial_ring(hydra.field, variables_subs, internal_ordering=:degrevlex);
    zero_vec = zero_matrix(P_subs, length(gens(P)) - length(variables_subs), 1);
    image = [vec(zero_vec[:, 1]); variables_subs];
    phi = hom(P, P_subs, image);
    polys_downsized_subs = map(phi, polys_downsized_subs);
    
    gb_subs = [polys_downsized_subs[2 + 1:hydra.rounds_head];
               polys_downsized_subs[2 + hydra.rounds_head + 1:2 * hydra.rounds_head + 2]];
    D_1 = 2^(2 * hydra.rounds_head - 2 - 1)
    D_2 = 2 * D_1

    variables_subs = gens(P_subs);
    B = recursive_basis(variables_subs);
    B_small = B[D_1 + 1:D_2];
    t_1 = time()
    mat_polys = map(mon -> divrem(variables_subs[2 * hydra.rounds_head - 2] * mon, gb_subs)[2], B_small)
    mat = matrix(map(poly -> map(mon -> coeff(poly, mon), B), mat_polys))
    t_2 = time()
    t_mat = t_2 - t_1
    if print_time
        println("Time needed for matrix generation: ", t_mat, "s")
    end
    
    Q, x = polynomial_ring(hydra.field)

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

K = GF(7741)
N_trials = 10
max_rounds = 15

println("Hydra Characteristic Polynomial")
println("Number of trials: ", N_trials)
println("Field: ", K)
println("Maximal number of rounds: ", max_rounds)
println("-----------------------------------------------------------")

# Trial run for precompilation
hydra_trial = Hydra_constructor(field=K, rounds_head=2)
_t_mat_trial, _t_det_trial, _t_rand_trial = bench_char_poly(hydra_trial; print_time=false)

for rounds_head in 3:max_rounds
    println("Rounds: ", rounds_head)
    times_mat = []
    times_det = []
    times_rand = []
    for _i in 1:N_trials
        hydra = Hydra_constructor(field=K, rounds_head=rounds_head)
        t_mat, t_det, t_rand = bench_char_poly(hydra)
        push!(times_mat, t_mat)
        push!(times_det, t_det)
        push!(times_rand, t_rand)
        println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
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
    println("-----------------------------------------------------------")
end
