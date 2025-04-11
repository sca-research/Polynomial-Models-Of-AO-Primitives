using Oscar
using Statistics
include("Hydra.jl")
include("Hydra_polynomial_model.jl")

# ----------------------------------------------------------------------
# Hydra parameters
nr_thrds = 48
info_level = 2
N_trials = 10
K = GF(7741)
max_rounds = 15
with_gb = true

sep_1 = repeat("+", 70)
sep_2 = repeat("-", 70)

println("F4 Step Degrees For Hydra")
println("Number of threads: ", nr_thrds)
println("Number of trials: ", N_trials)
println("Field: ", K)
println("Maximal number of rounds: ", max_rounds)
println("With Groebner Basis: ", with_gb)
println(sep_2)

# ----------------------------------------------------------------------
# Function Declarations
function bench_gb_computation(hydra,
                              with_gb;
                              nr_thrds=1, 
                              info_level=0)
    """
    Given a Hydra instance, it generates the
    polynomial model for two samples with a random nonce, and
    benchmarks its DRL Gröbner basis computation with F4.
    Optionally, a change of variables can be performed.

    INPUT:
    "hydra" -- A Hydra instance.
    "with_gb" -- Boolean to perform change of variables
                 to quadratic DRL Gröbner basis.
    "nr_thrds" -- Number of threads for Gröbner basis computation.
                  Default value is 1.
    "info_level" -- Integer, if greater than zero, then F4 progress is
                    printed in console.
                    Default value is 0.

    OUTPUT:
    The running time of F4.
    """
    m = 2
    polys = generate_Hydra_polynomials_m_samples(hydra=hydra, m=m, info_level=info_level);
    polys = transform_Hydra_polynomial_system(hydra, polys, m);

    if with_gb
        _, _, polys_downsized_subs = non_linear_variable_substitution_Hydra_polynomial_system(hydra, 
                                                                                              polys, 
                                                                                              m; 
                                                                                              transformed=true,
                                                                                              info_level=info_level);

        P = parent(polys_downsized_subs[1])
        variables_subs = map(i -> "x_subs_i" * string(i), 1:2 * hydra.rounds_head - 2)
        P_subs, variables_subs = polynomial_ring(hydra.field, variables_subs);
        induce(variables_subs, degrevlex(variables_subs));
        zero_vec = zero_matrix(P_subs, length(gens(P)) - length(variables_subs), 1);
        image = [vec(zero_vec[:, 1]); variables_subs];
        phi = hom(P, P_subs, image);
    
        polys_downsized_subs = map(phi, polys_downsized_subs);
        polys = polys_downsized_subs
    end

    t_1 = time()
    _gb = groebner_basis_f4(ideal(polys), 
                            nr_thrds=nr_thrds, 
                            info_level=info_level);
    t_2 = time()
    t_gb = t_2 - t_1
    return t_gb
end

# ----------------------------------------------------------------------
# Trial run for precompilation
hydra_trial = Hydra_constructor(field=K, rounds_head=2)
_t_trial = bench_gb_computation(hydra_trial, 
                                with_gb; 
                                nr_thrds=nr_thrds, 
                                info_level=0)

# ----------------------------------------------------------------------
# Execution of small scale experiments
for rounds_head in 3:max_rounds
    println("Rounds: ", rounds_head)
    times = []
    for _i in 1:N_trials
        hydra = Hydra_constructor(field=K, rounds_head=rounds_head)
        t_gb = bench_gb_computation(hydra,
                                    with_gb;
                                    nr_thrds=nr_thrds,
                                    info_level=info_level)
        push!(times, t_gb)
        println(sep_1)
    end
    println("Times: ", times)
    println("Mean: ", mean(times))
    println("Standard deviation: ", std(times))
    println(sep_2)
end
