using Oscar
include("Poseidon.jl")

# ----------------------------------------------------------------------
# Poseidon parameters
nr_thrds = 16

primes = [10007, # d = 3 
          10009, # d = 5
         ]
branches = [2, 4]
r_f_min = 2
r_f_max = 4
r_p_min = 1
r_p_max = 4

# ----------------------------------------------------------------------
# Generate log file
open("Poseidon_Compression_Second_Last_Guessing_Experiment.log", "w") do file
    text = "Prime" * "\t" 
    text *= "n" * "\t" 
    text *= "2 * r_f" * "\t\t" 
    text *= "r_p" * "\t" 
    text *= "Time GB" * "\t\t" 
    text *= "D"
    write(file, text * "\n")
    println(text)
end

# ----------------------------------------------------------------------
# Execution of small scale experiments
for p in primes
    field = GF(p)
    for n in branches
        for r_f in r_f_min:r_f_max
            for r_p in r_p_min:r_p_max
                # Generate random Poseidon instance
                poseidon = Poseidon_constructor(field=field, 
                                                n=n, 
                                                r_f=r_f, 
                                                r_p=r_p,
                                                info_level=0)
                # Generate Poseidon polynomials
                polynomials = generate_compression_polynomials(poseidon; 
                                                               info_level=0)
                P = parent(polynomials[1])
                variables = gens(P)
                # Guess the first n / 2 S-Box input in the second to last roud.
                guess = map(i -> variables[n * (2 * r_f + r_p - 3) + i] - rand(field), 1:Int64(n / 2))
                # Compute DRL Groebner basis.
                t = time()
                gb = groebner_basis_f4(ideal(polynomials) + ideal(guess), 
                                       nr_thrds=nr_thrds, 
                                       info_level=0)
                t = round(time() - t; digits=2)
                # Compute vector space dimension
                A, _ = quo(P, ideal(gb))
                D = factor(Int64(vector_space_dimension(A)))
                open("Poseidon_Compression_Second_Last_Guessing_Experiment.log", "a") do file
                    text = string(p) * "\t" 
                    text *= string(n) * "\t" 
                    text *= string(2 * r_f) * "\t\t" 
                    text *= string(r_p) * "\t" 
                    text *= string(t) * "\t\t" 
                    text *= string(D)
                    write(file, text * "\n")
                    println(text)
                end
            end
        end
    end
end
