from sage.all import *
load("GMiMC_erf_sponge.sage")

def GMiMC_erf_sponge_gc_in_range(field, n, rate, r_min, r_max, problem, info_level):
    result = []
    for r in range(r_min, r_max + 1):
        gmimc = GMiMC_erf_sponge(field=field,
                                 n=n,
                                 rate=rate,
                                 r=r,
                                 info_level=info_level)
        if problem == "pre":
            polynomials = gmimc.generate_preimage_polynomials(info_level=info_level)
            gc = gmimc.generic_coordinates(polynomials, problem="pre")
        elif problem == "CICO":
            polynomials = gmimc.generate_CICO_polynomials(info_level=info_level)
            gc = gmimc.generic_coordinates(polynomials, problem="CICO")
        result.append([r, gc])

    return result

def print_result(result):
    r_min = result[0][0]
    r_max = result[-1][0]
    success_for_all_rounds = True
    for res in result:
        if res[1] == False:
            success_for_all_rounds = False
            break
    if success_for_all_rounds:
        interval = "[" + str(r_min) + ", " + str(r_max) + "]"
        print("Generic coordinates for all round numbers in the interval:", interval)
    else:
        success_rounds = []
        failure_rounds = []
        for res in result:
            if res[1] == True:
                success_rounds.append(res[0])
            else:
                failure_rounds.append(res[0])
        print("Generic coordinates for round numbers:", success_rounds)
        print("No generic coordinates round numbers:", failure_rounds)

if __name__ == "__main__":
    info_level = 0

    q_1 = 1798650311944395247515796855756291049112378607019380099367795112915886931969
    q_2 = 2**127 + 45
    q_3 = 2**64 - 2**32 + 1
    params = [#[q, n, r_min, r_max]
               [q_1, 4, 130, 160],
               [q_2, 8, 160, 190],
               [q_3, 16, 310, 350],
             ]

    print(100 * "-")
    for param in params:
        field = GF(param[0])
        n = param[1]
        r_min = param[2]
        r_max = param[3]
        for rate in range(1, n):
            print("Field:", field)
            print("n:", n)
            print("Rate:", rate)
            problem = "pre"
            print("Problem:", problem)
            result = GMiMC_erf_sponge_gc_in_range(field, n, rate, r_min, r_max, problem, info_level)
            print_result(result)
            print(100 * "-")
            print("Field:", field)
            print("n:", n)
            print("Rate:", rate)
            problem = "CICO"
            print("Problem:", problem)
            result = GMiMC_erf_sponge_gc_in_range(field, n, rate, r_min, r_max, problem, info_level)
            print_result(result)
            print(100 * "-")
