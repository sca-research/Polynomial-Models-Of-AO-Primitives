from sage.all import *
load("GMiMC_univariate.sage")

def GMiMC_univariate_gc_in_range(field, n, r_min, r_max, mode, info_level):
    result = []
    for r in range(r_min, r_max + 1):
        gmimc = GMiMC(field=field,
                      n=n,
                      r=r,
                      mode=mode,
                      info_level=info_level)
        polynomials = gmimc.generate_polynomials(info_level=info_level)
        gc = gmimc.generic_coordinates(polynomials)
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

    r_min = 5
    r_max = 500

    q = 2**127 + 45
    params = [#[q, n, mode]
               [q, 3, "erf"],
               [q, 3, "crf"],
               [q, 4, "erf"],
               [q, 4, "crf"],
               [q, 5, "erf"],
               [q, 5, "crf"],
             ]

    print(100 * "-")
    for param in params:
        field = GF(param[0])
        n = param[1]
        mode = param[2]
        print("Field:", field)
        print("n:", n)
        print("Mode:", mode)
        result = GMiMC_univariate_gc_in_range(field, n, r_min, r_max, mode, info_level)
        print_result(result)
        print(100 * "-")
