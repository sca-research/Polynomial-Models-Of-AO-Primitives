from sage.all import *
load("GMiMC.sage")

def GMiMC_gc_in_range(field, n, r_min, r_max, mat, mode, info_level):
    result = []
    for r in range(r_min, r_max + 1):
        gmimc = GMiMC(field=field,
                      n=n,
                      r=r,
                      key_schedule_matrix=mat,
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
    r_max = 100

    q_1 = 2**64 - 2**32 + 1
    mat_1 = matrix.circulant([2, 1, 1])
    q_2 = 2**31 - 1
    mat_2 = matrix([[5, 7, 1, 3],
                    [4, 6, 1, 1],
                    [1, 3, 5, 7],
                    [1, 1, 4, 6],
                   ])
    params = [#[q, n, mat, mode]
               [q_1, 3, mat_1, "erf"],
               [q_1, 3, mat_1, "crf"],
               [q_2, 4, mat_2, "erf"],
               [q_2, 4, mat_2, "crf"],
             ]
    
    print(100 * "-")
    for param in params:
        field = GF(param[0])
        n = param[1]
        mat = param[2]
        mode = param[3]
        print("Field:", field)
        print("n:", n)
        print("Mode:", mode)
        print("Key schedule matrix:\n" + str(mat))
        result = GMiMC_gc_in_range(field, n, r_min, r_max, mat, mode, info_level)
        print_result(result)
        print(100 * "-")
