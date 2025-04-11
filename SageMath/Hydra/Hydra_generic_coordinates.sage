from sage.all import *
load("Hydra.sage")
load("Hydra_polynomial_model.sage")

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
    field = GF(2**127 + 45)
    matrix_body_E = matrix.circulant(vector(field, [3, 2, 1, 1]))
    matrix_head = matrix(field,
                         [[3, 1, 1, 1, 1, 1, 1, 1],
                          [7, 3, 1, 1, 1, 1, 1, 1],
                          [4, 1, 4, 1, 1, 1, 1, 1],
                          [3, 1, 1, 8, 1, 1, 1, 1],
                          [7, 1, 1, 1, 7, 1, 1, 1],
                          [8, 1, 1, 1, 1, 5, 1, 1],
                          [5, 1, 1, 1, 1, 1, 2, 1],
                          [4, 1, 1, 1, 1, 1, 1, 6],
                         ])

    min_samples = 2
    max_samples = 10

    min_head_rounds = 2
    max_head_rounds = 60

    print("Field:", field)
    print("Matrix external:\n" + str(matrix_body_E))
    print("Matrix head:\n" + str(matrix_head))

    for m in range(min_samples, max_samples + 1):
        print("Number of samples:", m)
        result = []
        for r_H in range(min_head_rounds, max_head_rounds + 1):
            hydra = Hydra(field=field, 
                          matrix_body_E=matrix_body_E,
                          matrix_head=matrix_head,
                          rounds_head=r_H,
                          info_level=info_level)
            gc = is_in_generic_coordinates_Hydra_polynomials_m_samples(hydra=hydra, m=m)
            result.append([r_H, gc])
        print_result(result)
        print(100 * "-")
