from sage.all import *

class Hades:
    """
    Implementation of Hades with multivariate key.
    """
    def __init__(self,
                 field=GF(101),
                 n=2,
                 d=None,
                 r_f=2,
                 r_p=2,
                 constants=None,
                 matrix_f=None,
                 matrix_p=None,
                 key_schedule_matrix=None,
                 info_level=1):
        """
        Initialization of the Hades instance.

        INPUT:
        "field" -- A finite field.
        "n" -- Number of blocks.
        "d" -- Exponent of power permutation.
        "r_f" -- Number of full rounds in the beginning and at the end.
        "r_p" -- Number of partial rounds.
        "constants" -- Round constants, expects a list of 2 * r_f + r_p vectors of length n.
        "matrix_f" -- Expects a n x n invertible matrix over the finite field.
        "matrix_p" -- Expects a n x n invertible matrix over the finite field.
        "key_schedule_matrix" -- Matrix for the key schedule.
                                 Expects a n x n matrix over the finite field.
        "info_level" -- Parameter to print information if > 0.
        """
        self.field = field
        self.q = self.field.order()
        self.n = n
        self.r_f = r_f
        self.r_p = r_p

        if d is None:
            counter = 2
            while gcd(counter, self.q - 1) != 1:
                counter += 1
            self.d = counter
        else:
            if gcd(d, self.q - 1) != 1:
                raise Exception("Given exponent does not induce a power permutation over field of size " + str(self.q) + ".")
            self.d = d
        
        if matrix_f is None:
            M = MatrixSpace(self.field, self.n, self.n)
            matrix_f = M.random_element()
            while matrix_f.det() == self.field(0):
                matrix_f = M.random_element()
        if self.field(matrix_f.det()) == self.field(0):
            raise Exception("Matrix for full rounds is not invertible.")
        self.matrix_f = matrix(self.field, matrix_f)
        self.matrix_f_inv = self.matrix_f.inverse()

        if matrix_p is None:
            matrix_p = copy(self.matrix_f)
        if self.field(matrix_p.det()) == self.field(0):
            raise Exception("Matrix for partial rounds is not invertible.")
        self.matrix_p = matrix(self.field, matrix_p)
        self.matrix_p_inv = self.matrix_p.inverse()

        if key_schedule_matrix is None:
            M = MatrixSpace(self.field, self.n, self.n)
            key_schedule_matrix = M.random_element()
        self.key_schedule_matrix = matrix(self.field, key_schedule_matrix)

        if constants is None:
            V = VectorSpace(self.field, self.n)
            constants = [V.random_element() for _ in range(0, 2 * self.r_f + self.r_p)]
        self.constants = [vector(self.field, const) for const in constants]
        
        if info_level > 0:
            print("Hades Parameters")
            print("Field:", self.field)
            print("n:", self.n)
            print("d:", self.d)
            print("r_f:", self.r_f)
            print("r_p:", self.r_p)
            print("Matrix full rounds:\n" + str(self.matrix_f))
            print("Matrix partial rounds:\n" + str(self.matrix_p))
            print("Key schedule matrix:\n" + str(self.key_schedule_matrix))
            print("Admissible key schedule matrix:", self.admissible_key_schedule_matrix())
            print("Constants:", self.constants)
    
    def admissible_key_schedule_matrix(self):
        """
        Verifies whether all powers of key schedule matrix have all entries non-zero.
        """
        tmp_mat = identity_matrix(self.field, self.n, self.n)
        for i in range(0, 2 * self.r_f + self.r_p - 1):
            tmp_mat = self.key_schedule_matrix * tmp_mat
            if self.field(0) in list(vector(tmp_mat)):
                return False
        return True

    def power_perm(self, x):
        """
        Power permutation.

        INPUT:
        "x" -- A field element.

        OUTPUT:
        d-th power.
        """
        return x**self.d
    
    def full_round(self, plain, key, constants):
        """
        Full Hades round.

        INPUT:
        "plain" -- A vector of field elements.
        "key" -- A vector of field elements.
        "constants" -- A vector of field elements.

        OUTPUT:
        Full keyed SPN permutation.
        """
        cipher = copy(plain)
        for i in range(0, self.n):
            cipher[i] = self.power_perm(cipher[i])
        cipher = self.matrix_f * cipher
        cipher = cipher + key
        cipher = cipher + constants
        return cipher

    def partial_round(self, plain, key, constants):
        """
        Full Hades round.

        INPUT:
        "plain" -- A vector of field elements.
        "key" -- A vector of field elements.
        "constants" -- A vector of field elements.

        OUTPUT:
        Partial keyed SPN permutation.
        """
        cipher = copy(plain)
        cipher[0] = self.power_perm(cipher[0])
        cipher = self.matrix_p * cipher
        cipher = cipher + key
        cipher = cipher + constants
        return cipher

    def encrypt(self, plain, key):
        """
        Hades encryption function.

        INPUT:
        "plain" -- A list or vector of n field elements.
        "key" -- A list or vector of n field elements.

        OUTPUT:
        A vector of n field elements.
        """
        plain = vector(self.field, plain)
        key = vector(self.field, key)
        cipher = plain + key
        for r in range(0, self.r_f):
            key = self.key_schedule_matrix * key
            cipher = self.full_round(cipher, 
                                     key,
                                     self.constants[r])
        for r in range(self.r_f, self.r_f + self.r_p):
            key = self.key_schedule_matrix * key
            cipher = self.partial_round(cipher, 
                                        key, 
                                        self.constants[r])
        for r in range(self.r_f + self.r_p, 2 * self.r_f + self.r_p):
            key = self.key_schedule_matrix * key
            cipher = self.full_round(cipher, 
                                     key,
                                     self.constants[r])
        return cipher

    def generate_variables(self):
        """
        Generates variables for the iterated Hades polynomial model.

        OUTPUT:
        A list of strings.
        """
        variables = []
        for i in range(0, 2 * self.r_f + self.r_p - 1):
            for j in range(0, self.n):
                variables.append("x_" + str(j + 1) + "_" + str(i + 1))
        for j in range(0, self.n):
            variables.append("y_" + str(j + 1))
        return variables
    
    def generate_polynomials(self, 
                             plain=None, 
                             cipher=None, 
                             order="degrevlex",
                             info_level=1):
        """
        Generates iterated Hades polynomial model.

        INPUT:
        "plain" -- A list or vector of n field elements.
                   If no value is provided a random vector is used.
        "cipher" -- A list or vector of n field elements.
                    If no value is provided a random key is used to encrypt plain.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        print_key = False
        V = VectorSpace(self.field, self.n)
        if plain is None:
            plain = V.random_element()
        else:
            plain = vector(plain)
        if cipher is None:
            print_key = True
            key = V.random_element()
            cipher = self.encrypt(plain, key)
        else:
            cipher = vector(cipher)
        
        if info_level > 0:
            print("Plain:", plain)
            if print_key:
                print("Key:", key)
            print("Cipher:", cipher)
            print("Order:", order)

        variables = self.generate_variables()
        P = PolynomialRing(self.field, variables, order=order)
        variables = list(P.gens())
        key = vector(P, variables[self.n * (2 * self.r_f + self.r_p - 1):])
        polynomials = []

        next_state = plain + key
        for r in range(0, self.r_f):
            current_state = next_state
            next_state = vector(P, variables[self.n * r:self.n * (r + 1)])
            key = self.key_schedule_matrix * key
            polys = self.full_round(current_state,
                                    key,
                                    self.constants[r]) - next_state
            polynomials = polynomials + list(polys)
        for r in range(self.r_f, self.r_f + self.r_p):
            current_state = next_state
            next_state = vector(P, variables[self.n * r:self.n * (r + 1)])
            key = self.key_schedule_matrix * key
            polys = self.partial_round(current_state,
                                       key,
                                       self.constants[r]) - next_state
            polynomials = polynomials + list(polys)
        for r in range(self.r_f + self.r_p, 2 * self.r_f + self.r_p ):
            current_state = next_state
            if r == 2 * self.r_f + self.r_p - 1:
                next_state = cipher
            else:
                next_state = vector(P, variables[self.n * r:self.n * (r + 1)])
            key = self.key_schedule_matrix * key
            polys = self.full_round(current_state,
                                    key,
                                    self.constants[r]) - next_state
            polynomials = polynomials + list(polys)

        return polynomials

    def transform_polynomials_erf(self,
                                  polynomials):
        """
        Inverts the matrix in every round of the Hades iterated polynomial model.

        INPUT: 
        "polynomials" -- Hades iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.
        
        OUPUT:
        A list of polynomials.
        """
        polynomials_transformed = []
        for r in range(0, 2 * self.r_f + self.r_p):
            polys = vector(polynomials[self.n * r:self.n * (r + 1)])
            if r >= self.r_f and r < self.r_f + self.r_p:
                polys = self.matrix_p_inv * polys
            else:
                polys = self.matrix_f_inv * polys
            polynomials_transformed = polynomials_transformed + list(polys)
        return polynomials_transformed
