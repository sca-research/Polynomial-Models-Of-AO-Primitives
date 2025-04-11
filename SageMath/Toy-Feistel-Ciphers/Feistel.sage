from sage.all import *

class Feistel:
    """
    Toy implementation of a Feistel cipher with expanding round function (erf) or contracting round function (crf).
    As degree increasing function the monomial function x |-> x^d is used, where d > 1.
    No key schedule is deployed in the cipher, i.e. the master key is added in every round.
    """
    
    def __init__(self,
                 field,
                 n=2,
                 d=2,
                 rounds=2,
                 constants=None,
                 mat=None,
                 mode="erf",
                 info_level=1):
        """
        Initialization of the Feistel instance.

        INPUT:
        "field" -- A finite field.
        "n" -- Number of branches, default value is 2.
        "d" -- Degree of power function in non-linear layer.
               Default value is 2.
        "rounds" -- Number of rounds, default value is 2.
        "constants" -- Affine round constants.
                       Expects a list r of vectors over the finit field. 
                       If not constants are provided, random ones are generated.
        "mat" -- Matrix of the affine layer.
                 If no matrix is provided, a shift permutation is used.
        "mode" -- Expects a string in ["erf", "crf"].
        "info_level" -- Integer, if > 0, then Feistel parameters are printed in console.
                        Default value is 1.
        """
        self.field = field
        self.n = n
        self.rounds = rounds
        self.d = d
        
        if mat is None:
            mat = matrix.zero(self.field, self.n)
            mat[0, self.n - 1] = 1
            for i in range(0, self.n - 1):
                mat[i + 1, i] = 1
        self.mat = matrix(self.field, mat)
        if self.mat.det() == self.field(0):
            raise Exception("Matrix is not invertible")
        self.mat_inverse = self.mat.inverse()
        
        if constants is None:
            V = VectorSpace(self.field, n)
            constants = [V.random_element() for j in range(0, self.rounds)]
        self.constants = [vector(self.field, const) for const in constants]

        if mode not in ["erf", "crf"]:
            print("GMiMC mode " + str(mode) + " not implemented.")
            return
        self.mode = mode
    
        if info_level > 0:
            print("Feistel Parameters")
            print("Field:", self.field)
            print("n:", self.n)
            print("Rounds:", self.rounds)
            print("Mode:", self.mode)
            print("Matrix:\n" + str(self.mat))
            print("Constants:", self.constants)

    def power_map(self, x):
        """
        Raises input to d-th power.

        INPUT:
        "x" -- A field element or polynomial.

        OUTPUT:
        d-th power of x.
        """
        return x**self.d
    
    def encryption_round_function_erf(self, x_in, key, r):
        """
        Expanding round function of the encryption.
        Applies the following Feistel round function:
        (x_1, ..., x_n), k |-> M (x_1 + x_n^d, ..., x_{n - 1} + x_n^d, x_n) + k + c.

        INPUT:
        "x_in" -- A vector of field elements or polynomials.
        "key" -- A vector of field elements or polynomials.
        "r" -- Current round.

        OUTPUT:
        A vector of field elements or polynomials.
        """
        x_out = copy(x_in)
        val = self.power_map(x_out[-1])
        for i in range(0, self.n - 1):
            x_out[i] += val
        x_out = self.mat * x_out
        x_out += self.constants[r]
        x_out += key
        return x_out

    def decryption_round_function_erf(self, x_in, key, r):
        """
        Expanding round function of the decryption.
        Applies the following Feistel round function:
        (x_1, ..., x_n), k |-> M^{-1} (x - k - c) |-> (x_1 - x_n^d, ..., x_{n - 1} - x_n^d, x_n).

        INPUT:
        "x_in" -- A vector of field elements.
        "key" -- A vector of field elements.
        "r" -- Current round.

        OUTPUT:
        A vector of field elements.
        """
        x_out = copy(x_in)
        x_out -= key
        x_out -= self.constants[r]
        x_out = self.mat_inverse * x_out
        val = self.power_map(x_out[-1])
        for i in range(0, self.n - 1):
            x_out[i] -= val
        return x_out
    
    def encryption_round_function_crf(self, x_in, key, r):
        """
        Contracting round function of the encryption.
        Applies the following Feistel round function:
        (x_1, ..., x_n), k |-> M (x_1, ..., x_{n - 1}, x_n + (x_1 + ... + x_{n - 1})^d) + k + c.

        INPUT:
        "x_in" -- A vector of field elements or polynomials.
        "key" -- A vector of field elements or polynomials.
        "r" -- Current round.

        OUTPUT:
        A vector of field elements or polynomials.
        """
        x_out = copy(x_in)
        val = self.field(0)
        for i in range(0, self.n - 1):
            val += x_out[i]
        val = self.power_map(val)
        x_out[self.n - 1] += val
        x_out = self.mat * x_out
        x_out += self.constants[r]
        x_out += key
        return x_out
    
    def decryption_round_function_crf(self, x_in, key, r):
        """
        Contracting Round function of the decryption.
        Applies the following Feistel round function:
        (x_1, ..., x_n), k |-> M^{-1} (x - k - c) |-> (x_1, ..., x_{n - 1}, x_n - (x_1 + ... + x_{n - 1})^d).

        INPUT:
        "x_in" -- A vector of field elements.
        "key" -- A vector of field elements.
        "r" -- Current round.

        OUTPUT:
        A vector of field elements.
        """
        x_out = copy(x_in)
        x_out -= key
        x_out -= self.constants[r]
        x_out = self.mat_inverse * x_out
        val = self.field(0)
        for i in range(0, self.n - 1):
            val += x_out[i]
        val = self.power_map(val)
        x_out[self.n - 1] -= val
        return x_out
    
    def encryption(self, plain, key):
        """
        Encrypts a plaintext.

        INPUT:
        "plain" -- A vector of field elements.
        "key" -- A vector of field elements.

        OUTPUT:
        A vector of field elements.
        """
        if self.mode == "erf":
            round_function = self.encryption_round_function_erf
        elif self.mode == "crf":
            round_function = self.encryption_round_function_crf

        cipher = vector(self.field, plain)
        key = vector(self.field, key)
        cipher += key
        for r in range(0, self.rounds):
            cipher = round_function(cipher, key, r)
        return cipher
    
    def decryption(self, cipher, key):
        """
        Decrypts a ciphertext.

        INPUT:
        "cipher" -- A vector of field elements.
        "key" -- A vector of field elements.

        OUTPUT:
        A vector of field elements.
        """
        if self.mode == "erf":
            round_function = self.decryption_round_function_erf
        elif self.mode == "crf":
            round_function = self.decryption_round_function_crf

        plain = vector(self.field, cipher)
        key = vector(self.field, key)
        for r in range(self.rounds - 1, -1, -1):
            plain = round_function(plain, key, r)
        plain -= key
        return plain
    
    def generate_variables(self):
        """
        Generates variables for iterated polynomial model.
        Intermediate state variables are denoted as x_j__i, where 1 <= j <= n and 1 <= i <= r - 1.
        Key variables are denotes as y_j, where 1 <= j <= n.

        OUTPUT:
        A list of strings
        """
        variables = []
        for i in range(0, self.rounds - 1):
            for j in range(0, self.n):
                variables.append("x_" + str(j + 1) + "__" + str(i + 1))
        for j in range(0, self.n):
            variables.append("y_" + str(j + 1))
        return variables
    
    def generate_polynomials(self, 
                             plain=None, 
                             cipher=None, 
                             termorder="degrevlex",
                             info_level=1):
        """
        Generates iterated polynomial model for Feistel cipher.
        
        INPUT:
        "plain" -- A vector of field elements.
                   If no vector is provided, a random one is generated.
        "ciphter" -- A vector of field elements.
                     If no one is provided a random key is generated, and the plaintext is encrypted under the key.
        "termorder" -- Term order of the polynomial ring.
                       Default order is degrevlex.
        "info_level" -- Integer, if > 0, then plain/ciphertext are printed to console.
                        Default value is 1.
        """
        if self.mode == "erf":
            round_function = self.encryption_round_function_erf
        elif self.mode == "crf":
            round_function = self.encryption_round_function_crf

        V = VectorSpace(self.field, self.n)
        print_key = False
        if plain is None:
            plain = V.random_element()
        plain = vector(self.field, plain)
        if cipher is None:
            print_key = True
            key = V.random_element()
            cipher = self.encryption(plain, key)
        cipher = vector(self.field, cipher)

        if info_level > 0:
            print("Plaintext:", plain)
            if print_key:
                print("Key:", key)
            print("Ciphertext:", cipher)
            print("Term order:", termorder)
        
        variables = self.generate_variables()
        P = PolynomialRing(self.field, variables, order=termorder)
        variables = P.gens()
        key_variables = vector(variables[self.n * (self.rounds - 1):])
        polynomials = []

        current_state = plain + key_variables
        next_state = vector(variables[:self.n])
        polys = round_function(current_state, key_variables, 0) - next_state
        polynomials = polynomials + list(polys)
        for r in range(1, self.rounds - 1):
            current_state = next_state
            next_state = vector(variables[self.n * r:self.n * (r + 1)])
            polys = round_function(current_state, key_variables, r) - next_state
            polynomials = polynomials + list(polys)
        current_state = next_state
        next_state = cipher
        polys = round_function(current_state, key_variables, self.rounds - 1) - next_state
        polynomials = polynomials + list(polys)

        return polynomials

    def transform_polynomials(self, polynomials):
        """
        Transform the iterated polynomial system so that only one non-linear polynomial is present per round.
        Also, eliminates <= r * (n - 1) linear variables in the non-linear polynomials.

        INPUT:
        "polynomials" -- List of polynomials, expects the same ordering as in the output of generate_polynomials.

        OUTPUT:
        A tuple (linear polynomials, non-linear polynomials).
        """
        transformed_polynomials = []
        for i in range(0, self.rounds):
            polys = vector(polynomials[self.n * i:self.n * (i + 1)])
            polys = self.mat_inverse * polys
            if self.mode == "erf":
                for j in range(1, self.n - 1):
                    polys[j] = polys[j] - polys[0]
            transformed_polynomials = transformed_polynomials + list(polys)

        lin_polys = [poly for poly in transformed_polynomials if poly.degree() == 1]
        M, v = Sequence(lin_polys).coefficients_monomials()
        M = M.echelon_form()
        lin_polys = list(vector(M * v))
        non_lin_polys = [poly.reduce(lin_polys) for poly in transformed_polynomials if poly.degree() > 1]

        return lin_polys, non_lin_polys
    
    def generic_coordinates(self, lin_polys):
        """
        Verifies whether transformed iterated Feistel polynomial system is in generic coordinates.

        INPUT:
        "lin_polys" -- Expects the linear polynomials returned from transform_polynomials.

        OUTPUT:
        True or False
        """
        generic_coord = False
        P = lin_polys[0].parent()
        variables = P.gens()
        zero_input = len(variables) * [0]

        lin_polys_top = [poly - poly(zero_input) for poly in lin_polys]
        if self.mode == "erf":
            lin_polys_top = lin_polys_top + [variables[self.n * (i + 1) - 1] for i in range(0, self.rounds)]
        elif self.mode == "crf":
            lin_polys_top = lin_polys_top + [sum(variables[self.n * i:self.n * (i + 1) - 1]) for i in range(0, self.rounds)]
        M, v = Sequence(lin_polys_top).coefficients_monomials()

        generic_coord = M.rank() == self.n * self.rounds

        return generic_coord
    
    def compute_Groebner_basis(self, polynomials, info_level=1):
        """
        Tries to computes a DRL Groebner basis for the Feistel via a linear change of coordinates.
        It is successful if the polynomial system is in generic coordinates.

        INPUT:
        "polynomials" -- List of polynomials, expects the same ordering as in the output of generate_polynomials.
        "info_level" -- Integer, if > 0, then checks success of substitution process and prints it to console.
                        Default value is 1.

        OUTPUT:
        A tuple (linear polynomials, substitution polynomials, non-linear DRL Groebner basis).
        """
        lin_polys, non_lin_polys = self.transform_polynomials(polynomials)
        if info_level > 0:
            print("Is in generic coordinates:", self.generic_coordinates(lin_polys))

        variables_subs = ["x_subs_" + str(i + 1) for i in range(0, self.rounds)]
        P = polynomials[0].parent()
        P_subs = PolynomialRing(self.field, list(P.gens()) + variables_subs, order=P.term_order())
        variables = [P_subs(var) for var in P.gens()]
        variables_subs = [P_subs(var) for var in variables_subs]
        lin_polys = [P_subs(poly) for poly in lin_polys]
        
        if self.mode == "erf":
            substitution_polys = [variables[self.n * (i + 1) - 1] - variables_subs[i] for i in range(0, self.rounds)]
        elif self.mode == "crf":
            substitution_polys = [sum(variables[self.n * i:self.n * (i + 1) - 1]) - variables_subs[i] for i in range(0, self.rounds)]
        substitution_polys = [poly.reduce(lin_polys) for poly in substitution_polys]
        M, v = Sequence(substitution_polys).coefficients_monomials()
        M = M.echelon_form()
        substitution_polys = list(vector(M * v))
        
        gb_subs = [P_subs(poly).reduce(substitution_polys) for poly in non_lin_polys]
        gb_subs = [poly / poly.lc() for poly in gb_subs]

        if info_level > 0:
            lms = sorted([poly.lm() for poly in gb_subs], reverse=True)
            variables_d = sorted([var**self.d for var in variables_subs], reverse=True)
            is_gb = lms == variables_d
            substituion_succcess = set.union(*[set(poly.variables()) for poly in gb_subs]) <= set(variables_subs)
            print("Ideal of leading terms equal to (x_subs_1^d, ..., x_subs_n^d):", is_gb)
            print("All terms of substituted Groebner basis contained in (x_subs_1, ..., x_subs_n):", substituion_succcess)
        
        return lin_polys, substitution_polys, gb_subs
