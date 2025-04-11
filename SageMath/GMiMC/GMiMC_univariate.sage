from sage.all import *

class GMiMC:
    """
    Implementation of GMiMC-erf and GMiMC-crf with univariate key.
    """
    def __init__(self,
                 field=GF(2**3),
                 n=5,
                 r=10,
                 constants=None,
                 mode="erf",
                 info_level=1):
        """
        Initialization of the GMiMC_erf instance with univariate key.

        INPUT:
        "field" -- A finite field.
        "n" -- Number of blocks.
        "r" -- Number of rounds.
        "constants" -- Round constants, expects a list of r field elements.
        "mode" -- Expects a string in ["erf", "crf"].
        "info_level" -- Parameter to print information if > 0.
        """
        self.field = field
        self.n = n
        self.r = r

        if constants is None:
            constants = [self.field.random_element() for i in range(0, self.r)]
        self.constants = [self.field(const) for const in constants]
        if len(self.constants) < self.r:
            print("Number of constants is less than number of rounds.")
            return
        
        if mode not in ["erf", "crf"]:
            print("GMiMC mode " + str(mode) + " not implemented.")
            return
        self.mode = mode

        if info_level > 0:
            print("GMiMC Parameters")
            print("Field:", self.field)
            print("n:", self.n)
            print("r:", self.r)
            print("Mode:", self.mode)
            print("Constants:", self.constants)

    def round_function_erf(self, x_in, key, constant):
        """
        GMiMC_erf round function.
        
        INPUT:
        "x_in" -- A vector of n elements.
        "key" -- A field element or variable.
        "constant" -- A field element.
        
        OUTPUT:
        A vector of n field elements.
        """
        x_out = copy(list(x_in))
        tmp = (x_in[0] + key + constant)**3
        for i in range(1, self.n):
            x_out[i] += tmp
        x_out = x_out[1:] + x_out[:1]
        return vector(x_out)
    
    def round_function_crf(self, x_in, key, constant):
        """
        GMiMC_crf round function.
        
        INPUT:
        "x_in" -- A vector of n elements.
        "key" -- A field element or variable.
        "constant" -- A field element.
        
        OUTPUT:
        A vector of n field elements.
        """
        x_out = copy(list(x_in))
        tmp = key + constant
        for j in range(1, self.n):
            tmp += x_out[j]
        tmp = tmp**3
        x_out[0] += tmp
        x_out = x_out[1:] + x_out[:1]
        return vector(x_out)
    
    def encrypt(self, plain, key):
        """
        GMiMC_erf encryption function.

        INPUT:
        "plain" -- A list or vector of n field elements.
        "key" -- A field element..

        OUTPUT:
        A list of n field elements.
        """
        if self.mode == "erf":
            round_function = self.round_function_erf
        elif self.mode == "crf":
            round_function = self.round_function_crf

        plain = vector(self.field, copy(plain))
        cipher = copy(plain)
        for i in range(0, self.r):
            cipher = round_function(cipher,
                                    key,
                                    self.constants[i])
        return cipher
    
    def generate_variables(self):
        """
        Generates variables for the iterated GMiMC_erf polynomial model.

        OUTPUT:
        A list of strings.
        """
        variables = []
        for i in range(0, self.r - 1):
            for j in range(0, self.n):
                variables.append("x_" + str(j + 1) + "_" + str(i + 1))
        variables.append("y")
        return variables

    def generate_polynomials(self,
                             plain=None,
                             cipher=None,
                             order="degrevlex",
                             info_level=1):
        """
        Generates iterated GMiMC_erf polynomial model.

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
        if self.mode == "erf":
            round_function = self.round_function_erf
        elif self.mode == "crf":
            round_function = self.round_function_crf

        print_key = False
        V = VectorSpace(self.field, self.n)
        if plain is None:
            plain = V.random_element()
        else:
            plain = vector(plain)
        if cipher is None:
            print_key = True
            key = self.field.random_element()
            cipher = self.encrypt(plain, key)
        else:
            cipher = vector(cipher)
        
        if info_level:
            print("Plain:", plain)
            if print_key:
                print("Key:", key)
            print("Cipher:", cipher)
            print("Order:", order)
        
        variables = self.generate_variables()
        P = PolynomialRing(self.field, variables, order=order)
        variables = list(P.gens())
        polynomials = []

        key = variables[-1]

        current_state = copy(plain)
        next_state = vector(P, variables[:self.n])
        polys = round_function(current_state,
                               key,
                               self.constants[0]) - next_state
        polynomials = polynomials + list(polys)
        for i in range(1, self.r - 1):
            current_state = next_state
            next_state = vector(P, variables[self.n * i:self.n * (i + 1)])
            polys = round_function(current_state,
                                   key,
                                   self.constants[i]) - next_state
            polynomials = polynomials + list(polys)
        current_state = next_state
        next_state = copy(cipher)
        polys = round_function(current_state,
                               key,
                               self.constants[self.r - 1]) - next_state
        polynomials = polynomials + list(polys)

        return polynomials

    def transform_polynomials_erf(self,
                                  polynomials):
        """
        Transforms GMiMC-erf iterated polynomial model so that every round contains only one non-linear polynomial.

        INPUT: 
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.
        
        OUPUT:
        A list of polynomials.
        """
        polynomials_transformed = copy(polynomials)
        for i in range(0, self.r):
            for j in range(1, self.n - 1):
                polynomials_transformed[i * self.n + j] -= polynomials_transformed[i * self.n]
        return polynomials_transformed
    
    def generic_coordinates(self,
                            polynomials):
        """
        Verfification whether transformed iterated GMiMC polynomial system is in generic coordinates.

        INPUT:
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.

        OUTPUT:
        True or False.
        """
        P = polynomials[0].parent()
        variables = P.gens()
        key = variables[-1]
        zero_input = len(variables) * [0]
        if self.mode == "erf":
            polynomials = self.transform_polynomials_erf(polynomials)
        lin_polys_top = [poly - poly(zero_input) for poly in polynomials if poly.degree() == 1]
        lin_polys_top.append(key)
        if self.mode == "erf":
            for i in range(1, self.r):
                lin_polys_top.append(variables[self.n * (i - 1)] + key)
        elif self.mode == "crf":
            for i in range(1, self.r):
                lin_polys_top.append(sum(variables[self.n * (i - 1) + 1:self.n * i]) + key)
        M, _ = Sequence(lin_polys_top).coefficients_monomials()

        generic_coord = M.rank() == self.n * (self.r - 1) + 1

        return generic_coord
    
    def compute_Groebner_basis(self,
                               polynomials,
                               info_level=1):
        """
        Computes a DRL Groebner basis for GMiMC_erf via an affine change of coordinates.

        INPUT: 
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generate_polynomials function.
        "info_level" -- Parameter to print information if > 0.
        
        OUTPUT:
        "lin_polys" -- Linear polynomials in Groebner basis.
        "substitution_polys" -- Polynomials used for substitution.
        "gb" -- DRL Groebner basis.
        """
        if self.mode == "erf":
            polynomials = self.transform_polynomials_erf(polynomials)

        variables_subs = ["x_subs_" + str(i + 1) for i in range(0, self.r)]
        P = polynomials[0].parent()
        P_subs = PolynomialRing(self.field, list(P.gens()) + variables_subs, order=P.term_order())
        variables = [P_subs(var) for var in P.gens()]
        key = variables[-1]
        variables_subs = [P_subs(var) for var in variables_subs]

        lin_polys = [P_subs(poly) for poly in polynomials if poly.degree() == 1]
        substitution_polys = []

        poly = key + self.constants[0] - variables_subs[0]
        substitution_polys.append(poly)
        if self.mode == "erf":
            for i in range(1, self.r):
                poly = variables[self.n * (i - 1)] + key  + self.constants[i] - variables_subs[i]
                substitution_polys.append(poly)
        elif self.mode == "crf":
            for i in range(1, self.r):
                poly = sum(variables[self.n * (i - 1) + 1:self.n * i]) + key + self.constants[i] - variables_subs[i]
                substitution_polys.append(poly)

        M, v = Sequence(lin_polys).coefficients_monomials()
        M = M.echelon_form()
        lin_polys = list(vector(M * v))

        substitution_polys = [poly.reduce(lin_polys) for poly in substitution_polys]
        M, v = Sequence(substitution_polys).coefficients_monomials()
        M = M.echelon_form()
        substitution_polys = list(vector(M * v))

        gb = [P_subs(poly).reduce(lin_polys).reduce(substitution_polys) for poly in polynomials if poly.degree() > 1]
        gb = [poly / poly.lc() for poly in gb]

        if info_level > 0:
            lms = set([poly.lm() for poly in gb])
            variables_cubed = set([var**3 for var in variables_subs[self.n - 1:]])
            is_gb = variables_cubed.issubset(lms)
            substituion_succcess = set.union(*[set(poly.variables()) for poly in gb]).issubset(set(variables_subs[self.n - 1:]))
            print("(x_subs_" + str(self.n) + "^3, ..., x_subs_n^3) contained in ideal of leading terms:", is_gb)
            print("All terms of substituted Groebner basis contained in (x_subs_" + str(self.n) + ", ..., x_subs_n):", substituion_succcess)

        return lin_polys, substitution_polys, gb
        