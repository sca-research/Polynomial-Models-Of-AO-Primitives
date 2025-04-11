from sage.all import *
class GMiMC_erf_sponge:
    """
    Implementation of GMiMC-Hash.
    """
    def __init__(self,
                 field=GF(2**3),
                 n=5,
                 rate=1,
                 capacity=None,
                 r=10,
                 constants=None,
                 info_level=1):
        """
        Initialization of the GMiMC_erf instance.

        INPUT:
        "field" -- A finite field.
        "n" -- Number of blocks.
        "rate" -- Rate of the sponge.
        "r" -- Number of rounds.
        "constants" -- Round constants, expects a list of r field elements.
        "info_level" -- Parameter to print information if > 0.
        """
        self.field = field
        self.n = n
        self.rate = rate
        if capacity is None:
            self.capacity = self.n - self.rate
        else:
            self.capacity = capacity
        if self.rate + self.capacity != self.n:
            print("Rate + Capacity != n.")
            return
        self.r = r

        if constants is None:
            constants = [self.field.random_element() for i in range(0, self.r)]
        self.constants = [self.field(const) for const in constants]
        if len(self.constants) < self.r:
            print("Number of constants is less than number of rounds.")
            return

        if info_level > 0:
            print("GMiMC Parameters")
            print("Field:", self.field)
            print("n:", self.n)
            print("Rate:", self.rate)
            print("Capacity:", self.capacity)
            print("r:", self.r)
            print("Constants:", self.constants)


    def round_function(self, x_in, constant):
        """
        GMiMC_erf round function.
        
        INPUT:
        "x_in" -- A vector of n elements.
        "constant" -- A field element.
        
        OUTPUT:
        A vector of n field elements.
        """
        x_out = copy(list(x_in))
        tmp = (x_in[0] + constant)**3
        for i in range(1, self.n):
            x_out[i] += tmp
        x_out = x_out[1:] + x_out[:1]
        return vector(x_out)
    
    def permutation(self, x_in):
        """
        GMiMC_erf permutation.

        INPUT:
        "x_in" -- A vector of n elements.

        OUTPUT:
        A vector of n field elements.
        """
        x_out = copy(x_in)
        for i in range(0, self.r):
            x_out = self.round_function(x_out,
                                        self.constants[i])
        return x_out
    
    def sponge(self, message, IV=None):
        """
        GMiMC_erf sponge function.

        INPUT:
        "message" -- A list of rate-many field elements.
        "IV" -- A list of capacity-many field elements.
                If no IV is provided zero is used.
        
        OUTPUT:
        A list of rate-many field elements.
        """
        if IV is None:
            IV = self.capacity * [self.field(0)]
        x_in = vector(self.field, message + IV)
        x_out = self.permutation(x_in)
        return list(x_out)[:self.rate]
    
    def generate_variables(self):
        """
        Generates variables for the iterated GMiMC_erf polynomial model.

        OUTPUT:
        A list of strings.
        """
        variables = []
        for j in range(0, self.rate):
            variables.append("x_in_" + str(j + 1))
        for i in range(0, self.r - 1):
            for j in range(0, self.n):
                variables.append("x_" + str(j + 1) + "_" + str(i + 1))
        for j in range(0, self.capacity):
            variables.append("x_out_" + str(j + 1))
        return variables

    def generate_preimage_polynomials(self,
                                      IV=None,
                                      hash=None,
                                      order="degrevlex",
                                      info_level=1):
        """
        Generates iterated GMiMC_erf preimage polynomial model.

        INPUT:
        "IV" -- A list of capacity-many field elements.
                If no IV is provided zero is used.
        "hash" -- A list or vector of rate-many field elements.
                  If no value is provided a random message is used to produce a hash.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        print_message = False
        if IV is None:
            IV = self.capacity * [self.field(0)]
        if hash is None:
            print_message = True
            message = [self.field.random_element() for i in range(0, self.rate)]
            hash = self.sponge(message, IV=IV)
        if info_level:
            print("IV:", IV)
            if print_message:
                print("Message:", message)
            print("Hash:", hash)
            print("Order:", order)
        
        variables = self.generate_variables()
        P = PolynomialRing(self.field, variables, order=order)
        variables = list(P.gens())
        polynomials = []

        input_variables = variables[:self.rate]
        output_variables = variables[-self.capacity:]
        variables = variables[self.rate:-self.capacity]

        current_state = vector(P, input_variables + IV)
        next_state = vector(P, variables[:self.n])
        polys = self.round_function(current_state,
                                    self.constants[0]) - next_state
        polynomials = polynomials + list(polys)
        for i in range(1, self.r - 1):
            current_state = next_state
            next_state = vector(P, variables[self.n * i:self.n * (i + 1)])
            polys = self.round_function(current_state,
                                        self.constants[i]) - next_state
            polynomials = polynomials + list(polys)
        current_state = next_state
        next_state = vector(P, hash + output_variables)
        polys = self.round_function(current_state,
                                    self.constants[self.r - 1]) - next_state
        polynomials = polynomials + list(polys)

        return polynomials

    def transform_preimage_polynomials(self,
                                       polynomials):
        """
        Transforms GMiMC iterated preimage polynomial model so that every round contains only one non-linear polynomial.

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

    def compute_preimage_Groebner_basis(self,
                                        polynomials,
                                        transformed=False,
                                        info_level=1):
        """
        Computes a DRL Groebner basis for GMiMC_erf preimage polynomial system via linear change of coordinates.

        INPUT: 
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.
        "transformed" -- Boolean whether polynomials are already transformed. Default value: False.
        "info_level" -- Parameter to print information if > 0.
        
        OUTPUT:
        "lin_polys" -- Linear polynomials in Groebner basis.
        "substitution_polys" -- Polynomials used for substitution.
        "gb" -- DRL Groebner basis.
        """
        if not transformed:
            polynomials = self.transform_preimage_polynomials(polynomials)

        variables_subs = ["x_subs_" + str(i + 1) for i in range(0, self.r)]
        P = polynomials[0].parent()
        P_subs = PolynomialRing(self.field, list(P.gens()) + variables_subs, order=P.term_order())
        input_variables = [P_subs(var) for var in P.gens()[:self.rate]]
        variables = [P_subs(var) for var in P.gens()[self.rate:-self.capacity]]
        variables_subs = [P_subs(var) for var in variables_subs]

        lin_polys = [P_subs(poly) for poly in polynomials if poly.degree() == 1]
        substitution_polys = []

        poly = input_variables[0] + self.constants[0] - variables_subs[0]
        substitution_polys.append(poly)
        for i in range(1, self.r):
            poly = variables[self.n * (i - 1)] + self.constants[i] - variables_subs[i]
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
            lms = sorted([poly.lm() for poly in gb], reverse=True)
            variables_cubed = sorted([var**3 for var in variables_subs], reverse=True)
            is_gb = lms == variables_cubed
            substituion_succcess = set.union(*[set(poly.variables()) for poly in gb]) <= set(variables_subs)
            print("Ideal of leading terms equal to (x_subs_1^3, ..., x_subs_n^3):", is_gb)
            print("All terms of substituted Groebner basis contained in (x_subs_1, ..., x_subs_n):", substituion_succcess)

        return lin_polys, substitution_polys, gb
    
    def generate_CICO_polynomials(self,
                                      IV_in=None,
                                      IV_out=None,
                                      order="degrevlex",
                                      info_level=1):
        """
        Generates iterated GMiMC_erf CICO polynomial model.

        INPUT:
        "IV_in" -- A list of capacity-many field elements.
                   If no IV is provided zero is used.
        "IV_out" -- A list of rate-many field elements.
                    If no IV is provided zero is used.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        if IV_in is None:
            IV_in = self.capacity * [self.field(0)]
        if IV_out is None:
            IV_out = self.rate * [self.field(0)]
        if info_level:
            print("IV_in:", IV_in)
            print("IV_out:", IV_out)
            print("Order:", order)
        
        variables = self.generate_variables()
        P = PolynomialRing(self.field, variables, order=order)
        variables = list(P.gens())
        polynomials = []

        input_variables = variables[:self.rate]
        output_variables = variables[-self.capacity:]
        variables = variables[self.rate:-self.capacity]

        current_state = vector(P, input_variables + IV_in)
        next_state = vector(P, variables[:self.n])
        polys = self.round_function(current_state,
                                    self.constants[0]) - next_state
        polynomials = polynomials + list(polys)
        for i in range(1, self.r - 1):
            current_state = next_state
            next_state = vector(P, variables[self.n * i:self.n * (i + 1)])
            polys = self.round_function(current_state,
                                        self.constants[i]) - next_state
            polynomials = polynomials + list(polys)
        current_state = next_state
        next_state = vector(P, output_variables + IV_out)
        polys = self.round_function(current_state,
                                    self.constants[self.r - 1]) - next_state
        polynomials = polynomials + list(polys)

        return polynomials

    def transform_CICO_polynomials(self,
                                   polynomials):
        """
        Transforms GMiMC iterated CICO polynomial model so that every round contains only one non-linear polynomial.

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

        for i in range(0, self.rate):
            lin_polys = [poly for poly in polynomials_transformed[-self.n * (i + 1):] if poly.degree() == 1] 
            M, v = Sequence(lin_polys).coefficients_monomials()
            M = M.echelon_form()
            lin_polys = list(vector(M * v))
            polynomials_transformed[-self.n * (i + 1)] = polynomials_transformed[-self.n * (i + 1)].reduce(lin_polys)
            polynomials_transformed[-self.n * (i + 1) + 1:] = lin_polys

        return polynomials_transformed
    
    def compute_CICO_Groebner_basis(self,
                                    polynomials,
                                    transformed=False,
                                    info_level=1):
        """
        Computes a DRL Groebner basis for GMiMC_erf preimage polynomial system via linear change of coordinates.

        INPUT: 
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.
        "transformed" -- Boolean whether polynomials are already transformed. Default value: False.
        "info_level" -- Parameter to print information if > 0.
        
        OUTPUT:
        "lin_polys" -- Linear polynomials in Groebner basis.
        "substitution_polys" -- Polynomials used for substitution.
        "gb" -- DRL Groebner basis.
        """
        if not transformed:
            polynomials = self.transform_CICO_polynomials(polynomials)

        variables_subs = ["x_subs_" + str(i + 1) for i in range(0, self.r - self.rate)]
        P = polynomials[0].parent()
        P_subs = PolynomialRing(self.field, list(P.gens()) + variables_subs, order=P.term_order())
        input_variables = [P_subs(var) for var in P.gens()[:self.rate]]
        variables = [P_subs(var) for var in P.gens()[self.rate:-self.capacity]]
        variables_subs = [P_subs(var) for var in variables_subs]

        lin_polys = [P_subs(poly) for poly in polynomials if poly.degree() == 1]
        substitution_polys = []

        poly = input_variables[0] + self.constants[0] - variables_subs[0]
        substitution_polys.append(poly)
        for i in range(1, self.r - self.rate):
            poly = variables[self.n * (i - 1)] + self.constants[i] - variables_subs[i]
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
            lms = sorted([poly.lm() for poly in gb], reverse=True)
            variables_cubed = sorted([var**3 for var in variables_subs], reverse=True)
            is_gb = lms == variables_cubed
            substituion_succcess = set.union(*[set(poly.variables()) for poly in gb]) <= set(variables_subs)
            print("Ideal of leading terms equal to (x_subs_1^3, ..., x_subs_n^3):", is_gb)
            print("All terms of substituted Groebner basis contained in (x_subs_1, ..., x_subs_n):", substituion_succcess)

        return lin_polys, substitution_polys, gb
    
    def generic_coordinates(self, 
                            polynomials,
                            problem="pre"):
        """
        Verfification whether transformed iterated GMiMC preimage polynomial system is in generic coordinates.

        INPUT:
        "polynomials" -- GMiMC iterated polynomial model.
                         Expects list to be in the order of the output of the generation function.
        "problem" -- Type of polynomial system.
                     Expects value in ["pre", "CICO"].

        OUTPUT:
        True or False.
        """
        P = polynomials[0].parent()
        input_variables = P.gens()[:self.rate]
        variables = P.gens()[self.rate:-self.capacity]
        zero_input = len(P.gens()) * [0]
        if problem == "pre":
            polynomials = self.transform_preimage_polynomials(polynomials)
        elif problem == "CICO":
            polynomials = self.transform_CICO_polynomials(polynomials)
        lin_polys_top = [poly - poly(zero_input) for poly in polynomials if poly.degree() == 1]
        lin_polys_top.append(input_variables[0])
        if problem == "pre":
            for i in range(1, self.r):
                lin_polys_top.append(variables[self.n * (i - 1)])
        elif problem == "CICO":
            for i in range(1, self.r - self.rate):
                lin_polys_top.append(variables[self.n * (i - 1)])
        M, _ = Sequence(lin_polys_top).coefficients_monomials()

        generic_coord = M.rank() == self.n * self.r

        return generic_coord
    