class Feistel_MiMC:
    """
    Implementation of MiMC.
    """
    
    def __init__(self,
                 field,
                 rounds=3,
                 constants=None,
                 info_level=1):
        """
        Initializatoin of the Feistel-MiMC instance.

        INPUT:
        "field" -- A finite field.
        "r" -- Number of rounds.
        "constants" -- Round constants, expects a list of r field elements.
        "info_level" -- Parameter to print information if > 0.
        """
        
        self.field = field
        self.q = self.field.order()
        
        self.branches = 2
        self.rounds = rounds
        
        if constants is None:
            constants = [self.field.random_element() for j in range(0, self.rounds)]
        self.constants = [self.field(const) for const in constants]
        
        if info_level > 0:
            print("Feistel-MiMC Parameters")
            print("Field:", self.field)
            print("r:", self.rounds)
            print("Constants:", self.constants)
    
    def power_perm(self, x):
        """
        MiMC power permutation.

        INPUT:
        "x" -- A field element.

        OUTPUT:
        Cubed field element.
        """
        return x**3
        
    def encryption_round_function(self, plain, key, round):
        """
        Feistel-MiMC encryption round function.

        INPUT:
        "plain" -- A two element vector of field elements.
        "key" -- A field elemnt.
        "round" -- Current round.

        OUTPUT:
        A two element vector of field elements.
        """
        L = plain[0]
        R = plain[1]
        return vector([R + self.power_perm(L + key + self.constants[round]), 
                       L])
    
    def decryption_round_function(self, cipher, key, round):
        """
        Feistel-MiMC decryption round function.

        INPUT:
        "cipher" -- A two element vector of field elements.
        "key" -- A field elemnt.
        "round" -- Current round.

        OUTPUT:
        A two element vector of field elements.
        """
        L = cipher[1]
        R = cipher[0]
        return vector([L, 
                       R - self.power_perm(L + key + self.constants[round])])
    
    def encryption(self, plain, key):
        """
        Feistel-MiMC encryption function

        INPUT:
        "plain" -- A two element vector of field elements.
        "key" -- A field element.

        OUTPUT:
        A two element vector of field elements.
        """
        cipher = vector(plain)
        for i in range(0, self.rounds):
            cipher = self.encryption_round_function(cipher, key, i)
        return vector(cipher) + vector([key, self.field(0)])
    
    def decryption(self, cipher, key):
        """
        Feistel-MiMC decryption function

        INPUT:
        "plain" -- A two element vector of field elements.
        "key" -- A field element.

        OUTPUT:
        A two element vector of field elements.
        """
        plain = vector(cipher) - vector([key, self.field(0)])
        for i in range(self.rounds - 1, -1, -1):
            plain = self.decryption_round_function(plain, key, i)
        return vector(plain)
    
    def generate_variables(self):
        """
        Generates variables for the iterated Feistel-MiMC polynomial model.

        OUTPUT:
        A list of strings.
        """
        variables = []
        for i in range(0, self.rounds - 1):
                variables.append(["x_L_" + str(i + 1), "x_R_" + str(i + 1)])
        variables.append(["y"])
        return variables
    
    def generate_polynomials(self, 
                             plain=None, 
                             cipher=None, 
                             order="degrevlex",
                             info_level=1):
        """
        Generates iterated Feistel-MiMC polynomial model.

        INPUT:
        "plain" -- A two element vector of field elements.
                   If no value is provided a random one is used.
        "cipher" -- A two element vector of field elements.
                    If no value is provided, the plaintext is
                    encrypted under a random key.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        print_key = False
        if plain is None:
            plain = [self.field.random_element(), self.field.random_element()]
        plain = vector(plain)
        if cipher is None:
            print_key = True
            key = self.field.random_element()
            cipher = self.encryption(plain, key)
        cipher = vector(cipher)
        
        if info_level:
            print("Plain:", plain)
            if print_key:
                print("Key:", key)
            print("Cipher:", cipher)
            print("Order:", order)

        variables = self.generate_variables()
        P = PolynomialRing(self.field, flatten(variables), order=order)
        variables = [vector([P(el) for el in var]) for var in variables]
        key_variable = variables[-1][0]

        polynomials = []
        polynomials.append(self.encryption_round_function(plain, 
                                                          key_variable, 
                                                          0) - variables[0])
        for r in range(1, self.rounds - 1):
            polynomials.append(self.encryption_round_function(variables[r - 1], 
                                                              key_variable, 
                                                              r) - variables[r])
        polynomials.append(self.encryption_round_function(variables[self.rounds - 2], 
                                                          key_variable, 
                                                          self.rounds - 1) + vector([key_variable, 
                                                                                     self.field(0)]) - cipher)
        
        return flatten([list(polys) for polys in polynomials])
        