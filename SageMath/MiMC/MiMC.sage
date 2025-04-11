class MiMC:
    """
    Implementation of MiMC.
    """
    
    def __init__(self,
                 field,
                 rounds=3,
                 constants=None,
                 info_level=1):
        """
        Initializatoin of the MiMC instance.

        INPUT:
        "field" -- A finite field.
        "r" -- Number of rounds.
        "constants" -- Round constants, expects a list of r field elements.
        "info_level" -- Parameter to print information if > 0.
        """
        self.field = field
        self.q = self.field.order()
        
        self.d = 3
        if gcd(self.d, self.q - 1) != 1:
            raise Exception("Cubing does not induce a permutation over finite field of size " + str(self.q) + ".")
        self.d_inv = self.d.xgcd(self.q - 1)[1]
        while self.d_inv < 0:
            self.d_inv += self.q - 1
        
        self.rounds = rounds
       
        if constants is None:
            constants = [self.field.random_element() for j in range(0, self.rounds)]
        self.constants = [self.field(const) for const in constants]

        if info_level > 0:
            print("MiMC Parameters")
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
        return x**self.d
    
    def inverse_power_perm(self, x):
        """
        Inverse MiMC power permutation.

        INPUT:
        "x" -- A field element.

        OUTPUT:
        Cubic root of field element.
        """
        return x**self.d_inv
    
    def encryption(self, plain, key):
        """
        MiMC encryption function.

        INPUT:
        "plain" -- A field element.
        "key" -- A field elemnt.

        OUTPUT:
        A field element.
        """
        cipher = copy(plain)
        for r in range(0, self.rounds):
            cipher = self.power_perm(cipher + key + self.constants[r])
        cipher += key
        return cipher
    
    def decryption(self, cipher, key):
        """
        MiMC decryption function.

        INPUT:
        "cipher" -- A field element.
        "key" -- A field elemnt.

        OUTPUT:
        A field element.
        """
        plain = copy(cipher)
        plain -= key
        for r in range(self.rounds - 1, -1, -1):
            plain = self.inverse_power_perm(plain) - key - self.constants[r]
        return plain
    
    def generate_variables(self, two_plain_text=False):
        """
        Generates variables for the iterated MiMC polynomial model.

        INPUT:
        "two_plain_text" -- Boolean to generate variables for two
                            plain/ciphertext samples.
                            Default is "False".
        OUTPUT:
        A list of strings.
        """
        variables = []
        if two_plain_text:
            for i in range(0, self.rounds - 1):
                variables.append("u_" + str(i + 1))
            for i in range(0, self.rounds - 1):
                variables.append("v_" + str(i + 1))
        else:
            for i in range(0, self.rounds - 1):
                variables.append("x_" + str(i + 1))
        variables.append("y")
        return variables
    
    def generate_polynomials(self, 
                             plain=None, 
                             cipher=None, 
                             order="degrevlex",
                             info_level=1):
        """
        Generates iterated MiMC polynomial model.

        INPUT:
        "plain" -- A field elemtn.
                   If no value is provided a random one is used.
        "cipher" -- A field elemtn.
                    If no value is provided a random one is used.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        print_key = False
        if plain is None:
            plain = self.field.random_element()
        else:
            plain = self.field(plain)
        if cipher is None:
            print_key = True
            key = self.field.random_element()
            cipher = self.encryption(plain, key)
        else:
            cipher = self.field(cipher)
        
        if info_level:
            print("Plain:", plain)
            if print_key:
                print("Key:", key)
            print("Cipher:", cipher)
            print("Order:", order)

        variables = self.generate_variables(two_plain_text=False)
        P = PolynomialRing(self.field, variables, order=order)
        variables = list(P.gens())
        key = variables[-1]
        polynomials = []

        polynomials.append(self.power_perm(plain + key + self.constants[0]) - variables[0])
        for r in range(1, self.rounds - 1):
            polynomials.append(self.power_perm(variables[r - 1] + key + self.constants[r]) - variables[r])
        polynomials.append(self.power_perm(variables[self.rounds - 2] + key + self.constants[self.rounds - 1]) + key - cipher)

        return polynomials
    
    def generate_two_plain_text_polynomials(self, 
                                            plain_1=None, 
                                            cipher_1=None, 
                                            plain_2=None, 
                                            cipher_2=None, 
                                            order="degrevlex",
                                            info_level=1):
        """
        Generates iterated MiMC polynomial model for two
        plain/ciphertext samples.

        INPUT:
        "plain" -- A field element.
                   If no value is provided a random one is used.
        "cipher" -- A field element.
                    If no value is provided, the plaintext is
                    encrypted under a random key.
        "order" -- Term order for the polynomial ring.
                   Default order is degrevlex.
        "info_level" -- Parameter to print information if > 0.

        OUTPUT:
        A list of polynomials.
        """
        print_key = False
        if plain_1 is None:
            plain_1 = self.field.random_element()
        else:
            plain_1 = self.field(plain_1)
        if plain_2 is None:
            plain_2 = self.field.random_element()
        else:
            plain_2 = self.field(plain_2)
        if (cipher_1 is None) or (cipher_2 is None):
            print_key = True
            key = self.field.random_element()
        if cipher_1 is None:
            cipher_1 = self.encryption(plain_1, key)
        else:
            cipher_1 = self.field(cipher_1)
        if cipher_2 is None:
            cipher_2 = self.encryption(plain_2, key)
        else:
            cipher_2 = self.field(cipher_2)
        
        if info_level > 0:
            if print_key:
                print("Key:", key)
            print("Plain 1:", plain_1)
            print("Cipher 1:", cipher_1)
            print("Plain 2:", plain_2)
            print("Cipher 2:", cipher_2)
            print("Order:", order)

        variables = self.generate_variables(two_plain_text=True)
        P = PolynomialRing(self.field, variables, order=order)
        variables = [P(var) for var in variables]
        key = variables[-1]
        polynomials = []

        # First sample
        polynomials.append(self.power_perm(plain_1 + key + self.constants[0]) - variables[0])
        for r in range(1, self.rounds - 1):
            polynomials.append(self.power_perm(variables[r - 1] + key + self.constants[r]) - variables[r])
        polynomials.append(self.power_perm(variables[self.rounds - 2] + key + self.constants[self.rounds - 1]) + key - cipher_1)

        # Second sample
        polynomials.append(self.power_perm(plain_2 + key + self.constants[0]) - variables[self.rounds - 1])
        for r in range(1, self.rounds - 1):
            polynomials.append(self.power_perm(variables[self.rounds - 1 + r - 1] + key + self.constants[r]) - variables[self.rounds - 1 + r])
        polynomials.append(self.power_perm(variables[self.rounds - 1 + self.rounds - 2] + key + self.constants[self.rounds - 1]) + key - cipher_2)

        return polynomials
