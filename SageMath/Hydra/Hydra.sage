class Hydra:
    """
    Implementation of Hydra.
    """
    def __init__(self,
                 field=GF(2**127 + 45),
                 rounds_body_E_1=2,
                 rounds_body_E_2=4,
                 rounds_body_I=42,
                 rounds_head=39,
                 d=None,
                 matrix_body_E=None,
                 matrix_body_I=None,
                 matrix_head=None,
                 constants_body=None,
                 constants_head=None,
                 initial_value=None,
                 info_level=1):
        """
        Initialization of Hydra instance.

        INPUT:
        "field" -- A prime field.
        "rounds_body_E_1" -- Number of external rounds at beginning of body.
        "rounds_body_E_2" -- Number of external rounds at end of body.
        "rounds_body_I" -- Number of internal rounds of body.
        "rounds_head" -- Number of rounds in the heads.
        "d" -- Exponent of power permutation.
               If no exponent is provided, then smallest expoenent that induces
               a permutation over the field is chosen.
        "matrix_body_E" -- 4x4 matrix for external rounds of body.
        "matrix_body_I" -- 4x4 matrix for internal rounds of body.
        "matrix_head" -- 8x8 matrix for the heads.
        "constants_body" -- Constants for the body, expects a list of 
                            "rounds_body_E_1 + rounds_body_I + rounds_body_E_2"
                            many 4 element lists of integers/field elements.
                            If no constants are provide, random constants are generated.
        "constants_head" -- Constants for the heads, expects a list of "rounds_head" many 
                            8 element lists of integers/field elements.
                            If no constants are provide, random constants are generated.
        "initial_value" -- Initial value of the body, expects a 4 element list/vector
                           of field elements.
                           If no initial value is provided, then [0, 0, 0, 0] is used.
        "info_level" -- Integer, if greater than zero, then Ciminion parameters are
                        printed in console.
        """
        self.field = field
        self.rounds_body_E_1 = rounds_body_E_1
        self.rounds_body_E_2 = rounds_body_E_2
        self.rounds_body_I = rounds_body_I
        self.rounds_head = rounds_head

        q = self.field.order()
        if d is None:
            d = 2
            while gcd(d, q - 1) != 1:
                d += 1
            self.d = d
        
        if matrix_body_E is None:
            self.matrix_body_E = matrix.circulant(vector(self.field, [3, 2, 1, 1]))
        else:
            self.matrix_body_E = matrix_body_E

        if matrix_body_I is None:
            self.matrix_body_I = matrix(self.field, 
                                        [[1, 1, 1, 1],
                                         [1, 4, 1, 1],
                                         [3, 1, 3, 1],
                                         [4, 1, 1, 2]])
        else: 
            self.matrix_body_E = matrix_body_E

        self.matrix_rolling_function = block_diagonal_matrix([self.matrix_body_I, self.matrix_body_I])

        if matrix_head is None:
            self.matrix_head = matrix(self.field,
                                      [[3, 1, 1, 1, 1, 1, 1, 1],
                                       [7, 3, 1, 1, 1, 1, 1, 1],
                                       [4, 1, 4, 1, 1, 1, 1, 1],
                                       [3, 1, 1, 8, 1, 1, 1, 1],
                                       [7, 1, 1, 1, 7, 1, 1, 1],
                                       [8, 1, 1, 1, 1, 5, 1, 1],
                                       [5, 1, 1, 1, 1, 1, 2, 1],
                                       [4, 1, 1, 1, 1, 1, 1, 6]])
        else:
            self.matrix_head = matrix_head

        if constants_body is None:
            self.constants_body = []
            for i in range(0, self.rounds_body_E_1 + self.rounds_body_I + self.rounds_body_E_2):
                con = [self.field.random_element() for j in range(0, 4)]
                self.constants_body.append(con)
        else:
            self.constants_body = constants_body
            
        if constants_head is None:
            self.constants_head = []
            for i in range(0, self.rounds_head):
                con = [self.field.random_element() for j in range(0, 8)]
                self.constants_head.append(con)
        else:
            self.constants_head = constants_head
        
        if initial_value is None:
            self.initial_value = vector(self.field, 4 * [0])
        else:
            self.initial_value = initial_value
        
        if info_level > 0:
            print("Hydra parameters")
            print("Field:", self.field)
            print("Rounds body E_1:", self.rounds_body_E_1)
            print("Rounds body E_2:", self.rounds_body_E_2)
            print("Rounds body I:", self.rounds_body_I)
            print("Rounds head:", self.rounds_head)
            print("d:", self.d)
            print("Matrix body E:\n" + str(self.matrix_body_E))
            print("Matrix body I:\n" + str(self.matrix_body_I))
            print("Matrix rolling function:\n" + str(self.matrix_rolling_function))
            print("Matrix head:\n" + str(self.matrix_head))
            print("Constants body:", self.constants_body)
            print("Constants head:", self.constants_head)
            print("Initial value:", self.initial_value)
        
    def body_round_function_external(self, v_in, d, mat, constants):
        """
        External round function of the body.

        INPUT:
        "v_in" -- 4 element input vector over field.
        "d" -- Exponent of power permutation.
        "mat" -- 4x4 matrix over field.
        "constants" -- 4 element list/vector of integers/field elements.

        OUTPUT:
        4 element vector over field.
        """
        v_out = vector(v_in.parent().base_ring(), 4 * [0])
        for i in range(0, 4):
            v_out[i] += v_in[i]**d
        
        v_out = mat * v_out + vector(self.field, constants)

        return v_out
    
    def body_round_function_internal(self, v_in, mat, constants):
        """
        Internal round function of the body.

        INPUT:
        "v_in" -- 4 element input vector over field.
        "mat" -- 4x4 matrix over field.
        "constants" -- 4 element list/vector of integers/field elements.

        OUTPUT:
        4 element vector over field.
        """
        R = v_in.parent().base_ring()
        v_out = vector(R, 4 * [0])
        s1 = R(0)
        s2 = R(0)

        for i in range(0, 4):
            s1 += (-1)**i * v_in[i]
            s2 += (-1)**floor(i / 2) * v_in[i]
        s = (s1**2 + s2)**2

        for i in range(0, 4):
            v_out[i] += v_in[i] + s
        
        v_out = mat * v_out + vector(self.field, constants)

        return v_out
    
    def body(self, nonce, key):
        """
        Body function of Hydra.

        INPUT:
        "nonce" -- Integer/field element.
        "key" -- 4 element list/vector of integers/field elements.

        OUTPUT:
        8 element vector of field elements.
        """
        key = vector(self.field, key)
        v_in = vector(self.field, 4 * [0])
        v_in += self.initial_value
        v_in[0] += self.field(nonce)

        v_in = self.matrix_body_E * (v_in + key)

        y = vector(self.field, 4 * [0])
        z = vector(self.field, 4 * [0])
        y += v_in

        for i in range(0, self.rounds_body_E_1):
            y = self.body_round_function_external(y, self.d, self.matrix_body_E, self.constants_body[i])
            z += y
        
        for i in range(0, self.rounds_body_I):
            y = self.body_round_function_internal(y, self.matrix_body_I, self.constants_body[self.rounds_body_E_1 + i])
            z += y

        for i in range(0, self.rounds_body_E_2 - 1):
            y = self.body_round_function_external(y, self.d, self.matrix_body_E, self.constants_body[self.rounds_body_E_1 + self.rounds_body_I + i])
            z += y
        y = self.body_round_function_external(y, self.d, self.matrix_body_E, self.constants_body[self.rounds_body_E_1 + self.rounds_body_I + self.rounds_body_E_2 - 1])
        y += key

        y_z = vector(self.field, 8 * [0])
        for i in range(0, 4):
            y_z[i] += y[i]
            y_z[i + 4] += z[i]
        
        return y_z

    def rolling_function(self, v_in, mat):
        """
        Rolling function of Hydra.

        INPUT:
        "v_in" -- 8 element list/vector of integers/field elements.
        "mat" -- 4x4 matrix over field.

        OUTPUT:
        8 element vector of field elements.
        """
        R = v_in.parent().base_ring()
        v_out = vector(R, 8 * [0])
        s1 = R(0)
        s2 = R(0)
        t1 = R(0)
        t2 = R(0)

        for i in range(0, 4):
            s1 += (-1)**i * v_in[i]
            s2 += (-1)**floor(i / 2) * v_in[i + 4]
            t1 += (-1)**i * v_in[i + 4]
            t2 += (-1)**floor(i / 2) * v_in[i]
        v = s1 * s2
        w = t1 * t2

        for i in range(0, 4):
            v_out[i] += v_in[i] + v
            v_out[i + 4] += v_in[i + 4] + w
        
        v_out = mat * v_out

        return v_out
    
    def head_round_function(self, v_in, key, mat, constants):
        """
        Round function of Hydra heads.

        INPUT:
        "v_in" -- 8 element list/vector of integers/field elements.
        "key" -- 8 element list/vector of integers/field elements.
        "mat" -- 8x8 matrix over field.
        "constants" -- 8 element list/vector of integers/field elements.

        OUTPUT:
        8 element vector over field.
        """
        R = v_in.parent().base_ring()
        v_out = vector(R, 8 * [0])
        s = R(0)

        for i in range(0, 8):
            s += (-1)**floor(i / 4) * v_in[i]
        s = s**2

        for i in range(0, 8):
            v_out[i] += v_in[i] + s

        v_out = mat * v_out + key + vector(self.field, constants)

        return v_out
    
    def head(self, y_z, key, samples=1):
        """
        Head function of Hydra.

        INPUT:
        "y_z" -- 8 element list/vector of integers/field elements.
        "key" -- 4 element list/vector of integers/field elements.
        "samples" -- Number of Hydra samples.
                     Default is set to 1.
        
        OUTPUT:
        List of field elements of length 8 * samples.
        """
        key = vector(self.field, key)
        large_key = vector(self.field, list(key) + list(self.matrix_body_E * key))

        v_in = vector(self.field, 8 * [0])
        v_in += y_z
        out = []
        
        for i in range(0, samples):
            v_out = vector(self.field, 8 * [0])
            v_out += v_in
            for j in range(0, self.rounds_head):
                v_out = self.head_round_function(v_out, large_key, self.matrix_head, self.constants_head[j])
            v_out += v_in
            out = out + list(v_out)
            v_in = self.rolling_function(v_in, self.matrix_rolling_function)

        return out
    
    def key_stream(self, nonce, key, samples=1):
        """
        Generates a Hydra key stream.

        INPUT:
        "nonce" -- An integer/field element.
        "key" -- 4 element list/vector of integers/field elements.
        "samples" -- Number of Hydra samples.
                     Default is set to 1.
        
        OUTPUT:
        List of field elements of length 8 * samples.
        """
        y_z = self.body(nonce, key)
        out = self.head(y_z, key, samples=samples)
        return out
