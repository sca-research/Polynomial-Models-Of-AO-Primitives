from sage.all import *
load("Feistel.sage")

field = GF(101)
n = 3
r = 4
mode = "erf"

feistel = Feistel(field=field, 
                  n=n, 
                  rounds=r,
                  mode=mode)
print(70 * "-")
V = VectorSpace(feistel.field, feistel.n)
plain = V.random_element()
key = V.random_element()
cipher = feistel.encryption(plain, key)
print("Plain:", plain)
print("Key:", key)
print("Cipher:", cipher)
print("Correct decryption:", plain == feistel.decryption(cipher, key))
print(70 * "-")
polynomials = feistel.generate_polynomials(plain=plain, cipher=cipher)

polys_lin, polys_subs, gb_subs = feistel.compute_Groebner_basis(polynomials)
print(70 * "-")
print("Is Groebner basis:", ideal(polys_lin + polys_subs + gb_subs).basis_is_groebner())
print("Variety:", ideal(polys_lin + polys_subs + gb_subs).variety())
