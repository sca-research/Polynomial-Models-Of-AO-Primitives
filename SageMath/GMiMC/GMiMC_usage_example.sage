from sage.all import *
load("GMiMC.sage")

field = GF(101)
n = 3
r = 10
mode = "erf"

gmimc = GMiMC(field=field, 
              n=n, 
              r=r,
              mode=mode)
print(70 * "-")
V = VectorSpace(gmimc.field, gmimc.n)
plain = V.random_element()
key = V.random_element()
cipher = gmimc.encrypt(plain, key)
print("Plain:", plain)
print("Key:", key)
print("Cipher:", cipher)
print(70 * "-")
polynomials = gmimc.generate_polynomials(plain=plain, cipher=cipher)

polys_lin, polys_subs, gb_subs = gmimc.compute_Groebner_basis(polynomials)
print(70 * "-")
print("Is Groebner basis:", ideal(polys_lin + polys_subs + gb_subs).basis_is_groebner())
