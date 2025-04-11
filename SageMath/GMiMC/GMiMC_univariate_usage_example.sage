from sage.all import *
load("GMiMC_univariate.sage")

field = GF(101)
n = 3
r = 7
mode = "erf"

gmimc = GMiMC(field=field, 
              n=n, 
              r=r,
              mode=mode)
print(70 * "-")
V = VectorSpace(gmimc.field, gmimc.n)
plain = V.random_element()
key = gmimc.field.random_element()
cipher = gmimc.encrypt(plain, key)
print("Plain:", plain)
print("Key:", key)
print("Cipher:", cipher)
print(70 * "-")
polynomials = gmimc.generate_polynomials(plain=plain, cipher=cipher)

polys_lin, polys_subs, gb_subs = gmimc.compute_Groebner_basis(polynomials)
print(70 * "-")
print("Variety:", ideal(polys_lin + polys_subs + gb_subs).variety())
