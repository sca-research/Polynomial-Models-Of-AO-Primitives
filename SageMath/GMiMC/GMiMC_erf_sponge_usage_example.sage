from sage.all import *
load("GMiMC_erf_sponge.sage")

field = GF(101)
n = 4
rate = 2
r = 12

gmimc = GMiMC_erf_sponge(field=field, 
                         n=n, 
                         rate=rate, 
                         r=r)
print(70 * "-")
print("Preimage polynomials")
polys = gmimc.generate_preimage_polynomials()
polys_lin, polys_subs, gb_subs = gmimc.compute_preimage_Groebner_basis(polys)
print(70 * "-")
print("Is Groebner basis:", ideal(polys_lin + polys_subs + gb_subs).basis_is_groebner())
print(70 * "-")
print("CICO polynomials")
polys = gmimc.generate_CICO_polynomials()
polys_lin, polys_subs, gb_subs = gmimc.compute_CICO_Groebner_basis(polys)
print(70 * "-")
print("Is Groebner basis:", ideal(polys_lin + polys_subs + gb_subs).basis_is_groebner())
