#  SageMath version 10.4
#  Curve-ID: brainpoolP256t1 from https://www.ietf.org/rfc/rfc5639.txt

# Toy example around brainpoolP256t1.
# Key takeaways :
#    - Arithmetic over GF(p) is faster than over GF(p**2)
#    - log method is much faster than discrete_log method in SageMath
#    - DLP over GF(p) is much faster than over GF(p**2)
#    - So working out equations of twists and solving DLP over GF(p) give a huge edge over the GF(p**2) naive computation


from random import randint
import time

# Parameters of the curve

p = 0xa9fb57dba1eea9bc3e660a909d838d726e3bf623d52620282013481d1f6e5377
a = 0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5374
b = 0x662c61c430d84ea4fe66a7733d0b76b7bf93ebc4af2f49256ae58101fee92b04
G = (0xa3e8eb3cc1cfe7b7732213b23a656149afa142c47aafbc2b79a191562e1305f4,0x2d996c823439c56d7f7b22e14644417e69bcb6de39d027001dabe8f35b25c9be)
n = 0xa9fb57dba1eea9bc3e660a909d838d718c397aa3b561a6f7901e0e82974856a7

# Define base field and curve

print("[+] Creating GF(p).")
t = time.time()
k = GF(p)
print(" |   Time needed to create GF(p) :", round(time.time()-t,3)," seconds.")

# Creating the curve and checking that the order of the generator G given as reference matches

print("[+] Creating the brainpoolP256t1 curve E and checking orders.")
E = EllipticCurve(k,[a,b])
G = E(G)

E1,E2 = E.twists()
if E1 != E:
    E1,E2 == E2,E1
print(" |   Done : ",G.order() == n)


# According to SageMath documentation :
# "The elliptic curve factorization method (ECM) is the fastest way to factor a known composite integer if one of the factors is
#     relatively small (up to approximately 80 bits / 25 decimal digits). "


print("[+] Checking prime factorization of the order of the curves.")
n1 = E1.order()
n2 = E2.order()
l1 = ecm.factor(n1)
l2 = ecm.factor(n2)

# brainpoolP256t1 has prime order but it's twist has 6 out of 7 relatively small prime factors.
# 2**42 makes Baby Step Giant Step possible.

print(" |   Prime facorization of E's order :", l1)
print(" |-  bit length of each prime in the factorization :", [a.bit_length() for a in l1])
print(" |   and of it's twist :", l2)
print(" |-  bit length of each prime in the factorization :",[a.bit_length() for a in l2])


# Work subgroup of order 4233394996199 in the twist

print("[+] Restricting to the subgroup of order 4233394996199 on the twist over GF(p).")
cofactor = 4233394996199
Gt = E2.gens()[0]
h = n2//cofactor
P = Gt*h
print(" |   Done : ",P.order() == cofactor)

# Generate a challenge DLP in the group generated by P and solve it using the log method instead of the (way) slower discrete_log method

print("[+] Generating a challenge for DLP on the twist over GF(p).")
key = randint(1,cofactor)
B = key * P
print(" |   Solving this DLP")
t = time.time()
dl = B.log(P)
print(" |   Time to solve the DLP :", round(time.time() -t,3) ,"seconds.")
print(" |-  Successfully solved the DLP :",  dl== key)

# Playing around with twists and base change

d = -1
print("[+] Let d = -1 in GF(p), d is a square mod p :",(legendre_symbol(d,p)==1))

# There is only one non trivial twist isomorphism class

print("[+] Creating the quadratic twist Ed of E by d over GF(p).")
Ed = EllipticCurve(k,[a*(d**2),b*(d**3)])
print(" |   Ed is isomorphic to the quadradtic twist of E over GF(p) :", Ed.is_isomorphic(E2))

# Convincing evidence that working over GF(p**2) is slow

print("[+] Creating GF(p**2)")
t = time.time()
k2 = GF(p**2)
print(" |   Time needed to create GF(p**2) :", round(time.time()-t,3)," seconds.")
print(" |   Ed is isomorphic to E over GF(p) :",Ed.is_isomorphic(E1))
Edk2 = Ed.base_extend(k2)
E1k2 = E1.base_extend(k2)
print(" |   Ed is isomorphic to E over GF(p**2) : ", Edk2.is_isomorphic(E1k2))

# Generate a challenge DLP in the group generated by P and solve it using the log method

print("[+] Restricting to the subgroup of order 4233394996199 of E over GF(p**2).")
G2t = E1k2.gens()[0]
h = G2t.order()//cofactor
P = G2t*h
print(" |   Done : ",P.order() == cofactor)
print("[+] Generating a challenge for DLP on E over GF(p**2).")
key = randint(1,cofactor)
B = key * P
print(" |   Solving this DLP")
t = time.time()
dl = B.log(P)
print(" |   Time to solve the DLP :", round(time.time() -t,3) ,"seconds.")
print(" |-  Successfully solved the DLP :", dl == key)









