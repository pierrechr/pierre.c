import os
import time
from Crypto.Util.number import inverse

###
# Code from CryptoHack
###

p = 2**192 - 237
a = -3
b = 1379137549983732744405137513333094987949371790433997718123
order = 6277101735386680763835789423072729104060819681027498877478



def dbl(P1):
    X1, Z1 = P1
    XX = X1**2 % p
    ZZ = Z1**2 % p
    A = 2 * ((X1 + Z1) ** 2 - XX - ZZ) % p
    aZZ = a * ZZ % p
    X3 = ((XX - aZZ) ** 2 - 2 * b * A * ZZ) % p
    Z3 = (A * (XX + aZZ) + 4 * b * ZZ**2) % p
    return (X3, Z3)


def diffadd(P1, P2, x0):
    X1, Z1 = P1
    X2, Z2 = P2
    X1Z2 = X1 * Z2 % p
    X2Z1 = X2 * Z1 % p
    Z1Z2 = Z1 * Z2 % p
    T = (X1Z2 + X2Z1) * (X1 * X2 + a * Z1Z2) % p
    Z3 = (X1Z2 - X2Z1) ** 2 % p
    X3 = (2 * T + 4 * b * Z1Z2**2 - x0 * Z3) % p

    return (X3, Z3)


def swap(bit, P1, P2):
    if bit == 1:
        P1, P2 = P2, P1
    return P1, P2


def scalarmult(scalar, x0):
    R0 = (x0, 1)
    R1 = dbl(R0)
    n = scalar.bit_length()
    pbit = 0
    for i in range(n - 2, -1, -1):
        bit = (scalar >> i) & 1
        pbit = pbit ^^ bit
        if pbit:
            R0, R1 = R1, R0
        R1 = diffadd(R0, R1, x0)
        R0 = dbl(R0)
        pbit = bit
    if bit:
        R0 = R1
    if R0[1] == 0:
        return "Infinity"
    return R0[0] * inverse(R0[1], p) % p

###
# Start of the solve
###


k = GF(p)
E = EllipticCurve(k,[a,b])
k2 = GF(p**2)
E2 = E.base_extend(k2)
# Trying to work in the twists
at = 2400511531176665852807485769429804650846392565173449649319
bt = 3079716538687524666927258369235666088983893946802511434274
Et = EllipticCurve(k,[at,bt])
#Et1,Et2 = E.twists()
r = k(at/a)
d = k(r.sqrt())
d**3 == k(bt/b)
sqd = k2(d).sqrt()
def phi(x):
    return k(x/d)
def phi_inv(x):
    return k(x*d)



x0 = 0
G0t = Et.lift_x(x0)
# G0 is a point on Et but not on E
#x0 = 0
#G0 = Et.lift_x(x0)
g0 = G0t.order()
print("ordre de G0 :",factor(g0))
x1 = 2
# G1 is a point on E
G1 = E.lift_x(x1)
g1 = G1.order()
print("ordre de G1 :",factor(g1))
# unhandleable cofactors for DLP
# working in subgroups with 'small' prime factors
cofact0 = 1749002286417992230333906793
cofact1 = 1763644255088983457164385909
G0_ = cofact0 * G0t
G1_ = cofact1 * G1
g0_ = G0_.order()
g1_ = G1_.order()
print("ordre de G0_ :",factor(g0_))
print("ordre de G1_ : ",factor(g1_))


for i in range(5):

    print("[+] Solving Challenge Number : ", i)
    #rand = int.from_bytes(os.urandom(24),"big")
    #privkey = min( rand %order,(order-rand)%order)
    pubkey0 ,  pubkey1 = (9901514941572862118381493448393005071048924854231695716, 4276122774780356582662611332798552909375719748724546125957)
    
    
    P0=Et.lift_x(phi_inv(pubkey0))*cofact0
    t0 = time.time()
    d0 = P0.log(G0_)
    print("(*) time for first DLP :", time.time()-t0)
    P1 = (E.lift_x(pubkey1))*cofact1
    t = time.time()
    d1 = P1.log(G1_)
    print("(*) time for second DLP :", time.time()-t)

    guess1 = CRT_list([d0%g0_ ,d1%g1_],[g0_,g1_]) %order
    guess2 = CRT_list([d0%g0_ ,-d1%g1_],[g0_,g1_]) %order
    guess3 = CRT_list([-d0%g0_ ,d1%g1_],[g0_,g1_]) %order
    guess4 = CRT_list([-d0%g0_ ,-d1%g1_],[g0_,g1_]) %order
    deltat =time.time()-t0
    print("(**) total time :",deltat )
    print("guess 1: ", guess1)
    print("guess 2: ", guess2)
    print("guess 3: ", guess3)
    print("guess 4: ", guess4)













