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
# Trying to work in the twists
at = 2400511531176665852807485769429804650846392565173449649319
bt = 3079716538687524666927258369235666088983893946802511434274
Et = EllipticCurve(k,[at,bt])
#Et1,Et2 = E.twists()
r = k(at/a)
d = k(r.sqrt())
# Checking the equation of the quadratic twist, actually useless.
d**3 == k(bt/b)
sqd = k2(d).sqrt()
def phi(x):
    return k(x/d)
def phi_inv(x):
    return k(x*d)




# G0 is a point on Et but not on E
x0 = 0
G0t = Et.lift_x(x0)
g0 = G0t.order()
print("G0's order :",factor(g0))

# G1 is a point on E
x1 = 2
G1 = E.lift_x(x1)
g1 = G1.order()
print("G1's order :",factor(g1))
# unhandleable cofactors for DLP
# working in subgroups with 'small' prime factors
cofact0 = 1749002286417992230333906793
cofact1 = 1763644255088983457164385909
G0_ = cofact0 * G0t
G1_ = cofact1 * G1
g0_ = G0_.order()
g1_ = G1_.order()
print("G0_'s order :",factor(g0_))
print("G1_'s order : ",factor(g1_))


for i in range(5):

    print("+ Solving Challenge Number : ", i)
    rand = int.from_bytes(os.urandom(24),"big")
    privkey = min( rand %order,(order-rand)%order)
    pubkey0, pubkey1 = scalarmult(privkey,x0),scalarmult(privkey,x1)

    P0=Et.lift_x(phi_inv(pubkey0))*cofact0
    t0 = time.time()
    d0 = P0.log(G0_)
    print("| time for first DLP :", time.time()-t0)
    P1 = (E.lift_x(pubkey1))*cofact1
    t = time.time()
    d1 = P1.log(G1_)
    print("| time for second DLP :", time.time()-t)
    for s1,s2 in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        if CRT_list([s1*d0 ,s2*d1],[g0_,g1_]) %order == privkey:
            print("|- \U0001F600 "+ " "+ "Challenge solved, privkey recovered in ", time.time() - t0," seconds")










