import math

P = Primes()

for i in range(2,11):
    p = P.unrank(i)
    m = (p-1)//2
    k = GF(p)
    kx.<x> = PolynomialRing(k)
    H  = sum([binomial(m,i)**2*x**i for i in range(m+1)])
    LH = set([H(a) for a in k if a not in [k(0),k(1)]])
    Lt = set()
    for a in k:
        for b in k:
            try:
                E = EllipticCurve([a,b])
            except:
                pass
            else:
                Lt.add(E.trace_of_frobenius())

    print("p = {},len(LH) ={},len(Lt) = {}, Hasse={} ".format(p,len(LH),len(Lt),2*math.trunc(float(sqrt(p)*2))+1))
    m = float((sqrt(p)-1)**2/4)
    M = float((sqrt(p)+1)**2/4)

    floor(M)-ceil(m)+1, ((p-1)//2)%2, LH


