P = Primes()

for i in range(2,20):
    p = P.unrank(i)
    m = (p-1)//2
    k = GF(p)
    Lt = set()
    kxy.<x,y> = PolynomialRing(k)
    S = set(k)
    S.remove(k(0))
    S.remove(k(1))
    for l in S:
            R = y**2-x*(x-1)*(x-l)
            ER = EllipticCurve(R)
            Lt.add(ER.order())
    m = ceil(float(-(sqrt(p)+1)**2/4))
    M = floor(float(-(sqrt(p)-1)**2/4))
    mE  = -4*M
    ME = -4*m
    print(mE,ME,Lt)
    assert all([ a in list(range(mE,ME+1)) for  a in Lt]) == True


