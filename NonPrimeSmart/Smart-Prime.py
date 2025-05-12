import time
db = HilbertClassPolynomialDatabase()

def get_p_D(pbit=250,delta=2**10,D=43):

    d = D.bit_length()
    r  = 2**((pbit-d)//2)
    inv4 = ZZ(pow(4,-1,D))

    for m in range(r,r+delta):
        p = D*m*(m+1)+inv4
        if is_prime(p):
            return p,D

    return None


def get_curve(p,D):
    k = GF(p)
    if D < 10**4:
        HD = db[D]
    else:
        K.<a> = NumberField(x**2+D)
        HD = K.hilbert_class_polynomial()

    HDq = HD.change_ring(k)
    j = HDq.any_root()
    assert j not in [0,1728]

    A = 3*j*(1728-j)
    B = 2*j*(1728-j)**2
    E = EllipticCurve(k,[A,B])
    P = E.random_point()

    if p*P != E.zero() :
        E = E.quadratic_twist()

    return E

def generate_challenge(E):
    p = E.base_ring().characteristic()
    P = E.random_point()
    sk = randint(1,p)
    return (P,sk, sk*P)

def SmartAttack(P,Q):
    E = P.curve()
    p = E.base_ring().characteristic()
    k = GF(p)
    Eqp = EllipticCurve(Qp(p,2), [ ZZ(a) + randint(0,p)*p for a in E.a_invariants() ])

    P_Qps = Eqp.lift_x(ZZ(P.x()), all=True)
    for P_Qp in P_Qps:
        if k(P_Qp.y()) == P.y():
            break

    Q_Qps = Eqp.lift_x(ZZ(Q.x()), all=True)
    for Q_Qp in Q_Qps:
        if k(Q_Qp.y()) == Q.y():
            break

    pP = p*P_Qp
    pQ = p*Q_Qp

    x_P,y_P = pP.xy()
    x_Q,y_Q = pQ.xy()

    phi_P = -(x_P/y_P)
    phi_Q = -(x_Q/y_Q)
    k = phi_Q/phi_P
    return ZZ(k)

print("Generating Parameters")
p,D = get_p_D()
E =  get_curve(p,D)
assert  E.order() == p
print("Generating Challenge DLP")
P, sk,Q = generate_challenge(E)
print("Secret key recovered with Smart's attack : ",SmartAttack(P,Q) == sk)












