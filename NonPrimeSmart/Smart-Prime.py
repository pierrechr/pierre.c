import time
import sympy.ntheory as snt
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial as HCP

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


def check_anomalous(p,D):
    k = GF(p)
    assert D < 10**4
    HD = db[D]
    HDq = HD.change_ring(k)
    j = HDq.any_root()
    assert j not in [0,1728]
    
    A = 3*j*(1728-j)
    B = 2*j*(1728-j)**2
    E = EllipticCurve(k,[A,B])
    P = E.random_point()
    
    if p*P != E.zero() :
        E = E.quadratic_twist()
        P = E.random_point()
        if p*P == E.zero() :
            return E
        
    return None

def generate_challenge(E):
    p = E.base_ring().characteristic()
    P = E.random_point()
    sk = randint(1,q)
    return (P,sk, sk*P)
    
def SmartAttack(P,Q):
    E = P.curve()
    p = E.base_ring().characteristic()
    k = GF(p)
    Eqp = EllipticCurve(Qp(p), [ ZZ(a) + randint(0,p)*p for a in E.a_invariants() ])

    P_Qps = Eqp.lift_x(ZZ(P.x()), all=True)
    for P_Qp in P_Qps:
        if k(P_Qp.y()) == P.x():
            break

    Q_Qps = Eqp.lift_x(ZZ(Q.x()), all=True)
    for Q_Qp in Q_Qps:
        if k(Q_Qp.y()) == Q.y():
            break

    p_times_P = p*P_Qp
    p_times_Q = p*Q_Qp

    x_P,y_P = p_times_P.xy()
    x_Q,y_Q = p_times_Q.xy()

    phi_P = -(x_P/y_P)
    phi_Q = -(x_Q/y_Q)
    k = phi_Q/phi_P
    return ZZ(k)


p,D = get_p_D()
E = check_anomalous(p,D)
assert  E != None
P, sk,Q = generate_challenge(E)
SmartAttack(P,Q) == sk










