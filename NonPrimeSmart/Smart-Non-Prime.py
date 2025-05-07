import time
from tqdm import tqdm
db = HilbertClassPolynomialDatabase()

def get_p_D(pbit,max_iter,max_D):
    """
    Given a bit length pbit and a range delta, returns a prime p and absolute dicriminant D such that
    4*p**2 = 1 + m*(m+1)*D
    or returns None if such a solution was not found in the given ranges.
    """
    for i in tqdm(range(1,max_D)):
        D = 4*i+3
        r  = 2**pbit
        inv4 = ZZ(pow(4,-1,D))
        for m in range(r,r+max_iter):
            q = D*m*(m+1)+ inv4
            if q.is_square():
                p = ZZ(sqrt(q))
                if is_prime(p):
                    return p,D
    return None


def check_q_D(q,D):
    p = ZZ(sqrt(q))
    k = GF(q)
    
    if D <10**4:
        HD = db[-D]
    else:
        K.<a> = NumberField(x**2+D)
        HD = K.hilbert_class_polynomial()

    Hq = HD.change_ring(k)
    j = Hq.any_root()
    if j not in [0,1728]:
                    A = 3*j*(1728-j)
                    B = 2*j*(1728-j)**2
                    E = EllipticCurve(k,[A,B])
    if j == 0:
                    print("cas j = 0")
                    E = EllipticCurve(k,[0,1])
    if  j == 1728:
                    print("cas j = 1728")
                    E = EllipticCurve(k,[1,0])

    P = E.random_point()
    if q*P != E.zero():
        E = E.quadratic_twist()
    P = E.random_point()
    return q*P == E.zero()


########################
# Generating a Curve and a challenge
########################


p = 79998601
q,D = p**2, 5107

assert check_q_D(q,D) == True

def get_curve(q,D):
    k= GF(q)
    if D < 10**4:
        HD = db[-D]
    else :
        K.<a> = NumberField(x**2+D)
        HD = K.hilbert_class_polynomial()
    HDq = HD.change_ring(k)
    roots = HDq.any_root()
    j = roots
    if j not in [0,1728]:
                    A = 3*j*(1728-j)
                    B = 2*j*(1728-j)**2
                    E = EllipticCurve(k,[A,B])
    if j == 0:
                    E = EllipticCurve(k,[0,1])
    if  j == 1728:
                    E = EllipticCurve(k,[1,0])
    if E.order() != q:
                    E = E.quadratic_twist()
    return E

t = time.time()
E = get_curve(q,D)
print("Curve generated in ",time.time()-t," sec")


def challenge_DLP(E):
    p = E.base_field().characteristic()
    secret_key = randint(1,q)
    for i in range(1,p):
        Lifts = E.lift_x(ZZ(i),all=True)
        if len(Lifts) !=0:
            P = Lifts[0]
            break
    Q = secret_key * P
    return P,secret_key,Q

P,sk,Q = challenge_DLP(E)


########################
#  Clean Solve
########################

def Smart_Attack(P,Q):
    E = P.curve()
    p = E.base_ring().characteristic()
    q = p**2
    prec = 10**2
    kE = E.base_field()
    modE = kE.modulus()
    Qur.<a> = Qp(p,prec=prec).extension(modE,prec=prec)
    assert Qur.ramification_index() == 1
    k = Qur.residue_field()
    mod = k.modulus()
    ak = k.gens()[0]
    akE = kE.gens()[0]
    assert mod == modE

    from sage.rings.padics.padic_generic import ResidueLiftingMap
    from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic

    f = ResidueLiftingMap._create_(kE,Qur)
    phi_kE_to_k =  FiniteFieldHomomorphism_generic(Hom(kE,k))
    phi_k_to_kE =  FiniteFieldHomomorphism_generic(Hom(k,kE))

    assert phi_kE_to_k(akE) == ak

    # randomize the lift to avoid the canonical lift
    EQq =EllipticCurve(Qur,[f(a).lift_to_precision(prec) +randint(0,p)*p for a in E.a_invariants()])

    # Beware that the lifting map f lifts with precision 1 by default.
    # One changes that to make use of the isomorphism E(Qq)/E_1(Qq) to E(Fq)
    P_Qurs = EQq.lift_x(f(P.x()).lift_to_precision(prec),all=True)
    Q_Qurs = EQq.lift_x(f(Q.x()).lift_to_precision(prec),all=True)


    J,K = P_Qurs
    assert ( phi_k_to_kE(k(J.y())) == P.y() ) or (phi_k_to_kE(k(K.y())) == P.y())
    # E(phi_k_to_kE(k(J.x())),phi_k_to_kE(k(J.y())))


    for P_Qur in P_Qurs:
        if phi_k_to_kE(k(P_Qur.y())) == P.y():
            break
    for Q_Qur in Q_Qurs:
        if phi_k_to_kE(k(Q_Qur.y())) == Q.y():
            break

    pP = q*P_Qur
    pQ = q*Q_Qur
    x_P,y_P = pP.xy()
    x_Q,y_Q = pQ.xy()

    phi_P = -(x_P/y_P)
    phi_Q = -(x_Q/y_Q)
    sk_guess = phi_Q/phi_P
    sk_guess = ZZ(sk_guess%p)


    P2 = p*P
    Q2 = (Q-sk_guess*P)

    P2_Qurs = EQq.lift_x(f(P2.x()).lift_to_precision(prec),all=True)
    Q2_Qurs = EQq.lift_x(f(Q2.x()).lift_to_precision(prec),all=True)


    for P2_Qur in P2_Qurs:
        if phi_k_to_kE(k(P2_Qur.y())) == P2.y():
            break
    for Q2_Qur in Q2_Qurs:
        if phi_k_to_kE(k(Q2_Qur.y())) == Q2.y():
            break

    pP2 = q*P2_Qur
    pQ2 = q*Q2_Qur
    x_P,y_P = pP2.xy()
    x_Q,y_Q = pQ2.xy()

    phi_P = -(x_P/y_P)
    phi_Q = -(x_Q/y_Q)
    sk_guess2 = phi_Q/phi_P
    sk_guess2 = ZZ(sk_guess2%p)

    return  sk_guess+p*sk_guess2

sk == Smart_Attack(P,Q)
