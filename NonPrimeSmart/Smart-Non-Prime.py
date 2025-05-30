from tqdm import tqdm
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
from sage.rings.padics.padic_generic import ResidueLiftingMap


###################
# Generating parameters
###################


def get_p_D_n(pbit,max_iter,max_D,n):
    for i in tqdm(range(1,max_D)):
        D = 4*i+3
        r  = 2**pbit
        inv4 = ZZ(pow(4,-1,D))
        for m in range(r,r+max_iter):
            q = D*m*(m+1)+ inv4
            p,b = q.nth_root(n,truncate_mode =1)
            if b == True:
                if is_prime(p):
                    return p,D
    return None


def get_curve(q,D,n):
    db = HilbertClassPolynomialDatabase()
    p = ZZ(q.nth_root(n))
    k = GF(q)
    if D < 10**4:
        HD = db[-D]
    else :
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
    return E



def challenge_DLP(E):
    kE = E.base_field()
    p = kE.characteristic()
    q =p ** kE.degree()
    secret_key = randint(1,q)
    for i in range(1,p):
        Lifts = E.lift_x(ZZ(i),all=True)
        if len(Lifts) !=0:
            P = Lifts[0]
            break
    Q = secret_key * P
    return P,secret_key,Q


###################
# Smart's Attack
###################


def Smart_Attack(P,Q):
    E = P.curve()
    assert E == Q.curve()
    kE = E.base_field()
    p = kE.characteristic()
    n = kE.degree()
    q = p**n
    prec = 10

    modE = kE.modulus()
    Qur.<a> = Qp(p,prec=prec).extension(modE,prec=prec)
    assert Qur.ramification_index() == 1
    k = Qur.residue_field()
    mod = k.modulus()
    ak = k.gens()[0]
    akE = kE.gens()[0]
    assert mod == modE

    f = ResidueLiftingMap._create_(kE,Qur)
    phi_kE_to_k =  FiniteFieldHomomorphism_generic(Hom(kE,k))
    phi_k_to_kE =  FiniteFieldHomomorphism_generic(Hom(k,kE))
    assert phi_kE_to_k(akE) == ak

    EQq =EllipticCurve(Qur,[f(a).lift_to_precision(prec) +randint(0,p)*p for a in E.a_invariants()])

    P_,Q_ = P,Q
    sk_guess_list = []

    for _ in range(n):
        P_Qurs = EQq.lift_x(f(P.x()).lift_to_precision(prec),all=True)

        Q_Qurs = EQq.lift_x(f(Q.x()).lift_to_precision(prec),all=True)

        for P_Qur in P_Qurs:
            if phi_k_to_kE(k(P_Qur.y())) == P.y():
                break
        for Q_Qur in Q_Qurs:
            if phi_k_to_kE(k(Q_Qur.y())) == Q.y():
                break
            
        qP = q*P_Qur
        qQ = q*Q_Qur
        x_P,y_P = qP.xy()
        x_Q,y_Q = qQ.xy()

        phi_P = -(x_P/y_P)
        phi_Q = -(x_Q/y_Q)
        sk_guess  = ZZ(phi_Q/phi_P %p)
        sk_guess_list.append(sk_guess)
        guess =  sum([a*p**i for i,a in enumerate(sk_guess_list)])
        if guess*P_ == Q_:
            break
        Q = Q-sk_guess*P
        P = p*P

    return guess

###################
# Examples
###################

print("* Tiny example with q=p**n, n=4")
print("|- Generating Parameters")
n = 4
p,D = get_p_D_n(3,2**10,2**25,n)
E = get_curve(p**n,D,n)
# p,D,n = 191 ,99763 , 4 is another example
assert E.order() == p**n
print("|- Generating Challenge DLP")
P,sk,Q = challenge_DLP(E)
print("|-Secret key recovered with Smart's attack :",Smart_Attack(P,Q) == sk)


print("* Tiny example with q=p**n, n=2")
print("|- Generating Parameters")
p, D, n = 79998601, 5107, 2
E = get_curve(p**n,D,n)
assert E.order() == p**n
print("|- Generating Challenge DLP")
P,sk,Q = challenge_DLP(E)
print("|- Secret key recovered with Smart's attack : ",Smart_Attack(P,Q) == sk)

