---
layout: post
title: "Non commuting endomorphisms"
date: 2026-03-28
---

<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

TEST
Andrew Sutherland [Problem Set 6](https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_ps6/) starts with a nice exercice.
Let 
$$p = 7,\quad  k = \mathbb{F}_{p^2}$ and $$ E/k : y^2 = x^3 + (1 + i)x$$ with $$i^2 = -1$$.

The goal is to give the explicit structure of 
$$\mathop{End}(E)$$
 as a quaternion algebra.

# Defining finite fields

First of all we will need to define the base field and its quadratic extension (see below why one needs this extension).
Thus one constructs those fields, but we have no guarantee that such a naive method allows to coerce elements from one to another.
That's why one constructs an isomorphism between fields of order $$p^2$$.


```
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic

k.<i> = GF(7**2, modulus = x**2 + 1)
k2 = GF(7**4)
kp = k2.subfields()[1][0]

phi_kp_to_k =  FiniteFieldHomomorphism_generic(Hom(kp, k))
phi_k_to_kp =  FiniteFieldHomomorphism_generic(Hom(k, kp))
```


## First question, warm up.


```
E0 = EllipticCurve(k, [1+i, 0])
E0.trace_of_frobenius()
```

This shows that the characteristic polynomial of the Frobenius is 
$$X^2 - 14 X + 49 = (X - 7)^2$$ so $$\pi_{E} = [7]$$.


## Frobenius, isn't it ?  

The easiest path for this question is to take random points $$P$$ on $$E$$ and try to apply $$pi_p$$, the $p$-th power Frobenius.
If one fails, it means that $$(x; y)$$ is on $$E$$ but not $$(x^7; y^7)$$ : $$pi_p$$ won't be an endomorphism of $E$.

```
while True:
    P = E0.random_point()
    x, y, _= P
    try:
        E0.lift_x(x**7)
    except:
        print("pi_p not an endomorphism of E")
        break
```


## Let's twist again 

Since the polynomial in x defining E has no constant term in short Weierstrass form, one has $$ j(E) = 1728 $$ and it is isomorphic (over the algebraic closure thought) to $$ E_2/k : y^2 = x^3 + x$$.
Actually, $$E$$ and $$E2$$ are quadratic twists, so isomorphic over a quadratic extension of the base field.

Thus, one has the following three curves

```
E0 = EllipticCurve(k, [1+i, 0])
E = EllipticCurve(kp, [1+phi_k_to_kp(i), 0])
E2 = EllipticCurve(k2, [1, 0])

E.j_invariant()
E2.j_invariant()
```

Something nice with $$E2$$ is that it DOES have an isogeny of degree $$7$$.
This allows to get an endomorphism $$\tau$$ of $$E$$ such that $$\tau^2 = [-7]$$

```
psi = E.base_extend(k2).isomorphism_to(E2)
piE2 = E2.frobenius_isogeny()

tau = psi**(-1)*piE2*psi

(tau**2).rational_maps()
E.multiplication_by_m(-7)
```

```
coef1 = tau.rational_maps()[0].subs(x=1)
assert coef1**49 == coef1
r1 = phi_kp_to_k(coef1)

coef2 = tau.rational_maps()[1].subs(y=1)
assert coef2**49 == coef2
r2 = phi_kp_to_k(coef2)


r1**8
r2**8
```
