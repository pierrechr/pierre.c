---
layout: post
title: "Non commuting endomorphisms"
date: 2026-03-28
---

<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

TEST
Andrew Sutherland [Problem Set 6](https://ocw.mit.edu/courses/18-783-elliptic-curves-spring-2021/resources/mit18_783s21_ps6/) starts with a nice exercice.
Let $p = 7$, $k = \mathbb{F}_{p^2}$ and let $E/k : y^2 = x^3 + (1 + i)x$ where $i^2 = -1$.

The goal is to give the explicit structure of $\mathop{End}(E)$ as a quaternion algebra.

# Defining finite fields

First of all we will need to define $k$ and its quadratic extension $k_2$ (see below why one needs this extension).
Thus one constructs those fields, but we have no guarantee that such a naive method allows to coerce elements from $k$ to $k_2$.
That's why one constructs an isomorphism between $k$ and the unique subfield of order $p^2$ of $k_2$.


```
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic

k.<i> = GF(7**2, modulus = x**2 + 1)
k2 = GF(7**4)
kp = k2.subfields()[1][0]

phi_k_to_kp =  FiniteFieldHomomorphism_generic(Hom(k, kp))
```


## First question, warm up.


```
E0 = EllipticCurve(k, [1+i, 0])
E0.trace_of_frobenius()

```

shows that the characteristic polynomial of the $p^2$-th power Frobenius is $X^2 - 14 X + 49 = (X - 7)^2$.
So $\pi_{E} = [7]$.


## Non ordinary curve

The next question ...

This is confirmed by and (reference)

```
E0.is_supersingular()
```


Then one can write the curves.
A first sanity check shows that $j(E) = 1728$, thus we introduce $E2 : y^2 = x^3 + x$.


```
E0 = EllipticCurve(k, [1+i, 0])
E = EllipticCurve(kp, [1+phi_k_to_kp(i), 0])
E2 = EllipticCurve(k2, [1, 0])

E.j_invariant()
E2.j_invariant()
```


```
psi = E.base_extend(k2).isomorphism_to(E2)
piE2 = E2.frobenius_isogeny()

tau = psi**(-1)*piE2*psi
(tau**2).rational_maps()
E.multiplication_by_m(-7)


coef1 = tau.rational_maps()[0].subs(x=1)
assert coef1**49 == coef1
r1 = phi_kp_to_k(coef1)

coef2 = tau.rational_maps()[1].subs(y=1)
assert coef2**49 == coef2
r2 = phi_kp_to_k(coef2)


r1**8
r2**8
```
