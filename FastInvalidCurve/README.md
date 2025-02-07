A Comprehensive Introduction To Invalid Curve Attack Using Twists Of Elliptic Curves.

***

## Key Takeaways

+ Invalid Curve Attacks are real threats, commonly exploited. The attack boils down to solve Discrete Log Problems (DLP) on Elliptic Curves (EC).
+ We comprehensively introduce the quadratic twist Ed of E, we describe morphisms from E to Ed.
+ Using the twist, we illustrate the attack **while keeping computations in GF(p)**, avoiding any need of construction of GF(p^2) neither solving DLP in the extended curve E over GF(p^2).
+ We **illustrate the edge gained by dodging GF(p^2)** on an example for the brainnpoolP256t1 curve and in a implementation of the Invalid Curve Attack on a toy curve.

## Materials 

+ _FastInvalidCurveAttack.pdf_ This paper has an expository goal.
+ _brainnpoolP256t1_example.py_ Toy example, it aims at convincing the reader that dodging GF(p^2) is a huge gain.
+ _InvalidCurveAttack.py_ Toy example of the attack.
