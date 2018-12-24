# IndexCalculusDiscreteLogarithm
Method of Index Calculus to find the Discrete Logarithm

Given an element `a` of $F_p$ and a generator `g`, find `X` such as $a=g \mod p$

Pari/GP implementation:
  - open gp
  - `\r [path]/main.gp`
  - `LogDiscret(a, g)` where `g` is `Mod(g, p)`
