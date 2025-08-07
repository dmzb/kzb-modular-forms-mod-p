// This code verifies some claims made in Example 6.3: q-expansions of mod 2 modular forms of level 5

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 6.3

p := 2;
N := 5;
F := GF(p);

// Compute bases of classical modular forms
M2 := ModularForms(N, 2);
M4 := ModularForms(N, 4);
B2 := Basis(M2);
B4 := Basis(M4);
prec := 1000; // set precision

// Get mod p q-expansions:
PP<q> := PowerSeriesRing(F);
x2 := PP!qExpansion(B2[1], prec);
f1 := PP!qExpansion(B4[1], prec);
f2 := PP!qExpansion(B4[2], prec);
f3 := PP!qExpansion(B4[3], prec);

f1+O(q^21);
f2+O(q^21);
f3+O(q^21);

// Observe that y4 is a square
f3+O(q^21);

// Construct a square root
coeffs := Coefficients(PP!qExpansion(B4[3], 2*prec));
y2 := &+[q^((i+1)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2));
y2+O(q^21);
y2^2-f3;

// More relations:
f1-x2^2;
f2-(y2^2+x2*y2);
