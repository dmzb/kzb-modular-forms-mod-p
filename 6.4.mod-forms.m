// This code verifies some claims made in Example 6.4: q-expansions of mod 3 modular forms of level 7

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 6.4

p := 3;
N := 7;
F := GF(p);

// Compute bases of classical modular forms
M2 := ModularForms(N, 2);
M4 := ModularForms(N, 4);
M6 := ModularForms(N, 6);
B2 := Basis(M2);
B4 := Basis(M4);
B6 := Basis(M6);
prec := 1000; // set precision

// Get mod p q-expansions:
PP<q> := PowerSeriesRing(F);
x2 := PP!qExpansion(B2[1], prec);
f1 := PP!qExpansion(B4[1], prec);
f2 := PP!qExpansion(B4[2], prec);
f3 := PP!qExpansion(B4[3], prec);
b61 := PP!qExpansion(B6[1], prec);
h1 := PP!qExpansion(B6[2], prec);
h2 := PP!qExpansion(B6[3], prec);
h3 := PP!qExpansion(B6[4], prec);
h4 := PP!qExpansion(B6[5], prec);

f1+O(q^21);
f2+O(q^21);
f3+O(q^21);

// Ethereal generator
y2 := f1;

// More relations:
f2-(x2*y2-y2^2);
f3-(x2^2+x2*y2+y2^2);
b61-y2^3;
h1-(2*x2^3+x2*y2^2+2*h4);
h2-(x2^2+x2*y2^2+y2^3+h4);
h3-(x2^3+2*y2^3);
