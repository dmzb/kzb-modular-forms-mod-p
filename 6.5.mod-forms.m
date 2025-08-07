// This code verifies some claims made in Example 6.5: q-expansions of mod 2 modular forms of level 13

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 6.5

p := 2;
N := 13;
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
b41 := PP!qExpansion(B4[1], prec);
f1 := PP!qExpansion(B4[2], prec);
f2 := PP!qExpansion(B4[3], prec);
f3 := PP!qExpansion(B4[4], prec);
f4 := PP!qExpansion(B4[5], prec);
b61 := PP!qExpansion(B6[1], prec);
h1 := PP!qExpansion(B6[2], prec);
h2 := PP!qExpansion(B6[3], prec);
h3 := PP!qExpansion(B6[4], prec);
h4 := PP!qExpansion(B6[5], prec);
h5 := PP!qExpansion(B6[6], prec);
h6 := PP!qExpansion(B6[7], prec);

f1+O(q^21);
f2+O(q^21);
f3+O(q^21);
f4+O(q^21);

// Construct ethereal generator
coeffs := Coefficients(PP!qExpansion(B4[3], 2*prec)+PP!qExpansion(B4[5], 2*prec));
y2 := &+[q^((i+1)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2));
f2+f4-y2^2;

// More relations:
b41-x2^2;
f1-(x2*y2+y2^2);
f2-(y2^2+f4);
b61-x2^3;
h1-(x2^2*y2+x2*y2^2);
h2-(x2*y2^2+x2*f4+h5+h6);
h3-(x2*f3+h5+h6);
h4-(x2*f4+h5+h6);
y2*f3-x2*f4;
y2*f4-(y2^3+x2*f3+x2*f4);
f3^2-x2*h6;
