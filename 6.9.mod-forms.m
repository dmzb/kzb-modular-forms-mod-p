// This code verifies some claims made in Example 6.9: q-expansions of mod 3 modular forms of level 91

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 6.9

p := 3;
N := 91;
F<a> := GF(p^2); // there are irrational newforms to look for
M2 := ModularForms(N, 2);
M4 := ModularForms(N, 4);
M6 := ModularForms(N, 6);
B2 := Basis(M2);
B4 := Basis(M4);
B6 := Basis(M6);
prec := 1000; // set precision

PP<q> := PowerSeriesRing(F);
x1 := PP!qExpansion(B2[1], prec);
x2 := PP!qExpansion(B2[2], prec);
x3 := PP!qExpansion(B2[3], prec);
x4 := PP!qExpansion(B2[4], prec);
x5 := PP!qExpansion(B2[5], prec);
x6 := PP!qExpansion(B2[6], prec);
x7 := PP!qExpansion(B2[7], prec);
x8 := PP!qExpansion(B2[8], prec);
x9 := PP!qExpansion(B2[9], prec);
x10 := PP!qExpansion(B2[10], prec);

// Construct cube roots
coeffs := Coefficients(PP!qExpansion(B6[22], 2*prec));
x11 := &+[coeffs[i]*q^((i+20)/3) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+2)/3));
x11^3-B6[22];

coeffs := Coefficients(PP!qExpansion(B6[4]+2*B6[7]+B6[16], 2*prec));
x12 := &+[coeffs[i]*q^((i+2)/3) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+2)/3));
x12^3-(B6[4]+2*B6[7]+B6[16]);

// Create oldforms
y7 := PP!qExpansion(Basis(ModularForms(7, 2))[1], prec) + 2*PP!qExpansion(Basis(ModularForms(7, 4))[1], prec);
B7 := [y7];

y13 := 2*PP!qExpansion(Basis(ModularForms(13, 2))[1], prec) + PP!qExpansion(Basis(ModularForms(13, 4))[1], prec);
B13 := [y13];

oldForms := [];
for b in B7 do
	Ub := Uoperator(b, 13);
	Append(~oldForms, b);
	Append(~oldForms, Ub);
	end for;
for b in B13 do
	Ub := Uoperator(b, 7);
	Append(~oldForms, b);
	Append(~oldForms, Ub);
	end for;

// Linear relation: oldForms[4] = 2*oldForms[1] + oldForms[2] + oldForms[3]
old1 := oldForms[1];
old2 := oldForms[2];
old3 := oldForms[3];
oldBasis := [old1, old2, old3];

// More relations:
x11-oldForms[4];
y7-(x2+x4+x5+2*x8+x10+2*x11);
y13-(x2+x4+x5+x8+x10+x11);
