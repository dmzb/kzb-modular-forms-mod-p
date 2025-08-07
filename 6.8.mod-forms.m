// This code verifies some claims made in Example 6.8: q-expansions of mod 2 modular forms of level 65

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 6.8

p := 2;
N := 65;
F<a> := GF(p^2); // there are irrational newforms to look for
M2 := ModularForms(N, 2);
M4 := ModularForms(N, 4);
B2 := Basis(M2);
B4 := Basis(M4);
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

// Notice that B4[24] is a square mod 2
B4[24]+O(q^200);

// Construct a square root
coeffs := Coefficients(PP!qExpansion(B4[24], 2*prec));
x9 := &+[q^((i+25)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2)); // note: need to start at the correct exponent
x9^2-B4[24];

// Here's another square
B4[9]+B4[13]+O(q^30);

// Construct a square root
coeffs := Coefficients(PP!qExpansion(B4[9]+B4[13], 2*prec));
x10 := &+[q^((i+7)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2));
x10^2-(B4[9]+B4[13]);

// Create oldforms
coeffs := Coefficients(PP!qExpansion(Basis(ModularForms(5, 4))[3], 2*prec)); 
y5 := &+[q^((i+1)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2));
B5 := [y5];

coeffs := Coefficients(PP!qExpansion(Basis(ModularForms(13, 4))[3]+Basis(ModularForms(13, 4))[5], 2*prec));
y13 := &+[q^((i+1)/2) : i in [1..#coeffs] | coeffs[i] ne 0] + O(q^((#coeffs+3)/2));
B13 := [y13];

oldForms := [];
for b in B5 do
	Ub := Uoperator(b, 13);
	Append(~oldForms, b);
	Append(~oldForms, Ub);
	end for;
for b in B13 do
	Ub := Uoperator(b, 5);
	Append(~oldForms, b);
	Append(~oldForms, Ub);
	end for;

// Linear relation: oldForms[4] = oldForms[1] + oldForms[2] + oldForms[3]
oldForms[4]-(oldForms[1]+oldForms[2]+oldForms[3]);
old1 := oldForms[1];
old2 := oldForms[2];
old3 := oldForms[3];
oldBasis := [old1, old2, old3];

// More relations:
x9-old2;
x6+x9-oldForms[4];
x9^2+x1*x9-B4[14];
y5-(x2+x3+x5+x6+x9);
