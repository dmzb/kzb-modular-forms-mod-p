// This code verifies some claims made in Example 7.1: log canonical ring of X_ns^+(3)

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 7.1

H := normalizerNonsplitCartan(3);
for tup in [[0, 0], [0, 1728], [2, 0]] do
	signatureH(H,tup[1],tup[2]);
	end for;

// Log canonical ring of X_H in tame characteristics (e.g. char p = 5)
p := 5;
g := 0; // genus of X is 0
delta := 1; // X has 1 cusp
B := 3*p; // degree bound from KZB
F := GF(p); // may need to search over extensions of GF(p)
sig := [1/2, 1/2, 1/2]; // signature of X as a stacky curve, from stacky structure code above

// Construct projective model of X (it's just P1 for this example)
P3<x,y,z> := PolynomialRing(F, 3);
P2 := ProjectiveSpace(P3);
X := Curve(P2, x+y+z);
// g := Genus(X); // should match g above

// Construct log canonical divisor
K := CanonicalDivisor(X);
pts := RationalPoints(X);
sPts := [Divisor(Random(pts)) : i in [1..#sig]];
if delta eq 0 then
	D := K;
	else D := K + &+[Divisor(Random(pts)) : i in [1..delta]]; // log divisor
	end if;

// This code will return a list of dimensions of Riemann--Roch spaces up to weight B:
//for k in [1..B] do
//	kD := k*D + &+[Floor(k*sig[i])*Divisor(Random(pts)) : i in [1..#sig]];
//	<k, Dimension(RiemannRochSpace(kD))>;
//	end for;

// Display number of generators in each degree up to degree bound
gens := getGeneratorsUpToDegree(<D,sig,sPts>,B);
[#gens[i] : i in [1..#gens]];

// Display degree of relations in minimal presentation of canonical ring
P := PolynomialRing(F, &cat([ [i : j in [1..#gens[i]]] : i in [1..#gens]] ));
    P := ProjectiveSpace(P);
      phi := map<X -> P | &cat(gens)>;
      C := Image(phi);
    [Degree(b) : b in MinimalBasis(C)];


// Log canonical ring of X_H in characteristic 2
p := 2;
g := 0; // genus of X is 0
delta := 1; // X has 1 cusp
// m := 1; // ramification jump at wild pts is 1
B := 3*p; // degree bound from KZB
F := GF(p); // may need to search over extensions of GF(p)
sig := [6/4]; // refined signature of stacky curve

// Construct projective model of X (it's just P1 for this example)
P3<x,y,z> := PolynomialRing(F, 3);
P2 := ProjectiveSpace(P3);
X := Curve(P2, x+y+z);
g := Genus(X); // should match g above

// Construct log canonical divisor
K := CanonicalDivisor(X);
pts := RationalPoints(X);
sPts := [Divisor(Random(pts)) : i in [1..#sig]];
if delta eq 0 then
	D := K;
	else D := K + &+[Divisor(Random(pts)) : i in [1..delta]]; // log divisor
	end if;

// This code will return a list of dimensions of Riemann--Roch spaces up to weight B:
//for k in [1..B] do
//	kD := k*D + &+[Floor(k*sig[i])*Divisor(Random(pts)) : i in [1..#sig]];
//	<k, Dimension(RiemannRochSpace(kD))>;
//	end for;

// Display number of generators in each degree up to degree bound
gens := getGeneratorsUpToDegree(<D,sig,sPts>,B);
[#gens[i] : i in [1..#gens]];

// Display degree of relations in minimal presentation of canonical ring
P := PolynomialRing(F, &cat([ [i : j in [1..#gens[i]]] : i in [1..#gens]] ));
    P := ProjectiveSpace(P);
      phi := map<X -> P | &cat(gens)>;
      C := Image(phi);
    [Degree(b) : b in MinimalBasis(C)];
