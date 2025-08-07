// This code verifies some claims made in Example 5.13: the canonical ring of a (2,3,7) curve and a wild Z/5Z root stack over it

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 5.13

// Stacky structures

// Log canonical ring of underlying (2,3,7) curve
p := 5;
g := 0; // genus of X is 0
delta := 0; // X has no log structure
B := 3*7; // degree bound from KZB
F := GF(p);
sig := [1/2, 2/3, 6/7]; // signature of X as a stacky curve, from stacky structure code above

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

// Log canonical ring of wild root stack over (2,3,7) curve
p := 5;
g := 0; // genus of X is 0
delta := 0; // X has no log structure
m := 1; // set ramification jump
B := 3*7; // degree bound from KZB
F := GF(p); // may need to search over extensions of GF(p)
sig := [1/2, 2/3, 6/7, (p-1)*(m+1)/p]; // refined signature of stacky curve

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
