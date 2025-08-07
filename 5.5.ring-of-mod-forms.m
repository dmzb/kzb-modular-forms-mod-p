// This code verifies some claims made in Example 5.5: the ring of modular forms of level 5 in char 0 and char 2

// LOAD FUNCTIONS

load "functions.m";

// BEGIN EXAMPLE 5.5

// Stacky structures

N := 5;
for tup in [[0, 0], [0, 1728], [2, 0], [3, 0]] do
	signatureH(borel(N),tup[1],tup[2]);
	end for;
// Output is <N, p, j, stackiness> where
//	N = level
//	p = characteristic
//	j = stacky j-invariant
//	stackiness = list of orders of nontrivial aut groups, with multiplicity
//		e.g. <* 2^^2, 3 *> means one two stacky points of order 2 and one of order 3

// Log canonical ring of X_0(5) in tame characteristics (e.g. char p = 7)
p := 7;
g := 0; // genus of X is 0
delta := 2; // X has 2 cusps
B := 3*p; // degree bound from KZB
F := GF(p); // may need to search over extensions of GF(p)
sig := [1/2, 1/2]; // signature of X as a stacky curve, from stacky structure code above

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

// Log canonical ring of X_0(5) in characteristic 2
p := 2;
g := 0; // genus of X is 0
delta := 2; // X has 2 cusps
m := 1; // ramification jump at wild pts is 1
B := 3*p; // degree bound from KZB
F := GF(p); // may need to search over extensions of GF(p)
sig := [(p-1)*(m+1)/p]; // refined signature of stacky curve

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
