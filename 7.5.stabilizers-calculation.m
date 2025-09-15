// This code verifies some claims made in the proof of Theorem 7.5 (calculation of stabilizers)

// DZB: improve this slightly before submitting.

P<t> := PolynomialRing(Rationals());

C3:=CyclicGroup(3);
 C4:=CyclicGroup(4);
AC3:=AutomorphismGroup(C3);
 phi:= hom< C4 -> AC3 | < C4.1,AC3.2 >>;                
 G3:=SemidirectProduct(C3,C4,phi);
CT3 := CharacterTable(G3);  
 chi3 := CT3[5];
 {H : sub in Subgroups(G3) |  (not IsAbelian(H)) and (not IsIrreducible(Restriction(chi3,H))) where H is sub`subgroup };
 {*<Order(g), t^2 - (Integers()!(chi3(g)))*t + 1 > : g in G3*};
 {H : sub in Subgroups(G3) |  (not IsCyclic(H)) and IsIrreducible(Restriction(chi3,H)) where H is sub`subgroup };

Q := SmallGroups(8)[4];
 AQ:=AutomorphismGroup(Q);
phi:= hom< C3 -> AQ | < C3.1,AQ.2*AQ.1 > >;
 G2:=SemidirectProduct(Q,C3,phi);
CT2 := CharacterTable(G2);  
 chi2 := CT2[4];
 {H : sub in Subgroups(G2) |  (not IsCyclic(H)) and (not IsIrreducible(Restriction(chi2,H))) where H is sub`subgroup };
 {*<Order(g), t^2 - (Integers()!(chi2(g)))*t + 1 > : g in G2*};
