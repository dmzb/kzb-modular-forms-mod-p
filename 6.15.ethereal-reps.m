// This code verifies some claims made in Example 6.15: mod 2 ethereal representations of level 65

// LOAD FUNCTIONS

load "functions.m";
 
// LOAD MOD FORMS FROM EX 6.8

load "6.8.mod-forms.m";

// BEGIN EXAMPLE 6.15

// Newforms
f1 := x2+x3+x5+x6;
f2 := x2+x3+x5+x9;
f3 := x2+x5+x6+x9;
f4 := x2+x5+x9;
E := [f1, f2, f3, f4];

// Check they're eigenforms:
for i in [1..#E] do
	for ell in [3, 7, 11, 17, 19] do
		<i, ell, heckeAction(E[i], ell)+O(q^100)>;
		end for;
	end for;

// Trace distributions
traceDist(f1);
traceDist(f2);
traceDist(f3);
traceDist(f4);
// Output is <#T0, #T1, #T2, p0, p1, p2> where
//	Tk = terms of f with coefficient k, for k = 0,1,2
//	pk = proportion of terms = #Tk/(#T0+#T1+#T2)

// Classical newforms
new1 := Newforms(M2)[1][1];
new2 := Newforms(M2)[2][1]; // irrational
new3 := Newforms(M2)[3][1]; // irrational

degs := [1];
Append(~degs, Degree(BaseRing(Parent(new1))));
Append(~degs, Degree(BaseRing(Parent(new2))));
Append(~degs, Degree(BaseRing(Parent(new3))));
m := LCM(degs);
F<a> := GF(p^m);
PP<q> := PowerSeriesRing(F);

// Obtaining mod 2 q-expansions
new1 := qExpansion(new1, prec);
BaseRing(Parent(new1));
// new1 has rational q-expansion, so we can reduce it mod 2
n1 := PP!new1;

new2 := qExpansion(new2, prec);
BaseRing(Parent(new2));
// new2 has irrational q-expansion, so we need the following workaround
K2<b2> := BaseRing(Parent(new2));
OK2 := MaximalOrder(K2);
PK2 := Factorization(p*OK2)[1][1];
FK2<a2>, redK2 := ResidueClassField(PK2);
Embed(FK2, F);
n2 := &+[redK2(Coefficients(new2)[i])*q^i : i in [1..#Coefficients(new2)]] + O(q^prec);

new3 := qExpansion(new3, prec);
BaseRing(Parent(new3));
// new3 has irrational q-expansion, so we need the following workaround
K3<b3> := BaseRing(Parent(new3));
OK3 := MaximalOrder(K3);
PK3 := Factorization(p*OK3)[1][1];
FK3<a3>, redK3 := ResidueClassField(PK3);
Embed(FK3, F);
n3 := &+[redK3(Coefficients(new3)[i])*q^i : i in [1..#Coefficients(new3)]] + O(q^prec);

// Relations
f1-n1;
f1-n2;
f1-n3;
