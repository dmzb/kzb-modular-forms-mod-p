// This code verifies some claims made in Example 6.16: mod 3 ethereal representations of level 7

// LOAD FUNCTIONS

load "functions.m";
 
// LOAD MOD FORMS FROM EX 6.4

load "6.4.mod-forms.m";

// BEGIN EXAMPLE 6.16

g := x2+2*y2;
traceDist(g);
// Output is <#T0, #T1, #T2, p0, p1, p2> where
//	Tk = terms of y2 with coefficient k, for k = 0,1,2
//	pk = proportion of terms = #Tk/(#T0+#T1+#T2)

// Trace distributions by conjugacy class for subgroups of GL(2,3)

p := 3;
F := GF(p);
G := GL(2, F);
subs := [H`subgroup : H in Subgroups(G)];

for H in subs do
	<Order(H), { Trace(h) : h in H }>;
	end for;

// Collect subgroups matching trace set of g above
goodsubs := [];
for H in subs do
	if { Trace(h) : h in H } eq { 0, 2 } then
		Append(~goodsubs, H);
		end if;
	end for;

for H in goodsubs do
	<Order(H), { Trace(h) : h in H }>;
	end for;

// Compute trace distributions
RR := RealField(5);
for H in goodsubs do
	trace0 := 0;
	trace1 := 0;
	trace2 := 0;
	cclasses := ConjugacyClasses(H);
	for c in cclasses do
		if Trace(c[3]) eq 0 then
			trace0 := trace0+c[2];
		elif Trace(c[3]) eq 1 then
			trace1 := trace1+c[2];
		elif Trace(c[3]) eq 2 then
			trace2 := trace2+c[2];
			end if;
		end for;
	<Order(H), trace0, trace1, trace2, RR!(trace0/Order(H)), RR!(trace1/Order(H)), RR!(trace2/Order(H))>;
	end for;
