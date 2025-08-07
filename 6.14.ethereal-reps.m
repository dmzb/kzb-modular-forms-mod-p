// This code verifies some claims made in Example 6.14: mod 2 ethereal representations of level 13

// LOAD FUNCTIONS

load "functions.m";
 
// LOAD MOD FORMS FROM EX 6.5

load "6.5.mod-forms.m";

// BEGIN EXAMPLE 6.14

traceDist(y2);
// Output is <#T0, #T1, #T2, p0, p1, p2> where
//	Tk = terms of y2 with coefficient k, for k = 0,1,2
//	pk = proportion of terms = #Tk/(#T0+#T1+#T2)
