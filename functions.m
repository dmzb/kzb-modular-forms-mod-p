// LIST OF FUNCTIONS FOR COMPUTATIONS OF STACKY STRUCTURES, CANONICAL RINGS, RINGS OF MODULAR FORMS, q-EXPANSIONS, ETC.

// The first few functions are for computing stacky structures of modular curves

// Input: an integer.
// Output: representatives of Aut E acting on E[N], 
// 	where E is the elliptic curve over K with j(E) = j
// 	and K is a field of characteristic p
function autReps(N,p,j)
  ZZ := Integers();
  GL2 := GL(2,Integers(N));
  SL2 := SL(2,Integers(N));

  i := SL2![0,-1,1,0];       
  wInit := SL2![-1,1 ,-1,0];

  if not p in [2,3] then
    if j eq 0 then
      auts := sub<GL2 | [wInit,-Identity(SL2)]>;
    else 
      auts := sub<GL2 | i>;
    end if;
  elif p eq 2 then
    w := Random({w : w in Conjugates(SL2,wInit) |      w^2*(w^2*i*w)*w eq i * (w^2*i*w)  });
    auts := sub<GL2 | [i,w]>;
  elif p eq 3 then
    w := Random({w : w in Conjugates(SL2,wInit) | i*w eq w^2*i});
    auts := sub<GL2 | [i,w]>;
  end if;
  return auts;   
end function;

// Compute H structures and stackyness above 0 = 12^3 in characteristic p
// Optional paramater auts for precomputed auts
function signatureH(H,p,j : auts := {})
    N := Characteristic(BaseRing(H));
      GL2 := GL(2,ResidueClassRing(N));
      generic := -Identity(H) in H select 2 else 1;
    mPerm,perm := CosetAction(GL2,H);
    auts := #(auts) eq 0 select autReps(N,p,j) else auts;
      orbits     := {{i^(mPerm(aut)) : aut in auts} : i in GSet(perm)}; 
      stackyness := {* #auts/(generic*#orb) : orb in orbits*};
    return <N,p,j, <{*n : n in stackyness | n gt 1*}>>;  
end function;

// Gamma0(N)
function borel(N)
  R := ResidueClassRing(N);
  G := GL(2,R);  
  gens := [G! [1,1,0,1]] 
      cat [G![a,0,0,1] : a in [1..N] | IsCoprime(a,N) ]
      cat [G![1,0,0,a] : a in [1..N] | IsCoprime(a,N) ];
  return sub<G | gens>;
end function;

// Gamma1(N)
function borel1(N)
  R := ResidueClassRing(N);
  G := GL(2,R);  
  gens := [G! [1,1,0,1]] 
      cat [G![1,0,0,a] : a in [1..N] | IsCoprime(a,N) ];
  return sub<G | gens>;
end function;

// Non-split Cartan subgroup
function nonsplitCartan(N)
	R := ResidueClassRing(N);
	del := [a : a in R | not IsSquare(a)][1];
	G := GL(2,R);
	gens := [G![a, del*b, b, a] : a,b in R | (a ne 0 or b ne 0) and IsUnit(a^2 - del*b^2)];
	return sub<G | gens>;
	end function;

// Normalizer of non-split Cartan subgroup
function normalizerNonsplitCartan(N)
  R := ResidueClassRing(N);
  del := [a : a in R | not IsSquare(a)][1];
  G := GL(2,R);
  gens := [G![0,-1,1,0]]
      cat [G![a,del*b,b,a] : a,b in [1..N] | (not a*b eq 0) and IsCoprime(a*b,N) ]
      cat [G![a,del*b,-b,-a] : a,b in [1..N] | (not a*b eq 0) and IsCoprime(a*b,N) ];
  return sub<G | gens>;
end function;

// The next few functions are for computing canonical rings

// Constructs stacky divisor
function nD(data,n)
  D := data[1];
  sig := data[2];
  sPts := data[3];
  return n*D + &+[Floor(n*sig[i])*sPts[i] : i in [1..#sig]];
end function;

// Multiplication functions
function rrMultiplication(D1,D2)
  R1,m1 := RiemannRochSpace(D1);
  R2,m2 := RiemannRochSpace(D2);
  R3,m3 := RiemannRochSpace(D1+D2);
  i := Inverse(m3);
  S := sub<R3|[i(m1(a)*m2(b))  : a in Basis(R1), b in Basis(R2)]>;
  return S, m3, Dimension(S), Dimension(R3);
end function;

// Get generators of canonical ring
// Input is [D, sig, sPts], and N (stopping degree)
// Output is [gens1,...,gensN]
function getGeneratorsUpToDegree(data,N)
  R,m := RiemannRochSpace(nD(data,1));
  gens := [[m(b) : b in Basis(R)] ];
  for i in [2..N] do
    genImages := &cat[
    [ &*[tup[j] : j in [1..#part]] : tup in CartesianProduct([gens[j] : j in part])]
      : part in Partitions(i) | not part eq [i] ];
    R,m := RiemannRochSpace(nD(data,i));
    mInv := Inverse(m);
    imgBas := Basis(sub<R|[mInv(f) : f in genImages]>);
    gens[i] := [m(b) : b in ExtendBasis(imgBas, R) | not b in imgBas ];    
  end for;
  return gens;
end function;

// The next two functions are for computing Hecke actions on q-expansions

PP<q> := PowerSeriesRing(GF(2));
function heckeAction(f, ell) // f should have a_1 = 1, or else need to tweak starting place using init
	exps := Exponents(f);
	if #exps eq 0 then
		init := 1;
	else
		init := exps[1];
		end if;
	coefs := Coefficients(f);
	newcoefs := [];
	for i in [1..#coefs] do
		j := i+init-1; // Tracks actual term in q-expansion
		k := ell*j;
		if k lt #coefs then
			newterm := coefs[k];
			if j mod ell eq 0 then
				m := Numerator(j/ell);
				newterm := newterm+ell*coefs[m];
				end if;
			Append(~newcoefs, newterm);
			end if;
		end for;
	Tf := &+[newcoefs[i]*q^i : i in [1..#newcoefs]] + O(q^(#newcoefs+1));
	return Tf;
	end function;

function Uoperator(f, ell)
	exps := Exponents(f);
	if #exps eq 0 then
		init := 1;
	else
		init := exps[1];
		end if;
	coefs := Coefficients(f);
	Uf := &+[coefs[i]*q^(ell*i) : i in [1..#coefs]] + O(q^(ell*#coefs+1));
	return Uf;
	end function;

// The last two functions are for analyzing q-expansions
function traceSets(f)
	exps := Exponents(f);
	if #exps eq 0 then
		init := 1;
	else
		init := exps[1];
		end if;
	coefs := Coefficients(f);
	trace0 := [];
	trace1 := [];
	trace2 := [];
	for i in [1..#coefs] do
		j := i+init-1; // Adjusts by first term in q-expansion of f
		if IsPrime(j) then
			if coefs[i] eq 0 then
				Append(~trace0, j);
			elif coefs[i] eq 1 then
				Append(~trace1, j);
			elif coefs[i] eq 2 then
				Append(~trace2, j);
				end if;
			end if;
		end for;
	return <trace0, trace1, trace2>;
	end function;

function traceDist(f)
	T := traceSets(f);
	tot := &+[#set : set in T];
	RR := RealField(5);
	return <#T[1], #T[2], #T[3], RR!(#T[1]/tot), RR!(#T[2]/tot), RR!(#T[3]/tot)>;
	end function;

// END FUNCTIONS
