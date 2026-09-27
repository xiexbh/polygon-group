(* Exact verification accompanying Section 6.
   Run: math -script section6_verification.wl
   Or evaluate Get["/absolute/path/section6_verification.wl"] in Mathematica.
   No optional packages, floating-point inputs, sampling, or external files.
   Every FAILED or TIMEOUT result is fatal; success ends with ALL CHECKS PASSED.
   Tested with Mathematica 14.0. Increase $checkTimeLimit if necessary.
*)
ClearAll["Global`*"];
$checkTimeLimit = 240;
passCount = 0;
startTime = AbsoluteTime[];
Print["Section 6 exact verification; ", $Version];
Print["Real parameter domain: 3/5 <= c < 1; t,h real."];

pass[tag_] := (passCount++; Print["PASS  ", tag]);
fail[tag_, result_] := (Print["FAILED  ", tag, ": ", InputForm[result]]; Exit[1]);
identity[tag_, assertion_, assumptions_:True] := Module[{answer},
  answer = TimeConstrained[FullSimplify[assertion, assumptions], $checkTimeLimit, $TimedOut];
  If[TrueQ[answer], pass[tag], fail[tag,answer]]];
falseTest[tag_, condition_, variables_List] := Module[{answer,started},
  started=AbsoluteTime[];
  answer = TimeConstrained[Reduce[condition, variables, Reals], $checkTimeLimit, $TimedOut];
  If[SameQ[answer,False], pass[tag], fail[tag,answer]]];
domain = 3/5 <= c < 1;
realDomain = domain && Element[{t,h},Reals];
omega=(1+I Sqrt[3])/2;
aa=Sqrt[3c](Sqrt[3]+I)/(2(1-c));
bb=(1+c omega)/(1-c);
dd=(1+2c)/(1-c);
ss=Sqrt[1-c^2];
bar[z_]:=ComplexExpand[Conjugate[z]];
sq[z_]:=Together[ComplexExpand[Re[z]]^2+ComplexExpand[Im[z]]^2];

(* A: normalization and analytic identities used in the disjointness proofs. *)
identity["cos(2 Pi/7)>3/5", Cos[2Pi/7]>3/5];
identity["B normalization", sq[bb]-sq[aa]==1,domain];
identity["F normalization", dd^2-sq[aa](2+c)==1,domain];
imat={{-bb Conjugate[omega],0,bar[aa]},{0,1,0},{-bar[aa],0,bar[bb]}};
jmat=DiagonalMatrix[{1,1,-1}];
rotmat={{c,-ss,0},{ss,c,0},{0,0,1}};
imat1=rotmat.imat.Transpose[rotmat];
identity["I0 preserves the Hermitian form",bar[Transpose[imat]].jmat.imat==jmat,domain];
identity["I0 has order three",MatrixPower[imat,3]==IdentityMatrix[3],domain];
identity["adjacent reflections commute",imat.imat1==imat1.imat,domain];
jreflection=DiagonalMatrix[{1,-1,1}];
sigmaLinear=rotmat.jreflection;
identity["conjugation reverses the reflection",bar[imat]==MatrixPower[imat,2],domain];
identity["J reverses the mirror index",jreflection.imat1.jreflection==Transpose[rotmat].imat.rotmat,domain];
identity["sigma0 sends I0 to inverse I1",sigmaLinear.bar[imat].sigmaLinear==MatrixPower[imat1,2],domain];
identity["sigma0 sends I1 to inverse I0",sigmaLinear.bar[imat1].sigmaLinear==MatrixPower[imat,2],domain];

qDomain=domain && q>0 && q^2==c^2+c+1;
falseTest["BB scalar inequalities",qDomain &&
 (2(q-1+c)/(q+1-c)-1<=0 || 2(q-1+c)/(q+1-c)-1<=2c^2-1),{c,q}];
ee=4c^2+5c+3-2q;
falseTest["BF radical inequality: positive factor",qDomain && ee<=0,{c,q}];
falseTest["BF radical inequality: squared comparison",qDomain &&
 ee^2(q+1-c)-16q^2<=0,{c,q}];
identity["rational upper bound for Q",
 (c+4/5)^2-(c^2+c+1)==3(c-3/5)/5];
identity["rational lower bound for Q",
 c^2+c+1-(c+2/3)^2==5/9-c/3];
identity["FF++ gap",(5c-2)/(2+c)-(2c^2-1)==2c(1-c)(c+3)/(2+c)];
identity["FF+- gap",(5c-2)/(2+c)-(1+2c)(2c^2-1)/(2+c)==
 (1-c)(4c^2+6c-1)/(2+c)];
r0={c+omega,ss}; r1={2c^2-1+omega c,(2c+omega)ss};
weighted=Exp[-I Pi/6]r0+Exp[I Pi/6]r1;
identity["adjacent FF weighted coefficient norm",Total[sq /@ weighted]==(1+2c)^2,domain];
falseTest["adjacent FF strict separation",domain && Sqrt[3](1+2c)<=3+c,{c}];
(* Real scalar products; x=cos(k theta), y=sin(k theta). *)
realAngles=domain&&Element[{x,y},Reals]&&x^2+y^2==1;
qB={bar[aa],0}; qBm={aa x,aa y};
qFp=aa{c+Conjugate[omega],ss};
qFkp=aa{c x-ss y+Conjugate[omega]x,ss x+c y+Conjugate[omega]y};
qFkm=aa{x+Conjugate[omega](c x-ss y),y+Conjugate[omega](ss x+c y)};
identity["FF++ normal scalar product",ComplexExpand[Re[qFp.bar[qFkp]]]==sq[aa](2+c)x,realAngles];
identity["FF+- normal scalar product",ComplexExpand[Re[qFp.bar[qFkm]]]==sq[aa](1+2c)x,realAngles];
identity["BB+- normal scalar product",ComplexExpand[Re[(qB bb).bar[qBm bar[bb]]]]==
 sq[aa](1+4c+c^2)x/(2(1-c)^2),realAngles];
identity["BF normal scalar product",ComplexExpand[Re[(qB bb).bar[qFkp]]]==
 sq[aa]((2+c)x+(1+2c)(c x-ss y))/(2(1-c)),realAngles];

(* B: rows in the scaled affine coordinates U=sqrt(c)u,V=sqrt(c)s v. *)
acoef=(3+I Sqrt[3])/2;
trig[k_Integer]:={ChebyshevT[Abs[k],c],If[k==0,0,Sign[k]ChebyshevU[Abs[k]-1,c]]};
row["B",k_Integer,1]:=Append[acoef trig[k],1+c omega];
row["B",k_Integer,-1]:=bar[row["B",k,1]];
row["F",k_Integer,1]:=Append[Conjugate[acoef](trig[k+1]+omega trig[k]),1+2c];
row["F",k_Integer,-1]:=Append[Conjugate[acoef](trig[k]+omega trig[k+1]),1+2c];
bp[k_]:=row["B",k,1]; bm[k_]:=row["B",k,-1];
fp[k_]:=row["F",k,1]; fm[k_]:=row["F",k,-1];
zeta[v_]:=(1+I v)/(1-I v);
phase[r_,rr_]:=Together[Inverse[{r[[1;;2]],rr[[1;;2]]}].
 {r[[3]]-(1-c)zeta[t],rr[[3]]-(1-c)zeta[h]}];
norm[uv_]:=Together[(1-c^2)(sq[uv[[1]]]-c)+sq[uv[[2]]]];
side[uv_,r_]:=Together[sq[r[[1;;2]].uv-r[[3]]]-(1-c)^2];
cb={bm[0],bp[1],bm[1],bp[-1],bm[-1],fp[0],fm[0],fp[-1],fm[-1]};
cf={bp[0],bm[0],bp[1],bm[1],fm[0],fm[1],fm[-1]};

(* Denominators are proved positive, not silently discarded. *)
denominatorCache={};
positiveNumerator[expression_,tag_] := Module[{rat,den},
 rat=Together[expression]; den=Factor[Denominator[rat]];
 If[!MemberQ[denominatorCache,den],
   falseTest["positive denominator / "<>tag,domain&&den<=0,{c,t,h}];
   AppendTo[denominatorCache,den]];
 Factor[Numerator[rat]]];

(* Nonempty disks, including all the listed pairs rather than one example. *)
witnesses=Join[({"B0 / "<>ToString[#],bp[0],cb[[#]],-c(1-c)^2}& /@Range[2,9]),
 Table[{"F0 / "<>ToString[j],fp[0],cf[[j]],
 If[j<=5,-c(1-c)^2,-c(1+c)(1-c)^3/(c^2+c+1)]},{j,7}]];
Do[
  det=Det[{item[[2,1;;2]],item[[3,1;;2]]}];
  falseTest["nonzero determinant / "<>item[[1]],domain&&sq[det]==0,{c}];
  uv=phase[item[[2]],item[[3]]];
  identity["negative witness / "<>item[[1]],(norm[uv]/.{t->0,h->0})==item[[4]],domain];
,{item,witnesses}];

(* C: simple complex slice and its circles. *)
identity["slice u=sqrt(c)",sq[aa Sqrt[c]-bb]==1,domain];
identity["second real slice root",sq[aa/Sqrt[c]-bb]==1,domain];
sliceUV={c,Sqrt[c(1-c)(1-c^2)](x+I y)};
rad=1/Sqrt[3c(1+c)]; center=(2c+1)/(2Sqrt[c(1+c)]);
identity["slice B1+ circle",side[sliceUV,bp[1]]==
 (1-c)^2 3c(1+c)((x-center)^2+(y+rad/2)^2-rad^2),
 domain&&Element[{x,y},Reals]];
identity["slice B1- circle",side[sliceUV,bm[1]]==
 (1-c)^2 3c(1+c)((x-center)^2+(y-rad/2)^2-rad^2),
 domain&&Element[{x,y},Reals]];
Do[identity["coincident slice equations / "<>ToString[j],
 side[sliceUV,pairs[[j,1]]]==side[sliceUV,pairs[[j,2]]],
 domain&&Element[{x,y},Reals]],
 {pairs,{{{bm[1],fp[0]},{bp[1],fm[0]},{bm[-1],fm[-1]},{bp[-1],fp[-1]}}}},{j,4}];
identity["orthogonal circles",center^2==1+3rad^2/4,domain];

(* D: two complete retained-region equivalences. *)
pBB=-3c-2Sqrt[3]c(t+h)+(4+3c)(t^2+h^2)-8c t h-
 2Sqrt[3]c t h(t+h)+(8+c)t^2 h^2;
pBF=-3c+2Sqrt[3]c h-4Sqrt[3]c t+(4+3c)h^2+(8+c)t^2-
 (8+4c)t h-2Sqrt[3]c h t^2+(4+3c)t^2 h^2;
Do[
  tag=If[j==1,"BB triangle","BF triangle"];
  rr=If[j==1,bp[1],fp[0]]; poly=If[j==1,pBB,pBF];
  uv=phase[bp[0],rr];
  identity[tag<>" norm formula",norm[uv]==(1-c)^2 poly/(3(1+t^2)(1+h^2)),realDomain];
  np=positiveNumerator[norm[uv],tag];
  hp=positiveNumerator[side[uv,#],tag]& /@cb;
  target=If[j==1,t>=0&&h>=0,t>=0&&h<=t];
  falseTest[tag<>" closed-ball region equivalence",
    domain&&np<=0&&Xor[And@@Thread[hp>=0],target],{c,t,h}];
,{j,2}];
identity["BB cutting line",(pBB/.t->0)==(4+3c)h^2-2Sqrt[3]c h-3c];
identity["BF first cutting line",(pBF/.t->0)==(4+3c)h^2+2Sqrt[3]c h-3c];
identity["BF diagonal cutting line",(pBF/.h->t)==(1+t^2)((4+3c)t^2-2Sqrt[3]c t-3c)];

(* E: no relatively open two-dimensional retained region for these pairs.
   Deleted zero functions are exactly the selected supporting equations.
   The opposite-index cases follow from the explicit isometries in the text.
*)
cases={{"B0+ B1-",bp[0],bm[1],cb},{"B0+ F0-",bp[0],fm[0],cb},
 {"F0+ B0-",fp[0],bm[0],cf},{"F0+ B1+",fp[0],bp[1],cf},
 {"F0+ F0-",fp[0],fm[0],cf},{"F0+ F1-",fp[0],fm[1],cf}};
Do[
 Print["START no-ridge test / ",item[[1]]];
 uv=phase[item[[2]],item[[3]]]; np=positiveNumerator[norm[uv],item[[1]]];
 others=Select[item[[4]],!SameQ[#,item[[2]]]&&!SameQ[#,item[[3]]]&];
 hp=positiveNumerator[side[uv,#],item[[1]]]& /@others;
 If[MemberQ[hp,0],fail["an unselected bisector vanishes identically / "<>item[[1]],hp]];
 (* Prove nontrivial restriction for every c, not only as a joint polynomial. *)
 Do[coeffs=Values[CoefficientRules[Expand[polynomial],{t,h}]];
    falseTest["nontrivial restriction / "<>item[[1]],domain&&And@@Thread[coeffs==0],{c}],
   {polynomial,hp}];
 falseTest["closed-ball strict test / "<>item[[1]],
   domain&&np<=0&&And@@Thread[hp>0],{c,t,h}];
,{item,cases}];
uv=phase[bp[0],bm[1]];
identity["mixed-pair product factor",side[uv,fp[0]]==(1-c)^2 8t h(1+t h)/((1+t^2)(1+h^2)),realDomain];
falseTest["mixed-pair strict B sides force t>0,h<0",
 domain&&positiveNumerator[norm[uv],"mixed signs"]<0&&
 positiveNumerator[side[uv,bm[0]],"mixed signs"]>0&&
 positiveNumerator[side[uv,bp[1]],"mixed signs"]>0&&Not[t>0&&h<0],{c,t,h}];
falseTest["mixed-pair 1+th is positive",
 domain&&positiveNumerator[norm[uv],"mixed factor"]<0&&t>0&&h<0&&1+t h<=0,{c,t,h}];
Print["ALL CHECKS PASSED (",passCount,"). Elapsed seconds: ",Round[AbsoluteTime[]-startTime]];
