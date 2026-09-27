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


(* A single supporting bisector has three real coordinates (t,x,y).
   These chart coordinates are auxiliary and are used only in this file.
   N(c) is affine in c, and N(0)>=0. Thus N(c)<=0,0<c<1 implies N(1)<=0.
   Testing N(1)<=0 proves a stronger, parameter-free implication. *)
z=x+I Sqrt[3]y;
phaseValue=1-zeta[Sqrt[3]t];
charts={
 {"B",{c+(1-c)phaseValue/acoef,
 (1-c)(c+(z-c phaseValue)/acoef)},
 {bm[0],bp[1],fp[0]},{bm[1],fm[0]}},
 {"F",{c+(1-c)z/acoef,
 (1-c)(c+phaseValue/Conjugate[acoef]-(c+omega)z/acoef)},
 {bp[0],bm[1],fm[0]},{bm[0],bp[1]}}
};
Do[
 uv=item[[2]];
 nn=Factor[Together[norm[uv]/(1-c)^2]];
 identity[item[[1]]<>" supporting equation",side[uv,If[item[[1]]=="B",bp[0],fp[0]]]==0,
 domain&&Element[{t,x,y},Reals]];
 identity[item[[1]]<>" norm is affine in c",nn==(1-c)(nn/.c->0)+c(nn/.c->1),Element[{c,t,x,y},Reals]];
 identity[item[[1]]<>" norm at c=0 is a sum of squares",
 (nn/.c->0)==(sq[uv[[1]]]+sq[uv[[2]]]/.c->0),Element[{t,x,y},Reals]];
 exprs=Join[{nn/.c->1},(Factor[Together[side[uv,#]/(1-c)^2]]&/@Join[item[[3]],item[[4]]])];
 Do[falseTest[item[[1]]<>" positive chart denominator",Denominator[Together[expr]]<=0,{t,x,y}],{expr,exprs}];
 nums=Factor[Numerator[Together[#]]]&/@exprs;
 identity[item[[1]]<>" local inequalities independent of c",D[nums,c]==ConstantArray[0,Length[nums]],Element[{c,t,x,y},Reals]];
 Do[falseTest[item[[1]]<>" redundant local side "<>ToString[j],
 nums[[1]]<=0&&And@@Thread[nums[[2;;4]]>=0]&&nums[[4+j]]<0,{x,y,t}],{j,2}];
,{item,charts}];

(* The two other F competitors are separated on the appropriate retained
   component by a proper Giraud disk; see the connectedness argument in text. *)
Do[
 uv=phase[fp[0],item[[2]]];
 np=positiveNumerator[norm[uv],"F separation"];
 hp=positiveNumerator[side[uv,item[[3]]],"F separation"];
 falseTest[item[[1]],domain&&np<=0&&hp>=0,{c,t,h}];
,{item,{{"F1- lies strictly beyond B1-",fm[1],bm[1]},
 {"F-1- lies strictly beyond B0+",fm[-1],bp[0]}}}];
identity["positive remote F witness at x1",side[{c,c(1-c)},fm[1]]==6c(1+2c)(1-c)^2,domain];
identity["positive opposite remote F witness at x1",side[{c,c(1-c)},fm[-1]]==6c(1+2c)(1-c)^2,domain];
Print["ALL SIDE CHECKS PASSED (",passCount,"). Elapsed seconds: ",Round[AbsoluteTime[]-startTime]];
