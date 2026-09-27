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


Do[
 uv=phase[bp[0],item[[2]]];
 np=positiveNumerator[norm[uv],"boundary cells"];
 other=Select[cb,!SameQ[#,item[[2]]]&];
 hp=positiveNumerator[side[uv,#],"boundary cells"]&/@other;
 falseTest[item[[1]]<>" no extra boundary inside triangle",
 domain&&np<=0&&item[[3]]&&Or@@Thread[hp<=0],{c,t,h}];
,{item,{{"BB",bp[1],t>0&&h>0},{"BF",fp[0],t>0&&h<t}}}];
Do[
 uv=phase[item[[2]],item[[3]]];
 np=positiveNumerator[norm[uv],"lower dimensional cells"];
 hp=positiveNumerator[side[uv,#],"lower dimensional cells"]&/@item[[4]];
 falseTest[item[[1]],domain&&np<=0&&And@@Thread[hp>=0]&&item[[5]],{c,t,h}];
,{item,{
 {"mixed B pair lies on an adjacent mirror",bp[0],bm[1],cb,
   side[phase[bp[0],bm[1]],bm[0]]!=0&&side[phase[bp[0],bm[1]],bp[1]]!=0},
 {"opposite F pair only at polygon vertex",fp[0],fm[0],cf,t!=0||h!=0},
 {"adjacent opposite F pair is empty",fp[0],fm[1],{bm[1]},True},
 {"B0 F0- lies on adjacent mirror",bp[0],fm[0],cb,
   side[phase[bp[0],fm[0]],bm[0]]!=0&&side[phase[bp[0],fm[0]],bp[1]]!=0},
 {"F0 B0- lies on a listed boundary",fp[0],bm[0],cf,
   side[phase[fp[0],bm[0]],bp[0]]!=0&&side[phase[fp[0],bm[0]],bm[1]]!=0},
 {"F0 B1+ lies on a listed boundary",fp[0],bp[1],cf,
   side[phase[fp[0],bp[1]],bp[0]]!=0&&side[phase[fp[0],bp[1]],bm[1]]!=0}
}}];
Print["ALL EXTRA CHECKS PASSED (",passCount,")."];
