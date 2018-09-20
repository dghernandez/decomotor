(* ::Package:: *)

(*from main_syn_v50_FIG1*)

idname = "test"; (*name for set of files*)
id=601;(*number for files*)
(*id=501, 503, 505, usual distr with two peaks*)
(*id=502, 504, 506,  normal distr with sig=0.5*)

(*LOADING FUNCTIONS*)
SetDirectory["~/Documents/decomotor-master/main/"];
(*SetDirectory["./"];*)

dirscr = "../scripts/";

<<uBIAmotor`;

(*Dumping RESULTS somewhere*)
exitfile = 1; (*set to 1 for output*)
dirres = "../resdump/res_synth/";

(*some functions and parameters for distributions*)
muJ2 = 0.5;
sigJ2 = 0.1;
muJ1 = 0.0;
sigJ1 = 0.5;
muh = -0.7;
sigh = 0.1;
fe[x_] := 1. If[x > 0, -1. x Log2[x], 0.];
Hp[p_] := Sum[fe[p[[k]]], {k, Length[p]}];

(*PARAMETERS*)
T = 18;
hcutoff=0.01;
Nmax=500;
choicedist=2;(*choosing distribution for interactions from 2 gauss(=2) or 1 gauss (=1)*)
emax=0.1;
ne = 20;
nmix = 3(*100*);

(*empirical scaling of m(nf,m(1)), to run less reshufflings*)
mc[x_,m1_:Exp[-2.67]]:=m1 x^(-0.23)Exp[-0.015 (Log[x]^2)];(*see mlevel_nf.nb*)
nf=Table[0.005 Sqrt[2]^(i-1),{i,27}];(*new threshold in nf*)

(*To calculate probabilities of all posible states, given interactions*)
e0 = Tuples[{0, 1}, T];
pJs[T_, h_, Jh_] := 
  Map[Exp, e0.h + 
    Map[Sum[Jh[[k, 2]] (Times @@ #[[Jh[[k, 1]]]]), {k, 1, 
        Length[Jh]}] &, e0]];

(*RUNNING*)
(*Parallel*)Do[
 
Ncon2=Round[1. alfa T /9];

dist2 = MixtureDistribution[{1, 1}, {NormalDistribution[-muJ2, sigJ2], NormalDistribution[muJ2, sigJ2]}];
dist1 = NormalDistribution[muJ1,sigJ1];
hh1 = RandomVariate[NormalDistribution[muh, sigh], T];

Jhl=RandomSample[Subsets[Range[T],{2}],Ncon2]
	~Join~RandomSample[Subsets[Range[T],{3}],Ncon2]
	~Join~RandomSample[Subsets[Range[T],{4}],Ncon2];
Jh = If[choicedist==2,
  Table[{Jhl[[k]], RandomVariate[dist2]}, {k, Length[Jhl]}],
  Table[{Jhl[[k]], RandomVariate[dist1]},{k, Length[Jhl]}]];

(*generating probability,and calculating entropy per spin*)
pp = pJs[T, hh1,Jh];(*//AbsoluteTiming*)
pp = pp/Total[pp];
sn = 1. Hp[pp]/T;
Mtot = 1600;
(*sampling*)
dat0 = RandomChoice[pp -> e0, Mtot];
 
 Do[(*over different M*)
  mat = dat0[[1 ;; M]];
  (*--------------------------------------*)
  
  p = 1. Total[mat]/M;
  (*h1 = h0goodpatt[mat,hcutoff,Nmax];*)(*{h,n,pM,z,+-1,patt}*)
  h1 = h0goodpatt2[mat,hcutoff,Nmax];(*{h,n,pM,z,+-1,patt}*)
  Nm1 = Length[h1];
  (*Print[{alfa,arch,M,Nm1}];*)
  {hef, Jvw} = couplings[h1, p, M];(*//AbsoluteTiming*)(*h from couplings and couplings Jvw*)(*Bottleneck ~25min *)
  {mi, h1ef} = margmagn[M, h1, hef, Jvw,emax,ne];(*marginal magnetizetion and e-correction to fields*)
  
  iopt = epsilonopt[h1[[All, 1]], h1ef, M,emax,ne];(*limit to taylor serie Mean|h0|~Mean|h(e)|*)
  
  (*mlevel = signlevellower2[mat,iopt,nmix,ne,hcutoff,Nmax];*)(*//AbsoluteTiming*)
  res = Sort[Join[h1, mi[[All, {iopt}]], h1ef[[All, {iopt}]], mi, h1ef, 2], #1[[7]] > #2[[7]] &];
  (*------Significance test------*)
  mlevel1 = signlevel[mat,iopt,nmix,ne,hcutoff,Nmax,emax];
  mlevel=mc[#,mlevel1[[1]] ]& /@ nf;
  
  win = Select[res, #[[7]] > mlevel[[27]] &];
  winx = Table[Select[win, #[[7]] > mlevel[[j]] &][[All, 6]], {j, 27}];
  tsm = Map[truesimdetected[Jhl, #] &, winx];
  
 (*tsm[[13]]: nf=0.32, tsm[[17]]: nf=1.28*)
 (*tsm[[i, {2,3}]]: {#det/#orig, #det/#guess}*)
  Print[{{alfa, arch, M, Nm1}}~Join~{tsm[[13,;;3]]}];

  If[exitfile==1,  
  mydata = {alfa, arch, M, muJ2, sigJ2, muh, sigh, sn, Nm1}~Join~Flatten[Table[tsm[[i, 2;;3(*3 ;; 6*)]],{i,27}]];
  
  $stream = OpenAppend[dirres <> idname <> ToString[T] <>"_Nmax"<>ToString[Nmax]<>"mJ05_mlevel_id"<>ToString[id]<>".dat", BinaryFormat -> True]; 
  Export[$stream, {mydata}, "Table"]; 
  WriteString[$stream, "\n"]; Close[$stream];
  ];(*end if exitfile*)

  , {M,{(*200,283,400,566,800,1131,1600*)1000}}]
 , {arch,Range[10(*400*)]},{alfa,{2(*,4*)}}];//AbsoluteTiming




