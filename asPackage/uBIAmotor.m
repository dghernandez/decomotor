(* ::Package:: *)

BeginPackage["uBIAmotor`"];

Unprotect["uBIAmotor`*"];
ClearAll["uBIAmotor`*", "uBIAmotor`Private`*"];

(******************************************************************************)

ftab4::usage = "ftab4[numberfile, vocal_feature, hl:1] extracts a matrix from files in list_files.txt, corresponding to spikes and behavior.";
mayzero::usage = "mayzero[tab] inverts 1->0 (from column 2 to last) if 1 is more likely.";
h0goodpatt2::usage = "h0goodpatt2[tab, hcutoff: 0.01, Nmax: 1000, sfiltroh: 1], function that taking data matrix, extracts all patterns to be analyzed with their occurrences, under some specific criteria and return {h0,nn,ppM,zz,+-1,patt}.";
couplings::usage = "couplings[h1, p] returns the couplings Jvw and h extra coming from couplings (sum_w Jvw)." ;
margmagn::usage = "margmagn[M, h1, hef, Jvw, emax: 1., ne: 20] .mi: marginal magnetization as a function of epsilon; h1ef: epsilon correction to h0 as a function of epsilon (to later decide eopt as Mean[|h0|]~Mean[|h1ef|]).";
epsilonopt::usage = "epsilonopt[h0v, h1ef, M, emax: 1, ne: 20] finds the limit(column of mi) of the taylor approximation in epsilon. Criteria is Mean|h0|~Mean|h(epsilon)|.";
signlevellower4::usage = "signlevellower4[mat, iopt, nmix: 100, ne: 20, hcutoff: 1, Nmax: 10000, emax: 1.] gets significance level from reshuffling data. Needs previous functions: h0goopatt, couplings, margmagn.";
truesimdetected::usage = "truesimdetected[orig_, guess_] counts the number of true patt detected.";



(******************************************************************************)

Begin["`Private`"];

(******************************************************************************)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)


ftab4[a_, pae_: 0, hl_: 1, dir_: "ownCloud/EMORY/bia_j/CODE/dataset/"] := Block[{tab, T0, T, lis, v1, st, dat, sdat, pi, mpi, bpi, M, 
   q10, q20, q30, q40, q50, q60, q70, q80, q90, q35, q65, q33, q66},
   T0 = 40; T = 21;
   lis = Import[dir <> "list_files.txt", "Table"] // Flatten;
   dat = Quiet[Import[dir <> ToString[lis[[a]]], "LabeledData"]];
   M = ("number_of_trials" /. dat)[[1, 1]] // Floor;
   sdat = "spiketrains" /. dat;
   (*sdat=sdat /. Transpose[Partition[{},0.`]]\[Rule] {{}};*)
   
   tab = ConstantArray[0, {M, T - 1}]; st = {};
   Quiet[Do[
      If[sdat[[3, 1, i, 1]] != {},
        
        v1 = Ceiling[sdat[[3, 1, i, 1, 1]]*(T - 1)/T0];
        st = st~Join~{v1};
        Do[
         tab[[i, v1[[j]] ]] = 1;
         , {j, Length[v1]}]
        
        ];
      , {i, M}];];(*//so far, extraction of spike activity*)
   
   (*pae: 0 pitch, 1 amplitude, 2 entropy*)(*hl: 1=1 high, 0=0 high*)
   
   (*for pae=-1, we don't include vocal feature, just marginal activity*)
   If[pae != -1,
    
    Which[pae == 0 || pae == 1 || pae == 2,(*median, q50*)     
      pi = ("acoustic_data" /. dat)[[pae + 1, 1, 1, 1, 1]]; mpi = Median[pi];
      bpi = Table[If[pi[[i]] > mpi, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
 
      , pae == 10 || pae == 11 || pae == 12,(*outside [q15,q85]*)
        pi = ("acoustic_data" /. dat)[[pae - 10 + 1, 1, 1, 1, 1]];
      {q20, q40, q50, q60, q80} = Quantile[pi, {0.20, 0.40, 0.5, 0.60, 0.80}];
      bpi = Table[If[pi[[i]] < q20 || pi[[i]] > q80, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 20 || pae == 21 || pae == 22,(*above q60*)
      pi = ("acoustic_data" /. dat)[[pae - 20 + 1, 1, 1, 1, 1]];
      {q20, q40, q50, q60, q80} = Quantile[pi, {0.20, 0.40, 0.5, 0.60, 0.80}];
      bpi = Table[If[pi[[i]] > q60, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 30 || pae == 31 || pae == 32,(*below q40*)
      pi = ("acoustic_data" /. dat)[[pae - 30 + 1, 1, 1, 1, 1]];
      {q20, q40, q50, q60, q80} = Quantile[pi, {0.20, 0.40, 0.5, 0.60, 0.80}];
      bpi = Table[If[pi[[i]] < q40, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 40 || pae == 41 || pae == 42,(*q0-q20*)
      pi = ("acoustic_data" /. dat)[[pae - 40 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] < q20, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 50 || pae == 51 || pae == 52,(*q10-q30*)
      pi = ("acoustic_data" /. dat)[[pae - 50 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q10 && pi[[i]] < q30, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 60 || pae == 61 || pae == 62,(*q20-q40*)
      pi = ("acoustic_data" /. dat)[[pae - 60 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q20 && pi[[i]] < q40, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 70 || pae == 71 || pae == 72,(*q30-q50*)
      pi = ("acoustic_data" /. dat)[[pae - 70 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q30 && pi[[i]] < q50, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 80 || pae == 81 || pae == 82,(*q40-q60*)
      pi = ("acoustic_data" /. dat)[[pae - 80 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q40 && pi[[i]] < q60, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 90 || pae == 91 || pae == 92,(*q50-q70*)
      pi = ("acoustic_data" /. dat)[[pae - 90 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q50 && pi[[i]] < q70, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 100 || pae == 101 || pae == 102,(*q60-q80*)
      pi = ("acoustic_data" /. dat)[[pae - 100 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q60 && pi[[i]] < q80, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 110 || pae == 111 || pae == 112,(*q70-q90*)
      pi = ("acoustic_data" /. dat)[[pae - 110 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q70 && pi[[i]] < q90, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 120 || pae == 121 || pae == 122,(*q80-q100*)
      pi = ("acoustic_data" /. dat)[[pae - 120 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q80, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 130 || pae == 131 || pae == 132,(*q0-q30*)
      pi = ("acoustic_data" /. dat)[[pae - 130 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] < q30, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 140 || pae == 141 || pae == 142,(*q35-q65*)
      pi = ("acoustic_data" /. dat)[[pae - 140 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q35, q40, q50, q60, q65, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q35 && pi[[i]] < q65, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 150 || pae == 151 || pae == 152,(*q70-q100*)
      pi = ("acoustic_data" /. dat)[[pae - 150 + 1, 1, 1, 1, 1]];
      {q10, q20, q30, q40, q50, q60, q70, q80, q90} = Quantile[pi, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}];
      bpi = Table[If[pi[[i]] > q70, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 160 || pae == 161 || pae == 162,(*q0-q33*)
      pi = ("acoustic_data" /. dat)[[pae - 160 + 1, 1, 1, 1, 1]];
      {q33, q66} = Quantile[pi, {0.3333, 0.6667}];
      bpi = Table[If[pi[[i]] < q33, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 170 || pae == 171 || pae == 172,(*q33-q66*)
      pi = ("acoustic_data" /. dat)[[pae - 170 + 1, 1, 1, 1, 1]];
      {q33, q66} = Quantile[pi, {0.3333, 0.6667}];
      bpi = Table[If[pi[[i]] > q33 && pi[[i]] < q66, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2];
      
      , pae == 180 || pae == 181 || pae == 182,(*q66-q100*)
      pi = ("acoustic_data" /. dat)[[pae - 180 + 1, 1, 1, 1, 1]];
      {q33, q66} = Quantile[pi, {0.3333, 0.6667}];
      bpi = Table[If[pi[[i]] > q66, hl, 1 - hl], {i, Length[pi]}];(*hl: 1=1 high, 0=0 high*)
      tab = Join[Transpose[{bpi}], tab, 2]
      (*end which*)];
    
    (*end if*)];
   
   ; tab];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

mayzero[tab_] := Block[{tout = tab[[All, 2 ;;]], pp, MM = Length@tab, TT},
  pp = 1. Total@tout/MM;
  TT = Length@pp;
  Do[If[pp[[j]] > 0.5, 
    tout[[All, j]] = ConstantArray[1, MM] - tout[[All, j]]], {j, 
    TT}];
  Join[Transpose@{tab[[All, 1]]}, tout, 2]
  ];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*function that taking data matrix, extract all patterns to be \
analyzed with their occurrences, under some specific criteria \
|h0|>0.01, and return {h0,nn,ppM,zz,+-1,patt} *)

h0goodpatt2[tab_, hcutoff_: 0.01, Nmax_: 1000, sfiltroh_: 1] := 
  Block[{res, sm, sets, tset, cc, M, T, omax, p, cc0, tset1},
   M = Length @tab;
   T = Dimensions[tab][[2]];
   sm = Max[Total[Transpose[tab]]];
   sets = {};
   p = 1. Total[tab]/M;
   
   (*finding patterns*)
   Do[If[tab[[j]] != {}, sets = sets~Join~Subsets[Flatten[Position[tab[[j]], 1]], {1, sm}];];
    , {j, 1, Length[tab]}];
   tset = Tally[sets];
   
   cc = hcutoff 2./(M + 1);
   omax = 3; 
   While[(Max @Map[Times @@ p[[#]] &, Subsets[Range[T], {omax}]]) > cc, omax++];
   omax = omax - 1;
   
   cc0 = Complement[Select[Subsets[Range[T], {1, omax}], (Times @@ p[[#]]) > cc &], tset[[All, 1]] ];
   tset1 = Join[tset, Transpose[{cc0, Table[0, {i, Length[cc0]}]}]];
   
   (*calculating h0*)
   
   res = Table[Block[{hh = 0., pp = 0., zz = 0., V, nv = 0},
      V = tset1[[i, 1]];
      nv = tset1[[i, 2]];
      pp = Times @@ p[[V]];
      zz = (nv - M pp)/Sqrt[M pp (1 - pp)];
      hh = M^2 ((nv/M - pp)^2 - 1./M pp (1 - pp))/4;
      {hh, nv, pp M, zz, If[nv > pp M, 1, -1], tset1[[i, 1]]}], {i, 
      Length[tset1]}];
   
   If[sfiltroh != 1, res = Select[res, Abs[#[[1]]] > 1. hcutoff &]; ];
   
   res = Sort[res, Abs[#1[[1]]] > Abs[#2[[1]]] &];
   (*res=Sort[res,#1[[1]]>#2[[1]]&];*)(* esta parece ser mejor, en combinacion con Nmax~10^3 *)
   
   If[Length[res] > Nmax, res = res[[1 ;; Nmax]]];; 
   res];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*returns h extra coming from couplings (sum_w Jvw) and the couplings Jvw*)
couplings[h1_, p_] := Block[{resJ, reshef, tv, tw, Iv, Uv, hh = 0., j1, Nm1, pv, pw, pu, pin},
   Nm1 = Length @h1;
   reshef = ConstantArray[0, Nm1];
   j1 = Reap[Do[
       	tv = h1[[v, 6]];
       	tw = h1[[w, 6]];
       	Iv = Intersection[tv, tw];
       
       	If[Iv != {},
        		pv = Times @@ p[[tv]];
        		pw = Times @@ p[[tw]];
        		Uv = Union[tv, tw]; pu = Times @@ p[[Uv]];
        		pin = Times @@ p[[Iv]];
        		
                hh = (1/16.) pu (1. - pin) (M^2 pu (1. - pin) - 2. M (h1[[v, 2]] - M pv) (h1[[w, 2]] - M pw));
               (*meti el M^2 dentro del parentesis*)
        		
                reshef[[v]] += hh;
        		reshef[[w]] += hh;
        		Sow[{v, w} -> hh]
        	];
       
       , {v, 1, Nm1 - 1}, {w, v + 1, Nm1}]][[2, 1]];
   
   resJ = SparseArray[j1~Join~{{Nm1, Nm1} -> 0}];
   ; {reshef, resJ}];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*mi: marginal magnetization as a function of epsilon;
h1ef: epsilon correction to h0 as a function of epsilon
(to later decide eopt as Mean[|h0|]~Mean[|h1ef|])*)

margmagn[M_, h1_, hef_, Jvw_, emax_: 1., ne_: 20] := 
 Block[{mi, h1ef, de, Nm1, e, JvwT},
  JvwT = Transpose[Jvw];
  
  Nm1 = Length @h1;
  de = emax/(M ne);
  mi = ConstantArray[0., {Nm1, ne + 1}];
  h1ef = ConstantArray[0., {Nm1, ne + 1}]; 
  Do[h1ef[[v, 1]] = h1[[v, 1]], {v, Nm1}];
  
  Do[
   e = de (t - 1.);
   Do[
    
    If[t == 1,
      	h1ef[[v, t + 1]] = h1[[v, 1]] + 2. (e + de) hef[[v]];
      	mi[[v, t + 1]] = Tanh[(e + de) h1ef[[v, t + 1]](*(h1[[v,1]]+2. (e+de)hef[[v]])*)];
      	,
      	h1ef[[v, t + 1]] = h1[[v, 1]] + 2. (e + de) hef[[v]] + 1. (e + de) (Jvw[[v]].mi[[All, t]] + JvwT[[v]].mi[[All, t]]);
      	mi[[v, t + 1]] = Tanh[(e + de) h1ef[[v, t + 1]]];
      ];
    , {v, 1, Nm1}];
   (*Print[t]*)
   , {t, 1, ne}];
  
  ; {mi, h1ef}];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*finds the limit(column of mi) of the taylor approximation in epsilon.
Criteria is Mean|h0|~Mean|h(epsilon)|*)
(*--h1[[All,1]]=h0v--*)

epsilonopt[h0v_, h1ef_, M_, emax_: 1, ne_: 20] := 
  Block[{iopt = 1, zero, taylor},
   zero = Mean @Map[Sqrt[#^2] &, h0v];
   taylor = Table[{emax (i - 1.)/ne, (Mean @Map[Sqrt[#^2] &, h1ef[[All, i]] - h0v])/zero}, {i, 1, ne + 1}];
   iopt = 1; While[taylor[[iopt, 2]] < 1 && iopt <= ne, iopt++];
   iopt = iopt - 1; 
   iopt];
   
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*needs previous functions: h0goopatt, couplings, margmagn*)
(*gets significance level from reshuffling data*)

signlevellower4[mat_, iopt_, nmix_: 100, ne_: 20, hcutoff_: 1, Nmax_: 10000, emax_: 1.] := 
 Block[{mlist, mc002 = 0., mc004 = 0., mc008 = 0., mc016 = 0., 
   mc032 = 0., mc064 = 0., mc128 = 0., mc256 = 0., mc512 = 0., 
   mc1024 = 0., mc2048 = 0., mc4096 = 0., M, T, p, matr, h1r, Nm1r, 
   hefr, Jvwr, mir, h1efr},
  T = Dimensions[mat][[2]];
  M = Length[mat];
  p = 1. Total[mat]/M;
  mlist = {};
  
  Do[
   matr = Transpose[Table[RandomSample[mat[[All, j]]], {j, T}]];
   (*h1r=h0goodpatt[matr,hcutoff,Nmax];*)(* changed on sep06 *)
   
   h1r = h0goodpatt2[matr, hcutoff, Nmax];(*{h,n,pM,z,+-1,patt}*)
   
   Nm1r = Length[h1r];
   {hefr, Jvwr} = couplings[h1r, p];(*h from couplings and couplings Jvw*)(* Bottleneck ~ 25min *)
   {mir, h1efr} = margmagn[M, h1r, hefr, Jvwr, emax];(*marginal magnetizetion and e-correction to fields*)
   
   mlist = mlist~Join~Sort[mir[[All, iopt]], Greater][[;; Min[nmix, Length@mir]]];
   , {ss, nmix}];
    
  (*mc002=RankedMax[mlist,Min[Length@mlist,Floor[0.02*nmix]]];
  mc004=RankedMax[mlist,Min[Length@mlist,Floor[0.04*nmix]]];
  mc008=RankedMax[mlist,Min[Length@mlist,Floor[0.08*nmix]]];*)
  
  mc016 = RankedMax[mlist, Min[Length@mlist, Floor[0.16* nmix]]];
  mc032 = RankedMax[mlist, Min[Length@mlist, Floor[0.32* nmix]]];
  mc064 = RankedMax[mlist, Min[Length@mlist, Floor[0.64* nmix]]];
  
  mc128 = RankedMax[mlist, Min[Length@mlist, Floor[1.28* nmix]]];
  mc256 = RankedMax[mlist, Min[Length@mlist, Floor[2.56* nmix]]];
  mc512 = RankedMax[mlist, Min[Length@mlist, Floor[5.12* nmix]]];
  (*mc1024=RankedMax[mlist,Min[Length@mlist,Floor[10.24* nmix]]];
  mc2048=RankedMax[mlist,Min[Length@mlist,Floor[20.48* nmix]]];
  mc4096=RankedMax[mlist,Min[Length@mlist,Floor[40.96* nmix]]];
  *)
  {(*mc002,mc004,mc008,*)mc016, mc032, mc064, mc128, mc256, mc512(*,mc1024,mc2048,mc4096*)}];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*needs previous functions: h0goopatt, couplings, margmagn*)
(*gets \
significance level from reshuffling data*)

signlevel[mat_, iopt_, nmix_: 5, ne_: 20, hcutoff_: 1, Nmax_: 10000, 
   emax_: 1.] := 
  Block[{mc1 = 0., mc2 = 0., mc3 = 0., mc5 = 0., mc10 = 0., M, T, p, 
    matr, h1r, Nm1r, hefr, Jvwr, mir, h1efr},
   T = Dimensions[mat][[2]];
   M = Length[mat];
   p = 1. Total[mat]/M;
   
   Do[
    matr = Transpose[Table[RandomSample[mat[[All, j]]], {j, T}]];
    h1r = h0goodpatt2[matr, hcutoff, Nmax];(*{h,n,pM,z,+-1,patt}*)
   
     Nm1r = Length[h1r];
    {hefr, Jvwr} = 
     couplings[h1r, p];(*h from couplings and couplings Jvw*)(* 
    Bottleneck ~
    25min *)
    {mir, h1efr} = 
     margmagn[M, h1r, hefr, Jvwr, emax, 
      ne];(*marginal magnetizetion and e-correction to fields*)
    
    mc1 += RankedMax[mir[[All, iopt]], 1];
    mc2 += RankedMax[mir[[All, iopt]], 2];
    mc3 += RankedMax[mir[[All, iopt]], 3];
    mc5 += RankedMax[mir[[All, iopt]], 5];
    mc10 += RankedMax[mir[[All, iopt]], 10];
    (*mc1p+=RankedMax[mir[[All,iopt]],Quotient[Nm1r,
    100]];*)
    , {ss, nmix}];
   
   mc1 = mc1/nmix;
   mc2 = mc2/nmix;
   mc3 = mc3/nmix;
   mc5 = mc5/nmix;
   mc10 = mc10/nmix;
   (*mc1p=mc1p/nmix;*)
   
   {mc1, mc2, mc3, mc5, mc10}];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)

(*Counting the number of true patt detected*)

truesimdetected[orig_, guess_] := 
  Block[{res, td = 0, no = Length[orig], nd = Length[guess], aux, 
    ptruelist = {}, at, as},
   If[nd == 0,
    res = {{0, 0, no}, 0, 1, {}},
    
    Map[If[MemberQ[orig, #], td += 1; ptruelist = ptruelist~Join~{#}] &, guess];
    res = {{td, nd, no}, 1. td/no, 1. td/nd, Sort @(Tally @ (Length /@ ptruelist))};
    ];
   res];

(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(*kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk*)
(******************************************************************************)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["uBIAmotor`*"], Head[#] === Symbol &]];

End[];
EndPackage[];

