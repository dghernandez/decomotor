(* ::Package:: *)

(*LOADING FUNCTIONS*)
SetDirectory["~/ownCloud/EMORY/bia_j/CODE/main/"]
(*biadir = "ownCloud/EMORY/";*)

dirscr = "../scripts/";(*biadir<>"bia_j/CODE/scripts/";*)

(*loads matrix behavior+spikes  ftab[a_,pea_:0,hl_:1,dir_:"ownCloud/EMORY/bia_j/Datasets_gt200trials/"]*) 
DumpGet[dirscr <> "ftab4.mx"](*including large deviations for vocal output*)
(*transform tab with mayority rule->0 mayzero[tab]*)
DumpGet[dirscr <> "mayzero.mx"]
(*from matrix it gets h0 for patterns that |h0|>1, ({h,n,pM,z,+-1,patt}) h0goodpatt2[tab_,hcutoff_:1.,Nmax_:10000]*)
DumpGet[dirscr <> "h0goodpatt2.mx"](*sin filtro final en h, pero si con filtro para patt que no ocurren, el cual dep de h*)
(*from h1 and p it gets hef=sum_v Jvw and Jvw as {hef,Jvw}=couplings[h1,p]*)
DumpGet[dirscr <> "couplings.mx"]
(*gets marginal magnetizetion and e-correction to fields h0, {mi,h1ef}=margmagn[M_,h1_,hef_,Jvw_,emax_:1.,ne_:20]*)
DumpGet[dirscr <> "margmagn.mx"]
(*gets the column of mi and h1ef up to where the taylor expansion holds, epsilonopt[h0v_,h1ef_,M_,emax_:1,ne_:20]*)
DumpGet[dirscr <> "epsilonopt.mx"]
(*gets significance level from reshuffling data, needs previous func, list of m(n_false), for n_f={0.16, 0.32, 0.64, 1.28, 2.56, 5.12} *)
DumpGet[dirscr <> "signlevellower4.mx"]

(*Dumping RESULTS somewhere*)
idname = "test"; (*name for files*)
dirres = "../resdump/";(*biadir <> "bia_j/CODE/resdump/";*)(*Nmax=500; hcutoof=0.01*)
(*to print file names later*)
beh = Association[0 -> "Pitch", 1 -> "Amplitude", 2 -> "Entropy"];
listf = Import[(*biadir <> "bia_j/CODE*)"../dataset/list_files.txt", "Table"] // Flatten;

(*PARAMETERS*)
numfiles=2;(*number of files in mat format to run the analysis, located in /dataset/ and listed in /dataset/list_files.txt*)
hcutoff=0.01;(*to filter high order words that do not occur*)
Nmax=500;(*to build final Ising model from words with highest fields*)

ne = 20;(*iterations in delta epsilon, for magnetization*)
nmix = 10;(*numbers of reshuffling, to get magnetization m(n_f) at a certain full false discovery rate n_false (at least 10, if using n_f=0.3)*)

(*RUNNING*)
Do[
 ParallelDo[(*to run over available cores*)

  (*getting matrix behavior+spikes from spikes time and continuos behavior data (1st column behavior, and the rest time bins)*)
  hl=1;
  mat = ftab4[a, pae, hl,(* biadir <> "bia_j/CODE*)"../dataset/"];(*coding above median behavior=hl=1 and 0 otherwise*)
  (*mat: this is a matrix (M x T) of binary variables*)

  (*number of samples*)
  M = Length[mat];
  (*number of time bins*)
  T = Dimensions[mat][[2]];
  
  (*mat=mayzero[mat];*)(*uncomment if there you want to transform your dataset taking more common symbol as x_i=0*)

  (*getting empirical probabilities for each bin, used by prior (simple) model*)
  p = 1. Total[mat]/M;
  
  (*filtering out patterns that do not occur and are of high order, and getting fields h_{\mu} for rest*)
  h1 = h0goodpatt2[mat,hcutoff,Nmax];(*{h,n,pM,z,+-1,patt}*)
  Nm1 = Length[h1];
  
  (*h_effective coming from couplings and couplings Jvw themselves*)(*Bottleneck in time ~mins*)
  {hef, Jvw} = couplings[h1, p];
  
  (*marginal magnetizetion and e-correction to fields*)
  {mi, h1ef} = margmagn[M, h1, hef, Jvw];
  
  (*limit to taylor serie Mean|h0|~Mean|h(e)|*)
  iopt = epsilonopt[h1[[All, 1]], h1ef, M];
  
  (*saving results in files*)
  (* /resdump/idname_filenumber_pae{0,1,2}_patt.dat: list of words analyzed*)
  (* /resdump/idname_filenumber_pae{0,1,2}.dat:      list of properties of words analyzed*)
  res = Sort[Join[h1, mi[[All, {iopt}]], h1ef[[All, {iopt}]], mi, h1ef,2], #1[[7]] > #2[[7]] &];
  Export[dirres <> idname <> "_" <> ToString[a] <> "_pae" <> ToString[pae] <>".dat", res[[All, Select[Range[6 + 2 + 2 (ne + 1) ], # != 6 &]]],"Table"];
  Export[dirres <> idname <> "_" <> ToString[a] <> "_pae" <> ToString[pae] <>"_patt.dat", res[[All, 6]], "Table"];
  
  
  (*--Significance test--: magnetizations m(n_f) thresholds for different false discovery rate n_false + some data properties*)
  If[pae==0,
  (*if considering several behaviors with same neural activity, run this only for one matrix (as it only depends on marginal probabilities)*)
  mlevel = signlevellower4[mat, iopt, nmix, ne, hcutoff, Nmax];(*{mc016,mc032,...,mc512}*)
  Export[dirres <> idname <> "_" <> ToString[a] <> "_pae" <> ToString[pae] <>"_info.dat",mlevel~Join~{ Mean[p], StandardDeviation[p], M, iopt}, "Table"];
  ];
  (* /resdump/idname_filenumber_pae{0,1,2}_info.dat: vector of 10 components 
         list of m(n_false), for n_f={0.16, 0.32, 0.64, 1.28, 2.56, 5.12}
         Mean probability of bins being active
         StandardDeviation of that probability
         number of samples
         number of column for limit of taylor serie*)
  npatt03=Select[res[[All,7]],#>mlevel[[2]]&]//Length;
  Print["File: ",listf[[a]],"\tBehavior encoded (hl="<>ToString[hl]<>"): ",beh[pae],"\tNum. words (w/o filtering): ",npatt03];

  , {a,Range[numfiles]}]
 (*running pae(argument ftab) for 0:pitch, 1:amplitude and 2:entropy*)
 , {pae,{0,1,2}}]; //AbsoluteTiming

 Exit[];






