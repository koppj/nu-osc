(* -------------------------- FILE: GlbPlots.m --------------------------- *)
(*                                                                         *)
(* GlbPlots: A Mathmematica package for producing some standard GLoBES     *)
(*           plots                                                         *)
(* Author: Joachim Kopp                                                    *)
(* Date:   2009                                                            *)
(*                                                                         *)
(* ----------------------------------------------------------------------- *)


(* ======================================================================= *)
BeginPackage["GlbPlots`"];
(* ======================================================================= *)

(* Declare Functions *)
UpperError::usage =
"UpperError[t, CL]: Accepts a list of the form {{\theta, chi^2}, ...}
and determines the maximum value of \theta for which chi^2 <= CL."

LowerError::usage =
"LowerError[t, CL]: Accepts a list of the form {{\theta, chi^2}, ...}
and determines the minimum value of \theta for which chi^2 <= CL."

Th23Dm31Plot::usage =
"Th23Dm31Plot[Filename]: Contour plot in the s22th23 vs. dm31
plane. Filename is the name of the input file (Format: s22th23, dm31,
chi_NH, chi_IH). The function accepts all options that ContourPlot does
accept."

Th13DeltaPlotPatrick::usage = "" (*FIXME *)

Th13DeltaPlot::usage =
"Th13DeltaPlot[Filename]: Bobble plot in the s22th13 vs. delta_CP
plane. Filename is the name of the input file (Format: s22th13, deltaCP,
chi_NH, chi_IH). The function accepts all options that ContourPlot does
accept."

Th13FracDeltaPlot::usage =
"Th13FracDeltaPlot[Filename]: Sensitivity plot in the theta-13 vs.
fraction of delta_CP plane. Filename is the name of the input file
(Format: s22th13, deltaCP, chi_NH, chi_IH). The function accepts
all options that plot does accept. Filename can also be a list of several
file names."

Th13PrecisionPlot::usage =
"Th13PrecisionPlot[Filename]: Contour plot in the true theta-13 vs.  fraction
of true delta_CP plane of the th13 precision, i.e. the relative magnitude of
the error in the theta13 measurement. Filename is the name of the input file
(Format: true s22th13, true deltaCP, test logs222th13, chi_NH, chi_IH). The
function accepts all options that plot does accept. Filename can also be a list
of several file names."

LSensPlot::usage =
"LSensPlot[Filename]: Plots a theta-13 sensitivity reach (th13, MH, or CPV
  sensitivity) as a function of the baseline. Filename is the name of the input
file (Format: L, s22th13, deltaCP, chi_NH, chi_IH). The function accepts all
options that plot does accept."

LTh13FracDeltaPlot::usage =
"LTh13FracDeltaPlot[Filename]: ContourPlot of the fraction of delta_CP
values for which an experiment is sensitive to th13, MH, CPV as a
function of baseline L and true th13. Filename is the name of the input
file (Format: L, s22th13, deltaCP, chi_NH, chi_IH)."

ChiSq[1, "1\[Sigma]"] = 1;
ChiSq[1, "90%"]       = 2.7055;
ChiSq[1, "2\[Sigma]"] = 4;
ChiSq[1, "95%"]       = 3.8414;
ChiSq[1, "99%"]       = 6.6348;
ChiSq[1, "3\[Sigma]"] = 9;
ChiSq[1, "4\[Sigma]"] = 16;
ChiSq[1, "5\[Sigma]"] = 25;

ChiSq[2, "1\[Sigma]"] = 2.2957;
ChiSq[2, "90%"]       = 4.6051;
ChiSq[2, "2\[Sigma]"] = 6.1800;
ChiSq[2, "95%"]       = 5.9914
ChiSq[2, "99%"]       = 9.2103;
ChiSq[2, "3\[Sigma]"] = 11.8291;
ChiSq[2, "4\[Sigma]"] = 19.3339;
ChiSq[2, "5\[Sigma]"] = 28.7437;


(* ----------------------------------------------------------------------- *)
Begin["`Private`"];
(* ----------------------------------------------------------------------- *)

Get["LogTicks`"];

(* ----------------------------------------------------------------------- *)
(*                       Data processing functions                         *)
(* ----------------------------------------------------------------------- *)

(* ----------------------------------------------------------------------- *)
UpperError[t_List, CL_] := Module[
(* ----------------------------------------------------------------------- *)
  { xmin, xmax, deltax, x0, tInterpol, s },

  Off[InterpolatingFunction::dmval, FindRoot::lstol];
  xmin = Min[t[[All, 1]]];
  xmax = Max[t[[All, 1]]];
  deltax = (xmax - xmin)/15.0;

  x0 = Select[t, #[[2]] < CL &];
  s = Piecewise[
        {{xmin,     Length[x0] <= 0},
         {xmax,     Last[x0] == Last[t]}},
        
        tInterpol = Interpolation[t, InterpolationOrder -> 1];
        x //. FindRoot[tInterpol[x] == CL, {x, Last[x0][[1]], Last[x0][[1]] + deltax}]
      ];
  On[InterpolatingFunction::dmval, FindRoot::lstol];

  Return[s];
];


(* ----------------------------------------------------------------------- *)
LowerError[t_List, CL_] := Module[
(* ----------------------------------------------------------------------- *)
  { xmin, xmax, deltax, x0, tInterpol, s },

  Off[InterpolatingFunction::dmval, FindRoot::lstol];
  xmin = Min[t[[All, 1]]];
  xmax = Max[t[[All, 1]]];
  deltax = (xmax - xmin)/15.0;

  x0 = Select[t, #[[2]] < CL &];
  s = Piecewise[
        {{xmax,     Length[x0] <= 0},
         {xmin,     First[x0] == First[t]}},
        
        tInterpol = Interpolation[t, InterpolationOrder -> 1];
        x //. FindRoot[tInterpol[x] == CL, {x, First[x0][[1]] - deltax, First[x0][[1]]}]
      ];
  On[InterpolatingFunction::dmval, FindRoot::lstol];

  Return[s];
];


(* ----------------------------------------------------------------------- *)
(*                         Plotting functions                              *)
(* ----------------------------------------------------------------------- *)

(* ======================================================================= *)
(* Th23Dm31Plot                                                            *)
(* ======================================================================= *)
Options[Th23Dm31Plot] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[2, "1\[Sigma]"], ChiSq[2, "2\[Sigma]"], ChiSq[2, "3\[Sigma]"]},
  ContourLines   -> True,
  ContourShading -> False,
  ContourStyle   -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
Options[Th23Dm31Plot] = Union[Flatten[{
  {Evaluate[FilterRules[Options[ContourPlot], Except[Options[Th23Dm31Plot]]]]},
  Options[Th23Dm31Plot]
}]];

Th23Dm31Plot[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, s22th23Min, s22th23Max, dm31Min, dm31Max
  },

  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {4}]}];
  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
          Global`BestFit = First[Sort[Data, OrderedQ[{#1[[3]], #2[[3]]}] &]][[1;;2]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "IH"|"Inverted",
          Global`BestFit = First[Sort[Data, OrderedQ[{#1[[4]], #2[[4]]}] &]][[1;;2]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[4]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "Tot"|"Both"|All|"All"|"NHIH",
          Global`BestFit = First[Sort[Data, OrderedQ[{Min[#1[[3]], #1[[4]]],
                                                      Min[#2[[3]], #2[[4]]]}] &]][[1;;2]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], Min[#[[3]], #[[4]]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];
  s22th23Min = Min[(#[[1]] &) /@ Data];
  s22th23Max = Max[(#[[1]] &) /@ Data];
  dm31Min    = Min[(#[[2]] &) /@ Data];
  dm31Max    = Max[(#[[2]] &) /@ Data];

  Return[ContourPlot[DataInterpol[s22th23, dm31],
    {s22th23, s22th23Min, s22th23Max}, {dm31, dm31Min, dm31Max}, 
    FrameLabel -> {"\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                   "2\!\(\*SubscriptBox[\"\[Theta]\", \"23\"]\)", 
                   "\!\(\*SubsuperscriptBox[\"\[CapitalDelta]m\", \"31\", \"2\"]\)"}, 
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[Th13DeltaPlot], Except[{opts}]]], {opts}],
                  Options[ContourPlot]]}]]
  ]];
];


(* ======================================================================= *)
(* Th13DeltaPlotPatrick FIXME                                              *)
(* ======================================================================= *)
Options[Th13DeltaPlotPatrick] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[2, "1\[Sigma]"], ChiSq[2, "2\[Sigma]"], ChiSq[2, "3\[Sigma]"]},
  ContourLines   -> True,
  ContourShading -> False,
  ContourStyle   -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  PlotRange      -> All,
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
Options[Th13DeltaPlotPatrick] = Union[Flatten[{
  {Evaluate[FilterRules[Options[ContourPlot], Except[Options[Th13DeltaPlotPatrick]]]]},
  Options[Th13DeltaPlotPatrick]
}]];

Th13DeltaPlotPatrick[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, logs22th13Min, logs22th13Max
  },

  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {3}]}];
  DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]]} &) /@ Data, 
                               InterpolationOrder -> 1];
  logs22th13Min = Min[(#[[1]] &) /@ Data];
  logs22th13Max = Max[(#[[1]] &) /@ Data];

  Return[ContourPlot[DataInterpol[logS22th13, DeltaCP],
    {logS22th13, logs22th13Min, logs22th13Max}, {DeltaCP, -180, 180}, 
    FrameLabel -> {"\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                   "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\)", 
                   "\!\(\*SubscriptBox[\"\[Delta]\", \"CP\"]\)"}, 
    FrameTicks -> {LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic, 
                   LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic},
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[Th13DeltaPlotPatrick], Except[{opts}]]], {opts}],
                  Options[ContourPlot]]}]]
  ]];
];



(* ======================================================================= *)
(* Th13DeltaPlot                                                           *)
(* ======================================================================= *)
Options[Th13DeltaPlot] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[2, "1\[Sigma]"], ChiSq[2, "2\[Sigma]"], ChiSq[2, "3\[Sigma]"]},
  ContourLines   -> True,
  ContourShading -> False,
  ContourStyle   -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  PlotRange      -> All,
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
Options[Th13DeltaPlot] = Union[Flatten[{
  {Evaluate[FilterRules[Options[ContourPlot], Except[Options[Th13DeltaPlot]]]]},
  Options[Th13DeltaPlot]
}]];

Th13DeltaPlot[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, logs22th13Min, logs22th13Max
  },

  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {4}]}];
  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "IH"|"Inverted",
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[4]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "Tot"|"Both"|All|"All"|"NHIH",
          DataInterpol = Interpolation[({#[[1]], #[[2]], Min[#[[3]], #[[4]]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];
  logs22th13Min = Min[(#[[1]] &) /@ Data];
  logs22th13Max = Max[(#[[1]] &) /@ Data];

  Return[ContourPlot[DataInterpol[logS22th13, Mod[DeltaCP, 360] Degree],
    {logS22th13, logs22th13Min, logs22th13Max}, {DeltaCP, -180, 180}, 
    FrameLabel -> {"\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                   "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\)", 
                   "\!\(\*SubscriptBox[\"\[Delta]\", \"CP\"]\)"}, 
    FrameTicks -> {LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic, 
                   LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic},
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[Th13DeltaPlot], Except[{opts}]]], {opts}],
                  Options[ContourPlot]]}]]
  ]];
];



(* ======================================================================= *)
(* Th13FracDeltaPlot                                                       *)
(* ======================================================================= *)
Options[Th13FracDeltaPlot] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[2, "1\[Sigma]"], ChiSq[2, "2\[Sigma]"], ChiSq[2, "3\[Sigma]"]},
  PlotStyle      -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
Options[Th13FracDeltaPlot] = Union[Flatten[{
  {Evaluate[FilterRules[Options[Plot], Except[Options[Th13FracDeltaPlot]]]]},
  Options[Th13FracDeltaPlot]
}]];

Th13FracDeltaPlot[Filenames_List, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, logs22th13Min, logs22th13Max, CLevs,
    HierarchyOptions, PlotFunction, MyPlot, i
  },

  If [Head[OptionValue[Hierarchy]] =!= List,
    HierarchyOptions = ConstantArray[OptionValue[Hierarchy], Length[Filenames]],
  (* Else *)
    HierarchyOptions = Flatten[ConstantArray[OptionValue[Hierarchy],
                         Ceiling[Length[Filenames]/Length[OptionValue[Hierarchy]]]]]
  ];

  For[i=1, i <= Length[Filenames], i++,
    Data[i] = Import[Filenames[[i]], "Table"];
    Data[i] = Cases[Data[i], {Repeated[_Real | _Integer, {4}]}];

    Switch[HierarchyOptions[[i]],
      "NH"|"Normal",
            DataInterpol[i] = Interpolation[({#[[1]], #[[2]], #[[3]]} &) /@ Data[i],
                                            InterpolationOrder -> 1],
      "IH"|"Inverted",
            DataInterpol[i] = Interpolation[({#[[1]], #[[2]], #[[4]]} &) /@ Data[i],
                                            InterpolationOrder -> 1],
      "Tot"|"Both"|All|"All"|"NHIH",
            DataInterpol[i] = Interpolation[({#[[1]], #[[2]], Min[#[[3]], #[[4]]]} &) /@ Data[i],
                                            InterpolationOrder -> 1],
      _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
    ];
  ];
  logs22th13Min = Min[(#[[1]] &) /@ Data[1]];
  logs22th13Max = Max[(#[[1]] &) /@ Data[1]];

  PlotFunction[i_][logs22th13_, chi_] := NIntegrate[
    HeavisideTheta[DataInterpol[i][logs22th13, delta] - chi],
    {delta, 0, 2*Pi}, Method -> "LocalAdaptive"] / (2*Pi);

  CLevs = OptionValue[Contours];
  Off[NIntegrate::inumr];
  MyPlot = Plot[
       Evaluate[Flatten[Table[PlotFunction[i][logS22th13, CLevs[[j]]],
                        {i, 1, Length[Filenames]}, {j, 1, Length[CLevs]}]]],
       {logS22th13, logs22th13Min, logs22th13Max}, PlotRange -> {0, 1},
       Axes -> False, Frame -> True,
       FrameLabel -> {"\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                      "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\)", 
                      "Fraction of \!\(\*SubscriptBox[\"\[Delta]\", \"CP\"]\)"}, 
       FrameTicks -> {LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic, 
                      LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic},
       BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
       AspectRatio -> 0.9, ImageSize -> 400,
       Evaluate[Sequence @@ Evaluate[
           {FilterRules[Union[Evaluate[FilterRules[Options[Th13FracDeltaPlot], Except[{opts}]]],
                              {opts}], Options[Plot]]}]]
    ];
  On[NIntegrate::inumr];
  Return[MyPlot];
];

Th13FracDeltaPlot[Filename_String, opts:OptionsPattern[]] :=
  Th13FracDeltaPlot[{Filename}, opts];



(* ======================================================================= *)
(* Th13PrecisionPlot                                                       *)
(* ======================================================================= *)
Options[Th13PrecisionPlot] = {
  Hierarchy    -> "Both",
  CLev         -> ChiSq[1, "3\[Sigma]"],
  FrameLabel   -> {"True \!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                   "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\)",
                   "Relative Error on \!\(\*SuperscriptBox[\(sin\), " <>
                   "\(2\)]\) 2\!\(\*SubscriptBox[\(\[Theta]\), \(13\)]\)"},
  PlotStyle    -> Directive[Thick, Blue],
  FillingStyle -> Lighter[Blue, 0.8],
  Epilog       -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};

Th13PrecisionPlot[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataBest, DataWorst, DataBestInterpol, DataWorstInterpol,
    UpperLimit, LowerLimit, ErrorList,
    logs22th13Min, logs22th13Max, prec, MyPlot, MyEpilog, OldOptions, i, j
  },

  Off[OptionValue::nodef,OptionValue::optnf];

  (* Load data file *)
  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {5}]}];
  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
      Data = ({#[[1]], #[[2]], #[[3]], #[[4]]} &) /@ Data,
    "IH"|"Inverted",
      Data = ({#[[1]], #[[2]], #[[3]], #[[5]]} &) /@ Data,
    "Tot"|"Both"|All|"All"|"NHIH",
      Data = ({#[[1]], #[[2]], #[[3]], Min[#[[4]], #[[5]]]} &) /@ Data,
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];

  (* Determine th13 range *)
  logs22th13Min = Min[Data[[All, 1]]];
  logs22th13Max = Max[Data[[All, 1]]];

  (* Determine allowed range for each true th13, true delta *)
  Data = Split[Data, #1[[1;;2]] == #2[[1;;2]] &];
  DiscoveryLimit = Select[Data, LowerError[#[[All,{3,4}]], OptionValue[CLev]] == logs22th13Min &];
  DiscoveryLimit = Max[DiscoveryLimit[[All, All, 1]]];
  UpperLimit = 10^UpperError[#[[All,{3,4}]], OptionValue[CLev]] & /@ Data;
  LowerLimit = 10^LowerError[#[[All,{3,4}]], OptionValue[CLev]] & /@ Data;

  (* Consistency check *)
  ErrorList = {};
  For[j=1, j <= Length[Data], j++,
    If[UpperLimit[[j]] >= 10^Max[Data[[j, All, 3]]] ||
       LowerLimit[[j]] <= 10^Min[Data[[j, All, 3]]],
      AppendTo[ErrorList, Data[[j, 1, 1]]]
    ];
  ];
  If[ErrorList != {},
    Print["WARNING: Sampling range insufficient to determine th13 precision at logs22th13 = ",
          Union[ErrorList]];
  ];
  Data = Table[{Data[[j,1,1]], Data[[j,1,2]],
           Abs[(UpperLimit[[j]] - LowerLimit[[j]]) / 10^Data[[j,1,1]]]},  {j, 1, Length[Data]}];
  (*    Data = { #[[1,1]], #[[1,2]], Abs[(10^UpperError[#[[All,{3,4}]], OptionValue[CLev]]
           - 10^LowerError[#[[All,{3,4}]], OptionValue[CLev]]) / (10^#[[1,1]])] } & /@ Data;*)


  (* Determine best/worst allowed range for every th13 *)
  Data = Split[Data, #1[[1]] == #2[[1]] &];
  DataBest  = { #[[1,1]], Min[#[[All,3]]] } & /@ Data;
  DataWorst = { #[[1,1]], Max[#[[All,3]]] } & /@ Data;
  DataBestInterpol  = Interpolation[DataBest,  InterpolationOrder -> 1];
  DataWorstInterpol = Interpolation[DataWorst, InterpolationOrder -> 1];
  DeltaS22th13Min = 0;
  DeltaS22th13Max = Max[DataWorst[[All,2]]];

  (*  MyEpilog = {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.03,0.02}], {-1,-1}],
    Arrow[{{DiscoveryLimit, 0.03*DeltaS22th13Max}, {DiscoveryLimit+0.2, 0.03*DeltaS22th13Max}}],
    Arrow[{{DiscoveryLimit, 0.97*DeltaS22th13Max}, {DiscoveryLimit+0.2, 0.97*DeltaS22th13Max}}],
    Text["Discovery Limit", {DiscoveryLimit, 0.95*DeltaS22th13Max}, {1, 1}, {0, 1}],
    Dashed, Line[{{DiscoveryLimit, DeltaS22th13Min-10}, {DiscoveryLimit, DeltaS22th13Max+10}}]};*)

  (* Plot *)
  MyPlot = Plot[{DataBestInterpol[logs22th13], DataWorstInterpol[logs22th13]},
    {logs22th13, logs22th13Min, logs22th13Max},
    FrameTicks -> {LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic,
                   LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]], Automatic},
    Filling -> {1 -> {{2}, OptionValue[FillingStyle]}},
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    AspectRatio -> 0.9, ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[Th13PrecisionPlot], Except[{opts}]]],
                           {opts}], Options[Plot]]}]]
  ];

  On[OptionValue::nodef,OptionValue::optnf];
  Return[MyPlot];
];



(* ======================================================================= *)
(* LSensPlot                                                               *)
(* ======================================================================= *)
Options[LSensPlot] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[1, "1\[Sigma]"], ChiSq[1, "2\[Sigma]"], ChiSq[1, "3\[Sigma]"]},
  PlotStyle      -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
(*Options[LSensPlot] = Union[Flatten[{
  {Evaluate[FilterRules[Options[Plot], Except[Options[LSensPlot]]]]},
  Options[LSensPlot]
}]];*)

LSensPlot[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataSplit, LMin, LMax, logs22th13Min, logs22th13Max, S, SInterpol,
    MyPlots, OldOptions, i
  },

  Off[OptionValue::nodef,OptionValue::optnf];

  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {5}]}];
  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
          DataSplit = Split[({#[[1]], #[[2]], #[[4]]} &) /@ Data,
                            #1[[1]] == #2[[1]] &],
    "IH"|"Inverted",
          DataSplit = Split[({#[[1]], #[[2]], #[[5]]} &) /@ Data,
                            #1[[1]] == #2[[1]] &],
    "Tot"|"Both"|All|"All"|"NHIH",
          DataSplit = Split[({#[[1]], #[[2]], Min[#[[4]], #[[5]]]} &) /@ Data,
                            #1[[1]] == #2[[1]] &],
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];

  LMin          = Min[(#[[1]] &) /@ Data];
  LMax          = Max[(#[[1]] &) /@ Data];
  logs22th13Min = Min[(#[[2]] &) /@ Data];
  logs22th13Max = Max[(#[[2]] &) /@ Data];
  
  MyPlots = {};
  CLevs = OptionValue[Contours];
  For[i=1, i <= Length[CLevs], i++,
    S = {#[[1, 1]], UpperError[Transpose[{#[[All, 2]], #[[All, 3]]}], CLevs[[i]]]}& /@ DataSplit;
    SInterpol = Interpolation[S, InterpolationOrder -> 1];

    OldOptions = Options[Plot];
    SetOptions[Plot, FrameLabel -> {"L [km",
                       "\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                       "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\) reach"}];
    SetOptions[Plot,
               FrameTicks -> {Automatic, LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]],
                     Automatic, LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]]}];
    AppendTo[MyPlots, Plot[SInterpol[L], {L, LMin, LMax},
      PlotRange -> {logs22th13Min, logs22th13Max},
      PlotStyle -> OptionValue[PlotStyle][[i]], Axes -> False, Frame -> True,
      BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
      AspectRatio -> 0.9, ImageSize -> 400,
      Evaluate[Sequence @@ Evaluate[
          {FilterRules[Union[Evaluate[FilterRules[Options[LSensPlot], Except[{opts}]]], {opts}],
                    Options[Plot]]}]]
    ]];
    (SetOptions[Plot, #] &) /@ OldOptions;
  ];
  On[OptionValue::nodef,OptionValue::optnf];

  Return[Show[Sequence @@ MyPlots]];
];



(* ======================================================================= *)
(* LTh13FracDeltaPlot                                                      *)
(* ======================================================================= *)
Options[LTh13FracDeltaPlot] = {
  Hierarchy            -> "Both",
  CLev                 -> ChiSq[1, "3\[Sigma]"],
  ColorFunction        -> (Lighter[RGBColor[1,.3,.3], Clip[#, {0,1}]] &),
  ColorFunctionScaling -> False,
  Contours             -> {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},
  FrameLabel           -> {"L [km]", "\!\(\*SuperscriptBox[\"sin\", \"2\"]\) " <>
                           "2\!\(\*SubscriptBox[\"\[Theta]\", \"13\"]\)"},
  Epilog               -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]},
  LTransform           -> (1*# &),
  InterpolationOrder   -> 1
};

LTh13FracDeltaPlot[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, LMin, LMax, logs22th13Min, logs22th13Max, S, SInterpol,
    MyPlots, OldOptions, i
  },

  Off[OptionValue::nodef,OptionValue::optnf];

  Data = Import[Filename, "Table"];
  Data = Cases[Data, {Repeated[_Real | _Integer, {5}]}];
  Data = Map[(MapAt[OptionValue[LTransform], #, 1] &), Data];

  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
      DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]], #[[4]]} &) /@ Data,
                                   InterpolationOrder -> OptionValue[InterpolationOrder]],
    "IH"|"Inverted",
      DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]], #[[5]]} &) /@ Data,
                                   InterpolationOrder -> OptionValue[InterpolationOrder]],
    "Tot"|"Both"|All|"All"|"NHIH",
      DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]], Min[#[[4]], #[[5]]]} &) /@ Data,
                                   InterpolationOrder -> OptionValue[InterpolationOrder]],
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];

  LMin          = Min[(#[[1]] &) /@ Data];
  LMax          = Max[(#[[1]] &) /@ Data];
  logs22th13Min = Min[(#[[2]] &) /@ Data];
  logs22th13Max = Max[(#[[2]] &) /@ Data];
  
  PlotFunction[L_, logs22th13_, chi_] := NIntegrate[
    HeavisideTheta[DataInterpol[L, logs22th13, delta] - chi],
    {delta, 0, 2*Pi}, Method -> "LocalAdaptive"] / (2*Pi);

  CLevs = OptionValue[Contours];
  Off[NIntegrate::inumr];
  MyPlot = ContourPlot[PlotFunction[L, logS22th13, OptionValue[CLev]],
    {L, LMin, LMax}, {logS22th13, logs22th13Min, logs22th13Max},
    FrameTicks -> {Automatic, LogTicks[Floor[logs22th13Min], Ceiling[logs22th13Max]],
                   Automatic, LogTicksNoLabel[Floor[logs22th13Min], Ceiling[logs22th13Max]]},
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    AspectRatio -> 0.9, ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[LTh13FracDeltaPlot], Except[{opts}]]],
                           {opts}], Options[ContourPlot]]}]]
  ];
  On[NIntegrate::inumr];
  On[OptionValue::nodef,OptionValue::optnf];

  Return[MyPlot];
];



End[]; (* Private *)
EndPackage[];


