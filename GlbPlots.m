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

GlbPlot2d::usage =
"GlbPlot2d[Filename]: Contour plot in a 2-parameter plane
Filename is the name of the input file (Format: param1, param2,
chi_NH, chi_IH). The function accepts all options that ContourPlot
accepts."


(* ----------------------------------------------------------------------- *)
Begin["`Private`"];
(* ----------------------------------------------------------------------- *)

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
(* GlbPlot2d                                                               *)
(* ======================================================================= *)
Options[GlbPlot2d] = {
  Hierarchy      -> "Both",
  Contours       -> {ChiSq[2, "1\[Sigma]"], ChiSq[2, "2\[Sigma]"], ChiSq[2, "3\[Sigma]"]},
  ContourLines   -> True,
  ContourShading -> False,
  ContourStyle   -> {Directive[Thick, Blue], Directive[Thick, Blue, Dashed],
                     Directive[Thick, Blue, Dotted]},
  P1Function     -> Identity,
  P2Function     -> Identity,
  Epilog         -> {Text["GLoBES "<>ToString[DateList[][[1]]], Scaled[{0.97,0.02}], {1,-1}]}
};
Options[GlbPlot2d] = Union[Flatten[{
  {Evaluate[FilterRules[Options[ContourPlot], Except[Options[GlbPlot2d]]]]},
  Options[GlbPlot2d]
}]];

GlbPlot2d[Filename_String, opts:OptionsPattern[]] := Module[
  {
    Data, DataInterpol, p1Min, p1Max, p2Min, p2Max, chi2Min
  },

  Data = Import[Filename, "Table"];
  Data = {OptionValue[P1Function][#[[1]]], OptionValue[P2Function][#[[2]]],
           #[[3]], #[[4]]} & /@ Cases[Data, {Repeated[_Real | _Integer, {4}]}];
  Switch[OptionValue[Hierarchy],
    "NH"|"Normal",
          Global`BestFit = First[Sort[Data, OrderedQ[{#1[[3]], #2[[3]]}] &]][[1;;2]];
          chi2Min = Min[Data[[All,3]]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[3]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "IH"|"Inverted",
          Global`BestFit = First[Sort[Data, OrderedQ[{#1[[4]], #2[[4]]}] &]][[1;;2]];
          chi2Min = Min[Data[[All,4]]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], #[[4]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    "Tot"|"Both"|All|"All"|"NHIH",
          Global`BestFit = First[Sort[Data, OrderedQ[{Min[#1[[3]], #1[[4]]],
                                                      Min[#2[[3]], #2[[4]]]}] &]][[1;;2]];
          chi2Min = Min[Flatten[Data[[All,{3,4}]]]];
          DataInterpol = Interpolation[({#[[1]], #[[2]], Min[#[[3]], #[[4]]]} &) /@ Data, 
                                       InterpolationOrder -> 1],
    _, Print["Invalid option: Hierarchy -> ", OptionValue[Hierarchy]]; Return[];
  ];
  p1Min = Min[(#[[1]] &) /@ Data];
  p1Max = Max[(#[[1]] &) /@ Data];
  p2Min = Min[(#[[2]] &) /@ Data];
  p2Max = Max[(#[[2]] &) /@ Data];

  Return[ContourPlot[DataInterpol[p1, p2] - chi2Min,
    {p1, p1Min, p1Max}, {p2, p2Min, p2Max}, 
    BaseStyle -> {FontFamily -> "Times", FontSize -> 14, TextAlignment -> Left},
    ImageSize -> 400,
    Evaluate[Sequence @@ Evaluate[
        {FilterRules[Union[Evaluate[FilterRules[Options[GlbPlot2d], Except[{opts}]]], {opts}],
                  Options[ContourPlot]]}]]
  ]];
];


End[]; (* Private *)
EndPackage[];


