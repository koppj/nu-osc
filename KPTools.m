(*BeginPackage["KPTools`"];*)

Clear[Global`SaveFigures];

(* Reduce history length to save memory *)
$HistoryLength = 5;

(* Useful packages *)
<< ErrorBarPlots`;
<< LogTicks`;
<< MyGrid`;

(* Commands for plot legends *)
LegendLine[Style_] := Graphics[{Style,Line[{{-1,0},{1,0}}]},AspectRatio->0.6,ImageSize->30];

LegendBox[Style_]  := Graphics[{Style,Rectangle[{0,0}]},AspectRatio->0.6,ImageSize->30];

LegendBox[Style_,EdgeStyle_] :=
  Graphics[{Style,EdgeForm[EdgeStyle],Rectangle[{0,0}]},AspectRatio->0.6,ImageSize->30];

LegendBound[Style_] := 
  Graphics[{Style,Arrowheads[.6], AbsoluteThickness[1.5],Line[{{0,-1},{0,1}}],
            Arrow[{{0,0},{-1,0}}]},AspectRatio->0.6,ImageSize->30];

Options[LegendDataPoint] = {
  PlotStyle -> Directive[Black,AbsoluteThickness[1]],
  Frame     -> False,
  Axes      -> False,
  ImageSize -> 30
};
Options[LegendDataPoint] = Union[Flatten[{
    {Evaluate[FilterRules[Options[ErrorListPlot], Except[Options[LegendDataPoint]]]]},
  Options[LegendDataPoint]
}]];
LegendDataPoint[opts:OptionsPattern[]] :=
  ErrorListPlot[{{{0,0},ErrorBar[1,1]}},
    Evaluate[Sequence @@ Evaluate[
      {FilterRules[Union[Evaluate[FilterRules[Options[LegendDataPoint], Except[{opts}]]], {opts}],
                  Options[ErrorListPlot]]}]] ];
(*LegendDataPoint[Style_, opts:OptionsPattern[]] :=
  ErrorListPlot[{{{0,0},ErrorBar[1,1]}}, Frame->False, ImageSize->30];*)

(* Common options for plotting functions *)
MyOptions = {
  BaseStyle   -> {FontFamily -> "Times", FontSize -> 18, TextAlignment -> Left},
  PlotRange   -> All,
  Axes        -> False,
  Frame       -> True,
  AspectRatio -> 0.9,
  ImageSize   -> 400};
Off[SetOptions::optnf];
SetOptions[Show,MyOptions];
SetOptions[Plot,MyOptions];
SetOptions[LogPlot,MyOptions];
SetOptions[LogLinearPlot,MyOptions];
SetOptions[LogLogPlot,MyOptions];
SetOptions[Plot3D,MyOptions];
SetOptions[ListPlot,MyOptions];
SetOptions[ListLogPlot,MyOptions];
SetOptions[ListLogLinearPlot,MyOptions];
SetOptions[ListLogLogPlot,MyOptions];
SetOptions[ListPlot3D,MyOptions];
SetOptions[Histogram,MyOptions];
SetOptions[ContourPlot,MyOptions];
SetOptions[ListContourPlot,MyOptions];
SetOptions[DensityPlot,MyOptions];
On[SetOptions::optnf];


(* Numerical constants *)
NN = {
   (* Astrophysics *)
   GN -> 6.708*^-39/1*^18,    (* eV^-2, Newton's constant *)
   Msun -> 1.989*^30 kg,
   MMW -> 7*^11 Msun,         (* Mass of Milky way, Wikipedia article "Milky Way" *)

   (* Particle physics *)
   \[Alpha] -> 1/137, GF -> 1.16637*^-5 / GeV^2,
   MZ -> 91.19*^9, me -> 511.*^3, mmu -> 105.658367MeV, mtau -> 1776.82MeV, mp -> 938*^6,
   Mu -> 931.494028 MeV,      (* Atomatic mass unit *)
   a0 -> 0.529177*^-10 meter, (* Bohr radius *)

   (* Unit conversion *)
   eV -> 1., keV -> 1*^3, MeV -> 1*^6, GeV -> 1*^9, TeV -> 1*^12,
   pb -> 2.5766*^-27, (* eV^-2 *)
   fb -> 1*^-3 pb,
   meter -> 5.076*^6, km -> 1000*meter,
   cm -> 5.076*^4, fm -> 1*^-15 meter,
   pc -> 30.857*^15 meter, kpc -> 1*^3 pc,
   sec -> 1.523*^15, hours -> 3600 sec, days -> 24 hours, yrs -> 365 days,
   kg -> 5.62*^35 eV,
   c -> 2.9979*^10 (* cm/s *)
};


(* Likelihood for Poisson distributed random variable *)
PoissonLikelihood[obs_, th_] := 2.0 * (th - obs + If[obs > 0, obs*Log[obs/th], 0])

(* Likelihood for Gauss distributed random variable *)
GaussLikelihood[obs_, th_, sigma_] := (th - obs)^2/sigma^2


(* Compute p value corresponding to given number of sigmas (2-sided) *)
PValue[sigmas_] := Module[{},
    Return[2 NIntegrate[1/(Sqrt[2 \[Pi]]) Exp[-x^2/2], {x, sigmas, \[Infinity]}]];
];

(* Compute chi^2 value for given number of DOF and given p value *)
\[Chi]2[dof_, pvalue_] := Module[{},
    Return[x /. FindRoot[1 - CDF[ChiSquareDistribution[dof]][x] == pvalue, {x, 10}]];
];

(* Conversion between th and sin^2 (2 th) *)
S2ToRad[x_] := ArcSin[Sqrt[x]]/2;
RadToS2[x_] := Sin[2 x]^2;



(*Begin["`Private`"]

End[ ] (* Private *)

EndPackage[]*)

