/* ----------------------------------------------------------------- */
/* Include file containing LAr near detector specific stuff          */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */
/*                  Baseline, Detector properties                    */
/* ----------------------------------------------------------------- */

/* Target mass */
$target_mass =          1.0 /* kt */

/* Baseline */
$profiletype =          1
$baseline =             1.0 /* km */

include "wbb_lar_common.inc"


/* ----------------------------------------------------------------- */
/*                                Rules                              */
/* ----------------------------------------------------------------- */

/* Neutrino mode */
rule(#nu_e_appearance)<
        @signal = 0.80@#nu_e_sig
        @background =  0.01@#nu_mu_cc  : 0.80@#nu_e_beam : 0.01@#nu_mu_nc
        @sys_on_function  = "chiZero"
        @sys_on_errors    = {}
        @sys_off_function = "chiNoSysSpectrum"
        @sys_off_errors   = {}
>

rule(#nu_mu_disappearance)<
        @signal = 0.85@#nu_mu_cc
        @background = 0.005@#nu_mu_nc
        @sys_on_function  = "chiZero"
        @sys_on_errors    = {}
        @sys_off_function = "chiNoSysSpectrum"
        @sys_off_errors   = {}
>


/* Antineutrino mode */
rule(#anu_e_appearance)<
        @signal = 0.8@#anu_e_sig
        @background = 0.01@#anu_mu_cc  :  0.8@#anu_e_bar_beam   :  0.8@#anu_e_beam : 0.01@#anu_mu_nc
        @sys_on_function  = "chiZero"
        @sys_on_errors    = {}
        @sys_off_function = "chiNoSysSpectrum"
        @sys_off_errors   = {}
>       
        
rule(#anu_mu_disappearance)<
        @signal = 0.85@#anu_mu_cc
        @background = 0.005@#anu_mu_nc : 0.85@#anu_mu_plus_cc
        @sys_on_function  = "chiZero"
        @sys_on_errors    = {}
        @sys_off_function = "chiNoSysSpectrum"
        @sys_off_errors   = {}
>


