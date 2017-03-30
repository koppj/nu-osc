/* ----------------------------------------------------------------- */
/* Include file containing LAr detector specific stuff for far-only  */
/* simulation                                                        */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */
/*          Baseline/Detector properties,  Channels                  */
/* ----------------------------------------------------------------- */

/* Fiducial target mass */
$target_mass =       MASS      /* kt */

/* Baseline */
$profiletype =          1
$baseline =          1300.0    /* km */

include "wbb_lar_common.inc"


/* ----------------------------------------------------------------- */
/*                                Rules                              */
/* ----------------------------------------------------------------- */

/* Neutrino mode */
rule(#nu_e_appearance)<
        @signal = 0.80@#nu_e_sig
        @signalerror = 0.01 : 0.0001
        @background =  0.01@#nu_mu_cc  : 0.01@#nu_mu_nc : 0.80@#nu_e_beam
        @backgrounderror = 0.1 : 0.0001
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>

rule(#nu_mu_disappearance)<
        @signal = 0.85@#nu_mu_cc
        @signalerror = 0.05 : 0.02
        @background = 0.005@#nu_mu_nc
        @backgrounderror = 0.1 : 0.02
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>


/* Antineutrino mode */
rule(#anu_e_appearance)< 
        @signal = 0.8@#anu_e_sig
        @signalerror = 0.01 : 0.0001
        @background = 0.01@#anu_mu_cc  :  0.8@#anu_e_bar_beam   :  0.8@#anu_e_beam : 0.01@#anu_mu_nc
        @backgrounderror = 0.1 : 0.0001
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>       
        
rule(#anu_mu_disappearance)<
        @signal = 0.85@#anu_mu_cc
        @signalerror = 0.05 : 0.02
        @background = 0.005@#anu_mu_nc : 0.85@#anu_mu_plus_cc
        @backgrounderror = 0.1 : 0.02
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>


