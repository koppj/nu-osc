/* ----------------------------------------------------------------- */
/* Include file containing LAr detector specific stuff for far-only  */
/* simulation                                                        */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */
/*          Baseline/Detector properties,  Channels                  */
/* ----------------------------------------------------------------- */

/* Fiducial target mass */
//$target_mass =       MASS      /* kt */

/* Baseline */
//$profiletype =          1
//$baseline =             1300.0 /* km */

//include "wbb_lar_common.inc"


/* ----------------------------------------------------------------- */
/*                                Rules                              */
/* ----------------------------------------------------------------- */

/* CC */


rule(#rule7)<
        @signal = 0.85@#ANTI_nu_mu_bar_disappearance_CC
        @signalerror = 0.05 : 0.02
        @background = 0.005@#ANTI_nu_mu_bar_disappearance_NC : 0.85@#ANTI_nu_mu_beam
        @backgrounderror = 0.1 : 0.02
	//@backgroundcenter = 1.0 : 0.0
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>



rule(#rule6)<
        @signal = 0.85@#nu_mu_disappearance_CC
        @signalerror = 0.1 : 0.0001
        @background = 0.005@#nu_mu_disappearance_NC
        @backgrounderror = 0.1 : 0.0001
        @backgroundcenter = 1.0 : 0.0
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>

rule(#rule4)<
        @signal = 0.85@#nu_e_appearance_CC
        @signalerror = 0.0001 : 0.0001
        @background =  0.0003@#nu_mu_disappearance_CC  : 0.85@#nu_e_beam
        // The following line includes NC (pi0) background
        // Coefficients for # of NC BG events = 1% # of \nu_\mu disapp. events:
        //   120 GeV beam: 0.009507
        //    60 GeV beam: 0.008013
//        @background =  0.0003@#nu_mu_disappearance_CC  : 0.85@#nu_e_beam : 0.009507@#nu_mu_disappearance_NC // FIXME
        @backgrounderror = 0.022 : 0.0001
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>



rule(#rule5)<
        @signal = 0.85@#ANTI_nu_e_bar_appearance_CC
        @signalerror = 0.001 : 0.0001
        @background = 0.0003@#ANTI_nu_mu_bar_disappearance_CC  :  0.85@#ANTI_nu_e_bar_beam   :  0.85@#ANTI_nu_e_beam
        // The following line is for inclusion of NC (pi0) background
        // Coefficients for # of NC BG events = 1% # of \bar\nu_\mu disapp. events:
        //   120 GeV beam: 0.007582
        //    60 GeV beam: 0.006436
//        @background = 0.0003@#ANTI_nu_mu_bar_disappearance_CC  :  0.85@#ANTI_nu_e_bar_beam   :  0.85@#ANTI_nu_e_beam : 0.007582@#ANTI_nu_mu_bar_disappearance_NC // FIXME
        @backgrounderror = 0.022 : 0.0001
        @sys_on_function = "chiSpectrumTilt"
        @sys_off_function = "chiNoSysSpectrum"
>

