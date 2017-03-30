$version="3.0.0"


/* ----------------------------------------------------------------- */
/*                       Baseline/Detector properties                */
/* ----------------------------------------------------------------- */

/* Binning */
$emin =                 0.5
$emax =                60.0
$binsize = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
             0.125, 0.125, 0.125, 0.125, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4}


$sampling_min =         0.5
$sampling_max =        60.0
$sampling_stepsize = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                       0.125, 0.125, 0.125, 0.125, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}


//$emin =                 0.5
//$emax =                20.0
//$bins =               156
//
//$sampling_min =        0.5
//$sampling_max =       20.0


/* Technical information */
$filter_state = 0
$filter_value = 0.125


/* Energy Resolution */
energy(#mucc)<
        @type = 1
        @sigma_e = {0.0,0.15,0.0}
>

energy(#muqe)<
        @type = 1
        @sigma_e = {0.0,0.05,0.0}
>

energy(#elec)<
        @type = 1
        @sigma_e = {0.0,0.15,0.0}
>

energy(#hadr)<
        @type = 1
        @sigma_e = {0.0,0.40,0.0}
>

energy(#nc)<>
include "nc_smear_nu_ereco.dat"



/* Cross sections */
cross(#CC)<
        @cross_file = "XCC.dat"
>

cross(#NC)<
        @cross_file = "XNC.dat"
>

cross(#QE)<
        @cross_file = "XQE.dat"
>


/* ----------------------------------------------------------------- */
/*                                Channels                           */
/* ----------------------------------------------------------------- */

/* Neutrino mode */
channel(#nu_mu_cc)<
        @channel =      #fnal_plus:     +:      m:      m:      #CC:    #mucc
>

channel(#nu_mu_qe)<
        @channel =      #fnal_plus:     +:      m:      m:      #QE:    #muqe
>

channel(#nu_mu_nc)<
        @channel =      #fnal_plus:     +:      m:  NOSC_m:     #NC:    #nc
>

channel(#nu_e_beam)<
        @channel =      #fnal_plus:     +:      e:      e:      #CC:    #elec
>

channel(#nu_e_bar_beam)<
        @channel =      #fnal_plus:     -:      e:      e:      #CC:    #elec
>

channel(#nu_e_sig)<
        @channel =      #fnal_plus:     +:      m:      e:      #CC:    #elec
>


/* Antineutrino mode */
channel(#anu_mu_cc)<
        @channel =      #fnal_minus:    -:      m:      m:      #CC:    #mucc
>

channel(#anu_mu_qe)<
        @channel =      #fnal_minus:     -:      m:     m:      #QE:    #muqe
>

channel(#anu_mu_nc)<
        @channel =      #fnal_minus:    -:      m: NOSC_m:      #NC:    #nc
>

channel(#anu_e_beam)<
        @channel =      #fnal_minus:    -:      e:      e:      #CC:    #elec
>

channel(#anu_e_bar_beam)<
        @channel =      #fnal_minus:    +:      e:      e:      #CC:    #elec
>

channel(#anu_e_sig)<
        @channel =      #fnal_minus:    -:      m:      e:      #CC:    #elec
>

channel(#anu_e_bar_sig)<
        @channel =      #fnal_minus:    +:      m:      e:      #CC:    #elec
>

channel(#anu_mu_plus_cc)<
        @channel =      #fnal_minus:    +:      m:      m:      #CC:    #mucc
>

