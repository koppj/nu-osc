/* ----------------------------------------------------------------- */
/* Input file containing flux definitions for 60 GeV WBB             */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* Neutrinos */
nuflux(#fnal_plus)<
        @flux_file="dusel060e250i002dr280dz1300km_flux.txt"
        @time = 125 * (BEAM_POWER/2.0)   /* x 1e20 pot, BEAM_POWER is in MW */
        @power = 1
        @norm = 5.19*7.504*10^12*1300*1300/120.0 / 5.202
>

/* Anti-neutrinos */
nuflux(#fnal_minus)<
        @flux_file="dusel060e250ni002dr280dz1300km_flux.txt"
        @time = 125 * (BEAM_POWER/2.0)   /* x 1e20 pot, BEAM_POWER is in MW */
        @power = 1
        @norm = 5.19*7.504*10^12*1300*1300/120.0 / 5.202
>


