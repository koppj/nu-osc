/* ----------------------------------------------------------------- */
/* Input file containing flux definitions for 120 GeV WBB            */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* Neutrinos */
nuflux(#fnal_plus)<
        @flux_file="dusel120e250i002dr280dz_tgtz30_1300km_0kmoa_flux.txt"
        @time = 72 * (BEAM_POWER/2.3)   /* x 1e20 pot, BEAM_POWER is in MW */
        @power = 1
        @norm = 5.19*7.504*10^12*1300*1300/120.0 / 5.202
>

/* Anti-neutrinos */
nuflux(#fnal_minus)<
        @flux_file="dusel120e250ni002dr280dz_tgtz30_1300km_0kmoa_flux.txt"
        @time = 72 * (BEAM_POWER/2.3)   /* x 1e20 pot, BEAM_POWER is in MW */
        @power = 1
        @norm = 5.19*7.504*10^12*1300*1300/120.0 / 5.202
>


