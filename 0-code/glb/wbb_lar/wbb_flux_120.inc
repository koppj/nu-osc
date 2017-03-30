/* ----------------------------------------------------------------- */
/* Input file containing flux definitions for 120 GeV WBB            */
/* see README file for copyright information and version history     */
/* ----------------------------------------------------------------- */

/* 
   NOTE: Conversion from pot to time according to LBNE baseline design for CDR:
      FNAL beam 7.3 x 10^20 POT per year (708kW, 1.98e7s per yr, 120GeV)
   (information from Lisa Whitehead, 13 Aug 2010)
*/

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


