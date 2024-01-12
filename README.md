# nu-osc – GLoBES-based sensitivity studies and global fits

This is a code framework for simulations of neutrino oscillations in the SM and beyond.
It is intended for sensitivity studies as well as global fits. Being highly modular, it
comes with a comprehensive set of experiments already implemented, either natively in GLoBES,
or via external codes.

The full version (branch `master`) has many dependencies (nusquids, ROOT, ...) and is therefore
difficult to install on a new system. There's a lightweight version (branch `lightweight`) that
only requires [GLoBES](https://www.mpi-hd.mpg.de/~globes/) and the
[GNU Scientific Library](https://www.gnu.org/software/gsl/). It excludes support for some experiments
and some BSM scenarios.

## Installation

good luck!

## Running the Code

Run `./nu --help` to get an overview of the command line options.

Most important are:
 - `-e <experiment alias or glb file>`: specifies which experiments to include in the fit,
 - `-a <action>`: this is typically `PARAM_SCAN`
 - `-p <param>,<min>,<max>,<n_steos>[,LOGSCALE]`: defines a parameter to scan over
 - `-m <param>`: specifies a parameter to marginalize over
 - `-s <param>,<min>,<max>,<n_steos>[,LOGSCALE]`: parameters to include in a pre-scan.
This is a rough and scan, with systematic uncertainties switched off, to narrow down the
local minima of the manifold that will be marginalized over in the full fit. This is useful
if there are many parameters specified via the `-m` option, so there's a risk the (local)
minimizer won't find the true global minimum
 - `-t <param>=<value>`: can be used to define the true values of parameters.

A typical command line for a run studying CP violation could be

`./nu -a PARAM_SCAN -e t2k-fd.glb -p "DELTA_0,0,6.28,32" -p "TH13,0,0.2,20" > output.txt`

Several additional examples for how to run the code can also be found in the `run-all` script.

## Adding Experiments

To add support for a new GLoBES (AEDL) experiment, do the following:
- put all the files relevant to the experiment in a suitable directory, ideally
  a subdirectory of glb/
- make sure GLoBES can find the files by adding their location to GLB_PATH.
  This happens in nu.cc:load_exps()
- if your experiment needs multiple .glb files (for instance for near/far detectors) or
  any kind of specialized initializtaion, define an alias for the experiment,
  and when that alias is passed with the -e option (e.g. ./nu -e DUNE ...),
  have the code do all the necessary initializations. This should go in
  nu.cc:load_exps(). See the existing examples in that function for
  inspiration.
- if your experiment needs a user-defined chi^2 function, add it to sys.c and
  make sure it gets loaded in main(), where also other user-defined chi^2 functions
  are loaded

Adding new codes that do not use GLoBES, but separate fitting codes, modifications are
necessary in several places. A good strategy is to pick an experiment that's been added
that way (e.g. MiniBooNE), and search the source files for occurrences of that experiment's
name.
