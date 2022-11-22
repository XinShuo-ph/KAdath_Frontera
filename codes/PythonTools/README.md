# PythonTools

## Author(s)    : Samuel D. Tootle,  L. Jens Papenfort
## Contributor(s) : Konrad Topolski

Maintainer(s):  
Samuel D. Tootle - tootle@itp.uni-frankfurt.de  
License      : GPLv3+ for all other code  
--------------------------------------------------------------------------

# Overview

A common question regarding the generation of initial conditions is, 
"what resolution should I use?" This question is
not trivially answered as it depends on the physics of interest and the configuration.
For symmetric configurations where the interest is primarily on (post-)merger dynamics,
coarse resolution will likely be sufficient as the remaining violations will have minimal
impact on the results. However, research focused on the inspiral dynamics and high resolution
simulations focused on informing analytical models will desire initial conditions where
the violations are minimal.

Here we include a set of tools that allow for ingestion of FUKA initial data into Python 
for analysis. The vast majority of these tools have been available, in some form, since FUKAv1 
as the "kadath_readers" that were previously stored in each initial data solver directory. 
The previous kadath_readers have been removed from the FUKAv(1/2) solvers in favor of the
new centralized location.  The key benefit is

* When building the solvers, one no longer needs to build the Python libraries or manually
disable them from being built within the CMakelist.txt file for systems that do not have 
Python install
* The readers, once built, are composed into a central Python library `fukaID_readers` which
enables easy and consistent access to readers of all the ID configurations currently supported.
* The Python libary is written such that additional functionality can either be added via
the C++ readers or by external Python scripts.
* The rewritten readers are now composed of common functions that populate `System_of_eqs` 
KADATH object and the python Dictionary.  This allows for consistency amongst the readers and
for changes within these global routines to be accessible to the readers.
* An additional functionality not previously available with the old `kadath_readers` is the
ingestion of the `info` file.

# Nomenclature and Units

The readers do not perform any unit conversions.  As such, the raw values obtained from the readers are in geometric units such that `G = c = 1` unless otherwise specifically stated


# Organization

`PythonTools` is laid out in the following way

1. `CMakeLists.txt` is the file needed by CMake to compile the codes.  It attempts to find
the host Python3 and Boost libraries.  If the environment varialbes `PYTHON_ROOT` and `BOOST_ROOT`
are not set, build errors may be encountered.
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains a reader for each initial data type which will be built into a 
shared library that can be read by Python.
    * `src/include` directory contains tools that are solely used by the python readers.
1. `lib` contains the new centralized Python libraries
    * `lib/fukaID_readers` is the python library directory which contains various `__init__.py`
    files that construct the library once the readers are built.
    * `lib/fuka_plot_tools` is a new library that includes some basic plot tools which enable
    the `plot_fukaid_1D.py` and `plot_fukaid_2D.py` scripts
1. `plot_fukaid_1D.py` and `plot_fukaid_2D.py` are plot scripts to enable easy access to
visualizing the interpolated solution of various quantities from the FUKA ID slice.
1. `test.py` is a basic script that shows the minimal code to initiate a reader using the example
initial data stored in `Example_id`
1. `Example_id` stores example initial data for playing around with
1. `Example_2D.sh` is a script that shows an advanced setup

# Getting Started

1. Use of the `test.py` file gives a very easy introduction to how to import and access
the Python readers for each ID type
2. For help using `plot_fukaid_1D.py` or `plot_fukaid_2D.py`, simply run these scripts with `--help`
Simple examples of both:
    * ```./plot_fukaid_1D.py --bh -f ./Example_id/converged_BH_TOTAL_BC.0.5.0.0.09.dat --vars cPsi --extent -3 3 --log```
    * ```./plot_fukaid_2D.py --bh -f ./Example_id/converged_BH_TOTAL_BC.0.5.0.0.09.dat --vars cPsi --extent -3 3 -3 3 --log```
Notes:
* By default a `pickle` file is generated once the solution as been interpolated based on the input
coordinates.  This can be disabled using the `--no_pickle` flag.
* Replotting using the pickle file is done by simply using the `pickle` file instead of 
the `dat` file
* The plot is saved to file by default and no pop-up of the plot is generated.  To see the plot,
use the `--pltshow` flag

# Quiver Plots

The `plot_fukaid_2D.py` is capable of using quiver plots as well.  If used in conjunction with a
FUKA DAT file, the quiver plot will by default use the same filename as that supplied with `-f`,
however, one can also reuse a pickle file for the 2D plot and a separate ID file for the quiver.
An example of what such a command could look like is:
```bash
./plot_fukaid_2D.py \
--bhns \
-f converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11-cPsi_2D.pickle \
--vars cPsi --extent -50 50 -50 50 \
--landscape --pltshow --log --cbar --vmin -10 --vmax -2 --npts 512 \
--qvar shift --qpts 64 --qscale 0.05 \
-qf Example_id/converged_BHNS_ECC_RED.togashi.35.0.6.0.52.3.6.q0.487603.0.1.11.dat
```


# Acknowledgements

The origin Python readers were written by L. Jens Papenfort which have since been rewritten by Samuel Tootle.  Konrad Topoloski has also contributed to the python readers by developing the original BHNS reader which has since been rewritten to be included in this library.

# Outstanding Tasks

The following are on the list of things to do based on expected level of effort:

1. Include optional excision filling
1. Include norm computation

# Contributing

1. Feedback and bugs are always a welcomed contribution.  Please report via the bitbucket issue
tracker
2. Contributions of additional python scripts or useful tools that would extend the current
capabilities can either be suggested via the issue tracker or be submitted as a pull request.