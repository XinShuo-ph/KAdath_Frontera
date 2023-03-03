# FUKAv2 - Binary Black Hole (BBH) Initial Data

# Overview

Here lies the binary black hole initial data solver and diagnostic codes.  The initial data is constructed using maximal
slicing with a flat spacial metric.  Therefore, the dimensionless spins of the black holes are limited to 
approximately `[-0.85, 0.85]`.  The v2 code is a vast improvement over the v1 code as it uses super-imposed BH solutions
to initialize the binary instead of building the binary from flat space.  When comparing to FUKAv1, generating equal-mass non-rotating
ID is roughly `>4x` faster due to the reduced number of stages.  When comparing to FUKAv1, generating arbitrary spin and unequal mass, the cost savings is roughly `(4 + N)x` faster 
where N is, in the case of FUKAv1, the number of iterative solutions needed to acheive a given mass ratio and spin.

Note:  

- When referring to `Mtot` below, we will be referring to the sum of the individual Christodoulou masses `Mtot := (MCH_MINUS + MCH_PLUS)`
where plus and minus simply refer to their location on the x-axis.
- `q <= 1` is built into the FUKAv2 codes

# Organization

1. `CMakeLists.txt` is the file needed by CMake to compile the codes
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains the relevant source files:
    - `solve.cpp`: the one and done solve code
    - `reader.cpp`: the reader can provide diagonstics from ID solutions that are computed from the ID

# Basic Usage

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bbh.info`

# Your first run!

1. Generate the initial config file by running `solve` for the first time
2. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bbh.info`
3. This will result in the generation of a pair of files containing the solution: `BBH_ECC_RED.10.0.0.1.q1.0.0.09.<info/dat>`

We can deconstruct the name to make it understandable:

- `BBH_ECC_RED.` denotes a converged BBH solution after the eccentricity reduction stage is completed which uses 
3.5th order PN estimates for the orbital frequency and radial infall velocity. This is meant to distinguish the solution 
from earlier stages which will be discussed later. This also distinquishes it from checkpoints that can be turned on which 
are saved to file during each iteration of the solver
- `10`: separation distance in geometric units!
- `0.0.`: In this case, both BHs are not spinning
- `1.q1`: the total mass and mass ratio are equal to one
- `0.0.09`: `nshells = 0` for `bh1`, `nshells = 0` for `bh2` where the resolution in each domain is `9` collocation points in the radial and theta direction with `8` points in the phi direction
- `.info`: The info file contains all the steering parameters, values used in creating the domain decomposition, stored fields, stages, settings, and controls
- `.dat`: The dat file is a binary file that contains the numerical space and the variable fields

The default configuration for BBH has a total mass of 1M and non-spinning.
Aside from the diagnostics observed during the solver stage, we can use the reader
to verify the ID.  This can be done by running:

`./bin/Release/reader BBH_ECC_RED.10.0.0.1.q1.0.0.09.info`

Which results in the following:

```
###################### BH_MINUS ######################
            Center_COM = (-5.00000, 0, 0)
            Coord R_IN = +0.21488
               Coord R = +0.40658
           Coord R_OUT = +1.39076
               Areal R = +1.00000
                 LAPSE = [+0.39165, +0.42362]
                   PSI = [+1.56539, +1.57119]
                  Mirr = +0.50000[+0.50000]
                   Mch = +0.50000[+0.50000]
                   Chi = +0.00000[+0.00000]
                     S = +1.03396e-12
                 Omega = +2.54278e-02
###################### BH_PLUS ######################
            Center_COM = (+5.00000, 0, 0)
            Coord R_IN = +0.21488
               Coord R = +0.40658
           Coord R_OUT = +1.39076
               Areal R = +1.00000
                 LAPSE = [+0.39165, +0.42362]
                   PSI = [+1.56539, +1.57119]
                  Mirr = +0.50000[+0.50000]
                   Mch = +0.50000[+0.50000]
                   Chi = +0.00000[+0.00000]
                     S = +1.03402e-12
                 Omega = +2.54278e-02
###################### Binary ######################
                   RES = [+9,+9,+8]
                     Q = +1.00000
            Separation = +10.00000[+10.00000]
         Orbital Omega = +2.79950e-02
            Komar mass = +9.88086e-01
            Adm   mass = +9.89694e-01, Diff: +1.62621e-03
            Total Mirr = +1.00000
             Total Mch = +1.00000
           Adm moment. = +9.63771e-01
        Binding energy = -1.03061e-02
               M * Ome = +2.79950e-02
              E_b / mu = -4.12245e-02
                  Axis = -0.00000
                    Px = -4.31131e-16
                    Py = +4.59808e-16
                    Pz = +0.00000e+00
                  COMx = -2.21951e-13, A-COMx = -1.67758e-11
                  COMy = -2.28587e-13, A-COMy = +1.68405e-11
                A-COMz = +0.00000e+00
```

The first two blocks contain information related to the component BHs.  These details are covered in the 
[BH README](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/BH/).
The only additional parameter is the `Center_COM`.  This is the coordinate center of each object when shifted by the
"center-of-mass" of the binary or, more specifically, the location of the axis of rotation for the binary that can approximate
a quasi-stationary solution.

The third block contains information specifically related to the binary
- `res`: resolution in each numerical domain `[r,theta,phi]`
- `q`: mass ratio such that `q <= 1`
- `Separation`: coordinate separation `d [d/Mtot]`
- `Orbital Omega`: orbital frequency of the binary
- Komar mass
- ADM quantities
- `E_b / mu`: Binding energy normalized by the reduced mass
- `COM<x/y>`: shift in the coordinates to find a helical killing vector
- `A-COM<x/y/z>`: A numerical calculation of the COM based on the analytical prescription from Osokine+ (REF)

Note: The `Diff` noted by the ADM mass is the symmetric difference between the ADM and Komar mass.

# Understanding the BBH INFO file

Note: It is always best practice to generate new ID using the `initial_bbh.info`.  Using old initial
data unless for very small changes is inefficient.

Using your favorite text editor, you can open up the `initial_bbh.info`.  We will go through the file,
but we'll discuss only the details relevant to the BBH case.  For details on all the parameters you can
see more in [Configurator](https://bitbucket.org/fukaws/fuka/src/fuka/include/Configurator/) README.

## BBH Fixing parameters

```
binary
{
    com 0
    comy 0
    distance 10
    ecc_omega 0.02799495935918607
    global_omega 0.02799495935918607
    outer_shells 0
    q 1
    qpig 12.566370614359172
    res 9
    rext 20
    bh1
    {
        chi 0
        dim 3
        fixed_lapse 0.29999999999999999
        mch 0.5
        mirr 0.5
        nshells 0
        omega 0
        qpig 12.566370614359172
        res 9
        rin 0.10000000000000001
        rmid 0.29999999999999999
        rout 1.5
        velx 0
        vely 0
    }
    bh2
    {
        ....
    }
}
```

The above includes parameters that can be fixed by the user as well as parameters that are automated in the background
and should not be changed.  The parameters for each BH are simply copied from the isolated solution which can be read
in detail in the [BH README](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/BH/) - 
the same fixing applies also in the BBH.

The fixing parameters most relevant to the binary are

- `distance`: this is in geometric units! It is important to pick something reasonable.  A general rule for a binary with
a few orbits is `distance = 8 * Mtot`, however this strongly depends on `q` and the spins of component objects
- `outer_shells`: This allows for additional shells to be placed near the compactified domain.  This can be helpful
for more accurate quasi-equilibrium ID at lower resolution, but otherwise can be ignored and left to zero
- `q`: this parameter is computed.  Changing it by hand does nothing
- `res`: global resolution of the binary!
- `global_omega`: the orbital frequency of the binary.  This will be discussed in detail in the relevant sections below

## Fields

```
fields
{
    conf on
    lapse on
    shift on
}
```

Fields documents the fields that are used in the solver and stored in the `dat` file.  Changing this has no impact.

## Stages

```
stages
{
    corot_equal off ; deprecated - v1 stage
    fixed_omega off ; deprecated - v1 stage
    pre off         ; deprecated - v1 stage
    total off       ; deprecated - v1 stage
    ecc_red on
    total_bc on
}
```

For the BBH, only the `TOTAL_BC` and `ECC_RED` stages are relevant.  Old stages available in the original v1 solver are since deprecated.

## Relevant Controls

```
sequence_controls
{
    checkpoint off
    corot_binary off
    fixed_bin_omega on
    fixed_lapse off
    sequences on
    update_initial off
    use_boosted_co on
    use_pn off
    resolve off
    centralized_cos on
    co_use_shells off
    initial_regrid off
}
```

- `checkpoint`: this will result in checkpoints being saved to file during each solving iteration - mainly helpful for high 
resolution binary ID where walltimes or server failures are a concern prior to a converged solution being obtained
- `corot_binary`: the objects are no longer fixed based on `chi` and instead provide a corotating ID solution
- `fixed_lapse`: toggling this control enables a fixed lapse on the horizon - not recommended
- `fixed_bin_omega`: toggling this allows a fixed omega to be used even during the `TOTAL_BC` stage.  This is primarily used in the background for generating
a solution.  Use of the `ECC_RED` is recommended instead
- `sequences`: this toggle is enabled by default and essentially tells the driver routine to start from scratch.  
If this is enabled when attempting to use a previous solution, the previous fields and numerical space (i.e. the `dat` file) is ignored
- `update_initial`: For those familiar with the v1 solvers, there previously existed an `initial` section that would track 
where your initial guess parameters started at and what you ended up with at the end.  This is a historical artifact and mostly ignored with the v2 codes
- `use_boosted_co`: toggle whether to boost the isolated solutions prior to import.  This is essential when starting from scratch hence it is enabled by default
- `use_pn`: toggle whether to always use 3.5PN estimates.  It is important to ensure this is off if the user wants to specify their own `adot` and `global_omega`
parameters by hand (e.g. for iterative eccentricity reduction)
- `resolve`: force all implicit compact object solutions to be resolved regardless of an existing previous solution
- `centralized_cos`: stores all implicit COs into `$HOME_KADATH/COs`
- `initial_regrid`: in the event one wants to regrid the solution prior to resolving a previous converged solution, this can be enabled.  This is helpful when
  - One wants to resolve a previous solution with additional shells added
  - One wants to resolve using the default regrid of the converged solution.  This will provide a slightly more optimized domain decomposition and could possibly add additional spherical shells by default.
- `co_use_shells`: toggle whether isolated object solvers use additional spherical shells `nshells` as noted in the binary config file.  If this is disabled, `nshells` is only used when constructing the binary space.  For BBH ID, it is recommended to leave this `off`

## Sequence Settings

```
sequence_settings
{
    solver_max_iterations 15
    solver_precision 1e-08
    initial_resolution 9
}
```

- `solver_max_iterations`: set the number of iterations not to exceed
- `solver_precision`: Determines what is the maximum precision allowed by the solver that determines whether or not the solution has converged
- `initial_resolution`: by default all solutions are ran at at a default resolution of `[9,9,8]` prior to regridding to a higher resolution.  In the event one wants to increase the default `initial_resolution`, it can be done here

# Your Second run!

Now that you've generated the simplest case and we have a better understanding of the config file, we can try something more interesting

1. Open the initial config file in your favorite editor
2. For BH1 set:
    - `mch 0.1` 
    - `chi -0.5`
3. For BH2 set:
    - `mch 0.9` 
    - `chi 0.85`
3. Run (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bh.info`

This time around we see the iterative `chi` increase being done for the primary BH, but overall the only changes 
observed are related to the isolated BH solvers (see the BH README for details), but the binary solver itself
is consistent.

This results in the converged dataset of `BBH_ECC_RED.10.-0.5.0.85.1.q0.1.2.0.09.info/dat`, however, the other implicit solutions have been saved as well.

We can of course verify that the ID matches our expectation using

`./bin/Release/reader BBH_ECC_RED.10.-0.5.0.85.1.q0.1.2.0.09.info`

```
###################### BH_MINUS ######################
            Center_COM = (-9.03205, 0, 0)
            Coord R_IN = +0.04067
               Coord R = +0.06982
                SHELL1 = +0.25145
                SHELL2 = +0.46651
           Coord R_OUT = +1.53569
               Areal R = +0.19319
                 LAPSE = [+0.32625, +0.37747]
                   PSI = [+1.64806, +1.68725]
                  Mirr = +0.09659[+0.09659]
                   Mch = +0.10000[+0.10000]
                   Chi = -0.50000[-0.50000]
                     S = -5.00000e-03
                 Omega = +1.22944e+00
###################### BH_PLUS ######################
            Center_COM = (+0.96795, 0, 0)
            Coord R_IN = +0.28755
               Coord R = +0.47065
           Coord R_OUT = +1.53569
               Areal R = +1.57270
                 LAPSE = [+0.27142, +0.31501]
                   PSI = [+1.77690, +1.86146]
                  Mirr = +0.78635[+0.78635]
                   Mch = +0.90000[+0.90000]
                   Chi = +0.85000[+0.85000]
                     S = +6.88500e-01
                 Omega = -3.14511e-01
###################### Binary ######################
                   RES = [+9,+9,+8]
                     Q = +9.00000
            Separation = +10.00000[+10.00000]
         Orbital Omega = +2.76358e-02
            Komar mass = +1.00257e+00
            Adm   mass = +1.00232e+00, Diff: +2.49557e-04
            Total Mirr = +0.88294
             Total Mch = +1.00000
           Adm moment. = +1.01572e+00
        Binding energy = +2.31537e-03
               M * Ome = +2.76358e-02
              E_b / mu = +2.57263e-02
                    Px = -5.56455e-16
                    Py = -4.41494e-16
                    Pz = +0.00000e+00
                  COMx = -4.03205e+00, A-COMx = -4.01449e+00
                  COMy = +7.34448e-03, A-COMy = -9.13143e-03
                A-COMz = +0.00000e+00
```

# How BBH ID is Generated

## Initial Setup

To generate the initial setup for BBH ID, we first need to make some guesses based on the input 
Christodoulou masses and the coordinate separation

  1. `COM` - An estimate of the center-of-mass: this is purely Newtonian
  2. `global_omega` - An estimate of the orbital frequency: 3.5th PN estimate

Once these estimates are computed an interface code is ran
which 

- solves each BH configuration in isolation 
(See the [BH README](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/BH/) for more details).  
- obtains boosted isolated solution using the estimated `global_omega`

At this point, the binary numerical space and fields are constructed and the isolated solutions are interpolated onto 
the new grid using the idea of super-imposed solutions.  Specifically:

- a decay parameter `w` is chosen such that `w = distance / 2 =: decay_limit`
- the fields are interpolated such that the solutions decay exponentially away from each object as 
`decay_rate = exp(-(r_BH / w)^4)` where `r_BH` is the coordinate distance the respective BH
- The resulting value at a given point is then simply the sum of the background with the deviations from the background from the isolated solutions

For example, if we wanted to compute the initial guess for the lapse at a given point `x`, it would be

`lapse(x) = 1. + decay_rate_BH1 * (lapse_BH1 - 1.) + decay_rate_BH2 * (lapse_BH2 - 1.)`

This is then repeated for all fields in all numerical domains, with the compactified domain set to the asymptotic 
values of the fields, i.e. `lapse = psi = 1`, `shift = 0`.

## TOTAL_BC Stage

The `TOTAL_BC` stage solves the full XCTS system of equations consistently and all at once.  The only thing that distinguishes 
`TOTAL_BC` from `ECC_RED` is that, by default, `TOTAL_BC` fixes the variable `global_omega` using the quasi-equilibrium 
assumption of `Madm == Mkomar` where as `ECC_RED` uses user-defined values or 3.5PN estimates of the `global_omega` and `adot`.  
However, when running a new initial data sequence, the `global_omega` is initially fixed in the `TOTAL_BC` stage (you can see 
the control enabled within the `initial_bbh.info` file) using the initial 3.5PN estimate.

There is of course a very simple reason for this - otherwise the solution would diverge.  Put simply, in ID generation of 
binary objects the shift is very sensitive in such systems and, as such, can cause wild changes in the solution before 
converge is obtained, if at all.  Additionally, the quasi-equilibrium equations are evaluated at infinity on the surface of
 the compatified domain which is incredibly course.  To get around these challenges, `global_omega` is fixed for one solving 
 stage to allow all the fields to converge to a reasonable initial solution.

After this initial solution is obtained, the `fixed_omega` control is deactivated and the `TOTAL_BC` stage is reran with 
quasi-equilibrium.  This is also done to ensure an accurate `COM` is found before obtaining 3.5PN estimates of the 
`global_omega` and `adot` within the `ECC_RED` stage.

## ECC_RED Stage

In this stage `global_omega` and `adot` are fixed to generate a less eccentric binary than when using a quasi-equilibrium 
approximation. After an initial 3.5PN estimate
binary has been created, two new parameters appear called `ecc_omega` and `adot`.

- In the event `ecc_omega` or `adot` are not present in the config file, 3.5PN estimates will always be used
- In the event the control `use_pn` is set to `on`, the values for `ecc_omega` and `adot` will **ALWAYS** be overwritten
- In the event `ecc_omega` and `adot` are both set and `use_pn` is set to `off`, these parameters will only be used in the 
`ecc_red` stage regardless of whether previous stages are ran that change `global_omega`
  - For an example, say one wanted to generate initial data using configurations and eccentricity reduction parameters from
  a published dataset in the [SXS database](https://data.black-holes.org/waveforms/catalog.html), one could set `adot` and
  `ecc_omega` manually in the `initial_bbh.info` file and these values would be used in the `ECC_RED` stage. The
  `global_omega` parameter will also be updated to reflect this once a converged solution is obtained.

### Automated resolution increase

Once the above procedures are completed for the last stage and the input `res` is higher than the `initial_resolution`, 
the low resolution solution is interpolated onto the higher resolution numerical grid and is used as the initial guess before 
running the last stage again.  No additional iterative procedures are required at this point.