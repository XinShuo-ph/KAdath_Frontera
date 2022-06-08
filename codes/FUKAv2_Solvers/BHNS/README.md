# FUKAv2 - Neutron Star - Black Hole (BHNS)

# Overview

A considerable amount of the effort that went into the FUKAv2 solvers for the isolated objects (BH, NS) along with
the BBH and BNS solvers built up to constructing this solver.  The BHNS has the benefit of suffering from the sensitivity
of introducing a NS to a binary setup along with the resolution issues inherant to BH ID.  It has presented quite a challenge
even with a fairly abundant basis of literature on the topic, but, in the end, has become a work horse within my group.

With that said, if you have taken the time to work through the BBH and BNS v2 solvers, using the BHNS solver will feel quite familiar.

Note:  When referring to `Mtot` below, we will be referring to the sum of the ADM mass of the TOV solution as measured
at infinite separation with that of the Christodoulou mass of the companion BH - `Mtot := (MADM_MINUS + MCH_PLUS)`
where plus and minus simply refer to their location on the x-axis.

# Organization

1. `CMakeLists.txt` is the file needed by CMake to compile the codes
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains the relevant source files:
    - `solve.cpp`: the one and done solve code
    - `reader.cpp`: the reader can provide diagonstics from ID solutions that are computed from the ID
    - `kadath_readers.cpp`: Python libraries to allow for additional analysis of the initial data without needing to evolve it!

# Basic Usage

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bhns.info`

# Your first run!

1. Generate the initial config file by running `solve` for the first time
2. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bhns.info`
3. This will result in the generation of a pair of files containing the solution: `converged_BHNS_ECC_RED.togashi.35.0.0.2.8.q1.0.0.09.<info/dat>`

Note: In the event you have learned about generating FUKAv2 ID in the recommended order as discussed in the 
[FUKAv2 README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/) and
you have not disabled `centralized_cos`; you can look into the solver output to find that the solution of the 1.4M NS generated
in the previous runs has been reused.

We can deconstruct the name to make it understandable:

- `converged_BHNS_ECC_RED.` denotes a converged BHNS solution after the eccentricity reduction stage is completed which uses 
3.5th order PN estimates for the orbital frequency and radial infall velocity. This is meant to distinguish the solution 
from earlier stages which will be discussed later. This also distinquishes it from checkpoints that can be turned on which 
are saved to file during each iteration of the solver
- `togashi`: the leading name of the eosfile
- `35`: separation distance in geometric units!
- `0.0.`: In this case, both NSs are not spinning
- `2.8.q1`: the total mass `2.8` and mass ratio `1`
- `0.0.09`: `nshells = 0` for `ns1`, `nshells = 0` for `bh2` where the resolution in each domain is `9` collocation points in the radial and theta direction with `8` points in the phi direction
- `info`: The info file contains all the steering parameters, values used in creating the domain decomposition, stored fields, stages, settings, and controls
- `dat`: The dat file is a binary file that contains the numerical space and the variable fields

The default configuration for BHNS has a total mass of 2.8M and non-spinning.
Aside from the diagnostics observed during the solver stage, we can use the reader
to verify the ID.  This can be done by running:

`./bin/Release/reader converged_BHNS_ECC_RED.togashi.35.0.0.2.8.q1.0.0.09.info`

Which results in the following:

```
###################### Neutron Star ######################
            Center_COM = (-17.55727, 0, 0)
            Coord R_IN = +3.03469
               Coord R = [+5.76967,+6.26370] ([+8.52180,+9.25148] km)
           Coord R_OUT = +9.10407
               Areal R = +7.79621
                 NS Mb = +1.55255 (+0.46443,+1.08812,)
     Isolated ADM Mass = +1.40000
      Quasi-local Madm = +1.37888 Diff:+0.01508
         Quasi-local S = +0.00000
                   Chi = +0.00000 [+0.00000]
                 Omega = +0.00005
       Central Density = +1.37054e-03
        Central log(h) = +2.30573e-01
      Central Pressure = +2.32291e-04
    Central dlog(h)/dx = -1.39736e-15
Central Euler Constant = +7.63530e-01
     Integrated log(h) = +186.01316

###################### Black Hole ######################
            Center_COM = (+17.44273, 0, 0)
            Coord R_IN = +0.58273
               Coord R = +1.14670 [+1.69368km]
           Coord R_OUT = +9.10407
               Areal R = +2.80000
                LAPSE = [+0.40653, +0.43436]
                  PSI = [+1.56000, +1.56521]
                  Mirr = +1.40000
                   Mch = +1.40000 [+1.40000]
                   Chi = -0.00009 [+0.00000]
                     S = -0.00018
                 Omega = +0.00679

###################### Binary ######################
                     Q = +1.00000
            Separation = +35.00000 [+25.00000]
         Orbital Omega = +0.00732
            Komar mass = +2.61945
              Adm mass = +2.84343, Diff: +0.08200
            Total Mass = +2.80000 [+2.80000]
           Adm moment. = +8.15384
        Binding energy = +0.04343
            Minf * Ome = +0.02051
            E_b / Minf = +0.01551
                    Px = -1.40775e-15
                    Py = -1.49499e-15
                    Pz = +0.00000e+00
                  COMx = -0.05727, A-COMx = -0.42814
                  COMy = -0.18749, A-COMy = -0.15731
                A-COMz = +0.00000
```

The first two blocks contain information related to the component objects.
These details are covered in the 
[NS README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/NS/) and the 
[BH README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/BH/).
The only additional parameter is the `Center_COM`.  
This is the coordinate center of each object when shifted by the
"center-of-mass" of the binary or, more specifically, the location of the axis of rotation 
for the binary that can approximate a quasi-stationary solution.

The third block contains information specifically related to the binary
- `res`: resolution in each numerical domain `[r,theta,phi]`
- `q`: mass ratio such that `q <= 1`
- `Separation`: coordinate separation `d [d/Mtot] (d in km)`
- `Orbital Omega`: orbital frequency of the binary
- Komar mass
- ADM quantities
- `E_b / mu`: Binding energy normalized by the reduced mass
- `COM<x/y>`: shift in the coordinates to find a helical killing vector
- `A-COM<x/y/z>`: A numerical calculation of the COM based on the analytical prescription from Osokine+ (REF)

Note: The `Diff` noted by the ADM mass is the symmetric difference between the ADM and Komar mass.

# Understanding the BHNS INFO file

Note: It is always best practice to generate new ID using the `initial_bhns.info`.  Using old initial
data unless for very small changes is inefficient.

Using your favorite text editor, you can open up the `initial_bhns.info`.  We will go through the file,
but we'll discuss only the details relevant to the BHNS case.  For details on all the parameters you can
read more in the [Configurator README](https://bitbucket.org/fukaws/fuka/src/fukav2//include/Configurator/).

## BHNS Fixing parameters

```
binary
{
    com 0
    comy 0
    distance 35
    ecc_omega 0.0073243699259383068
    global_omega 0.0073243699259383068
    outer_shells 0
    q 1
    qpig 12.566370614359172
    res 9
    rext 70
    ns1
    {
        chi 0
        dim 3
        hc 1.26
        madm 1.3999999999999999
        mb 1.55
        n_inner_shells 0
        nc 0.0013699999999999999
        nshells 0
        omega 0
        ql_madm 1.3999999999999999
        qpig 12.566370614359172
        res 9
        rin 3.1000000000000001
        rmid 6.2000000000000002
        rout 9.3000000000000007
        eosfile togashi.lorene
        eostype Cold_Table
        h_cut 0
        interpolation_pts 2000
    }
    bh2
    {
        chi 0
        dim 3
        fixed_lapse 0.29999999999999999
        mch 1.3999999999999999
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
}
```

The above includes parameters that can be fixed by the user as well as parameters that are automated in the background
and should not be changed.  The parameters for each NS are simply copied from the isolated solution which can be read
in detail in the [NS README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/NS/) - 
the same fixing applies also in the BHNS.  The relevant parameters to discuss are

- `res` The resolution shown for the individual compact objects is the highest resolution the *isolated* dataset will be ran at.  This can be important for TOV solutions as the total baryonic mass is sensitive to the resolution.  `res 11` is recommended for production runs
- `n_inner_shells` these are additional spherical shells that can be placed in the stellar interior (only for the binary).  These are not recommended to be used as the solution is very sensitive to the domain decomposition.
- `nshells` describes the number of Additional spherical shells outside of the stellar or excision surface used in the binary ID solver and, optionally, used in the isolated ID solvers when the control `co_use_shells` is set to `on`.

The fixing parameters most relevant to the binary are

- `adot`: This is the radial infall velocity parameter when performing eccentricity reduction.  This will be discussed more
in the eccentricity reduction section below
- `distance`: this is in geometric units! It is important to pick something reasonable.  A general rule for a binary with
a few orbits is `distance = 9 * Mtot`, however this strongly depends on `q` and the spins of component objects
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
    logh on  ; log specific enthalpy
    phi on   ; fluid velocity potential
}
```

Fields documents the fields that are used in the solver and stored in the `dat` file.  Changing this has no impact.

## Stages

```
stages
{
    total off       ; v1 stage, deprecated
    total_bc on     ; hydrostatic-equilibrium stage
    ecc_red on      ; eccentricity reduction stage - hydro-rescaling with fixed `global_omega`
}
```

## Relevant Controls

```
sequence_controls
{
    checkpoint off
    corot_binary off
    fixed_bin_omega off
    fixed_lapse off
    sequences on
    update_initial off
    use_boosted_co on
    use_pn off
    resolve off
    centralized_cos on
}
```

- `checkpoint`: this will result in checkpoints being saved to file during each solving iteration - mainly helpful for high 
resolution binary ID where walltimes or server failures are a concern prior to a converged solution being obtained
- `corot_binary`: the objects are no longer fixed based on `chi` and instead provide a corotating ID solution
- `fixed_lapse`: toggling this control enables a fixed lapse on the horizon - not recommended
- `fixed_bin_omega`: toggling this allows a fixed omega to be used even when the `TOTAL` stage is activated.  This is primarily used in the background for generating
a solution from scratch.  Use of the `ECC_RED` or `TOTAL_BC` stage is recommended when doing this intentionally
- `sequences`: this toggle is enabled by default and essentially tells the driver routine to start from scratch.  
If this is enabled when attempting to use a previous solution, the previous fields and numerical space (i.e. the `dat` file) is ignored
- `update_initial`: For those familiar with the v1 solvers, there previously existed an `initial` section that would track 
where your initial guess parameters started at and what you ended up with at the end.  This is a historical artifact and mostly ignored with the v2 codes
- `use_boosted_co`: toggle whether to boost the isolated solutions prior to import.  This is essential when starting from scratch hence it is enabled by default
- `use_pn`: toggle whether to always use 3.5PN estimates.  It is important to ensure this is off if the user wants to specify their own `adot` and `global_omega`
parameters by hand (e.g. for iterative eccentricity reduction)
- `resolve`: force all implicit compact object solutions to be resolved regardless of an existing previous solution
- `centralized_cos`: stores all implicit COs into `$HOME_KADATH/COs`
- `co_use_shells`: toggle whether isolated object solvers use additional spherical shells `nshells` as noted in the binary config file.  If this is disabled, `nshells` is only used when constructing the binary space.

## Sequence Settings

```
sequence_settings
{
    solver_max_iterations 15
    solver_precision 1e-08
    initial_resolution 9
}
```

- `solver_max_iterations`: set the number of iterations not to be exceeded
- `solver_precision`: Determines what is the maximum precision allowed by the solver that determines whether or not the solution has converged
- `initial_resolution`: the resolution to use when solving the binary initially before regridding to a higher
resolution

# Your Second run!

Now that you've generated the simplest case and we have a better understanding of the config file, we can try something more interesting

1. Open the initial config file in your favorite editor
1. Set `distance 35`
1. (optional) set `res 11`
1. For `ns1` set:
    - `madm 1.18` 
    - `chi 0`
    - `res 11`
1. For `bh2` set:
    - `mch 2.42` 
    - `chi 0.52`
1. Run (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bhns.info`

This time around we see the iterative `chi` increase being done for the NS, but overall the only changes 
observed are related to the isolated NS solver (see the [NS README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/NS/) for details).  The binary solver itself is consistent.

This results in the converged dataset of `converged_BHNS_ECC_RED.togashi.28.0.52.0.3.6.q0.487603.0.0.11.info/dat`, however, the other implicit solutions have been saved as well

1. `converged_BHNS_TOTAL_BC_FIXED_OMEGA.`: is the initial solution after the import of the two isolated solutions have been solved in the binary space for a fixed COM and orbital frequency.  
The hydro fields are simply rescaled to enforce the
specified baryonic mass, but the fluid is not in hydrostatic equilibrium
2. `converged_BHNS_TOTAL_BC.`: this is the solution in complete hydrostatic equilibrium with the ADM linear momenta and the orbital frequency being fixed by the force-balance equation for the NS and varying the COM
4. `converged_BHNS_ECC_RED.`:  The final solution is one where the orbital frequency and radial infall velocity is fixed to either 3.5th order PN estimates based on the COM obtained in the `TOTAL_BC` stage or the values for `adot` and `ecc_omega` are used in the case of iterative eccentricity reduction

We can of course verify that the ID matches our expectation using

`./bin/Release/reader converged_BHNS_ECC_RED.togashi.35.0.0.52.3.6.q0.487603.0.0.11.info`

```
###################### Neutron Star ######################
            Center_COM = (-23.57751, 0, 0)
            Coord R_IN = +2.97850
               Coord R = [+5.96311,+6.22628] ([+8.80751,+9.19621] km)
           Coord R_OUT = +9.97836
               Areal R = +7.74500
                 NS Mb = +1.28557 (+0.36481,+0.92076,)
     Isolated ADM Mass = +1.18000
      Quasi-local Madm = +1.15768 Diff:+0.01892
         Quasi-local S = +0.00000
                   Chi = +0.00000 [+0.00000]
                 Omega = +0.00007
       Central Density = +1.22225e-03
        Central log(h) = +1.84227e-01
      Central Pressure = +1.58366e-04
    Central dlog(h)/dx = +2.49019e-15
Central Euler Constant = +7.58231e-01
     Integrated log(h) = +141.96974

###################### Black Hole ######################
            Center_COM = (+11.42249, 0, 0)
            Coord R_IN = +0.94427
               Coord R = +1.79591 [+2.65255km]
           Coord R_OUT = +9.97836
               Areal R = +4.66020
                LAPSE = [+0.37110, +0.39523]
                  PSI = [+1.60048, +1.62197]
                  Mirr = +2.33010
                   Mch = +2.42016 [+2.42000]
                   Chi = +0.52039 [+0.52000]
                     S = +3.04801
                 Omega = -0.04992

###################### Binary ######################
                     Q = +0.48760
            Separation = +35.00000 [+14.46188]
         Orbital Omega = +0.00807
            Komar mass = +3.55526
              Adm mass = +3.56874, Diff: +0.00378
            Total Mass = +3.60016 [+3.60000]
           Adm moment. = +13.78446
        Binding energy = -0.03142
            Minf * Ome = +0.02906
            E_b / Minf = -0.00873
                    Px = -3.48680e-15
                    Py = -2.14706e-15
                    Pz = +0.00000e+00
                  COMx = -6.07751, A-COMx = -6.05926
                  COMy = -0.03323, A-COMy = -0.10260
                A-COMz = +0.00000
```

# How NS ID is Generated

## Initial Setup

To generate the initial setup for the binary ID, we first need to make some guesses based on the input 
ADM masses and the coordinate separation

  1. `COM` - An estimate of the center-of-mass: this is purely Newtonian
  2. `global_omega` - An estimate of the orbital frequency: 3.5th PN estimate

Once these estimates are computed an interface code is ran
which 

- solves each component configuration in isolation.  
- obtains boosted isolated solutions using the estimated `global_omega`

At this point, the binary numerical space and fields are constructed and the isolated 
solutions are interpolated onto 
the new grid using the idea of superimposed solutions.  Specifically:

- a decay parameter `decay_limit := w` is chosen such that `w = distance / 2`
- the fields are interpolated such that the solutions decay exponetially away from each object as, e.g. 
`decay_rate = exp(-(r_NS / w)^4)`, where `r_NS` is the coordinate distance from the NS
- The resulting value at a given point is then simply the sum of the background with the deviations from the background from the isolated solutions

For example, if we wanted to compute the initial guess for the lapse at a given point `x`, it would be

`lapse(x) = 1. + decay_rate_NS * (lapse_NS(x) - 1.) + decay_rate_BH * (lapse_BH(x) - 1.)`

This is then repeated for all fields in all numerical domains, with the compactified domain set the asymptotic 
values of the fields, i.e. `lapse = psi = 1`, `log(h) = shift = 0`.

## Initial Solution

Once the initial guess has been setup for the BHNS, the first solver stage solves the full XCTS 
system of equations consistently and all at once, however, with regards to NS source, 
the matter is simply rescaled to achieve the desired baryonic mass.  The reason for this is two fold:

  1. The spacetime fields need to settle to a more accurate estimate before resolving the matter consistently
  2. The fluid velocity potential field `phi` needs to be initialized to the binary configuration

The output file from this stage, `converged_BHNS_TOTAL_BC_FIXED_OMEGA.`, can readily be 
discarded once a solution from a later stage is obtained.

## Hydrostatic Equilibrium

Now with a more stable background, the matter becomes a variable field and is allowed to be solved for
consistently in order to obtain a solution in hydrostatic equilibrium. Now
`global_omega` is determined by- and the ADM linear momenta are minimized by- using the force balance equation
and varying the COM.

Note: in the event one wants to later increase the resolution or make iterative changes 
(e.g. make small changes to the MADM, MB, CHI of one or both stars),
it can only be done using the solution from this stage, `converged_BHNS_TOTAL.*<info/dat>`, 
as the initial starting point.  Using the solutions from the hydro-rescaling stage
will more often than not cause the solution to diverge.  Therefore, these solutions can be useful
to retain.

### Automated resolution increase

In the event the binary resolution was set to something higher than the sequence_setting `initial_resolution`, 
the automated increase resolution will take place here and resolve the binary in hydro-static equilibrium.  
This is very important as increasing the resolution and repeating a matter-rescaling stage will result 
in a very inconsistent description of the fluid as it will include numerical errors from the interpolated solution at lower resolution.

## Eccentricity Reduction

Using either 3.5PN estimates or, in the case of iterative eccentricity reduction, 
`ecc_omega` and `adot` from the config file,
a final stage of matter rescaling is done based on the changes introduced by 
`ecc_omega` and `adot`.  This is the recommended
solution to use for evolutions and it is stored with a filename title of `converged_BHNS_ECC_RED.*<info/dat>`