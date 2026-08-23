# Hybrid Monte Carlo / Molecular Dynamics

This is a Hybrid Monte Carlo / Molecular Dynamics (MC/MD) model that is part of the **MOLDSTRUCT** toolbox. It can be used to simulate photon-matter interaction for smaller systems like single proteins.
The code is based on a modified version of GROMACS ([webpage](https://www.gromacs.org/)).
With functionality close to the normal GROMACS but with some additional parameters to control ionization. Running the code requires that atomic data (more details later) is supplied for the appropriate energy.
The model is developed by the Biophysics group at Uppsala University and is published here [link to article](PUT URL HERE WHEN PUBLISHED).

# Manual

This is a brief manual that will cover how to install, use the model and evaluate output.
Basic knowledge of GROMACS is assumed, check out [GROMACS webpage](https://www.gromacs.org/) for more information.

## Contents

- Installation
- List of input parameters
- Supplying atomic data
- Running a simulation
- Output
- Example(s)
- Limitations

## Installation

The installation process is the exact same as for a normal GROMACS installation.

1. Download the zip of the repository and extract it.
2. Place the `gromacs-4.5.4-MCMD` folder somewhere appropriate.
   This is only the installation files, you will choose later on where you want to install the software.
3. Go into `gromacs-4.5.4-MCMD` and create a new directory called `build` (or whatever you want).
4. Go into the newly created directoy.
5. Run `cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/install/location`, where `/path/to/install/location` is replaced to where the software will be installed.
6. Run `make install`
7. ==Congratualtions!== If everything worked out, the software should now be installed at your specified location.

If running into problems, try to search the issue as there are many GROMACS reasources out there that can help with the install process.

## List of input parameters

The model's parameters are set in the `.mdp` file like any other MD parameter.
They are all named `mcmd-*` and take effect only when `mdrun` is given
`-ionize`; without that flag the run is unmodified GROMACS 4.5.4 whatever the
`.mdp` says.

**Note the time unit.** `mcmd-pulse-peak-time` and `mcmd-pulse-fwhm` are in
**femtoseconds**, unlike `dt` and `tinit`, which stay in GROMACS' picoseconds.

```
mcmd-charge-transfer         (default = 1)    Classical over-the-barrier charge
                                              transfer between any pair of atoms.
mcmd-charge-transfer-idle    (default = 2000) Steps without a transfer before the
                                              pass goes idle and is only sampled
                                              periodically. 0 uses the default,
                                              negative never goes idle.
mcmd-charge-transfer-recheck (default = 100)  How often to sample the pass once
                                              idle, in steps. 0 uses the default.
mcmd-autostop                (default = 0)    Stop once E_kin/E_tot exceeds the
                                              threshold below.
mcmd-autostop-threshold      (default = 0.99) That threshold.
mcmd-initial-charges         (default = 0)    Read starting charges from
                                              "charges.txt" in the working
                                              directory: two columns, atom index
                                              and charge.
mcmd-detailed-output         (default = 0)    Write the per-step analysis output
                                              (see below). Off by default.
mcmd-charge-output-stride    (default = 50)   Steps between frames in
                                              charges_over_time.bin. The text
                                              outputs stay per-step.
mcmd-collisional-ionization  (default = 0)    NOT IMPLEMENTED. Requires
                                              collisional data.

mcmd-pulse-peak-time         (default = 0)    Centre of the Gaussian pulse [fs]
mcmd-pulse-fwhm              (default = 0)    Pulse duration, full width at half
                                              maximum [fs]
mcmd-pulse-photons           (default = 0)    Total photons in the pulse
mcmd-pulse-focal-diameter    (default = 0)    Focal spot diameter [nm]
mcmd-pulse-photon-energy     (default = 0)    Photon energy [eV]
```

`mcmd-charge-transfer-idle` and `mcmd-charge-transfer-recheck` are performance
knobs: they change how often the charge-transfer pass runs, not the physics it
computes. The defaults are tuned; see "Charge transfer cost" below before
changing them.

To allow for custom parameters in GROMACS the `.tpr` format version was raised from 73 to 74 at the same time. 
This means a stock GROMACS tool would read a MolDStruct `.tpr` and
silently misparse everything after `userint4`. So when doing a MolDStruct run, use the `grompp` from the MolDStruct install to generate the `.tpr`


### mcmd-detailed-output

This one gates **all** per-step output, not just the transition log:

| file | contents |
|---|---|
| `mean_charge_vs_time.txt` | time [ps], mean charge |
| `pulse_profile.txt` | time [ps], pulse intensity |
| `electron_data.txt` | time [ps], electron density, temperature, R_g/R_g(0) |
| `charges_over_time.bin` | every atom's charge, every step |
| `electronic_transition_log.txt` | one entry per electronic transition |

With it off you get only the end-of-run `charges.txt` and
`procces_statistics.txt`, so turn it on for any run you intend to analyse.
These files can become VERY large as it output data each step.

It is off by default because it is not free. `charges_over_time.bin` is one
`real` per atom per frame and dominates everything else, so it is written
every `mcmd-charge-output-stride` steps (default 50) rather than every step:

| system | per 5000 steps, stride 1 | stride 50 |
|---|---|---|
| 1960 atoms | 40 MB | 1.2 MB |
| 10000 atoms | 200 MB | 4.4 MB |
| 50000 atoms | 1000 MB | 21 MB |

The three text files stay per-step, since together they are under 100 bytes a
step and they are what the mean-charge and pulse curves are plotted from.
**This means `charges_over_time.bin` no longer has one frame per row of
`mean_charge_vs_time.txt`**: frame *k* is step *k* x stride, so take every
stride'th row of the text file to get the matching times. Set the stride to 1
for the old frame-per-step behaviour.

`electronic_transition_log.txt` scales with events rather than steps, and the
events are heavily front-loaded - on a 25000-step test run, 5915 of its 6224
entries fell in the first 5000 steps, because photoionization tracks the pulse
and the Auger cascades follow within a few hundred steps.

All five files are written through buffered handles and flushed every 1000
steps, so they can be tailed during a run without costing a file open per step.



### Reproducible runs

The Monte Carlo draws are seeded from the clock and the process id, so every
run is a different realisation, which is what you want for ensemble statistics.
Setting the environment variable `GMX_MCMD_SEED` pins the seed instead, making
a run exactly reproducible:

```
GMX_MCMD_SEED=12345 path/to/gromacs/bin/mdrun -deffnm explode -nt 1 -ionize
```

This is intended for debugging and regression testing. Do not set it for
production ensembles, or every member will be the same trajectory.

### Continuing a run

Restarting from a GROMACS checkpoint (`mdrun -cpi`) is supported. At the end of
every run the module writes `simulation_output/mcionize_state.bin`, holding the
electronic state, the free-electron plasma and the reference radius of gyration
the sample expansion is measured against. That file is read back on restart, so
the electronic dynamics continue rather than resetting to the ground state.

The file must be present to restart; the reference radius of gyration in
particular cannot be recovered from an already-expanded structure, and silently
re-deriving it would rescale the plasma volume.


## Regression test

`test/regression.sh` runs the shipped example for a few hundred steps with the
Monte Carlo seed pinned and checks three things:

1. the four event counts match `test/reference_<steps>.txt` exactly
2. net charge is conserved: `sum(charges) == photoionizations + Auger`
3. a run without `-ionize` writes no ionization output

```
test/regression.sh                     # against the stored reference
test/regression.sh --steps 1000        # longer run (needs its own reference)
test/regression.sh --update            # regenerate the reference
test/regression.sh --bindir /path/bin  # test an installed build
```

Exit status is 0 on pass, 1 on failure, and a failing run keeps its working
directory for inspection.

Check (1) is the point of the test. The module is stochastic, but
`GMX_MCMD_SEED` makes it deterministic, so any change that alters the physics
moves those numbers. It is only meaningful against a reference generated by
the same build on the same machine -- a different compiler or optimisation
level can shift the trajectory without anything being wrong. Regenerate with
`--update` when you change the physics deliberately, and say so in the commit.

Check (2) holds on any build: every electron removed from an atom is either a
photoionization or an Auger event, and charge transfer only moves charge
between atoms, so a failure there is a real bug rather than a platform
difference. It is what caught the charge-transfer accumulator bug.

## Supplying atomic data

To run the model we must supply it with atomic data, this includes energy levels and transition rates between states.
This data can be generated in any way you see fit as there are multiple softwares with this capabiities (CRETIN).
Formatting this data is probably the most cumbersome part of getting started with the model, but here I will do my best to guide you.
I suggest you write a script that transforms the data you generate onto this format such that itt can be recreated as some of these parameters are depedant on the choice of energy.

For each atomic species present in the system we require 2 files. They are

- `energy_levels_X.txt`
- `rate_transitions_X.txt`

where `X` is replaced with a the periodic table symbol in upper case,
Ex.

- `H` - Hydrogen
- `C` - Carbon
- `O` - Oxygen
- `FE` - Iron.

You only need to supply the species that are in your system.
The files must be placed in a folder called `Atomic_data` which is in the same location from where a simulation is run from.

Now I will go through the files and the format they must be in.
We start by defining some notation, the model uses K,L,M shells to specify electronic states. So we specify states with 3 integers. The groundstate of hydrogen would for example be `1 0 0` and so on.
I will denote a general initial state as `a b c` and a general final state as `a' b' c'`.

### `energy_levels_X.txt`

This file simply contains the energy level for different states. Each state has its own row and the format is the state `a b c` followed by the energy, thus `a b c E` with a single space as delimiter.
So for hydrogen two rows of the file might look like this:

```
0 0 0 0.0
1 0 0 13.5984
```

### `rate_transitions_X.txt`

This file contains all possible transitions for an atomic species. Each initial state has its own row and has the format

`initial state ; final state transition_rate transition-type`,

or

`a b c ; a' b' c' transition_rate transition-type ;`.

`transition-type` is 0, 1 or 2, corresponding to Auger-Meitner decay, fluorescence and
photoionization respectively. Only type 2 has its rate scaled by the pulse profile;
Auger and fluorescence rates are used as given.
Here a `;` is used to separate the initial and final states. In the case of multiple possible final states, we add one more final state, rate, and type after the first one like

`a b c ; a' b' c' transition_rate' transition-type' ; a'' b'' c'' transition_rate'' transition-type'' ;`.

This pattern continues for more final states.
We look at hydrogen again. These are the first two rows of the `rate_transitions_H.txt`
shipped with the example, generated for 1000 eV:

```
1 0 0 ;0 0 0 1.1527211049935397e-23 2 ;0 1 0 4.903126534799968e-19 1 ;0 0 1 9.318532169694317e-20 1 ;
0 1 0 ;0 0 0 3.5273211902944603e-25 2 ;1 0 0 0.0005225303312029761 0 ;0 0 1 7.547893250471743e-19 1 ;
```

Note that every final state carries its own rate _and_ its own type.

And that is that!
Remember that some of these parameters are dependant on the photon energy.
How we supply these files to the simulation is covered in the next section.

## Running a simulation

To run a simulation we follow the exact same steps as for a normal GROMACS simulation up we call `mdrun`.
We need to run it with the `-ionize` flag, and make sure that we have a the folder `Atomic_data` containing the `energy_levels_X.txt` and `rate_transitions_X.txt`.
The folder must be present in the same directory from where we run the simulation.

## Output

### MD output

All standard GROMACS output like the .trr and .edr files are still given as output.

### Additional output

Along with the normal MD output that a GROMACS sim would give, when `mcmd-detailed-output` is set to 1,
then additional output will be given in `/path_to_sim-directory/simulation_output` in the form of 4 `.txt` files.
These files can become VERY large as it output data each step.

### `electron_data.txt`

Contains information about electrons and other ionization goodies for each timestep.
Written as plain whitespace-separated columns, so it loads directly with e.g. `numpy.loadtxt`.

Column 1: time in picoseconds.

Column 2: electron density in electrons/cm^3.

Column 3: electron temperature in eV.

Column 4: volume expansion (current radius of gyration / initial radius of gyration)

Earlier versions of this manual listed a Debye length column. It was never
computed by the code and has been removed from the documentation.

### `mean_charge_vs_time.txt`

Contains information about the mean charge of the systems for each timestep.

Column 1: time in picoseconds.

Column 2: mean charge of the system in elementary charge units e.

### `pulse_profile.txt`

Contains information about the intensity of the laser pulse at each timestep.

Column 1: time in picoseconds.

Column 2: X-ray flunece [photons/nm^2].

### `electronic_transitions_log.txt`

This file logs all electronic transitions and at what time they occur in ps, as well as the current Monte Carlo timestep.

All observables are in the same units as the other files.

## Example

`example/` contains a complete worked run: hen egg-white lysozyme (PDB 1AKI,
1960 atoms) in vacuum, exploded by a 1e11-photon, 2.94 fs pulse at 1000 eV.

```bash
cd example
# edit run_example.sh to point at your gromacs bin directory, then
bash run_example.sh
```

which is just:

```bash
pdb2gmx -f 1aki.pdb -ff "charmm27" -water none
grompp  -f exp.mdp -c conf.gro -p topol.top -o explode.tpr
mdrun   -deffnm explode -v -nt 1 -ionize
```

Add `-pd` and raise `-nt` to use more cores (`-nt 8 -pd`) - `-pd` is required
whenever `-nt` > 1, because domain decomposition renumbers atoms into
rank-local indices and the charge-transfer pass needs a global search over all
positions. It is refused with a clear error rather than producing wrong
answers. Results are independent of the number of ranks.

The `-ionize` flag must be given, or none of the ionization code runs and the
simulation behaves as unmodified GROMACS 4.5.4.

The atomic data in `example/Atomic_data` is generated for 1000 eV, matching
`mcmd-pulse-photon-energy` in `exp.mdp`. Changing the photon energy requires
regenerating it (see "Supplying atomic data"); the other parameters can be
played with freely.

30000 steps at dt = 1 as covers 30 fs, with the pulse peaking at t = 0. On a
completed run every atom ends up at least singly charged:

```
min 1   mean 3.23   max 14   neutral atoms: 0/1960

```

The exact numbers vary run to run, because the Monte Carlo draws are seeded
from the clock - set `GMX_MCMD_SEED` to pin them.

Output lands in `simulation_output/` (see above); `charges.txt` holds the
final per-atom charges and `procces_statistics.txt` the event counts.

## Limitations

### Atomic species

In the way that the code is written it is limited to selection of atomic species.
At the moment the list consists of

- Hydrogen
- Carbon
- Nitrogen
- Oxygen
- Flourine
- Sodium
- Magnesium
- Silicon
- Phosphorus
- Sulfur
- Chlorine
- Calcium
- Iron
- Nickel

Note that an element appearing in this list only means the code knows its mass and
ground-state configuration. You still have to supply the corresponding
`energy_levels_X.txt` and `rate_transitions_X.txt` for every species in your system.

With a bit of programming knowledge this can easily be extended: the element
table lives in `gromacs-4.5.4-MCMD/src/kernel/atomic_data.c`, in the
`mcmd_elements[]` array near the top. Add a row with the rounded atomic mass,
the periodic-table symbol in upper case and the K/L/M ground-state occupations,
and increase `MCMD_NUM_ELEMENTS` in `atomic_data.h` to match. You then also
need to supply the corresponding data files.

### Parallelisation

The simulation runs on multiple cores with `-nt N -pd` (see above). The
molecular dynamics itself parallelises well; the ionization module runs on the
master rank, so speedup saturates once that serial section dominates. Expect a
useful gain up to roughly 8 ranks, less for small systems.

Domain decomposition is not supported, and neither is any parallelisation over
GPUs.

For ensemble statistics it is still usually more efficient to run many
independent single-core simulations at once than to parallelise one.

### Charge transfer cost

Charge transfer considers every pair of atoms within a 1 nm search radius. The
search uses a cell list, so its cost grows roughly linearly with system size
rather than quadratically, but it still runs on the master rank and is the
largest single piece of the ionization module. If it matters less to you than
the runtime, `mcmd-charge-transfer = 0` turns it off; measured on a 10000-atom protein it
was about 40% of the wall clock.

Several quantities the module needs are worked out once at setup rather than
re-derived every step: each atom's element (masses do not change), a
direct-index table from a K/L/M state to its entry in the rate data (replacing
a scan over every state of the element -- 324 of them for sulphur), and
scratch buffers for the Monte Carlo kernel, which previously allocated and
freed twice per atom per event. Each atom's state energy and outermost
occupied shell are cached and updated only when its configuration changes,
instead of being recomputed for the whole system at the start of every
transfer pass. Together these were worth about 8% of total run time on the
25000-step test, reproducing the previous results bit for bit.

The 1 nm radius is only a pre-filter; whether a transfer actually happens is
decided by the over-the-barrier criterion, which in practice is satisfied only
well inside a nanometer. If a run ever uses a pair close to that radius, a
warning is printed at the end telling you to raise `MCMD_CT_CUTOFF_NM` in
`src/kernel/mcionize.c`.

Charge transfer also dies out well before a long run does, so once nothing has
happened for `mcmd-charge-transfer-idle` steps (default 2000) the pass stops running every step
and is sampled every `mcmd-charge-transfer-recheck` steps instead (default 100). It is never
switched off outright, because Auger decay can create a new donor at any time
and the periodic check picks that up. On a 25000-step run this cut about 19%
off the wall clock while reproducing the un-skipped run bit for bit, on three
separate seeds. Setting `mcmd-charge-transfer-idle` to a negative value disables the skip and
runs the pass on every step; the end-of-run message reports what fraction of
steps it actually ran on.

Be careful lowering `mcmd-charge-transfer-idle`. The default window is generous on purpose: late
transfers are sparse rather than absent, and a window shorter than the gaps
between them clips real transfers. On the 25000-step test, `mcmd-charge-transfer-idle = 2000`
reproduces an unskipped run bit for bit, but `mcmd-charge-transfer-idle = 200`
with `mcmd-charge-transfer-recheck = 10` does not -- it changes the transfer and Auger counts, because
the run goes idle during a lull, misses a transfer, and diverges from there.
The saving over the default was under 2%, so there is nothing to gain by
tightening it.

Note that this is keyed on whether transfers are still happening, not on how
long ago the pulse was. Transfers continue long after the pulse ends: measured
on that run, 7.5% of them happened more than five pulse widths after the peak,
and 0.7% more than fifty, by which point the X-ray intensity was two hundred
orders of magnitude below its maximum.

### Systems blowing up (Too much!)

For high ionization we get huge forces, this can make the numerical integration unstable.
If you suspect this check the kinetic and potential energy of the system. As long as they look relativly smooth it should be okay.
The work around is usually to lower the stepsize. We recommend a step size of 1 as.

### Compatibility

Has only been tested on Linux systems.


### Renamed parameters
Previous versions used the default GROMACS userints and userreals. Current version uses custom parmeter names. These old ones are kept here as a reference.
These options used to be the generic `userint*`/`userreal*` slots. `grompp`
rejects an `.mdp` that still uses the old names and prints the replacement for
each one, because two of them also changed meaning:

| old | new |
|---|---|
| `userint1` | removed - `-ionize` alone enables the altered force field |
| `userint2` | `mcmd-charge-transfer` |
| `userint3` | `mcmd-autostop` |
| `userint4` | removed - use `mdrun -cpi`, see "Continuing a run" |
| `userint5` | `mcmd-detailed-output` |
| `userint6` | `mcmd-collisional-ionization` |
| `userint7` | `mcmd-charge-transfer-idle` |
| `userint8` | `mcmd-charge-transfer-recheck` |
| `userint9` | `mcmd-initial-charges` |
| `userreal1` | `mcmd-pulse-peak-time` - **now fs**, multiply by 1000 |
| `userreal2` | `mcmd-pulse-photons` |
| `userreal3` | `mcmd-pulse-fwhm` - **now the FWHM in fs**, multiply by 2354.82 |
| `userreal4` | `mcmd-pulse-focal-diameter` |
| `userreal5` | `mcmd-pulse-photon-energy` |
| `userreal6` | `mcmd-autostop-threshold` |