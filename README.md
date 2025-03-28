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
3. Place the `gromacs-4.5.4-MCMD` folder somewhere appropriate.
   This is only the installation files, you will choose later on where you want to install the software.
4. Go into `gromacs-4.5.4-MCMD` and create a new directory called `build` (or whatever you want).
5. Go into the newly created directoy. 
6. Run `cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/install/location`, where `/path/to/install/location` is replaced to where the software will be installed.
7. Run `make install`
8. ==Congratualtions!== If everything worked out, the software should now be installed at your specified location.

If running into problems, try to search the issue as there are many GROMACS reasources out there that can help with the install process.


## List of input parameters
The parameters of the model can be set like any other MD parameters in the `.mdp` file and can be accesed through the userints and userreals.
Userints are used to enable/disable certain features while the userreals are used for the simulation parameters.
```
userint1 - (default = 1) Alter forcefield. If set to zero the everything should run as unmodifed gromacs 4.5.4 (hopefully)
userint2 - (default = 1) Do charge transfer. Enables the charge transfer module 
userint3 - (default = 0) Autostop simulations after energy threshold is reached E_kin/E_tot > threshold. (Threshold set by userreal6)
userint4 - (default = 0) Read electronic states from file.  Useful for continued sims (Experimental)
userint5 - (default = 0) Enable logging of electronic dynamics. Big performance drop due to I/O.
userint6 - (default = 0) Enable collisional ionization. (Not implemented yet). Requires collisional data.
userint9 - (default = 0) Read charges from a file. File must be names "charges.txt" and be placed in the same directory from where mdrun is called. File consists of two rows, atom index and charge.
    
userreal1 - Peak of the gaussian [ps]
userreal2 - (default = 0.0) Total number of photons in the pulse
userreal3 - (default = 0.0) Width of the peak (sigma values of gaussian) [ps]
userreal4 - (default = 0.0) Diameter of the focal spot [nm]
userreal5 - (default = 0.0) Photon energy [eV]
userreal6 - (default = 0.99) Threshold for energy autostop (userint3)
```

If you want to use even more parameters for whatever reason, everything up to userint8/userreal8 is implemented and available.

## Supplying atomic data
To run the model we must supply it with atomic data, this includes energy levels, collisional parameters, transition rates between states, and statistical weights.
This data can be generated in any way you see fit as there are multiple softwares with this capabiities. 
Formatting this data is probably the most cumbersome part of getting started with the model, but here I will do my best to guide you. 
I suggest you write a script that transforms the data you generate onto this format such that itt can be recreated as some of these parameters are depedant on the choice of energy.

For each atomic species present in the system we require 2 files. They are
- `energy_levels_X.txt`
- `rate_transitions_to_gromacs_X.txt`

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

### `rate_transitions_to_gromacs_X.txt`
This file contains all possible transitions for an atomic species. Each initial state has its own row and has the format 
`initial state ; final state transition_rate`, or 
`a b c ; a' b' c' transition_rate ;`.
Here a `;` is used to separate the initial and final states. In the case of multiple possible final states, we add one more final state and rate after the first one like 
`a b c ; a' b' c' transition_rate' ; a'' b'' c'' transition_rate'' ;`.
This pattern continues for more final states.
We look at hydrogen again, the transitions `1 0 0` -> `0 0 0` and `0 1 0` -> `0 0 0`/`1 0 0`  could look like this:
```
1 0 0 ;0 0 0 6.230054227114906e-23 ;
0 1 0 ;0 0 0 1.8994473708741955e-24 ;1 0 0 0.008360485299247617 ;
```

And that is that!
Remember that some of these parameters are dependant on the photon energy.
How we supply these files to the simulation is covered in the next section.


## Running a simulation
To run a simulation we follow the exact same steps as for a normal GROMACS simulation up we call `mdrun`.
We need to run it with the `-ionize` flag, and make sure that we have a the folder `Atomic_data` containing the `energy_levels_X.txt` and `rate_transitions_to_gromacs_X.txt`.
The folder must be present in the same directory from where we run the simulation.


## Output
### MD output
All standard GROMACS output like the .trr and .edr files are still given as output.

### Additional output
Along with the normal MD output that a GROMACS sim would give, when the userint5 is set to 1, 
then additional output will be given in `/path_to_sim-directory/simulation_output` in the form of 4 `.txt` files.
If you do not need any of the data, it is highly recommended to turn off this output as it cuts into the performance.

### `electron_data.txt`

Contains information about electrons and other ionization goodies for each timestep.

Column 1: time in picoseconds.

Column 2: electron density in electrons/cm^3.

Column 3: electron temperature in eV.

Column 4: debye length in nm.

Column 5: volume expansion (current radius of gyration / initial radius of gyration)


### `mean_charge_vs_time.txt`

Contains information about the mean charge of the systems for each timestep.

Column 1: time in picoseconds.

Column 2: mean charge of the system in elementary charge units e.


### `pulse_profile.txt`

Contains information about the intensity of the laser pulse at each timestep.

Column 1: time in picoseconds.

Column 2: X-ray intensity [unit?].


### `electronic_transitions_log.txt`

This file logs all electronic transitions and at what time they occur in ps, as well as the current Monte Carlo timestep.

All observables are in the same units as the other files.

## Example
To help with getting started here is a simple example of running a simple explosion simulaiton, check the `example` folder.
In this example we will blow up a lysozyme protein, the structure is given in `1aki.pdb`.
To run the sim we must first make sure that we have all required atomic data, in this case the data is already provided in the `Atomic_data` folder this is generated for 600 eV photon energy. We also make sure that there exists an folder for the output `simulation_output`.

With the prerequisite data present, the rest of the procedure is very similar to a normal Gromacs run.

#### Generate topology
To generate a topology we can run
```
path/to/gromacs/bin/pdb2gmx -f 1aki.pdb -ff "charmm27" -water "SPC"
```

Feel free to use any other suitable forcefield, in case of a system with water, try to avoid using models with dummy-particles.

#### Configure
We call the preprocessor on our parameter file, at the end of it you can see the new parameters, feel free to play around with them to see how it influences the system. However, changing the photon energy would require new atomic data.

```
path/to/gromacs/bin/grompp -f exp.mdp -c conf.gro -p topol.top -o explode.tpr 
```
#### Run!
Finally we can run the simulation by calling the following command:
```
path/to/gromacs/bin/mdrun -deffnm explode -v -nt 1 -ionize
```
Important to note that it only works for serial simulations `-nt 1` and that the `-ionize` flag must be given to enable ionization.

There is also a simple bash-script provided for running all commands.

## Limitations
### Atomic species
In the way that the code is written it is limited to selection of atomic species.
At the moment the list consists of
- Hydrogen
- Carbon
- Nitrogen
- Oxygen
- Flourine
- Magnesium
- Phosphorus
- Sulfur
- Chlorine
- Calcium
- Iron

With a bit of programming knowledge this can easily be extended, look in the `md.c` file in `gromacs-4.5.4-MCMD/src/kernel`.
Almost at the top of the file you will find a list of masses, where you must add the mass of whatever you want to add. 
Then accordingly update the `Elements` array and the `elementConfigs` array.

### Only serial simulations
Currently the simulation can only utilize one core. One can of course run many simulations in parallel, as one might want enough simulations to evaluate statistics.


### Systems blowing up (Too much!)
For high ionization we get huge forces, this can make the numerical integration unstable. 
If you suspect this check the kinetic and potential energy of the system. As long as they look relativly smooth it should be okay.
The work around is usually to lower the stepsize. The simulations in the publication are done at 1as timestep.

