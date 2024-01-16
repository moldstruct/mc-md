# Hybrid Monte Carlo / Molecular Dynamics
This is a Hybrid Monte Carlo / Molecular Dynamics (MC/MD) model that is part of the **MOLDSTRUCT** toolbox. It can be used to simulate photon-matter interaction for smaller systems like single proteins. '
The code is based on a modified version of GROMACS [webpage](https://www.gromacs.org/).
With functionality close to the normal GROMACS but with some additional parameters to control ionization. Running the code requires that atomic data (more details later) is supplied for the appropriate energy.
The model is developed by the Biophysics group at Uppsala University and is published here [link to article](PUT URL HERE WHEN PUBLISHED).

# Todo 
- Write something about atomic data
- Guide to running sims
- example
   - Write guide
   - upload exampel files


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
8. Congratualtions! If everything worked out, the software should now be installed at your specified location.

If running into problems, try to search the issue as there are many GROMACS reasources out there that can help with the install process.


## List of input parameters
The parameters of the model can be set like any other MD parameters in the `.mdp`-file and can be accesed through the userints and userreals.
Userints are used to enable/disable certain features while the userreals are used for the simulation parameters.
```
userint1 - (default = 1)Alter forcefield. If set to zero the everything should run as unmodifed gromacs 4.5.4 (hopefully)
userint2 - (default = 1)Do charge transfer. Enables the charge transfer module 
userint3 - (default = 0)Use screened potential. Enables Debye shielding (Currently not implemented)
userint4 - (default = 0)Read electronic states from file. Reads electronic states from file. Useful for continued sims (Not well tested, only use if you can confirm validity)
userint5 - (default = 0)Enable logging of electronic dynamics. Big performance drop due to I/O. Do not use on distributed systems.
    
userreal1 - Peak of the gaussian [ps]
userreal2 - Total number of photons in the pulse
userreal3 - Width of the peak (sigma values of gaussian) [ps]
userreal4 - Diameter of the focal spot [nm]
userreal5 - Photon energy [eV]
```

If you want to use even more parameters for whatever reason, everything up to userint8/userreal8 is implemented and available.

## Supplying atomic data
To run the model we must supply it with atomic data, this includes 


## Running a simulation

## Output
### MD output

### Additional output
Along with the normal MD output that a GROMACS sim would give, when the userint5 is set to 1, 
then additional output will be given in /path_to_sim-directory/simulation_output in the form of 4 .txt files.
If you do not need any of the data, it is highly recommended to turn off this output as it cuts into the performance.

### electron_data.txt 

Contains information about electrons and other ionization goodies for each timestep.

Column 1: time in picoseconds.

Column 2: electron density in electrons/cm^3.

Column 3: electron temperature in eV.

Column 4: debye length in nm.

Column 5: volume expansion (current radius of gyration / initial radius of gyration)


### mean_charge_vs_time.txt 

Contains information about the mean charge of the systems for each timestep.

Column 1: time in picoseconds.

Column 2: mean charge of the system in elementary charge units e.


### pulse_profile.txt

Contains information about the intensity of the laser pulse at each timestep.

Column 1: time in picoseconds.

Column 2: X-ray intensity [unit?].


### electronic_transitions_log.txt

This file logs all electronic transitions and at what time they occur in ps, as well as the current Monte Carlo timestep.

All observables are in the same units as the other files.

## Examples
### Example 1

## Limitations
### Systems blowing up (Too much!)
For high ionization we get huge forces, this can make the numerical integration unstable. 
If you suspect this check the kinetic and potential energy of the system. As long as they look relativly smooth it should be okay.
The work around is usually to lower the stepsize. The simulations in the publication are done at 1as timestep.

### Atomic species
In the way that the code is written it currently can only handle a few atomic species.
At the moment the list consists of
- Hydrogen
- Carbon
- Nitrogen
- Oxygen
- Phosphorus
- Iron.

With a bit of programming knowledge this can easily be extended, look in the `md.c` file in `gromacs-4.5.4-MCMD/src/kernel`.



