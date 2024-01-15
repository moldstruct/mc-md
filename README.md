# mc-md
Hybrid Monte Carlo / Molecular Dynamics

## Important parameters
    userint1 - (default = 1)Alter forcefield. If set to zero the everything should run as unmodifed gromacs 4.5.4 (hopefully)
    userint2 - (default = 1)Do charge transfer. Enables the charge transfer module 
    userint3 - (default = 1)Use screened potential. Enables Debye shielding 
    userint4 - (default = 0)Read electronic states from file. Reads electronic states from file. Useful for continued sims 
    userint5 - (default = 0)Enable logging of electronic dynamics. Big performance drop due to I/O. Do not use on distributed systems.
    
    The userreal are mostly the FEL parameters
    userreal1 - Peak of the gaussian pulse in ps
    userreal2 - Total number of photons in the pulse
    userreal3 - Width of the peak in ps. This is the sigma value of the guassian.
    userreal4 - Diameter of the focal spot (nm) 
    userreal5 - Photon energy [eV]

If you want to use even more parameters for whatever reason, everything up to userint8/userreal8 is implemented and available.


## Extra output
Along with the normal MD output that a GROMACS sim would give, when the userint5 is set to 1, 
then additional output will be given in /path_to_sim-directory/simulation_output in the form of 4 .txt files.

### electron_data.txt  @@@ TODO: SPECIFY UNITS

Contains information about electrons and other ionization goodies for each timestep.

Column 1: time in picoseconds.

Column 2: electron density in (electrons per nm^3?).

Column 3: electron temperature in _.

Column 4: debye length in (nm?).

Column 5: rg_factor???.


### mean_charge_vs_time.txt 

Contains information about the mean charge of the systems for each timestep.

Column 1: time in picoseconds.

Column 2: mean charge of the system in elementary charge units e.


### pusle_profile.txt

Contains information about the intensity of the laser pulse at each timestep.

Column 1: time in picoseconds.

Column 2: X-ray intensity


### eletronic_transitions_log.txt

This file logs all electornic transitions and at what time they occur in ps, as well as the current Monte Carlo timestep.

All observables are in the same units as the other files.


