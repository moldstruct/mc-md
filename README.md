# mc-md
Hybrid Monte Carlo / Molecular Dynamics


Important parameters
    userint1 - Alter forcefield. If set to zero the everything should run as unmodifed gromacs 4.5.4 (hopefully) (default = 1)
    userint2 - Do charge transfer. Enables the charge transfer module (default = 1)
    userint3 - Use screened potential. Enables Debye shielding (default = 1)
    userint4 - Read electronic states from file. Reads electronic states from file. Useful for continued sims (default = 0)
    userint5 - Enable logging of electronic dynamics,writes a bunch of useful information, big performance drop due to I/O. (default = 0)

    The userreal are mostly the FEL parameters
    userreal1 - Peak of the gaussian pulse in ps
    userreal2 - Total number of photons in the pulse
    userreal3 - Width of the peak in ps. This is the sigma value of the guassian.
    userreal4 - Diameter of the focal spot (nm) 