# Path to gromacs bin (Change to install location)
path='/path/to/gromacs/bin'

$path/pdb2gmx -f 1aki.pdb -ff "charmm27" -water none

# Configure and run explosion sim
$path/grompp -f exp.mdp -c conf.gro -p topol.top -o explode.tpr
# 4 threads (-pd is required with more than one); about 10 minutes for the 30 fs
$path/mdrun -deffnm explode -v -nt 4 -pd -ionize



