# Path to gromacs bin (Change to install location)
path='/home/tomas/programs/gromacs-mc/bin'

mkdir simulation_output

$path/pdb2gmx -f 1aki.pdb -ff "charmm27" -water "SPCE"

# Configure and run explosion sim
$path/grompp -f exp.mdp -c conf.gro -p topol.top -o explode.tpr
$path/mdrun -deffnm explode -v -nt 1 -ionize



