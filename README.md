# MolDStruct MC/MD

Hybrid Monte Carlo / molecular dynamics for X-ray induced ionization and Coulomb explosion of proteins and other
small systems, part of the **MolDStruct** toolbox from the Biophysics group at Uppsala University.

The model is a modified GROMACS 4.5.4. With `mdrun -ionize`, every atom carries an electronic configuration (K, L, M
shell occupations) that changes by photoionization, Auger decay, fluorescence and charge transfer, sampled by kinetic
Monte Carlo every MD step. Ionized atoms lose their bonds. The released electrons either escape or join a trapped
electron cloud, which is solved self-consistently as a Poisson–Boltzmann density on a 3D grid. The cloud screens the
ions, ionizes them further by electron impact, and evaporates as the sample expands.

**The full documentation is the reference manual, [`manual/manual.pdf`](manual/manual.pdf).** It covers the physics,
the modified force field, every option and output file, the atomic data, the numerics and the validation. This README
is the short version.

Publication: [link to article](PUT URL HERE WHEN PUBLISHED). Please also cite GROMACS 4.5 (Pronk et al.,
Bioinformatics 29, 845 (2013)), and XATOM (Jurek et al., J. Appl. Cryst. 49, 1048 (2016)) if you use data generated
with XATOM.

## Installation

The build is the standard GROMACS 4.5.4 CMake build. It needs FFTW 3 in single precision.

```bash
cd gromacs-4.5.4-MCMD
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make -j8 install
```

Write run inputs with this build's `grompp`. The `.tpr` carries the `mcmd-*` options (format 80), which stock GROMACS
cannot read.

## Quick start: the example

`example/` holds lysozyme (1AKI, 1960 atoms) in vacuum, hit by 10¹¹ photons at 1000 eV (2.94 fs FWHM, 100 nm focus).
No atomic data files are needed: they are computed at start-up.

```bash
cd example
pdb2gmx -f 1aki.pdb -ff charmm27 -water none
grompp  -f exp.mdp -c conf.gro -p topol.top -o explode.tpr
mdrun   -deffnm explode -v -nt 4 -pd -ionize
```

The run covers 30 fs at a 2 as step and takes about 10 minutes on 4 threads. The seed is new every run
(`mcmd-seed = -1`), so the numbers vary by about a percent. Roughly:

- every atom ends up ionized, with a mean charge of ~3.5 per atom;
- there are ~4000 photoionizations, ~800 Auger decays, ~2100 electron-impact ionizations and ~820 charge transfers;
- at 30 fs ~1350 electrons have escaped and ~5600 are still trapped, at ~18 eV;
- the radius of gyration has grown 3.6-fold, to ~5.7 nm.

At start-up `mdrun` prints a **run summary**, and writes the same block to `md.log`. It lists the system, the pulse,
the seed, the atomic data and every process with its setting. Settings outside the validated model are listed as
warnings, for example an MD step above 2 as, non-default solver settings, recombination, or harmonic bonds that
cannot break.

## Your own system

- Equilibrate without `-ionize`, then run in vacuum with `pbc = no`, `dt = 0.000002` (2 as), `morse = yes` (so that
  bonds can break), and the pulse options below.
- No input files are needed for the atomic data: by default they are computed for the photon energy of the pulse
  (H to Ar except He; see below).
- Run `mdrun -ionize`; use `-nt N -pd` for several threads (domain decomposition is not supported).
- With a fixed `mcmd-seed` a run is reproducible on the same number of ranks. `mdrun -cpi` continues a run from a
  checkpoint.

## Main options

All options are `mcmd-*` entries in the `.mdp` file. The pulse options are in **fs**, unlike `dt`. The defaults are
the validated model; the full list is in the manual, chapter 10.

| option                                            | default    | meaning                                                                               |
| ------------------------------------------------- | ---------- | ------------------------------------------------------------------------------------- |
| `mcmd-pulse-peak-time`, `mcmd-pulse-fwhm`         | 0, 0       | pulse centre and FWHM (fs)                                                            |
| `mcmd-pulse-photons`, `mcmd-pulse-focal-diameter` | 0, 0       | photons in the pulse, focal spot (nm)                                                 |
| `mcmd-pulse-photon-energy`                        | 0          | photon energy (eV); atomic data files must be generated for it                        |
| `mcmd-atomic-data`                                | hydrogenic | built-in atomic data (H to Ar but He); `files`: read `./Atomic_data`                  |
| `mcmd-electrons`                                  | pbgrid     | trapped-electron cloud; `none`: every electron leaves (the original MolDStruct model) |
| `mcmd-collisions`, `mcmd-coll-data`               | 1, lotz    | electron-impact ionization; Lotz built in, or `files`                                 |
| `mcmd-charge-transfer`                            | 1          | over-the-barrier charge transfer                                                      |
| `mcmd-pb-smoothing`                               | 0.25       | electron size s_e (nm); the ion kinetic energy depends on it                          |
| `mcmd-pb-long-range`                              | auto       | all pairs up to 10 000 atoms, 2 nm cutoff + grid long range above                     |
| `mcmd-field-ionization`                           | 0          | over-the-barrier ionization in the cloud's field                                      |
| `mcmd-seed`                                       | −1         | random seed; −1 draws a new one each run (printed)                                    |
| `mcmd-detailed-output`                            | 1          | time-resolved output, including `transitions.bin`                                     |

Old `userint`/`userreal` options and the `GMX_MCMD_*` environment variables of development versions are refused, with
the replacement named (manual, section 10.8 and chapter 11).

## Atomic data

By default (`mcmd-atomic-data = hydrogenic`) the energies, rates and cross sections of every shell configuration are
computed at start-up from a screened-hydrogenic shell model at the photon energy of the pulse, for H to Ar except He.
The run writes them to `simulation_output/Atomic_data/`.

It is also possible to supply your own input data with `mcmd-atomic-data = files`: they are then read from
`Atomic_data/<table>_<SYMBOL>.txt` in the run directory. Such data can for example be generated with XATOM. The format
is specified in [`ATOMIC_DATA.md`](ATOMIC_DATA.md) and in chapter 8 of the manual. Only the photoionization cross
sections depend on the photon energy.

- `tools/make_atomic_data_xatom.py OUTDIR ENERGY` generates the data with XATOM (open source).
- `tools/make_atomic_data_cretin.py MODELS OUTDIR ENERGY` generates it from CRETIN atomic models.
- `tools/check_atomic_data.py Atomic_data --photon-energy E` checks that a directory satisfies the format.

The electron processes need no input data by default either. Electron-impact ionization uses Lotz's formula
(`mcmd-coll-data = lotz`), with the binding energies from the total energies of the configurations, and the
statistical weights follow from the shell occupations, `g = C(2,K) C(8,L) C(18,M)`. Other electron-impact data
(`mcmd-coll-data = files`) and the data for collisional (de-)excitation (`mcmd-coll-excitation = 1`, off by default)
are read from `Atomic_data/` and need `mcmd-atomic-data = files`. Three-body recombination (off by default) follows
from the ionization rates by detailed balance and needs no data of its own.

Elements up to the M shell are supported: H to Ar except He, with built-in data. Ca, Fe, Ni and I need `files`.
Helium is left out on purpose: a mass of ~4 is far more likely a heavy hydrogen (`pdb2gmx -heavyh`), and is refused.

## Output

Standard GROMACS output, plus `simulation_output/`:

- `charges.txt`: final per-atom charges; `procces_statistics.txt`: event counts;
- `transitions.bin`: every electronic event; `charges_over_time.bin`: all charges every 500 steps;
- `mean_charge_vs_time.txt`, `pulse_profile.txt`, `electron_data.txt`;
- the electron cloud: `pb_log.txt`, `screening_log.txt`, `collisional_log.txt`, `energy_ledger.txt`, `pbg_box.txt`;
- `mcionize_state.bin`: the restart state.

`tools/mcmd_output.py` reads them all. `mcmd_output.py plot [dir]` draws an overview of a run, `--gif` adds an
animation, and `summary`, `transitions` and `charges` convert the bin files into text files (manual, chapter 9).

## Regression test

`test/regression.sh` runs the example for 300 steps with a fixed seed, in two cases: without the electron cloud,
and with the full model. It checks the event counts against stored references and charge conservation, and, reading
the tables the run wrote with `mcmd-atomic-data = files`, that they reproduce the run exactly and that the right files
are required. It takes under a minute. Run `--update` after a deliberate change of the physics.

## Limitations

- Shell-averaged configurations up to the M shell.
- One temperature and a classical mean field for the trapped electrons.
- Three-body recombination and collisional excitation are available but not fully validated.
- Particle decomposition only, no GPUs.
- Tested on Linux.

See the manual, appendix C, for the validation and the limits of the model.
