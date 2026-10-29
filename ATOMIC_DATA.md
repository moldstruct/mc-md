# Atomic data format (`Atomic_data/`)

By default (`mcmd-atomic-data = hydrogenic`, since tpr version 80) the ionization module computes the atomic physics
at start-up and needs no input file: a screened-hydrogenic shell model (Scott & Hansen, HEDP 6, 39 (2010), screening
constants of Marchand, Caillé & Lee, JQSRT 43, 149 (1990), energies scaled to the NIST ionization energies) at the
photon energy of the pulse, for every `K L M` state of H to Ar except He (an atom of mass ~4 is refused: far more likely a heavy hydrogen). It is the model of the CRETIN screened-hydrogenic
data the published results used (rates within ~10 %). It writes what it computed to `simulation_output/Atomic_data/`
in the format below; run with `mcmd-atomic-data = files` on those files, the run is identical. Lotz's formula
supplies the electron-impact ionization; `mcmd-coll-data = files` and collisional excitation need files.

With `mcmd-atomic-data = files` the module reads the atomic physics of every element from plain-text tables in
`./Atomic_data/` of the run directory. This document specifies those tables completely, so that they can be produced
by any atomic-structure code (XATOM and CRETIN models have generators, see the last section).
`tools/check_atomic_data.py` checks a directory against everything below.

The reader is `src/kernel/atomic_data.c`; how each quantity enters the physics is in `src/kernel/mcionize.c`.
Written 2026-10-28 against source 56f1b52a.

## Atomic states

A state is the **shell occupation** of an atom, three integers `K L M`: the number of electrons in the n = 1, 2, 3
shells (capacities 2, 8, 18). There are no subshells, no terms or levels within a configuration, and no n ≥ 4: a
code that resolves subshells or levels must sum or average over them per shell configuration. Excited
configurations are separate states (`1 5 1` is carbon with a K hole and an electron promoted to M).

Elements are identified by their rounded atomic mass, through the table `mcmd_elements[]` in `atomic_data.c`
(mass, symbol, ground state): H 1, C 12, N 14, O 16, F 19, Na 23, Mg 24, Si 28, P 31, S 32, Cl 35, Ca 40, Fe 56,
Ni 59, I 127 (deuterium uses hydrogen's files). A file is `<table>_<SYMBOL>.txt` with the symbol as written there
(`energy_levels_C.txt`, `rate_transitions_MG.txt`, `rate_transitions_CL.txt`). A new element needs a row there.

## Which files are needed (`mcmd-atomic-data = files`)

| file | needed when | what for |
|---|---|---|
| `energy_levels_X.txt` | always | outermost binding energy per state |
| `rate_transitions_X.txt` | always | photoionization, Auger, fluorescence |
| `total_energies_X.txt` | `mcmd-charge-transfer = 1` or `mcmd-collisions = 1` (both default) | electron energies, charge-transfer energy matching, collisional transition energies |
| `collisional_ionization_table_X.txt` **or** `collisional_parameters_X.txt` | only with `mcmd-coll-data = files`; the table is used where present. With the default `mcmd-coll-data = lotz` the ionization comes from Lotz's formula and neither is read | electron-impact ionization (and, by detailed balance, three-body recombination) |
| `statistical_weight_X.txt` | optional (`mcmd-collisions = 1`); without it the weights are computed, `C(2,K) C(8,L) C(18,M)` | detailed balance, Lotz channels |
| `collisional_excitation_table_X.txt` **or** `collisional_excitation_X.txt` | optional; used with `mcmd-coll-excitation = 1` (off by default); the table is used where present | collisional excitation / de-excitation |

Only `rate_transitions_X.txt` depends on the photon energy (the photoionization cross-sections): regenerate it for
every `mcmd-pulse-photon-energy`. The other tables are properties of the atom.

## General rules

- One record per line, whitespace-separated; `;` separates the initial state from the final states in the
  transition tables. No comments, no header.
- **Lines must be shorter than 1024 characters** (the reader's `fgets` buffer): a longer line is split and misread.
- Energies are in **eV**, rates in **s⁻¹**, cross-sections in **cm²**.
- The value `-1` is the reader's "not found" marker and must not appear as an energy.
- Where a state appears on more than one line of a table, the first line wins.

## `energy_levels_X.txt`

```
K L M E                (carbon)
2 4 0 11.2603
1 1 0 86.2745
0 2 0 100000000000000.0
```

`E` (eV) is the binding energy of the **outermost** electron of the state: the energy to remove one electron from
its highest occupied shell. **Autoionizing states** (excitation energy above the ionization potential of their own
charge state, so no positive binding energy) get `1e14`. Used for the charge-transfer donor energy, and as the
fallback for the photoelectron energy when `total_energies` is absent. Every state of `rate_transitions` should
have an entry (`0 0 0` gets `0.0`).

## `total_energies_X.txt`

```
K L M T                (carbon)
2 4 0 1022.1926
1 4 0 723.6903
0 0 0 0.0
```

`T` (eV) is the **total binding energy** of the state, the energy to remove all its electrons:
`T = Σ ionization potentials up to this charge state − the state's excitation energy`. Defined for autoionizing
states too. Must cover every state reachable through the transition tables (a missing one is fatal in charge
transfer and at start-up with collisions). For every non-autoionizing state it must reproduce `energy_levels`:
`E(c) = T(c) − T(c minus one electron from its outermost shell)`.

Used for:
- the kinetic energy of emitted electrons (which decides escape or trapping):
  photoelectron `E_ph − (T_from − T_to)`, Auger electron `T_to − T_from` (both must be > 0);
- the binding energy of any shell for charge transfer: `ε(c, s) = T(c + e_s) − T(c)`;
- collisional transition energies `ΔE = T_i − T_f`.

## `rate_transitions_X.txt`

One line per initial state, then one `;`-terminated token per final state with **five** fields:

```
K L M ;K' L' M' rate type ;K'' L'' M'' rate type ;...        (oxygen, 600 eV)
1 1 0 ;1 0 0 7.416614e-21 2 ;0 1 0 0.0 2 ;2 0 0 7.291753e+11 1 ;
```

| type | process | `rate` | electrons |
|---|---|---|---|
| 0 | Auger–Meitner decay | rate, s⁻¹ | final = initial − 1 |
| 1 | fluorescence (radiative decay) | rate, s⁻¹ | final = initial (one electron moves inward) |
| 2 | photoionization at the photon energy | cross-section, cm² | final = initial − 1 |

- **The type column is mandatory.** It cannot be inferred from the state change: many Auger and photoionization
  transitions have the same occupation change, and a guessed type puts Auger rates under the pulse.
- Photoionization rate in the run = `σ × I(t)`, `I(t)` the pulse's photon flux in photons cm⁻² s⁻¹ (Gaussian in
  time; peak `N_ph / (π (d/2)²) / (σ_t √(2π))` with the focal diameter `d` in cm).
- **Every state reachable through any transition must have a line of its own**, also when it has no transitions
  (`K L M ;`); reaching a state without a line stops the run. Only the bare nucleus `0 0 0` may be absent.
- The element's ground state must have a line.
- A zero rate is an inactive entry and never fires (the CRETIN-derived files end with `0 0 0 ; 0 0 1 0.000 0;`).
- Two entries to the same final state are two channels; their rates add.

## `statistical_weight_X.txt` (optional)

The weight of a shell configuration counts the ways to place the electrons among the 2n² spin-orbitals of each
shell, `g = C(2,K) C(8,L) C(18,M)`, whatever atomic model produced the rest of the data; summing the subshell or
level weights of a configuration gives the same number. Since 2026-10-29 mdrun computes it when the file is absent.
A file, where present, is used as given, with a note on stderr if it departs from the formula.

```
K L M g                (carbon)
2 4 0 70
1 1 0 16
```

The degeneracy of the configuration: `g = C(2,K) C(8,L) C(18,M)` (the CRETIN files print it to 6 significant digits).
Where the file is given, it must cover every state in the collisional tables.

## `collisional_parameters_X.txt` (electron-impact ionization)

```
K L M ; K' L' M' c0 c1 c2 c3 s0 ; ...        (carbon)
2 4 0 ; 2 3 0 5.3345 -8.52561 3.1911 0 7.30793 ; 1 4 0 3.25123 -5.2269 1.97567 0 6.80217 ;
```

Each final state has **one electron fewer** (ionization only; this file does *not* list recombination or excitation).
The five numbers are CRETIN's Sampson-fit parameters. With `ΔE = T_i − T_f` (from `total_energies`, not from this
file) and `b = ΔE / T_e`:

```
q(b)  = c0 + c1 √(π b) erfcx(√b) + (c2 b + c3) e^b E1(b)
k_ion = 1.0899e-6 · s0 · q(b) / (ΔE √T_e)          cm³ s⁻¹ (T_e, ΔE in eV)
rate  = n_e k_ion e^(−b)                          s⁻¹
```

Three-body recombination (`mcmd-coll-recombination = 1`, off by default) is the inverse by detailed balance:
`rate = n_e² k_ion (g_i/g_f) · 1.656e-22 · T_e^(−3/2)`, attached to the final state. Every state used must have a
`rate_transitions` line, a total energy and a weight, and `ΔE > 0`, or the run stops at start-up.

## Built-in electron-impact ionization: `mcmd-coll-data = lotz` (default)

By default (since 2026-10-28, tpr version 79) no collisional data file is needed: every state of `rate_transitions`
loses one electron from each occupied shell, wherever that state exists, with Lotz's Maxwellian rate coefficient
(1967/68, b = 0)

```
<σv> = 6.7e-7 · a · q · T_e^(−3/2) · E1(P/T_e) / (P/T_e)      cm³ s⁻¹,  a = 4.5 (×1e-14 cm² eV²)
```

q the occupation of the shell, P = ΔE = T_i − T_f from `total_energies` (T_e, P in eV). Recombination follows by
detailed balance as for the other sources. For H C N O S this gives exactly the transitions of the CRETIN fit files,
with coefficients 13–25 % above them (10–300 eV); on lysozyme at 600 eV 10⁷ the explosion changes by −2.5 % in ion
KE and within the seed noise otherwise. `mcmd-coll-data = files` uses the data files below instead (a table where
present, else the fit file); tprs older than version 79 run as `files`. Collisional excitation (off by default) always
comes from files.

## `collisional_ionization_table_X.txt`, `collisional_excitation_table_X.txt` (tabulated, generator-independent)

The alternative to the two fit files below (2026-10-28): rate coefficients on a temperature grid, which any atomic
code, or a formula (Lotz, BEB, …), can supply. Where a table is present it is used and the fit file is not read.

```
T_e t1 t2 ... tn
K L M K' L' M' k1 k2 ... kn
...
T_e 0.1 0.10592537251772889 0.11220184543019636 ...                       (carbon)
2 4 0 2 3 0 2.6588838065e-57 1.4886556270e-54 ...
```

- First line: `T_e` and the grid temperatures in eV, positive and increasing; write them at full precision (the code
  uses them in ΔE/T). The grid should cover the cloud's temperatures (~1 eV to ~1 keV); outside it the coefficients
  are held at the edge, and the run counts the steps where that happened.
- Then one transition per line (no `;`, any length): initial and final occupations and the full Maxwellian rate
  coefficient ⟨σv⟩ (cm³ s⁻¹, Boltzmann factor included) at each grid temperature. Ionization table: one electron fewer
  in the final state; excitation table: same electron count, lower → upper.
- Any exponent is allowed (`3.1e-2650`): at low T_e the coefficients carry exp(−ΔE/T_e), far below the smallest
  double. A `0` means "too small to write"; the code then extrapolates from the neighbouring values.
- **How it is used:** with ΔE = T_i − T_f from `total_energies`, the code forms g = ln⟨σv⟩ + ΔE/T_e (the log of the
  coefficient without its Boltzmann factor, smooth in T_e) and interpolates it linearly in ln T_e. Ionization:
  rate = n_e e^g e^(−ΔE/T_e); recombination n_e² e^g (g_i/g_f) 1.656e-22 T_e^(−3/2); excitation n_e e^g e^(−ΔE/T_e);
  de-excitation n_e e^g (g_lower/g_upper): the same as for the fits, with e^g in place of the fit's coefficient.
- **Consistency:** the table's threshold must be the total-energy ΔE; if it is not, g is not smooth at low T_e, which
  the checker reports.
- `tools/make_collisional_tables.py <Atomic_data>` writes these tables from the CRETIN fit files (40 points per decade,
  0.1 eV – 100 keV: the interpolation is within 1·10⁻⁴ of the fits from 0.5 eV to 30 keV).

## `collisional_excitation_X.txt` (optional)

```
K L M ; K' L' M' c0 c1 c2 c3 a0 ; ...        (oxygen)
0 5 3 ; 0 4 4 0.290372 -0.65319 6.60995 -8.79012 4.7719 ;
```

Lower → upper state, **same electron count**, `ΔE = T_lower − T_upper > 0`. CRETIN `colex2` fit of the collision
strength, integrated over a Maxwellian (`A = a0 + 1`):

```
q(b)  = c0 e^b E1(b) + c1 + c3 b / A + (c2 − c3 b) b e^(Ab) E1(Ab)
k_up  = 8.01e-8 q(b) / √T_e                       cm³ s⁻¹;  rate = n_e k_up e^(−b)
k_down = k_up (g_lower / g_upper)                 (detailed balance)
```

## Producing the tables with another code

| table | from an atomic-structure code | notes |
|---|---|---|
| `energy_levels`, `total_energies` | configuration energies (XATOM: from the orbital eigenvalues, see below; FLYCHK: level energies) | sum or average subshell / level structure per K L M configuration; mark autoionizing states `1e14` in `energy_levels` |
| `statistical_weight` | combinatorics | `C(2,K) C(8,L) C(18,M)` |
| `rate_transitions` | photoionization cross-sections at the photon energy, Auger and fluorescence rates (XATOM and FLYCHK provide all three) | cm² for photo, s⁻¹ for the others; the type column; a line for every reachable state |
| collisional ionization | nothing: built in (`mcmd-coll-data = lotz`, the default); for other data, `collisional_ionization_table` (⟨σv⟩(T_e) from FLYCHK, FAC, BEB, …) with `mcmd-coll-data = files` | threshold = the total-energy ΔE |
| `collisional_excitation_table` (only if excitation is switched on) | ⟨σv⟩(T_e) per transition (FAC distorted-wave, FLYCHK, or van Regemorter with f from the fluorescence rates in `rate_transitions`, which cover every CRETIN excitation) | the CRETIN fit files (`collisional_parameters`, `collisional_excitation`) are only needed where no table is given |

Two generators write every table: `tools/make_atomic_data_xatom.py OUTDIR ENERGY` (XATOM, open source) and `tools/make_atomic_data_cretin.py MODELS OUTDIR ENERGY` (CRETIN `.z01` models).

**XATOM** (open source, Hartree–Fock–Slater) gives all of `rate_transitions`; what worked (2026-10-28):
- States: every K L M up to the ground-state occupation (photoionization, Auger, fluorescence and charge transfer
  never leave that set); each computed as its aufbau subshell configuration (2s before 2p, 3s before 3p), the
  subshell final states summed onto K L M.  Only elements whose ground state lies in n ≤ 3 (H–Cl): the module counts
  4s electrons in M.
- `xatom -Z <Z> -config <subshells> -pcs -PE <E> -decay -dat <file>` per state: `P` and `V` (core and valence
  photoionization, a₀² → cm²) are type 2, `A` type 0, `F` type 1 (rates in atomic units → s⁻¹); `E` (photoexcitation)
  has no channel here.  One-electron states with `-exact`.
- Energies from the orbital eigenvalues, `T(c) = T(c minus its outermost electron) − ε_outermost(c)`: XATOM's total
  energies are tens of eV off (N⁵⁺ 1s²: 1267 eV against 1218), its eigenvalues are what it prices transitions with.
  Drop channels these energies close (a few near-threshold photoionizations) and intra-shell ones (2p → 2s).
- Its photoionization matches Scofield's tables to < 0.3 % (the CRETIN models' does not: K +7…+23 %, valence
  3–12× low at 600 eV), and its S K width and fluorescence yield match Krause.  Weak point: neutral C and N
  outermost binding energies ~20 % low (open 2p shells).

Checking a directory:

```
tools/check_atomic_data.py Atomic_data --photon-energy 600 [--coll-data lotz|files] [--no-ct] [--no-collisions]
```

It reports what the reader would stop on (or misread) as errors and what is legal but suspicious as warnings, and
exits 1 on errors. On the current 600 eV tables (H C N O S) it passes; the 2024 tables of the published version fail
(no type column, no total energies).

## History

- Tables before 2026-08 had four fields per transition (no type) and no `total_energies`; their fluorescence rates
  were 10⁻¹⁵–10⁻⁹ of the physical values (effectively no fluorescence) and their energy levels covered only the
  states reachable by photoionization and Auger from the ground state. The reader refuses them.
- `collisional_parameters` used to mix ionization, recombination and excitation entries; recombination is now
  derived by detailed balance and excitation has its own file.
