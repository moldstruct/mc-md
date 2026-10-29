#!/usr/bin/env python3
"""Generate MolDStruct MC/MD Atomic_data files from CRETIN .z01 atomic models.

CRETIN is not open source; make_atomic_data_xatom.py makes the same files from XATOM (open source).

For each model it writes the files the ionization module reads:

    energy_levels_X.txt    binding energy of each configuration's outermost
                           electron, or 1e14 where the level lies above the
                           ionization limit of its own stage (autoionizing)
    total_energies_X.txt   total electronic binding energy of each
                           configuration, defined for autoionizing levels too
    rate_transitions_X.txt photoionization, fluorescence and Auger rates out of
                           each configuration
    collisional_parameters_X.txt
                           electron-impact ionization fit coefficients out of
                           each configuration ('sampson ionize' block)
    statistical_weight_X.txt
                           statistical weight g of each configuration
    collisional_excitation_X.txt
                           electron-impact excitation fit coefficients out of
                           each configuration ('colex2' block)

Collisional ionization (the last two files; skip with --no-collisional).
Each line of collisional_parameters_X.txt is

    K L M ; Kf Lf Mf c0 c1 c2 c3 s0 ; ...

one entry per ionizing transition, taken unchanged from the model.  The
coefficients define the cross section

    sigma(U) = pi a0^2 (Ry/dE)^2 s0 (1/U) [c0 + c1 U^-1/2 + c2 U^-1 + c3 ln U]

with U = E/dE the incident energy in units of the transition energy (c0 + c1 +
c2 = 0 for every transition, so sigma vanishes at threshold; c3 = 0 in all
current models).  Averaged over a Maxwellian of temperature T [eV]:

    k(T) = 1.0899e-6 s0 exp(-b) q / (dE sqrt(T))   cm^3/s,   b = dE/T
    q    = c0 + c1 sqrt(pi b) erfc(sqrt b) e^b + c2 b E1(b) e^b + c3 E1(b) e^b

(1.0899e-6 = sqrt(8/pi m_e) pi a0^2 Ry^2.)  dE is not stored: it is the
difference of the two configurations' total energies in total_energies_X.txt,
which reproduces the model's own transition energy.  Recombination is not
tabulated; it follows from detailed balance with the statistical weights.
This form reproduces Voronov (1997) for hydrogen to ~10% from 2 eV to 1 keV and
C/N/O ions to within ~1.0-1.6x at 30-300 eV; neutral atoms and He-like ions are
2-6x high at T <~ 30 eV (the models lump the L shell at the outer-shell
ionization energy).  `--check` prints this comparison.

Collisional excitation (collisional_excitation_X.txt; also skipped with
--no-collisional).  Each line is

    K L M ; Ku Lu Mu c0 c1 c2 c3 a0 ; ...

one entry per excitation from the lower configuration, taken unchanged from the
model's 'colex2' block (an optional sixth coefficient c4 must be 0; it is in
all current models).  The collision strength and its Maxwellian average:

    W(x) = c0 ln x + c1 + c2/(x+a0) + c3/(x+a0)^2,     x = E/dE
    k(T) = 8.01e-8 q(b) exp(-b) / sqrt(T)   cm^3/s,    b = dE/T
    q(b) = c0 E1(b) e^b + c1 + c3 b/(a0+1) + (c2 - c3 b) b E1((a0+1)b) e^((a0+1)b)

(q follows from integrating W over the Maxwellian, and matches numerical
integration; `--check`.)  k is per ion in the lower level, with no
statistical-weight factor, as for ionization: it lies 1.1-2.8x above van
Regemorter for every transition, independent of g.  dE is the difference of the
two configurations' total energies.  De-excitation follows from detailed
balance with the statistical weights.  `--check` prints both comparisons.

Photoionization cross sections depend on the photon energy, so
rate_transitions_X.txt has to be regenerated when that changes.  The two energy
files do not depend on it.

Standalone: the only dependency is numpy (scipy for --check and --tables).

Fluorescence: the Einstein A of each bound-bound transition in the model
(phxs_rate_cgs), type 1.  Photoexcitation is not written: it is off resonance
at the photon energies used, and the module has no channel for it.

Format and use (2026-10-28): the files are specified in
ATOMIC_DATA.md (top of this repository) and checked by tools/check_atomic_data.py (--validate runs
it on the output).  Collisional ionization in mdrun defaults to Lotz's formula (mcmd-coll-data = lotz), which needs
only total_energies and statistical_weight; collisional_parameters_X.txt (or its tabulated form,
collisional_ionization_table_X.txt, --tables) is read with mcmd-coll-data = files.  collisional_excitation_X.txt is
read when mcmd-coll-excitation = 1 (off by default).  statistical_weight_X.txt is always written: detailed balance
needs it whenever collisions are on.

Usage:
    make_atomic_data_cretin.py MODELS OUTDIR ENERGY [--no-collisional] [--tables] [--validate]
    make_atomic_data_cretin.py --check MODELS          (rate check, needs scipy)

    MODELS  a .z01 model file, or a directory of them
    OUTDIR  directory to write the files into (created if needed)
    ENERGY  photon energy in eV
"""

import argparse
import os
import subprocess
import re
import sys
from collections import defaultdict

import numpy as np


TOOLS = os.path.dirname(os.path.abspath(__file__))   # check_atomic_data.py, make_collisional_tables.py


# ---------------------------------------------------------------- #
# .z01 model parsing
# ---------------------------------------------------------------- #

class Transition(object):
    def __init__(self, t):
        self.iso1 = int(t[1])
        self.i1 = int(t[2])
        self.iso2 = int(t[3])
        self.i2 = int(t[4])
        self.de = None
        self.g1 = None
        self.g2 = None


class Phis(Transition):
    """Photoionization."""
    def __init__(self, line):
        t = line.strip().split()
        Transition.__init__(self, t)
        self.a0 = float(t[5])
        self.a1 = float(t[6])
        self.a2 = float(t[7])
        self.a3 = float(t[8])
        self.n = float(t[9])
        self.de = float(t[10])
        self.emin = float(t[11])
        self.emax = float(t[12])


class Phxs(Transition):
    """Bound-bound oscillator strength; drives fluorescence and photoexcitation."""
    def __init__(self, line):
        t = line.strip().split()
        Transition.__init__(self, t)
        self.f = float(t[5])
        self.wavelength = float(t[6])
        self.width = float(t[7]) if len(t) > 7 else 0.0


class Augxs(Transition):
    """Auger."""
    def __init__(self, line):
        t = line.strip().split()
        Transition.__init__(self, t)
        self.A = float(t[5])


class Sampson(Transition):
    """Electron-impact ionization ('sampson ionize'), see the module docstring."""
    def __init__(self, line):
        t = line.strip().split()
        Transition.__init__(self, t)
        self.c0 = float(t[5])
        self.c1 = float(t[6])
        self.c2 = float(t[7])
        self.c3 = float(t[8])
        self.s0 = float(t[9])


class Colex(Transition):
    """Electron-impact excitation ('colex2'), see the module docstring."""
    def __init__(self, line):
        t = line.strip().split()
        Transition.__init__(self, t)
        self.c0 = float(t[5])
        self.c1 = float(t[6])
        self.c2 = float(t[7])
        self.c3 = float(t[8])
        self.a0 = float(t[9])
        self.c4 = float(t[10]) if len(t) > 10 else 0.0


class Level(object):
    def __init__(self, line, iso_name=None, iso_energy=None):
        t = line.strip().split()
        self.iso = int(t[1])          # electron count of this ionization stage
        self.i = int(t[2])
        self.name = t[3]
        self.g = float(t[4])
        self.energy = float(t[5])     # excitation above the stage ground level
        self.config = [int(x) for x in t[6:-1]]
        self.nmax = int(t[-1])
        self.iso_name = iso_name
        self.iso_energy = iso_energy  # ionization potential of this stage


class LevelDB(object):
    def __init__(self, path):
        self.path = path
        self.levels = []
        self._g = {}
        self._iso_energy = {}
        self._level_energy = {}
        self.element, sections = self._read(path)
        self._model_data(sections["model"])
        self.phxs = self._transitions(sections, "phxs", Phxs)
        self.ai = self._transitions(sections, "augxs", Augxs)
        self.phis = self._transitions(sections, "phis", Phis)
        self.colis = self._transitions(sections, "sampson", Sampson) or []
        self.colex = self._transitions(sections, "colex2", Colex) or []

    @staticmethod
    def _read(path):
        """The element from a 'c atom X' header line (None without one) and the data blocks by name: each block
        from its 'data <name> ...' line up to 'end data', comment lines (starting with c) removed."""
        with open(path) as fh:
            content = fh.read()
        match = re.search(r"^c[ \t]+atom[ \t]+(\S+)", content, re.MULTILINE)
        element = match.group(1) if match else None
        content = re.sub(r"^c.*\n?", "", content, flags=re.MULTILINE)
        sections = {}
        for m in re.finditer(r"^[ \t]*(data[ \t]+(\S+).*?)^[ \t]*end data", content, re.MULTILINE | re.DOTALL):
            if m.group(2) in sections:
                raise ValueError("%s: two '%s' blocks" % (path, m.group(2)))
            sections[m.group(2)] = m.group(1).strip()
        if "model" not in sections:
            raise ValueError("%s: no 'data model' block" % path)
        return element, sections

    def _model_data(self, section):
        for chunk in re.split("^.*enot|enot", section)[1:]:
            lines = [l.strip() for l in chunk.split("\n") if l.strip()]
            enot = lines[0].split()
            iso_name, iso_energy = str(enot[1]), float(enot[2])
            for line in lines[1:]:
                lv = Level(line, iso_name, iso_energy)
                self.levels.append(lv)
                self._g[(lv.iso, lv.i)] = lv.g
                self._iso_energy[(lv.iso, lv.i)] = lv.iso_energy
                self._level_energy[(lv.iso, lv.i)] = lv.energy

    def _transitions(self, sections, name, cls):
        try:
            section = sections[name]
        except KeyError:
            return None
        out = []
        for line in section.split("\n")[1:]:
            line = line.strip()
            if not line or line.startswith("enot") or line.startswith("c"):
                continue
            t = cls(line)
            t.g1 = self._g[(t.iso1, t.i1)]
            t.g2 = self._g[(t.iso2, t.i2)]
            t.de = self.transition_energy(t.iso1, t.i1, t.iso2, t.i2)
            out.append(t)
        return out

    def transition_energy(self, iso1, i1, iso2, i2):
        return (self._level_energy[(iso2, i2)]
                + self._iso_energy[(iso1, i1)]
                - self._level_energy[(iso1, i1)])

    def level(self, iso, i):
        return next((l for l in self.levels if l.iso == iso and l.i == i), None)

    def phxs_of(self, iso1, i1, iso2, i2):
        return next((t for t in self.phxs if t.iso1 == iso1 and t.i1 == i1
                     and t.iso2 == iso2 and t.i2 == i2), None)

    def augx_of(self, iso1, i1, iso2, i2):
        return next((t for t in self.ai if t.iso1 == iso1 and t.i1 == i1
                     and t.iso2 == iso2 and t.i2 == i2), None)


# ---------------------------------------------------------------- #
# Cross sections
# ---------------------------------------------------------------- #

C_CM_S = 29979245800        # speed of light [cm/s]


def phis_crossection(energy, t):
    """Photoionization, from the fitted log-polynomial in the model."""
    if not (t.emin <= energy <= t.emax):
        return 0.0
    b = energy / t.de
    pblog = (t.a0 + t.a1 * np.log(b) + t.a2 * np.log(b) ** 2
             + t.a3 * np.log(b) ** 3)
    return 1e-18 * t.n * np.exp(pblog) * (13.606 / t.de)


# Gaussian (cgs) constants for the fluorescence rate.
E_ESU = 4.80320471e-10      # electron charge [esu]
M_G = 9.1093837e-28         # electron mass [g]


def phxs_rate_cgs(db, iso1, i1, iso2, i2):
    """Fluorescence: Einstein A of the upper level in Gaussian units,

        A_ul = 8 pi^2 e^2 nu^2 / (m_e c^3) * (g_l / g_u) * f_lu    [1/s]

    (C K alpha ~6e11/s.)"""
    t = db.phxs_of(iso1, i1, iso2, i2)
    v = C_CM_S / (t.wavelength * 1e-8)              # wavelength [A] -> cm
    gi = db.level(iso1, i1).g
    gf = db.level(iso2, i2).g
    return (8 * np.pi ** 2 * E_ESU ** 2 * v ** 2 / (M_G * C_CM_S ** 3)) * (gi / gf) * t.f


def auger_crossection(db, iso1, i1, iso2, i2):
    return db.augx_of(iso1, i1, iso2, i2).A


# ---------------------------------------------------------------- #
# Writers
# ---------------------------------------------------------------- #

def _config3(level):
    cfg = list(level.config)
    while len(cfg) < 3:
        cfg.append(0)
    return cfg


def write_energy_levels(db, path):
    """Binding energy of each configuration's outermost electron.

    A level whose excitation exceeds the ionization potential of its own stage
    is autoionizing; there is no bound outermost electron to quote, so those
    are written as the 1e14 placeholder rather than omitted.
    """
    with open(path, "w") as fh:
        for lv in db.levels:
            bd = lv.iso_energy - lv.energy
            if bd < 0:
                bd = 1e14
            for x in _config3(lv):
                fh.write(f"{x} ")
            fh.write(f"{bd}\n")


def write_total_energies(db, path):
    """Total electronic binding energy of each configuration.

        T(config) = sum(enot(1..N)) - excitation,   N = electron count

    Unlike the outermost binding energy this is defined for autoionizing
    levels too, which is what lets charge transfer price an electron entering
    any shell of an acceptor rather than only its outermost one.
    """
    ip = {}
    for lv in db.levels:
        ip.setdefault(lv.iso, lv.iso_energy)

    tg = {0: 0.0}
    for stage in sorted(q for q in ip if q > 0):
        if stage - 1 not in tg:
            raise ValueError("%s: ionization stage %d has no stage %d below "
                             "it, so total energies cannot be accumulated"
                             % (db.path, stage, stage - 1))
        tg[stage] = tg[stage - 1] + ip[stage]

    with open(path, "w") as fh:
        for lv in db.levels:
            cfg = _config3(lv)
            n = sum(cfg)
            if n != lv.iso:
                continue      # the bare-nucleus line carries no occupations
            for x in cfg:
                fh.write(f"{x} ")
            fh.write(f"{tg[n] - lv.energy}\n")


def collect_transitions(db, energy):
    """Transitions out of each configuration, with the module's type codes
    (0 Auger, 1 fluorescence, 2 photoionization); no photoexcitation (see the
    module docstring)."""
    PHIS, FLUOR, AUGXS = 2, 1, 0
    out = []

    # Photoionization
    for t in db.phis:
        out.append((_config3(db.level(t.iso1, t.i1)),
                    _config3(db.level(t.iso2, t.i2)),
                    phis_crossection(energy, t), PHIS))

    # Fluorescence: the inverse of the tabulated bound-bound transition
    for t in db.phxs:
        out.append((_config3(db.level(t.iso2, t.i2)),
                    _config3(db.level(t.iso1, t.i1)),
                    phxs_rate_cgs(db, t.iso1, t.i1, t.iso2, t.i2), FLUOR))

    # Auger
    if db.ai:
        for t in db.ai:
            out.append((_config3(db.level(t.iso1, t.i1)),
                        _config3(db.level(t.iso2, t.i2)),
                        auger_crossection(db, t.iso1, t.i1, t.iso2, t.i2), AUGXS))

    grouped = defaultdict(lambda: ([], [], []))
    for initial, final, sigma, kind in out:
        grouped[tuple(initial)][0].append(final)
        grouped[tuple(initial)][1].append(sigma)
        grouped[tuple(initial)][2].append(kind)
    return [(list(k), f, s, t) for k, (f, s, t) in grouped.items()]


def write_transitions(db, path, energy):
    with open(path, "w") as fh:
        for initial, finals, rates, kinds in collect_transitions(db, energy):
            for x in initial:
                fh.write(f"{x} ")
            fh.write(";")
            for k in range(len(finals)):
                for x in finals[k]:
                    fh.write(f"{x} ")
                fh.write(f"{rates[k]} ")
                fh.write(f"{kinds[k]} ;")
            fh.write("\n")
        fh.write("0 0 0 ; 0 0 1  0.000 0; \n")


def write_statistical_weights(db, path):
    """Statistical weight g of each configuration (one level per configuration
    in these models; a configuration listed twice with different g is an
    error rather than a silent choice)."""
    g = {}
    for lv in db.levels:
        key = tuple(_config3(lv))
        if key in g and g[key] != lv.g:
            raise ValueError("%s: configuration %s has two weights (%g, %g)"
                             % (db.path, key, g[key], lv.g))
        g[key] = lv.g
    with open(path, "w") as fh:
        for key in sorted(g, key=lambda k: (-sum(k), k)):
            fh.write("%d %d %d %g\n" % (key[0], key[1], key[2], g[key]))
    return len(g)


def write_collisional_parameters(db, path):
    """Electron-impact ionization coefficients, grouped by initial configuration."""
    out = defaultdict(list)
    for t in db.colis:
        ci = tuple(_config3(db.level(t.iso1, t.i1)))
        cf = tuple(_config3(db.level(t.iso2, t.i2)))
        if sum(ci) - sum(cf) != 1:
            raise ValueError("%s: collisional transition %s -> %s is not a "
                             "single ionization" % (db.path, ci, cf))
        out[ci].append((cf, (t.c0, t.c1, t.c2, t.c3, t.s0)))
    with open(path, "w") as fh:
        for ci in sorted(out, key=lambda k: (-sum(k), k)):
            fh.write("%d %d %d ;" % ci)
            for cf, p in out[ci]:
                fh.write(" %d %d %d " % cf + " ".join("%.6g" % x for x in p) + " ;")
            fh.write("\n")
    return len(db.colis)


def write_collisional_excitation(db, path):
    """Electron-impact excitation coefficients, grouped by lower configuration."""
    out = defaultdict(list)
    for t in db.colex:
        if t.c4 != 0:
            raise ValueError("%s: colex2 transition with c4 != 0 (not supported by the module)" % db.path)
        ci = tuple(_config3(db.level(t.iso1, t.i1)))
        cf = tuple(_config3(db.level(t.iso2, t.i2)))
        if sum(ci) != sum(cf) or t.iso1 != t.iso2:
            raise ValueError("%s: excitation %s -> %s changes the electron count" % (db.path, ci, cf))
        out[ci].append((cf, (t.c0, t.c1, t.c2, t.c3, t.a0)))
    with open(path, "w") as fh:
        for ci in sorted(out, key=lambda k: (-sum(k), k)):
            fh.write("%d %d %d ;" % ci)
            for cf, p in out[ci]:
                fh.write(" %d %d %d " % cf + " ".join("%.6g" % x for x in p) + " ;")
            fh.write("\n")
    return len(db.colex)


# ---------------------------------------------------------------- #
# Rate check (--check): closed form vs numerical integration, and vs the
# Voronov (1997) fits for the ground level of every ion stage of H, C, N, O.
# ---------------------------------------------------------------- #

RY_EV, A0_CM = 13.6057, 0.529177e-8
SAMPSON_PREF = (np.sqrt(8 * 1.602177e-12 / (np.pi * 9.10938e-28))
                * np.pi * A0_CM ** 2 * RY_EV ** 2)          # 1.0899e-6

# Voronov 1997, ADNDT 65, 1: (Z, electrons) -> dE [eV], P, A [cm^3/s], X, K
# (copied from Verner's cfit.f, https://www.pa.uky.edu/~verner/dima/col/cfit.f)
VORONOV = {
    (1, 1): (13.6, 0, 2.91e-08, 0.2320, 0.39),
    (6, 6): (11.3, 0, 6.85e-08, 0.1930, 0.25), (6, 5): (24.4, 1, 1.86e-08, 0.2860, 0.24),
    (6, 4): (47.9, 1, 6.35e-09, 0.4270, 0.21), (6, 3): (64.5, 1, 1.50e-09, 0.4160, 0.13),
    (6, 2): (392.1, 1, 2.99e-10, 0.6660, 0.02), (6, 1): (490.0, 1, 1.23e-10, 0.6200, 0.16),
    (7, 7): (14.5, 0, 4.82e-08, 0.0652, 0.42), (7, 6): (29.6, 0, 2.98e-08, 0.3100, 0.30),
    (7, 5): (47.5, 1, 8.10e-09, 0.3500, 0.24), (7, 4): (77.5, 1, 3.71e-09, 0.5490, 0.18),
    (7, 3): (97.9, 0, 1.51e-09, 0.0167, 0.74), (7, 2): (552.1, 0, 3.71e-10, 0.5460, 0.29),
    (7, 1): (667.0, 1, 7.77e-11, 0.6240, 0.16),
    (8, 8): (13.6, 0, 3.59e-08, 0.0730, 0.34), (8, 7): (35.1, 1, 1.39e-08, 0.2120, 0.22),
    (8, 6): (54.9, 1, 9.31e-09, 0.2700, 0.27), (8, 5): (77.4, 0, 1.02e-08, 0.6140, 0.27),
    (8, 4): (113.9, 1, 2.19e-09, 0.6300, 0.17), (8, 3): (138.1, 0, 1.95e-09, 0.3600, 0.54),
    (8, 2): (739.3, 0, 2.12e-10, 0.3960, 0.35), (8, 1): (871.4, 1, 5.21e-11, 0.6290, 0.16),
}


def sampson_rate(t, de, T):
    """Maxwellian rate coefficient [cm^3/s] of one transition at T [eV]."""
    from scipy.special import erfcx, exp1
    b = de / T
    q = (t.c0 + t.c1 * np.sqrt(np.pi * b) * erfcx(np.sqrt(b))
         + (t.c2 * b + t.c3) * exp1(b) * np.exp(b))
    return SAMPSON_PREF * t.s0 * np.exp(-b) * q / (de * np.sqrt(T))


def sampson_rate_numeric(t, de, T):
    U = np.logspace(0, 5, 200001)
    E = U * de
    sig = (np.pi * A0_CM ** 2 * (RY_EV / de) ** 2 * t.s0 / U
           * (t.c0 + t.c1 / np.sqrt(U) + t.c2 / U + t.c3 * np.log(U)))
    v = np.sqrt(2 * E * 1.602177e-12 / 9.10938e-28)
    f = 2 * np.sqrt(E / np.pi) * T ** -1.5 * np.exp(-E / T)
    return np.trapz(sig * v * f, E)


def e1x(x):
    """E1(x) e^x without overflow (asymptotic series for x > 50)."""
    from scipy.special import exp1
    if x > 50:
        return (1 - 1 / x + 2 / x ** 2 - 6 / x ** 3 + 24 / x ** 4) / x
    return exp1(x) * np.exp(x)


def colex_rate(t, de, T):
    """Maxwellian excitation rate coefficient [cm^3/s], per ion in the lower level."""
    b, A = de / T, t.a0 + 1
    q = (t.c0 * e1x(b) + t.c1 + t.c3 * b / A
         + (t.c2 - t.c3 * b) * b * e1x(A * b))
    return 8.01e-8 * q * np.exp(-b) / np.sqrt(T)


def colex_rate_numeric(t, de, T):
    from scipy.integrate import quad
    b, a0 = de / T, t.a0
    W = lambda x: t.c0 * np.log(x) + t.c1 + t.c2 / (x + a0) + t.c3 / (x + a0) ** 2
    q = b * quad(lambda x: W(x) * np.exp(-b * (x - 1)), 1, np.inf, limit=200)[0]
    return 8.01e-8 * q * np.exp(-b) / np.sqrt(T)


def van_regemorter_rate(f, de, T):
    """van Regemorter (1962), gbar = max(0.2, sqrt(3)/(2 pi) E1(b) e^b), per ion."""
    from scipy.special import exp1
    b = de / T
    gbar = max(0.2, np.sqrt(3) / (2 * np.pi) * exp1(b) * np.exp(b))
    return 3.16e-7 * f * (RY_EV / de) * np.sqrt(RY_EV / T) * np.exp(-b) * gbar


def voronov_rate(Z, N, T):
    dE, P, A, X, K = VORONOV[(Z, N)]
    U = dE / T
    return A * (1 + P * np.sqrt(U)) / (X + U) * U ** K * np.exp(-U)


def check(models):
    temps = (10, 30, 100, 300)
    zof = {"H": 1, "C": 6, "N": 7, "O": 8}
    for path in models:
        db = LevelDB(path)
        sym = symbol_for(db, path)
        if not db.colis:
            print("%s: no 'sampson ionize' block" % sym)
            continue
        worst = max(abs(t.c0 + t.c1 + t.c2) for t in db.colis)
        t0 = db.colis[0]
        de0 = db.transition_energy(t0.iso1, t0.i1, t0.iso2, t0.i2)
        num = [sampson_rate_numeric(t0, de0, T) / sampson_rate(t0, de0, T) for T in temps]
        print("%s: %d transitions, max |c0+c1+c2| = %.1e, numeric/closed form %s"
              % (sym, len(db.colis), worst, " ".join("%.4f" % x for x in num)))
        if db.colex:
            ex = db.colex
            num = max(abs(colex_rate_numeric(t, t.de, T) / colex_rate(t, t.de, T) - 1)
                      for t in ex[::max(1, len(ex) // 20)] for T in temps)
            vr = [colex_rate(t, t.de, 50.0) / van_regemorter_rate(db.phxs_of(t.iso1, t.i1, t.iso2, t.i2).f, t.de, 50.0)
                  for t in ex if db.phxs_of(t.iso1, t.i1, t.iso2, t.i2) is not None]
            print("    %d excitations: max |numeric/closed form - 1| = %.1e; model/van Regemorter at 50 eV: "
                  "median %.2f, 10-90%% %.2f-%.2f" % (len(ex), num, np.median(vr),
                                                     np.percentile(vr, 10), np.percentile(vr, 90)))
        Z = zof.get(sym)
        if Z is None:
            continue
        print("    electrons  model/Voronov at T = " + " ".join("%5d" % T for T in temps) + " eV")
        for N in range(Z, 0, -1):
            tr = [t for t in db.colis if t.iso1 == N and t.i1 == 1]
            if not tr or (Z, N) not in VORONOV:
                continue
            r = [sum(sampson_rate(t, db.transition_energy(t.iso1, t.i1, t.iso2, t.i2), T)
                     for t in tr) / voronov_rate(Z, N, T) for T in temps]
            print("    %9d  %s" % (N, " " * 20 + " ".join("%5.2f" % x for x in r)))
    return 0


# ---------------------------------------------------------------- #

SYMBOLS = ("H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN GA GE AS SE BR KR "
           "RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE").split()


def symbol_for(db, path):
    """Element symbol as the module expects it: periodic table, upper case.  From the 'c atom X' header, else from
    the nuclear charge, the electron count of the model's neutral stage (not from the file name: C_pqn.z01,
    models.z02)."""
    if db.element:
        return db.element.upper()
    z = max(lv.iso for lv in db.levels)
    if not 1 <= z <= len(SYMBOLS):
        raise ValueError("%s: no 'c atom X' header, and Z = %d is outside the symbol table" % (path, z))
    return SYMBOLS[z - 1]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("models", help="a .z01 model file, or a directory of them")
    ap.add_argument("outdir", nargs="?", help="directory to write the data files into")
    ap.add_argument("energy", nargs="?", type=float, help="photon energy in eV")
    ap.add_argument("--no-collisional", action="store_true",
                    help="skip collisional_parameters_X.txt and collisional_excitation_X.txt (enough for "
                         "mcmd-coll-data = lotz with excitation off; statistical_weight_X.txt is still written)")
    ap.add_argument("--tables", action="store_true",
                    help="also write the tabulated collisional_ionization_table_X.txt and "
                         "collisional_excitation_table_X.txt (tools/make_collisional_tables.py, needs scipy)")
    ap.add_argument("--validate", action="store_true",
                    help="check the output against ATOMIC_DATA.md (tools/check_atomic_data.py)")
    ap.add_argument("--tools", default=TOOLS, help="directory of check_atomic_data.py and make_collisional_tables.py (default: this script's)")
    ap.add_argument("--check", action="store_true",
                    help="only check the collisional rate formula (needs scipy)")
    args = ap.parse_args(argv)
    for flag, tool in (("tables", "make_collisional_tables.py"), ("validate", "check_atomic_data.py")):
        if getattr(args, flag) and not os.path.exists(os.path.join(args.tools, tool)):
            ap.error(f"--{flag} needs {tool}, not in {args.tools} (give its directory with --tools)")
    if not args.check and (args.outdir is None or args.energy is None):
        ap.error("OUTDIR and ENERGY are required unless --check is given")

    if os.path.isdir(args.models):
        paths = sorted(os.path.join(args.models, n)
                       for n in os.listdir(args.models) if n.endswith(".z01"))
        if not paths:
            ap.error("no .z01 models in %s" % args.models)
    elif os.path.isfile(args.models):
        paths = [args.models]
    else:
        ap.error("no such file or directory: %s" % args.models)

    if args.check:
        return check(paths)

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    for path in paths:
        db = LevelDB(path)
        sym = symbol_for(db, path)
        write_energy_levels(db, os.path.join(args.outdir, "energy_levels_%s.txt" % sym))
        write_total_energies(db, os.path.join(args.outdir, "total_energies_%s.txt" % sym))
        write_transitions(db, os.path.join(args.outdir, "rate_transitions_%s.txt" % sym), args.energy)
        ng = write_statistical_weights(db, os.path.join(args.outdir, "statistical_weight_%s.txt" % sym))
        extra = ", %d weights" % ng
        if not args.no_collisional:
            nc = write_collisional_parameters(db, os.path.join(args.outdir, "collisional_parameters_%s.txt" % sym))
            ne = write_collisional_excitation(db, os.path.join(args.outdir, "collisional_excitation_%s.txt" % sym))
            extra += ", %d collisional ionizations, %d excitations" % (nc, ne)
        print("%-3s %4d levels%s -> %s" % (sym, len(db.levels), extra, args.outdir))

    print("\n%d model(s) at %g eV." % (len(paths), args.energy))
    print("The energy, weight and collisional files are photon-energy "
          "independent; rate_transitions_X.txt is not.")
    status = 0
    if args.tables:
        if args.no_collisional:
            print("--tables: nothing to tabulate with --no-collisional")
        else:
            status |= subprocess.call([sys.executable, os.path.join(args.tools, "make_collisional_tables.py"),
                                       args.outdir])
    if args.validate:
        status |= subprocess.call([sys.executable, os.path.join(args.tools, "check_atomic_data.py"), args.outdir,
                                   "--photon-energy", str(args.energy)])
    return status


if __name__ == "__main__":
    sys.exit(main())
