#!/usr/bin/env python3
"""Check an Atomic_data directory against the format in ATOMIC_DATA.md, before a run.

    check_atomic_data.py <Atomic_data dir> [--elements H,C,N,O,S] [--photon-energy 600]
                         [--no-ct] [--no-collisions]

Errors are what mdrun would stop on (or silently misread); warnings are what is legal but suspicious.  Exit status
1 if there are errors.  The checks follow the reader (src/kernel/atomic_data.c) and its use in mcionize.c:
  - files present for the features used (charge transfer: total_energies; mcmd-coll-data = files:
    collisional_parameters); statistical_weight is optional (mdrun computes C(2,K) C(8,L) C(18,M) without it);
  - every line parses, lines are shorter than the reader's 1024-character buffer;
  - transitions: five fields per final state, type 0/1/2, electron count consistent with the type, rates finite and
    in a plausible range for their units (photo: cross-section, cm^2; Auger and fluorescence: s^-1);
  - every state the data can reach has its own rate line (only the bare nucleus 0 0 0 may be missing), the
    element's ground state is there, occupations within K <= 2, L <= 8, M <= 18;
  - total energies: present for every reachable state, consistent with energy_levels (outermost binding energy
    = T(c) - T(c minus its outermost electron)), giving positive photo- and Auger-electron energies;
  - collisional ionization: one electron fewer in the final state, all states known, dE > 0; excitation: same
    electron count, dE > 0; statistical weights, where given, compared with the shell combinatorics.
"""
import argparse
import math
import os
import sys
from math import comb

LINEBUF = 1023                     # fgets(line, 1024): longer lines are split and misread
CAP = (2, 8, 18)                   # K, L, M capacities
GROUND = {"H": (1, 0, 0), "C": (2, 4, 0), "N": (2, 5, 0), "O": (2, 6, 0), "F": (2, 7, 0), "MG": (2, 8, 2),
          "P": (2, 8, 5), "S": (2, 8, 6), "CL": (2, 8, 7), "CA": (2, 8, 10), "FE": (2, 8, 16), "NI": (2, 8, 18),
          "SI": (2, 8, 4), "NA": (2, 8, 1), "I": (2, 8, 18), "LI": (2, 1, 0), "BE": (2, 2, 0),
          "B": (2, 3, 0), "NE": (2, 8, 0), "AL": (2, 8, 3), "AR": (2, 8, 8)}
TYPES = {0: "Auger", 1: "fluorescence", 2: "photoionization"}
SENTINEL = 1e14                    # energy_levels: autoionizing level (no positive binding energy)
NOT_FOUND = -1.0                   # the reader's 'missing' value: an energy may not be exactly -1


class Report:
    def __init__(self):
        self.err, self.warn = [], []

    def e(self, msg):
        self.err.append(msg)

    def w(self, msg):
        self.warn.append(msg)


def nel(s):
    return sum(s)


def state(tok):
    return tuple(int(x) for x in tok[:3])


def read_table(path, rep, name):
    """'K L M value' per line -> {state: value}."""
    out = {}
    for n, line in enumerate(open(path), 1):
        if not line.strip():
            continue
        f = line.split()
        if len(f) != 4:
            rep.e(f"{name}:{n}: expected 'K L M value', got {line.strip()!r}")
            continue
        s, v = state(f), float(f[3])
        if s in out:
            rep.w(f"{name}:{n}: state {s} listed twice (the first wins)")
            continue
        if v == NOT_FOUND:
            rep.e(f"{name}:{n}: value -1 is the reader's 'not found' marker")
        if not math.isfinite(v):
            rep.e(f"{name}:{n}: value not finite")
        out[s] = v
    return out


def read_lines(path, rep, name, nfield):
    """'K L M ; K' L' M' f1 .. ; ...' per line -> list of (line number, initial, [(final, fields)])."""
    rows = []
    for n, raw in enumerate(open(path), 1):
        if len(raw.rstrip("\n")) > LINEBUF:
            rep.e(f"{name}:{n}: line has {len(raw.rstrip())} characters; the reader's buffer holds {LINEBUF}")
        line = raw.strip()
        if not line:
            continue
        toks = [t.split() for t in line.split(";")]
        if not toks[0] or len(toks[0]) != 3:
            rep.e(f"{name}:{n}: no initial state 'K L M ;'")
            continue
        ini, fins = state(toks[0]), []
        for t in toks[1:]:
            if not t:
                continue                       # the trailing ';'
            if len(t) != 3 + nfield:
                rep.e(f"{name}:{n}: final state token {' '.join(t)!r} has {len(t)} fields, expected {3 + nfield}"
                      + (" (rate files need 'K L M rate type': the type column)" if nfield == 2 else ""))
                continue
            fins.append((state(t), [float(x) for x in t[3:]]))
        rows.append((n, ini, fins))
    return rows


def check_occ(s, rep, where):
    if any(x < 0 or x > c for x, c in zip(s, CAP)):
        rep.e(f"{where}: occupation {s} outside 0..{CAP}")


def ln_text(tok):
    """ln of a non-negative number in text, any exponent ('3.1e-2650'), as mdrun reads it; -inf for 0."""
    t = tok.lower()
    m, _, x = t.partition("e")
    v = float(m)
    if v < 0 or not math.isfinite(v):
        raise ValueError(tok)
    return -math.inf if v == 0 else math.log(v) + (int(x) if x else 0) * math.log(10)


def check_table(path, name, excite, have, T, G, rep):
    """collisional_ionization_table / collisional_excitation_table (ATOMIC_DATA.md)."""
    lines = [l for l in open(path) if l.strip()]
    head = lines[0].split()
    if not head or head[0] != "T_e":
        rep.e(f"{name}: the first line must be 'T_e t1 t2 ... tn' (eV)")
        return
    grid = [float(x) for x in head[1:]]
    if len(grid) < 2 or any(t <= 0 for t in grid) or any(b <= a_ for a_, b in zip(grid, grid[1:])):
        rep.e(f"{name}: the T_e grid must have >= 2 positive, increasing temperatures")
        return
    lnT = [math.log(t) for t in grid]
    if grid[0] > 1.0 or grid[-1] < 1e3:
        rep.w(f"{name}: grid {grid[0]:g}..{grid[-1]:g} eV does not cover 1 eV .. 1 keV (outside, mdrun holds the edge)")
    rough = 0
    for n, line in enumerate(lines[1:], 2):
        f = line.split()
        where = f"{name}:{n}"
        if len(f) != 6 + len(grid):
            rep.e(f"{where}: expected 6 occupations and {len(grid)} coefficients, got {len(f)} fields")
            continue
        ini, fs = tuple(int(x) for x in f[:3]), tuple(int(x) for x in f[3:6])
        try:
            lk = [ln_text(x) for x in f[6:]]
        except ValueError:
            rep.e(f"{where}: coefficients must be numbers >= 0 (cm^3/s)")
            continue
        dn = nel(fs) - nel(ini)
        if dn != (0 if excite else -1):
            rep.e(f"{where}: {ini} -> {fs}: {'excitation must keep the electron count' if excite else 'ionization must remove one electron'} (change {dn:+d})")
        for s in (ini, fs):
            if s not in have:
                rep.e(f"{where}: state {s} has no rate_transitions line (fatal at start)")
            if s not in T and s != (0, 0, 0):
                rep.e(f"{where}: state {s} has no total energy (fatal at start)")
            if s not in G:
                rep.e(f"{where}: state {s} has no statistical weight (fatal at start)")
        if ini not in T or (fs not in T and fs != (0, 0, 0)):
            continue
        dE = T[ini] - (0.0 if fs == (0, 0, 0) else T[fs])
        if dE <= 0:
            rep.e(f"{where}: {ini} -> {fs}: dE = T_i - T_f = {dE:.3f} eV must be > 0 (fatal at start)")
            continue
        # mdrun interpolates g = ln k + dE/T linearly in ln T: g must be smooth.  A table whose threshold is not
        # the total-energy dE leaves a residual (dE_table - dE)/T that blows up at low T.
        g = [(x + dE / t, lt) for x, t, lt in zip(lk, grid, lnT) if math.isfinite(x)]
        if len(g) >= 3:
            slope = max(abs((g[i + 1][0] - g[i][0]) / (g[i + 1][1] - g[i][1])) for i in range(len(g) - 1))
            if slope > 30:
                rough += 1
                if rough <= 3:
                    rep.w(f"{where}: {ini} -> {fs}: ln<sigma v> + dE/T changes by up to {slope:.0f} per unit ln T: "
                          f"the table's threshold disagrees with dE = {dE:.2f} eV from the total energies?")
        if not any(math.isfinite(x) for x in lk):
            rep.w(f"{where}: {ini} -> {fs}: all coefficients 0 (never fires)")
    if rough > 3:
        rep.w(f"{name}: {rough} transitions in all with a rough ln<sigma v> + dE/T")


def check_element(d, el, a, rep):
    p = lambda f: os.path.join(d, f"{f}_{el}.txt")  # noqa: E731
    need = ["energy_levels", "rate_transitions"] + (["total_energies"] if a.ct or a.coll else []) + \
           (["collisional_parameters"] if a.coll and a.coll_data == "files"
            and not os.path.exists(p("collisional_ionization_table")) else [])
    for f in need:
        if not os.path.exists(p(f)):
            rep.e(f"{f}_{el}.txt missing" + (" (needed for charge transfer / collisions)" if f == "total_energies" else ""))
    if not os.path.exists(p("energy_levels")) or not os.path.exists(p("rate_transitions")):
        return
    E = read_table(p("energy_levels"), rep, f"energy_levels_{el}")
    T = read_table(p("total_energies"), rep, f"total_energies_{el}") if os.path.exists(p("total_energies")) else None
    G = read_table(p("statistical_weight"), rep, f"statistical_weight_{el}") \
        if os.path.exists(p("statistical_weight")) else None
    G_file = G is not None

    rows = read_lines(p("rate_transitions"), rep, f"rate_transitions_{el}", 2)
    have = {}
    for n, ini, fins in rows:
        if ini in have:
            rep.w(f"rate_transitions_{el}:{n}: initial state {ini} also on line {have[ini]} (the first line wins)")
        else:
            have[ini] = n
    if G is None:                       # what mdrun computes without the file
        G = {s: comb(2, s[0]) * comb(8, s[1]) * comb(18, s[2]) for s in have}
    reach = set(have)
    count = {0: 0, 1: 0, 2: 0}
    zero = []
    for n, ini, fins in rows:
        where = f"rate_transitions_{el}:{n}"
        check_occ(ini, rep, where)
        seen = set()
        for fs, (rate, typ) in fins:
            check_occ(fs, rep, where)
            t = int(typ)
            if rate == 0:
                # an inactive entry (the generator closes every file with '0 0 0 ; 0 0 1 0.000 0'): never fires
                zero.append(f"{ini} -> {fs}")
                continue
            reach.add(fs)
            if typ != t or t not in TYPES:
                rep.e(f"{where}: {ini} -> {fs}: type {typ} is not 0 (Auger), 1 (fluorescence) or 2 (photo)")
                continue
            count[t] += 1
            if (fs, t) in seen:
                rep.w(f"{where}: {ini} -> {fs} type {t} listed twice (both are used: the rates add)")
            seen.add((fs, t))
            dn = nel(fs) - nel(ini)
            if t in (0, 2) and dn != -1:
                rep.e(f"{where}: {ini} -> {fs}: {TYPES[t]} must remove one electron (change {dn:+d})")
            if t == 1 and dn != 0:
                rep.e(f"{where}: {ini} -> {fs}: fluorescence must keep the electron count (change {dn:+d})")
            if not math.isfinite(rate) or rate < 0:
                rep.e(f"{where}: {ini} -> {fs}: rate {rate} not finite and >= 0")
            elif t == 2 and rate > 1e-10:
                rep.e(f"{where}: {ini} -> {fs}: photoionization {rate:.3g} cannot be a cross-section in cm^2 "
                      f"(looks like a rate in 1/s: wrong type column?)")
            elif t in (0, 1) and 0 < rate < 1e-10:
                rep.e(f"{where}: {ini} -> {fs}: {TYPES[t]} {rate:.3g} cannot be a rate in 1/s "
                      f"(looks like a photoionization cross-section in cm^2: wrong type column?)")
            elif rate > 0:
                lo, hi, unit = {2: (1e-26, 1e-15, "cm^2"), 0: (1e9, 1e18, "1/s"), 1: (1e5, 1e17, "1/s")}[t]
                if not lo <= rate <= hi:
                    rep.w(f"{where}: {ini} -> {fs}: {TYPES[t]} {rate:.3g} outside {lo:g}..{hi:g} {unit}"
                          + (" (a fluorescence rate this small is the old-data symptom)" if t == 1 and rate < lo else ""))
    if zero:
        rep.w(f"rate_transitions_{el}: {len(zero)} zero-rate (inactive) entries, e.g. {zero[0]}")
    missing = sorted(s for s in reach - set(have) if s != (0, 0, 0))
    for s in missing:
        rep.e(f"rate_transitions_{el}: state {s} is reachable but has no line of its own (fatal when an atom gets "
              f"there; give it a line, even one without transitions: '{s[0]} {s[1]} {s[2]} ;')")
    g0 = GROUND.get(el)
    if g0 and g0 not in have:
        rep.e(f"rate_transitions_{el}: the ground state {g0} has no line")
    if g0 is None:
        rep.w(f"{el} is not in the element table of atomic_data.c (mcmd_elements[])")

    for s in sorted(have):
        if s not in E and s != (0, 0, 0):
            rep.w(f"energy_levels_{el}: no energy for state {s} (used by charge transfer and the escape fallback)")
    nsent = sum(1 for v in E.values() if v >= SENTINEL)
    if any(v < 0 for v in E.values()):
        rep.e(f"energy_levels_{el}: negative binding energies (autoionizing levels take {SENTINEL:g})")

    if T is not None:
        for s in sorted(reach):
            if s not in T and s != (0, 0, 0):
                rep.e(f"total_energies_{el}: no total energy for reachable state {s}")
        if (0, 0, 0) in T and T[(0, 0, 0)] != 0:
            rep.w(f"total_energies_{el}: T(0 0 0) = {T[(0, 0, 0)]}, expected 0 (bare nucleus)")
        bad = 0
        for s, e in E.items():
            if e >= SENTINEL or nel(s) == 0 or s not in T:
                continue
            sh = max(i for i in range(3) if s[i] > 0)            # outermost occupied shell
            s1 = tuple(x - (i == sh) for i, x in enumerate(s))
            t1 = 0.0 if s1 == (0, 0, 0) else T.get(s1)
            if t1 is not None and abs((T[s] - t1) - e) > 1e-3 * max(1.0, abs(e)):
                bad += 1
                if bad <= 3:
                    rep.w(f"energy_levels_{el}: E{s} = {e} but T{s} - T{s1} = {T[s] - t1:.4f}")
        if bad > 3:
            rep.w(f"energy_levels_{el}: {bad} states in all where E != T(c) - T(c minus outermost electron)")
        neg = 0
        for n, ini, fins in rows:
            for fs, (rate, typ) in fins:
                if ini not in T or (fs not in T and fs != (0, 0, 0)) or rate == 0:
                    continue
                tf = 0.0 if fs == (0, 0, 0) else T[fs]
                ek = (a.photon_energy - (T[ini] - tf)) if int(typ) == 2 else (tf - T[ini]) if int(typ) == 0 else None
                if ek is not None and ek <= 0 and (int(typ) == 0 or a.photon_energy):
                    neg += 1
                    if neg <= 3:
                        rep.e(f"rate_transitions_{el}:{n}: {TYPES[int(typ)]} {ini} -> {fs} would emit an electron of "
                              f"{ek:.2f} eV (from the total energies)")
        if neg > 3:
            rep.e(f"rate_transitions_{el}: {neg} transitions in all with a non-positive electron energy")

    if a.coll and T is not None and G is not None:
        for fname, excite in (("collisional_ionization_table", False), ("collisional_excitation_table", True)):
            if not excite and a.coll_data == "lotz":
                continue                                   # ionization from Lotz's formula: the file is not read
            if os.path.exists(p(fname)):
                check_table(p(fname), f"{fname}_{el}", excite, have, T, G, rep)
        for fname, excite in (("collisional_parameters", False), ("collisional_excitation", True)):
            table = "collisional_ionization_table" if not excite else "collisional_excitation_table"
            if not os.path.exists(p(fname)) or os.path.exists(p(table)) or (not excite and a.coll_data == "lotz"):
                continue                                   # a table is used instead, or Lotz's formula
            for n, ini, fins in read_lines(p(fname), rep, f"{fname}_{el}", 5):
                where = f"{fname}_{el}:{n}"
                for fs, par in fins:
                    dn = nel(fs) - nel(ini)
                    if not excite and dn != -1:
                        rep.e(f"{where}: {ini} -> {fs}: collisional ionization must remove one electron (change "
                              f"{dn:+d}); recombination is derived by detailed balance, excitation goes to "
                              f"collisional_excitation_{el}.txt")
                    if excite and dn != 0:
                        rep.e(f"{where}: {ini} -> {fs}: excitation must keep the electron count (change {dn:+d})")
                    for s in (ini, fs):
                        if s not in have:
                            rep.e(f"{where}: state {s} has no rate_transitions line (fatal at start)")
                        if s not in T and s != (0, 0, 0):
                            rep.e(f"{where}: state {s} has no total energy (fatal at start)")
                        if s not in G:
                            rep.e(f"{where}: state {s} has no statistical weight (fatal at start)")
                    if ini in T and (fs in T or fs == (0, 0, 0)):
                        dE = T[ini] - (0.0 if fs == (0, 0, 0) else T[fs])
                        if dE <= 0:
                            rep.e(f"{where}: {ini} -> {fs}: dE = T_i - T_f = {dE:.3f} eV must be > 0 (fatal at start)")
                    if not all(math.isfinite(x) for x in par):
                        rep.e(f"{where}: {ini} -> {fs}: parameters not finite")
                    if not excite and par[4] <= 0:
                        rep.w(f"{where}: {ini} -> {fs}: scale s0 = {par[4]} (the rate is proportional to it)")
        cb = lambda s: comb(2, s[0]) * comb(8, s[1]) * comb(18, s[2])  # noqa: E731
        off = [s for s, g in G.items() if G_file and abs(g - cb(s)) > 1e-5 * cb(s)]     # 6-digit printing is fine
        if off:
            rep.w(f"statistical_weight_{el}: {len(off)} weights differ from C(2,K) C(8,L) C(18,M), e.g. {off[:3]} "
                  f"(mdrun uses them as given)")
    return dict(states=len(have), photo=count[2], auger=count[0], fluor=count[1], sentinel=nsent)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", nargs="?", default="Atomic_data")
    ap.add_argument("--elements", help="comma list (default: every element with an energy_levels file)")
    ap.add_argument("--photon-energy", type=float, default=0.0, help="eV: also check photoelectron energies > 0")
    ap.add_argument("--no-ct", dest="ct", action="store_false", help="charge transfer off (total energies optional)")
    ap.add_argument("--coll-data", choices=("lotz", "files"), default="lotz",
                    help="as mcmd-coll-data: lotz (default; no ionization file needed) or files")
    ap.add_argument("--no-collisions", dest="coll", action="store_false",
                    help="mcmd-collisions = 0 (collisional files not needed)")
    a = ap.parse_args()
    els = a.elements.upper().split(",") if a.elements else sorted(
        f[len("energy_levels_"):-4] for f in os.listdir(a.dir) if f.startswith("energy_levels_") and f.endswith(".txt"))
    nerr = 0
    for el in els:
        rep = Report()
        info = check_element(a.dir, el, a, rep)
        head = f"{el}: " + (f"{info['states']} states, transitions photo {info['photo']} / Auger {info['auger']} / "
                            f"fluorescence {info['fluor']}, {info['sentinel']} autoionizing levels" if info else "")
        print(f"{head} -> {'OK' if not rep.err else f'{len(rep.err)} ERRORS'}"
              + (f", {len(rep.warn)} warnings" if rep.warn else ""))
        for m in rep.err[:20]:
            print("  ERROR  ", m)
        if len(rep.err) > 20:
            print(f"  ... {len(rep.err) - 20} more errors")
        for m in rep.warn[:10]:
            print("  warning", m)
        if len(rep.warn) > 10:
            print(f"  ... {len(rep.warn) - 10} more warnings")
        nerr += len(rep.err)
    sys.exit(1 if nerr else 0)


if __name__ == "__main__":
    main()
