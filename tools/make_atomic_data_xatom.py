#!/usr/bin/env python3
"""Generate MolDStruct MC/MD Atomic_data files with XATOM (2026-10-28).

The XATOM counterpart of make_atomic_data_cretin.py: same files, same format
(ATOMIC_DATA.md at the top of this repository), from Hartree-Fock-Slater
calculations instead of CRETIN models.  XATOM is open source (GPLv3, https://www.desy.de/~xraypac/xatom.html);
the xatom executable must be on PATH or given with --xatom.  Writes

    energy_levels_X.txt, total_energies_X.txt, rate_transitions_X.txt,
    statistical_weight_X.txt

and no collisional files: mdrun's default mcmd-coll-data = lotz needs only the
total energies and weights.

Elements: those of the module's table (mcmd_elements[] in atomic_data.c) up to
Z = 18: H to Ar except He, with the same ground states.  Ca, Fe and Ni
are refused: the module counts their 4s electrons in M, so a K L M state is no
single configuration XATOM could compute (and I is a placeholder there).

States: every K L M with K, L, M at most the ground-state occupation (C 2 4 0,
S 2 8 6, ...).  Photoionization, Auger and fluorescence
never leave that box, and charge transfer only fills shells up to the ground
occupation, so it is everything a run can reach.  (The CRETIN models also hold
states with electrons promoted outwards, which nothing in the module produces.)

Subshells: XATOM resolves 1s 2s 2p 3s 3p; the module has only K L M.  Each
K L M state is computed as its aufbau member (2s filled before 2p, 3s before
3p: L = 5 is 2s2 2p3), and every XATOM final state is summed onto its K L M
(L = 2s + 2p, M = 3s + 3p).  So a state's energy and rates are those of its
lowest-energy subshell arrangement.

Energies: from the orbital eigenvalues, not XATOM's total energies.  The HFS
total energy is off by tens of eV (N5+ 1s2: 1267 eV against 1218), while the
eigenvalues are what XATOM itself prices its transitions with and are close to
the measured ionization energies (S 3p 10.29 eV, measured 10.36).  So

    T(c) = T(c minus one electron from its outermost shell) - eps_top(c)

with eps_top the eigenvalue of the highest occupied subshell of the outermost
shell (2p if occupied, else 2s, ...), from T(0 0 0) = 0; energy_levels is then
-eps_top exactly.  A channel these energies close (photoelectron E_ph - dT <= 0,
Auger or fluorescence dT <= 0) is dropped and counted, as is any XATOM channel
that leaves the K L M state unchanged (2p -> 2s within L).

Conversions: P and V (core and valence
photoionization) cross sections a0^2 -> cm^2, type 2; A rates 1/t_au -> 1/s,
type 0; F rates the same, type 1.  E (photoexcitation) is not written: the
module has no channel for it.  One-electron states use XATOM's exact -Z/r
potential (-exact).

Usage:
    make_atomic_data_xatom.py OUTDIR ENERGY [--elements H,C,N,O,S] [--jobs 12] [--keep-raw DIR] [--validate]
                              [--xatom PATH]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from math import comb
from multiprocessing import Pool

HARTREE_EV = 27.211386245988
A0_CM2 = 0.529177210903e-8 ** 2           # a0^2 [cm^2]
AU_RATE_S = 1 / 2.4188843265857e-17       # 1/t_au [1/s]
TYPE = {"P": 2, "V": 2, "A": 0, "F": 1}

# symbol: (Z, ground K L M), as mcmd_elements[] in src/kernel/atomic_data.c
ELEMENTS = {"H": (1, (1, 0, 0)), "LI": (3, (2, 1, 0)), "BE": (4, (2, 2, 0)),
            "B": (5, (2, 3, 0)), "C": (6, (2, 4, 0)), "N": (7, (2, 5, 0)), "O": (8, (2, 6, 0)),
            "F": (9, (2, 7, 0)), "NE": (10, (2, 8, 0)), "NA": (11, (2, 8, 1)), "MG": (12, (2, 8, 2)),
            "AL": (13, (2, 8, 3)), "SI": (14, (2, 8, 4)), "P": (15, (2, 8, 5)), "S": (16, (2, 8, 6)),
            "CL": (17, (2, 8, 7)), "AR": (18, (2, 8, 8))}
UNSUPPORTED = {"HE": "not in the module table (a mass ~4 is refused)", "CA": "4s counted in M", "FE": "4s counted in M", "NI": "4s counted in M",
               "I": "placeholder in the module"}
TOOLS = os.path.dirname(os.path.abspath(__file__))     # check_atomic_data.py


def subshells(klm, ground):
    """XATOM -config string of the aufbau member of K L M."""
    k, l, m = klm
    s = [f"1s{k}"]
    if ground[1]:
        s += [f"2s{min(l, 2)}", f"2p{l - min(l, 2)}"]
    if ground[2]:
        s += [f"3s{min(m, 2)}", f"3p{m - min(m, 2)}"]
    return "_".join(s)


def to_klm(occ):
    """Subshell occupations (1s [2s 2p [3s 3p]]) -> K L M."""
    occ = list(occ)
    occ += [0] * (5 - len(occ))
    return occ[0], occ[1] + occ[2], occ[3] + occ[4]


def run_xatom(job):
    el, klm, energy, raw, exe = job
    Z, ground = ELEMENTS[el]
    cfg = subshells(klm, ground)
    with tempfile.TemporaryDirectory() as tmp:
        dat, out = os.path.join(tmp, "x.dat"), os.path.join(tmp, "x.out")
        cmd = [exe, "-Z", str(Z), "-config", cfg, "-pcs", "-PE", str(energy), "-decay",
               "-dat", dat, "-o", out] + (["-exact"] if sum(klm) == 1 else [])
        r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
        txt = open(out).read() if os.path.exists(out) else ""
        orb = dict(re.findall(r"^#\s+(\d[spd])\s+\d+\s+(-?[\d.]+)\s+-?[\d.]+\s+-?[\d.]+\s*$", txt, re.M))
        occ = [x[:2] for x in cfg.split("_") if int(x[2:]) > 0]          # occupied subshells, in order
        top = [x for x in occ if x[0] == max(y[0] for y in occ)][-1]     # last one of the outermost shell
        if r.returncode or top not in orb:
            raise RuntimeError(f"xatom failed for {el} {cfg}: rc {r.returncode}\n{r.stdout[-500:]}{r.stderr[-500:]}")
        lines = open(dat).read().splitlines() if os.path.exists(dat) else []
        if raw:
            name = f"{el}_{cfg}"
            open(os.path.join(raw, name + ".out"), "w").write(txt)
            open(os.path.join(raw, name + ".dat"), "w").write("\n".join(lines) + "\n")
    trans = []
    for ln in lines:
        head, _, tail = ln.partition("#")
        f = head.split()
        if not f or f[0] not in TYPE:
            continue
        ini, fin = tail.split("->")
        assert to_klm(map(int, ini.split())) == klm, (el, cfg, ln)
        scale = A0_CM2 if TYPE[f[0]] == 2 else AU_RATE_S
        trans.append((to_klm(map(int, fin.split())), float(f[3]) * scale, TYPE[f[0]]))
    return el, klm, -float(orb[top]) * HARTREE_EV, trans


def states(ground):
    return [(k, l, m) for k in range(ground[0] + 1) for l in range(ground[1] + 1)
            for m in range(ground[2] + 1) if k + l + m > 0]


def outermost(klm):
    return max(i for i in range(3) if klm[i] > 0)


def write_element(outdir, el, res, energy):
    T = {(0, 0, 0): 0.0}
    for s in sorted(res, key=sum):                      # fewer electrons first
        o = outermost(s)
        T[s] = T[tuple(x - (i == o) for i, x in enumerate(s))] + res[s][0]
    ground_T = defaultdict(float)                       # most bound state per electron count
    for s, t in T.items():
        ground_T[sum(s)] = max(ground_T[sum(s)], t)
    order = sorted(T, key=lambda s: (-sum(s), [-x for x in s]))
    with open(os.path.join(outdir, f"total_energies_{el}.txt"), "w") as fh:
        for s in order:
            fh.write(f"{s[0]} {s[1]} {s[2]} {T[s]!r}\n")
    with open(os.path.join(outdir, f"energy_levels_{el}.txt"), "w") as fh:
        for s in order:
            if s == (0, 0, 0):
                e = 0.0
            elif T[s] < ground_T[sum(s) - 1]:           # autoionizing
                e = 1e14
            else:
                o = outermost(s)
                e = T[s] - T[tuple(x - (i == o) for i, x in enumerate(s))]
            fh.write(f"{s[0]} {s[1]} {s[2]} {e!r}\n")
    with open(os.path.join(outdir, f"statistical_weight_{el}.txt"), "w") as fh:
        for s in order:
            if s != (0, 0, 0):
                fh.write(f"{s[0]} {s[1]} {s[2]} {comb(2, s[0]) * comb(8, s[1]) * comb(18, s[2])}\n")
    counts, dropped, longest = defaultdict(int), defaultdict(list), 0
    with open(os.path.join(outdir, f"rate_transitions_{el}.txt"), "w") as fh:
        for s in order:
            if s == (0, 0, 0):
                continue
            merged = defaultdict(float)                 # subshell channels onto K L M
            for fin, rate, typ in res[s][1]:
                dT = T[s] - T[fin]
                if fin == s:
                    dropped["same state"].append(rate)
                elif (typ == 2 and energy - dT <= 0) or (typ != 2 and dT >= 0):
                    dropped[("photo", "Auger", "fluorescence")[[2, 0, 1].index(typ)] + " closed"].append(rate)
                else:
                    merged[(fin, typ)] += rate
            line = f"{s[0]} {s[1]} {s[2]} ;" + "".join(
                f"{f[0]} {f[1]} {f[2]} {r:.10e} {t} ;" for (f, t), r in sorted(merged.items()))
            for (_, t) in merged:
                counts[t] += 1
            longest = max(longest, len(line))
            fh.write(line + "\n")
    if longest >= 1024:
        raise RuntimeError(f"{el}: rate_transitions line of {longest} characters (reader limit 1023)")
    return len(order) - 1, counts, dropped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("energy", type=float, help="photon energy in eV")
    ap.add_argument("--elements", default="H,C,N,O,S")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--keep-raw", metavar="DIR", help="also keep every XATOM .out and .dat here")
    ap.add_argument("--validate", action="store_true", help="run tools/check_atomic_data.py on the output")
    ap.add_argument("--tools", default=TOOLS, help="directory of check_atomic_data.py (default: this script's)")
    ap.add_argument("--xatom", default="xatom", help="the XATOM executable (default: xatom on PATH)")
    a = ap.parse_args(argv)
    if shutil.which(a.xatom) is None:
        ap.error(f"XATOM executable '{a.xatom}' not found (put xatom on PATH or give it with --xatom)")
    if a.validate and not os.path.exists(os.path.join(a.tools, "check_atomic_data.py")):
        ap.error(f"--validate needs check_atomic_data.py, not in {a.tools} (give its directory with --tools)")
    els = a.elements.upper().split(",")
    bad = [f"{e} ({UNSUPPORTED.get(e, 'not in the module table')})" for e in els if e not in ELEMENTS]
    if bad:
        ap.error("cannot generate " + ", ".join(bad) + "; supported: " + " ".join(ELEMENTS))
    os.makedirs(a.outdir, exist_ok=True)
    if a.keep_raw:
        os.makedirs(a.keep_raw, exist_ok=True)
    jobs = [(el, s, a.energy, a.keep_raw, a.xatom) for el in els for s in states(ELEMENTS[el][1])]
    res = defaultdict(dict)
    with Pool(a.jobs) as pool:
        for el, s, t, tr in pool.imap_unordered(run_xatom, jobs):
            res[el][s] = (t, tr)
    for el in els:
        n, c, d = write_element(a.outdir, el, res[el], a.energy)
        print(f"{el}: {n} states, photo {c[2]} / Auger {c[0]} / fluorescence {c[1]} channels -> {a.outdir}"
              + "".join(f"; dropped {len(v)} {k} (largest {max(v):.2e})" for k, v in sorted(d.items())))
    print(f"XATOM at {a.energy:g} eV; only rate_transitions_X.txt depends on the photon energy.")
    if a.validate:
        return subprocess.call([sys.executable, os.path.join(a.tools, "check_atomic_data.py"), a.outdir,
                                "--photon-energy", str(a.energy)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
