#!/usr/bin/env python3
"""Tabulate the collisional rate coefficients of an Atomic_data directory (2026-10-28).

Reads the CRETIN fit files (collisional_parameters_X.txt: Sampson ionization, collisional_excitation_X.txt: colex2
excitation) and total_energies_X.txt, and writes the generator-independent tables mdrun reads in their place
(collisional_ionization_table_X.txt, collisional_excitation_table_X.txt; format in ATOMIC_DATA.md): the full
Maxwellian rate coefficient <sigma v> (cm^3/s) of every transition on a logarithmic T_e grid, with the same formulas
as mcionize.c, written with any exponent needed (the low-T values are far below the smallest double).

    make_collisional_tables.py <Atomic_data> [--out DIR] [--tmin 0.1] [--tmax 1e5] [--per-decade 40] [--check]

--lotz: the ionization table from Lotz's formula instead of the fits (no CRETIN data used): every state of
rate_transitions_X.txt loses one electron from each occupied shell, wherever that state exists, with
<sigma v> = 6.7e-7 a q / T^1.5 E1(P/T) / (P/T) cm^3/s (Lotz 1967/68, b = 0; a = 4.5 in 1e-14 cm^2 eV^2), q the shell's
occupation and P = dE from the total energies.  Writes no excitation table (the excitation fit, if any, stays).

--check: the largest relative error of the coefficient interpolated as mdrun does (linear in ln T of
ln <sigma v> + dE/T) against the fit, on a fine grid over 0.5 eV .. 30 keV, per element and file.
"""
import argparse
import math
import os

import numpy as np
from scipy.special import erfcx, exp1

COLL_K, COLEX_K = 1.0899e-6, 8.01e-8          # mcionize.c MCMD_COLL_K, MCMD_COLEX_K


def e1x(b):
    """e^b E1(b), stable for large b (continued fraction / asymptotics beyond 50)."""
    b = np.asarray(b, float)
    out = np.empty_like(b)
    small = b < 50
    out[small] = np.exp(b[small]) * exp1(b[small])
    x = b[~small]
    out[~small] = (1 / x) * (1 - 1 / x + 2 / x ** 2 - 6 / x ** 3 + 24 / x ** 4 - 120 / x ** 5)
    return out


def read_fit(path):
    rows = []
    for line in open(path):
        toks = [t.split() for t in line.strip().rstrip(";").split(";")]
        if not toks or len(toks[0]) != 3:
            continue
        ini = tuple(int(x) for x in toks[0])
        for t in toks[1:]:
            if len(t) == 8:
                rows.append((ini, tuple(int(x) for x in t[:3]), [float(x) for x in t[3:]]))
    return rows


def read_T(path):
    return {tuple(int(x) for x in l.split()[:3]): float(l.split()[3]) for l in open(path) if l.strip()}


def ln_kred_ion(par, dE, T):
    """ln of the ionization coefficient without its Boltzmann factor (Sampson), as mcionize.c."""
    c0, c1, c2, c3, s0 = par
    b = dE / T
    q = c0 + c1 * np.sqrt(np.pi * b) * erfcx(np.sqrt(b)) + (c2 * b + c3) * e1x(b)
    k = COLL_K * s0 * q / (dE * np.sqrt(T))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.isfinite(q) & (q > 0) & (k > 0), np.log(np.where(k > 0, k, 1)), -np.inf)


def ln_kred_exc(par, dE, T):
    """ln of the excitation coefficient without its Boltzmann factor (colex2), as mcionize.c."""
    c0, c1, c2, c3, a0 = par
    b, A = dE / T, a0 + 1.0
    q = c0 * e1x(b) + c1 + c3 * b / A + (c2 - c3 * b) * b * e1x(A * b)
    k = COLEX_K * q / np.sqrt(T)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.isfinite(q) & (q > 0), np.log(np.where(k > 0, k, 1)), -np.inf)


def fmt_ln(lk):
    """A number given by its natural log, as text with any exponent ('0' for -inf)."""
    if not np.isfinite(lk):
        return "0"
    l10 = lk / math.log(10)
    e = math.floor(l10)
    m = 10 ** (l10 - e)
    if m >= 9.9999999995:
        m, e = 1.0, e + 1
    return f"{m:.10f}e{e:d}"


def interp_like_mdrun(lnT, g, lt):
    """mdrun's evaluation: linear in ln T, held at the grid edges."""
    return np.interp(lt, lnT, g)


def ln_lotz(q, P, T, a_lotz):
    """ln of Lotz's Maxwellian ionization coefficient (cm^3/s), stable for P >> T: E1(x)/x = e^-x e1x(x)/x."""
    x = P / T
    return math.log(6.7e-7 * a_lotz * q) - 1.5 * np.log(T) - x + np.log(e1x(x) / x)


def write_lotz(a, el, TE, T, out, n):
    states = set()
    for line in open(os.path.join(a.dir, f"rate_transitions_{el}.txt")):
        t = line.split(";")[0].split()
        if len(t) == 3:
            states.add(tuple(int(x) for x in t))
    rows = []
    for s in sorted(states):
        if s not in TE:
            continue
        for i in range(3):
            if s[i] == 0:
                continue
            f = tuple(x - (j == i) for j, x in enumerate(s))
            if f != (0, 0, 0) and (f not in states or f not in TE):
                continue
            P = TE[s] - (0.0 if f == (0, 0, 0) else TE[f])
            if P > 0:
                rows.append((s, f, s[i], P))
    fit = os.path.join(a.dir, f"collisional_parameters_{el}.txt")
    nfit = len(read_fit(fit)) if os.path.exists(fit) else 0
    with open(os.path.join(out, f"collisional_ionization_table_{el}.txt"), "w") as fo:
        fo.write("T_e " + " ".join(repr(float(t)) for t in T) + "\n")
        for s, f, q, P in rows:
            fo.write(" ".join(map(str, s + f)) + " " + " ".join(fmt_ln(x) for x in ln_lotz(q, P, T, a.lotz_a)) + "\n")
    print(f"{el} collisional_ionization_table (Lotz, a = {a.lotz_a}): {len(rows)} transitions "
          f"(the CRETIN fit file has {nfit}), {n} temperatures")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--out")
    ap.add_argument("--tmin", type=float, default=0.1)
    ap.add_argument("--tmax", type=float, default=1e5)
    ap.add_argument("--per-decade", type=int, default=40)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--lotz", action="store_true", help="ionization from Lotz's formula (see above)")
    ap.add_argument("--lotz-a", type=float, default=4.5, help="Lotz a, 1e-14 cm^2 eV^2 (default 4.5)")
    a = ap.parse_args()
    out = a.out or a.dir
    os.makedirs(out, exist_ok=True)
    n = int(round(math.log10(a.tmax / a.tmin) * a.per_decade)) + 1
    T = np.logspace(math.log10(a.tmin), math.log10(a.tmax), n)
    lnT = np.log(T)
    Tf = np.logspace(math.log10(0.5), math.log10(3e4), 2000)          # check grid
    els = sorted(f[len("total_energies_"):-4] for f in os.listdir(a.dir) if f.startswith("total_energies_"))
    for el in els:
        TE = read_T(os.path.join(a.dir, f"total_energies_{el}.txt"))
        if a.lotz:
            write_lotz(a, el, TE, T, out, n)
            continue
        for fit, tab, fn in (("collisional_parameters", "collisional_ionization_table", ln_kred_ion),
                             ("collisional_excitation", "collisional_excitation_table", ln_kred_exc)):
            p = os.path.join(a.dir, f"{fit}_{el}.txt")
            if not os.path.exists(p):
                continue
            rows = read_fit(p)
            worst = 0.0
            with open(os.path.join(out, f"{tab}_{el}.txt"), "w") as fo:
                fo.write("T_e " + " ".join(repr(float(t)) for t in T) + "\n")   # full precision: mdrun uses these T in dE/T
                for ini, fin, par in rows:
                    dE = TE[ini] - (0.0 if fin == (0, 0, 0) else TE[fin])
                    lkr = fn(par, dE, T)
                    lk = lkr - dE / T                                   # ln <sigma v>
                    fo.write(" ".join(map(str, ini + fin)) + " " + " ".join(fmt_ln(x) for x in lk) + "\n")
                    if a.check and np.isfinite(lkr).all():
                        ref = fn(par, dE, Tf)
                        ok = np.isfinite(ref)
                        got = interp_like_mdrun(lnT, lkr, np.log(Tf))
                        if ok.any():
                            worst = max(worst, float(np.max(np.abs(np.expm1(got[ok] - ref[ok])))))
            print(f"{el} {tab}: {len(rows)} transitions, {n} temperatures {a.tmin:g}..{a.tmax:g} eV"
                  + (f"; interpolation vs fit, max relative error {worst:.2e} (0.5 eV .. 30 keV)" if a.check else ""))


if __name__ == "__main__":
    main()
