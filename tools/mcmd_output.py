#!/usr/bin/env python3
"""Read the binary outputs of the MolDStruct ionization module (simulation_output/), as numpy arrays or as text.

    mcmd_output.py transitions <simulation_output> [--type photo,auger] [--atom 12] [--tmin 20 --tmax 40]
        one line per event:  t[fs]  atom  element  type  K L M -> K L M  [donor]
    mcmd_output.py summary <simulation_output>
        events per type and element, and per type in 5 fs bins
    mcmd_output.py charges <simulation_output> [--every 10]
        mean charge per element over time, from charges_over_time.bin
    mcmd_output.py plot <simulation_output> [-o overview.png] [--gif [--every N]] [--extract]
        one overview figure of everything the run wrote (panels for files that are missing are left out);
        default <run dir>/overview.png.  --gif: also explosion.gif, three views along the principal axes of the
        initial heavy-atom gyration tensor, atoms coloured by charge (reads the run's .trr; no GROMACS needed).
        --extract: also <run>/extracted/ with transitions.txt and charges_over_time.txt

<simulation_output> may be omitted (default: the current directory) or be the run directory holding it.

In Python:
    import mcmd_output as mo
    ev = mo.read_transitions("run/simulation_output")      # structured array: t, atom, partner, type, before, after
    el = mo.elements("run/simulation_output")              # element per atom (from masses.bin)
    q  = mo.read_charges("run/simulation_output")          # (frames, atoms) float32

transitions.bin: 16-byte header ("MCMDTRN1", int32 version, int32 record size), then 20-byte little-endian records:
float32 t [fs], int32 atom (0-based), int32 partner (charge transfer: the donor; else -1), uint8 type, uint8 K L M
before, uint8 K L M after, 1 byte padding.  Atoms print 1-based (as in the .gro/.pdb).
"""
import argparse
import os
import sys

import numpy as np

TYPES = ["auger", "fluorescence", "photo", "coll_ion", "recombination", "excitation", "deexcitation",
         "charge_transfer"]
REC = np.dtype([("t", "<f4"), ("atom", "<i4"), ("partner", "<i4"), ("type", "u1"), ("before", "u1", 3),
                ("after", "u1", 3), ("pad", "u1")])
# standard atomic weights; an atom is the nearest within 0.5 (as mdrun does: Ar 39.948 vs Ca 40.078), deuterium is H
STD_MASS = {"H": 1.008, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
            "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
            "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Mn": 54.938, "Fe": 55.845,
            "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Se": 78.971, "I": 126.904}


def _path(p, name):
    return os.path.join(p, name) if os.path.isdir(p) else p


def _read_text_log(f):
    """The text electronic_transition_log.txt of builds before 2026-10-27, as the same structured array (charge
    transfers there name no donor: partner -1)."""
    import re
    rows, t = [], 0.0
    for line in open(f):
        if line.startswith("Making"):
            t = float(line.split()[6]) * 1e3
        elif line.startswith("Transition of type"):
            w = line.split()
            if len(w) >= 16:
                rows.append((t, int(w[5]) - 1, -1, int(w[3]), tuple(map(int, w[9:12])), tuple(map(int, w[13:16])), 0))
            else:   # main mc-md / older builds log the type only: atom -1, states unknown (0 0 0)
                rows.append((t, -1, -1, int(w[3]), (0, 0, 0), (0, 0, 0), 0))
        elif line.startswith("Charge transfer"):
            tc = float(re.search(r"t: ([0-9.eE+-]+)", line).group(1)) * 1e3
            b = re.search(r"Previous state: \[(\d+), (\d+), (\d+)\]", line)
            n = re.search(r"New state: \[(\d+), (\d+), (\d+)\]", line)
            if b and n:
                rows.append((tc, -1, -1, 7, tuple(map(int, b.groups())), tuple(map(int, n.groups())), 0))
    return np.array(rows, REC)


def read_transitions(path):
    """Events of transitions.bin (a file or a simulation_output directory) as a structured array; falls back to
    the older text electronic_transition_log.txt (charge transfers there have atom -1: not logged)."""
    f = _path(path, "transitions.bin")
    if not os.path.exists(f) and os.path.isdir(path) and os.path.exists(os.path.join(path,
                                                                                 "electronic_transition_log.txt")):
        return _read_text_log(os.path.join(path, "electronic_transition_log.txt"))
    with open(f, "rb") as fh:
        magic = fh.read(8)
        if magic != b"MCMDTRN1":
            raise ValueError(f"{f}: not a transitions.bin (magic {magic!r})")
        version, recsize = np.frombuffer(fh.read(8), "<i4")
        if recsize != REC.itemsize:
            raise ValueError(f"{f}: record size {recsize}, this reader knows {REC.itemsize} (version {version})")
        data = fh.read()
    n = len(data) // REC.itemsize          # a run killed mid-write may leave a partial record
    return np.frombuffer(data[:n * REC.itemsize], REC)


def masses(path):
    return np.fromfile(_path(path, "masses.bin"), "<f4")


def element_of(m):
    if 1.5 <= m < 2.5:
        return "H"
    el, d = min(((e, abs(m - s)) for e, s in STD_MASS.items()), key=lambda t: t[1])
    return el if d < 0.5 else f"m{m:.0f}"


def elements(path):
    return np.array([element_of(float(m)) for m in masses(path)])


def read_charges(path):
    """charges_over_time.bin as (frames, atoms) float32 (one frame every mcmd-charge-output-stride steps)."""
    n = len(masses(path))
    q = np.fromfile(_path(path, "charges_over_time.bin"), "<f4")
    return q[:len(q) // n * n].reshape(-1, n)


def read_trr(path):
    """Frames of a GROMACS .trr (XDR, single or double precision) without GROMACS: yields dicts with step,
    t [ps], box (3x3 nm), x and v (natoms x 3, nm and nm/ps; None if the frame has none)."""
    import struct
    with open(path, "rb") as fh:
        data = fh.read()
    o, n = 0, len(data)
    while o + 8 <= n:
        magic, = struct.unpack_from(">i", data, o)
        if magic != 1993:
            raise ValueError(f"{path}: bad frame magic {magic} at byte {o}")
        slen, sl = struct.unpack_from(">ii", data, o + 4)
        o += 12 + (sl + 3) // 4 * 4
        (ir, e, box, vir, pres, top, sym, xs, vs, fs, nat, step, nre) = struct.unpack_from(">13i", data, o)
        o += 52
        rs = box // 9 if box else (xs // (3 * nat) if xs else 4)
        rf = ">d" if rs == 8 else ">f"
        t, = struct.unpack_from(rf, data, o)
        o += 2 * rs                        # t, lambda
        o += ir + e
        bx = None
        if box:
            bx = np.frombuffer(data, ">f8" if rs == 8 else ">f4", 9, o).reshape(3, 3).astype(float)
        o += box + vir + pres + top + sym
        arr = {}
        for key, size in (("x", xs), ("v", vs), ("f", fs)):
            arr[key] = (np.frombuffer(data, ">f8" if rs == 8 else ">f4", nat * 3, o).reshape(nat, 3).astype(float)
                        if size else None)
            o += size
        yield {"step": step, "t": t, "box": bx, "x": arr["x"], "v": arr["v"]}


def _run_dir(d):
    return os.path.dirname(os.path.abspath(d.rstrip("/")))


def _find_trr(d):
    run = _run_dir(d)
    cands = sorted(f for f in os.listdir(run) if f.endswith(".trr"))
    if not cands:
        return None
    return os.path.join(run, "e.trr" if "e.trr" in cands else cands[0])


def _charges_at(d, times_fs, t_end):
    """Per-atom charges at the given times (nearest frame of charges_over_time.bin, whose frames span 0..t_end)."""
    q = read_charges(d)
    tq = np.linspace(0, t_end, len(q)) if len(q) > 1 else np.zeros(1)
    return q[np.abs(tq[None, :] - np.asarray(times_fs)[:, None]).argmin(1)]


def cmd_transitions(a):
    ev = read_transitions(a.dir)
    el = elements(a.dir)
    keep = np.ones(len(ev), bool)
    if a.type:
        want = [TYPES.index(t.strip()) for t in a.type.split(",")]
        keep &= np.isin(ev["type"], want)
    if a.atom:
        keep &= ev["atom"] == a.atom - 1
    if a.tmin is not None:
        keep &= ev["t"] >= a.tmin
    if a.tmax is not None:
        keep &= ev["t"] <= a.tmax
    print("# t[fs]    atom el  type              K L M -> K L M  donor")
    out = sys.stdout
    for e in ev[keep]:
        b, f = e["before"], e["after"]
        donor = f"  {e['partner'] + 1} {el[e['partner']]}" if e["partner"] >= 0 else ""
        out.write(f"{e['t']:9.4f} {e['atom'] + 1:6d} {el[e['atom']]:2s}  {TYPES[e['type']]:16s}  "
                  f"{b[0]} {b[1]} {b[2]} -> {f[0]} {f[1]} {f[2]}{donor}\n")


def cmd_summary(a):
    ev = read_transitions(a.dir)
    el = elements(a.dir)
    els = sorted(set(el), key=lambda e: (e != "H", e))
    print(f"{len(ev)} events, {len(el)} atoms, t {ev['t'].min():.2f}..{ev['t'].max():.2f} fs\n")
    print(f"{'type':16s} {'total':>7s} " + " ".join(f"{e:>7s}" for e in els) + "   per atom")
    for k, name in enumerate(TYPES):
        m = ev["type"] == k
        if m.any():
            print(f"{name:16s} {m.sum():7d} " + " ".join(f"{(m & (el[ev['atom']] == e)).sum():7d}" for e in els)
                  + f"   {m.sum() / len(el):.3f}")
    bins = np.arange(0, np.ceil(ev["t"].max() / 5) * 5 + 5, 5)
    print(f"\nevents per 5 fs bin (bin start in fs)\n{'type':16s} " + " ".join(f"{b:5.0f}" for b in bins[:-1]))
    for k, name in enumerate(TYPES):
        m = ev["type"] == k
        if m.any():
            print(f"{name:16s} " + " ".join(f"{c:5d}" for c in np.histogram(ev['t'][m], bins)[0]))


def cmd_charges(a):
    q = read_charges(a.dir)
    el = elements(a.dir)
    els = sorted(set(el), key=lambda e: (e != "H", e))
    print("# frame  mean " + " ".join(f"{e:>6s}" for e in els) + "   (one frame per mcmd-charge-output-stride steps)")
    for i in range(0, len(q), a.every):
        print(f"{i:7d} {q[i].mean():5.3f} " + " ".join(f"{q[i][el == e].mean():6.3f}" for e in els))


# ---------------------------------------------------------------- overview figure
INK, MUTED, GRIDC, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#efece6"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#c2417f", "#8a6a00", "#6a51a3", "#008b8b", "#7a7a7a"]   # fixed order


def _table(path, names=None):
    """Columns of a '#'-headed text log as a dict (header names, or the names given) or None."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    head = [l for l in open(path) if l.startswith("#")][:1]
    a = np.loadtxt(path, comments="#", ndmin=2)
    if a.size == 0:
        return None
    if names is None:
        names = head[0][1:].split()[:a.shape[1]] if head else [f"c{i}" for i in range(a.shape[1])]
    return {n: a[:, i] for i, n in enumerate(names)}


def _time_per_step(d):
    """MD time step [fs] from energy_ledger.txt (step and t columns), for logs that count steps."""
    led = _table(os.path.join(d, "energy_ledger.txt"))
    if led is None or "step" not in led:
        return None
    k = np.argmax(led["step"] > 0)
    return led["t[ps]"][k] * 1e3 / led["step"][k] if led["step"][k] > 0 else None


def cmd_plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = a.dir.rstrip("/")
    run = os.path.dirname(os.path.abspath(d))
    out = a.o or os.path.join(run, "overview.png")
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "legend.frameon": False,
                         "lines.linewidth": 1.4})
    el = elements(d) if os.path.exists(os.path.join(d, "masses.bin")) else None
    els = sorted(set(el), key=lambda e: (e != "H", e)) if el is not None else []
    ecol = {e: CAT[i % len(CAT)] for i, e in enumerate(els)}
    pp = np.loadtxt(os.path.join(d, "pulse_profile.txt"), ndmin=2) if os.path.exists(
        os.path.join(d, "pulse_profile.txt")) else None
    mq = np.loadtxt(os.path.join(d, "mean_charge_vs_time.txt"), ndmin=2) if os.path.exists(
        os.path.join(d, "mean_charge_vs_time.txt")) else None
    t_end = mq[-1, 0] * 1e3 if mq is not None else None
    scr, col, pb, led = (_table(os.path.join(d, f)) for f in
                         ("screening_log.txt", "collisional_log.txt", "pb_log.txt", "energy_ledger.txt"))
    box = _table(os.path.join(d, "pbg_box.txt"), ["t[ps]", "R[nm]", "box[R]", "dx[nm]", "outside", "rebuilt", "N"])
    dts = _time_per_step(d)
    ev = read_transitions(d) if any(os.path.exists(os.path.join(d, f)) for f in
                                    ("transitions.bin", "electronic_transition_log.txt")) else None
    panels = []

    def panel(fn):
        panels.append(fn)
        return fn

    if pp is not None:
        @panel
        def _(ax):
            ax.plot(pp[:, 0] * 1e3, pp[:, 1], color=CAT[0])
            ax.set_title("pulse intensity (photons / (cm² s))")

    if mq is not None:
        @panel
        def _(ax):
            ax.plot(mq[:, 0] * 1e3, mq[:, 1], color=INK, label="all atoms")
            if el is not None and os.path.exists(os.path.join(d, "charges_over_time.bin")):
                q = read_charges(d)
                tq = np.linspace(0, t_end, len(q)) if len(q) > 1 else np.zeros(1)
                for e in els:
                    ax.plot(tq, q[:, el == e].mean(1), color=ecol[e], lw=1.1, label=e)
            ax.legend(fontsize=7, ncol=2)
            ax.set_title("mean ion charge")

    if os.path.exists(os.path.join(d, "charges.txt")) and el is not None:
        @panel
        def _(ax):
            qf = np.loadtxt(os.path.join(d, "charges.txt"), ndmin=2)[:, 1]
            qmax = int(np.ceil(qf.max())) + 1
            for e in els:
                h = np.bincount(np.rint(qf[el == e]).astype(int), minlength=qmax)[:qmax] / max(1, (el == e).sum())
                ax.plot(np.arange(qmax), h, color=ecol[e], marker="o", ms=4, lw=1.1, label=e)
            ax.legend(fontsize=7, ncol=2)
            ax.set_xlabel("final charge"); ax.set_title("final charge distribution (fraction of each element)")
            ax.set_xticks(range(0, qmax, 2 if qmax > 12 else 1))

    if ev is not None and len(ev):
        tb = np.arange(0, (t_end or ev["t"].max()) + 0.5, 0.5)
        present = [k for k in range(len(TYPES)) if (ev["type"] == k).any()]

        @panel
        def _(ax):
            for k in present:
                h = np.histogram(ev["t"][ev["type"] == k], tb)[0] / 0.5
                ax.plot(tb[:-1] + 0.25, h, color=CAT[k], label=TYPES[k].replace("_", " "))
            ax.legend(fontsize=7)
            ax.set_title("events per fs, by type")

        @panel
        def _(ax):
            n = len(el) if el is not None else 1
            for k in present:
                m = ev["type"] == k
                ax.plot(np.sort(ev["t"][m]), np.arange(1, m.sum() + 1) / n, color=CAT[k],
                        label=f"{TYPES[k].replace('_', ' ')} {m.sum() / n:.2f}")
            ax.legend(fontsize=7)
            ax.set_title("events per atom, cumulative (legend: total)")

    if scr is not None and dts:
        ts = scr["step"] * dts

        @panel
        def _(ax):
            ax.plot(ts, scr["escaped"], color=CAT[1], label="escaped")
            ax.plot(ts, scr["trapped"], color=CAT[0], label="trapped")
            if "evaporated" in scr:
                ax.plot(ts, scr["evaporated"], color=CAT[3], lw=1, label="evaporated (cumulative)")
            ax.legend(fontsize=7)
            ax.set_title("electrons")

        @panel
        def _(ax):
            ax.plot(ts, scr["T_e[eV]"], color=CAT[1], label="T$_e$")
            ax.plot(ts, scr["mean_E_trapped[eV]"], color=CAT[0], label="mean trapped-electron energy")
            ax.legend(fontsize=7)
            ax.set_title("trapped electrons: temperature, energy (eV)")

        @panel
        def _(ax):
            ax.plot(ts, scr["R[nm]"], color=INK, label="sample R")
            if box is not None:
                ax.plot(box["t[ps]"] * 1e3, box["R[nm]"] * box["box[R]"], color=CAT[0], lw=1,
                        label="PB grid half-width")
            ax.legend(fontsize=7)
            ax.set_title("size (nm)")

    if col is not None:
        @panel
        def _(ax):
            ax.plot(col["t[ps]"] * 1e3, col["n_coll_ion"], color=CAT[3], label="collisional ionizations")
            ax.plot(col["t[ps]"] * 1e3, col["n_coll_rec"], color=CAT[4], label="recombinations")
            ax.legend(fontsize=7)
            ax.set_title("electron-impact processes, cumulative")

    if led is not None:
        @panel
        def _(ax):
            for k, (c, lab) in enumerate((("J_ii[eV]", "ion-ion jumps of charge changes"),
                                          ("J_ie[eV]", "ion-electron jumps"), ("W[eV]", "electron-field work on ions"),
                                          ("E_kin_e[eV]", "trapped electrons' kinetic energy"))):
                if c in led:
                    ax.plot(led["t[ps]"] * 1e3, led[c] / 1e3, color=CAT[k], label=lab)
            ax.legend(fontsize=7)
            ax.set_title("energy ledger, cumulative (keV)")

    if pb is not None:
        @panel
        def _(ax):
            ax.plot(pb["t[ps]"] * 1e3, pb["wall[ms]"] / 1e3, color=INK, lw=1)
            ax.set_title(f"PB solve: wall time per solve (s); {len(pb['t[ps]'])} solves, "
                         f"{pb['wall[ms]'].sum() / 6e4:.1f} min")

        if box is not None:
            @panel
            def _(ax):
                ax.plot(box["t[ps]"] * 1e3, box["N"], color=INK)
                ax.set_title(f"PB grid: cells per side (dx {box['dx[nm]'][0]:.3f} → {box['dx[nm]'][-1]:.3f} nm)")

    if not panels:
        sys.exit(f"nothing to plot in {d}")
    nc = 3
    nr = (len(panels) + nc - 1) // nc
    fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.8 * nr), squeeze=False)
    axs = axs.ravel()
    on = None
    if pp is not None and pp[:, 1].max() > 0:
        half = pp[pp[:, 1] >= 0.5 * pp[:, 1].max(), 0] * 1e3
        on = (half.min(), half.max())
    for ax, fn in zip(axs, panels):
        fn(ax)
        ax.title.set_fontsize(8.5); ax.title.set_color(INK)
        if "final charge" not in ax.get_title():
            ax.set_xlabel("time (fs)")
            if on:
                ax.axvspan(*on, color=BAND, zorder=0, lw=0)
            if t_end:
                ax.set_xlim(0, t_end)
    for ax in axs[len(panels):]:
        ax.axis("off")
    stats = ""
    if ev is not None and el is not None:
        stats = f"; {len(el)} atoms, {len(ev)} events"
    fig.suptitle(f"{run}{stats}" + (f"; final mean charge {mq[-1, 1]:.2f}" if mq is not None else "") +
                 ("  (band: pulse FWHM)" if on else ""), x=0.01, ha="left", fontsize=9, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(out)
    if a.gif:
        make_gif(d, os.path.join(os.path.dirname(out), "explosion.gif"), every=a.every)
    if a.extract:
        extract(d, os.path.join(run, "extracted"))


def make_gif(d, out, every=1, fps=10):
    """explosion.gif: three views along the principal axes of the initial heavy-atom gyration tensor (longest
    first), atoms coloured by charge, H smaller, zooming out with the expansion (never back in)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    trr = _find_trr(d)
    if trr is None:
        print(f"no .trr next to {d}: no gif"); return
    frames = [f for i, f in enumerate(read_trr(trr)) if i % every == 0 and f["x"] is not None]
    t = np.array([f["t"] for f in frames]) * 1e3
    X = np.array([f["x"] for f in frames])
    m, el = masses(d), elements(d)
    h = el == "H"
    mq = np.loadtxt(os.path.join(d, "mean_charge_vs_time.txt"), ndmin=2)
    Q = _charges_at(d, t, mq[-1, 0] * 1e3)
    w = m * ~h
    c0 = (X[0] * w[:, None]).sum(0) / w.sum()
    G = ((X[0] - c0).T * w) @ (X[0] - c0) / w.sum()
    ev, V = np.linalg.eigh(G)
    V = V[:, ::-1]                                       # longest axis first
    rg = np.sqrt(ev[::-1])
    P = np.array([(x - (x * w[:, None]).sum(0) / w.sum()) @ V for x in X])
    lims = np.maximum.accumulate([np.percentile(np.abs(p[~h]), 99.5) * 1.15 for p in P])
    order = np.argsort(~h)                               # H first, drawn below
    size = np.where(h, 0.8, 3.0) * min(1.0, (2000 / len(m)) ** 0.5)
    size = np.maximum(size[order], 0.15)
    vmax = max(1.0, float(np.percentile(Q[-1], 99.5)))
    pp = np.loadtxt(os.path.join(d, "pulse_profile.txt"), ndmin=2) if os.path.exists(
        os.path.join(d, "pulse_profile.txt")) else None
    on = None
    if pp is not None and pp[:, 1].max() > 0:
        half = pp[pp[:, 1] >= 0.5 * pp[:, 1].max(), 0] * 1e3
        on = (half.min(), half.max())
    views = [((1, 2), "view along axis 1 (longest)"), ((0, 2), "view along axis 2"), ((0, 1), "view along axis 3")]
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.9), dpi=90)
    fig.patch.set_facecolor("white")
    sc, bars = [], []
    for ax, ((i, j), title) in zip(axs, views):
        sc.append(ax.scatter(P[0][order, i], P[0][order, j], s=size, c=Q[0][order], cmap="plasma", vmin=0,
                             vmax=vmax, lw=0, rasterized=True))
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color("#c8c6c0")
        ax.set_title(title, fontsize=9, color=MUTED)
        bars.append((ax.plot([], [], color=INK, lw=2)[0], ax.text(0, 0, "", ha="center", fontsize=8, color=INK)))
    cb = fig.colorbar(sc[0], cax=fig.add_axes([0.3, 0.12, 0.4, 0.03]), orientation="horizontal")
    cb.set_label("ion charge", fontsize=9)
    head = fig.suptitle("", x=0.01, ha="left", fontsize=10, color=INK)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.2, wspace=0.04)
    name = "/".join(_run_dir(d).split(os.sep)[-2:])
    axes_txt = f"R_g axes {rg[0]:.2f} / {rg[1]:.2f} / {rg[2]:.2f} nm"

    def draw(k):
        lim = lims[k]
        bar = next(v for v in (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500) if v >= 0.3 * lim)
        for s_, ((i, j), _), ax, (ln, tx) in zip(sc, views, axs, bars):
            s_.set_offsets(np.c_[P[k][order, i], P[k][order, j]])
            s_.set_array(Q[k][order])
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ln.set_data([-0.9 * lim, -0.9 * lim + bar], [-0.9 * lim] * 2)
            tx.set_position((-0.9 * lim + bar / 2, -0.85 * lim)); tx.set_text(f"{bar:g} nm")
        pulse = "  pulse on (FWHM)" if on and on[0] <= t[k] <= on[1] else ""
        head.set_text(f"{name}: {len(m)} atoms ({axes_txt})    t = {t[k]:5.1f} fs    "
                      f"mean charge {Q[k].mean():.2f}{pulse}")
        return sc

    FuncAnimation(fig, draw, frames=len(t), blit=False).save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(out)


def extract(d, outdir):
    """Human-readable copies: transitions.txt and charges_over_time.txt (time, then one column per atom)."""
    n = len(masses(d))
    el = elements(d)
    os.makedirs(outdir, exist_ok=True)
    mq = np.loadtxt(os.path.join(d, "mean_charge_vs_time.txt"), ndmin=2)
    t_end = mq[-1, 0] * 1e3
    if any(os.path.exists(os.path.join(d, f)) for f in ("transitions.bin", "electronic_transition_log.txt")):
        ev = read_transitions(d)
        with open(os.path.join(outdir, "transitions.txt"), "w") as fo:
            fo.write("# t[fs] atom(1-based) element type K L M -> K L M  donor(1-based, charge transfer)\n")
            for e in ev:
                b, f = e["before"], e["after"]
                donor = f"  {e['partner'] + 1} {el[e['partner']]}" if e["partner"] >= 0 else ""
                at = f"{e['atom'] + 1:6d} {el[e['atom']]:2s}" if e["atom"] >= 0 else "     ? ? "
                fo.write(f"{e['t']:9.4f} {at}  {TYPES[e['type']]:16s}  {b[0]} {b[1]} {b[2]} -> {f[0]} {f[1]} {f[2]}"
                         f"{donor}\n")
        print(os.path.join(outdir, "transitions.txt"))
    q = read_charges(d)
    tq = np.linspace(0, t_end, len(q)) if len(q) > 1 else np.zeros(1)
    np.savetxt(os.path.join(outdir, "charges_over_time.txt"), np.c_[tq, q], fmt="%.4g",
               header="t[fs] then the charge of atoms 1.." + str(n) + " (one row per charges_over_time.bin frame)")
    print(os.path.join(outdir, "charges_over_time.txt"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("transitions"); p.add_argument("dir", nargs="?", default=".", help="simulation_output or the run directory holding it (default: .)")
    p.add_argument("--type", help="comma list of: " + ", ".join(TYPES))
    p.add_argument("--atom", type=int, help="1-based atom number")
    p.add_argument("--tmin", type=float); p.add_argument("--tmax", type=float)
    p.set_defaults(f=cmd_transitions)
    p = sub.add_parser("summary"); p.add_argument("dir", nargs="?", default=".", help="simulation_output or the run directory holding it (default: .)"); p.set_defaults(f=cmd_summary)
    p = sub.add_parser("charges"); p.add_argument("dir", nargs="?", default=".", help="simulation_output or the run directory holding it (default: .)"); p.add_argument("--every", type=int, default=10)
    p.set_defaults(f=cmd_charges)
    p = sub.add_parser("plot"); p.add_argument("dir", nargs="?", default=".", help="simulation_output or the run directory holding it (default: .)"); p.add_argument("-o", help="output png (default <run>/overview.png)")
    p.add_argument("--gif", action="store_true", help="also explosion.gif (three principal-axis views, by charge)")
    p.add_argument("--every", type=int, default=1, help="gif: every Nth trr frame")
    p.add_argument("--extract", action="store_true",
                   help="also <run>/extracted/: transitions.txt, charges_over_time.txt")
    p.set_defaults(f=cmd_plot)
    a = ap.parse_args()
    # the run directory works too: use its simulation_output
    if os.path.isdir(os.path.join(a.dir, "simulation_output")):
        a.dir = os.path.join(a.dir, "simulation_output")
    try:
        a.f(a)
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    main()
