/*
 * mcionize.c
 *
 * The MolDStruct hybrid Monte Carlo / MD ionization module.  See mcionize.h.
 *
 * This file is part of the MolDStruct modifications to GROMACS 4.5.4 and is
 * distributed under the same licence as GROMACS.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fftw3.h>
#if defined(GMX_THREADS)
#include <pthread.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>

#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && \
      !defined __CYGWIN__ && !defined __CYGWIN32__)
/* getpid(), used to decorrelate the RNG seed between jobs */
#include <sys/types.h>
#include <unistd.h>
#endif

#include "atomic_data.h"
#include "gmx_fatal.h"
#include "mcionize.h"
#include "mvdata.h"
#include "network.h"
#include "partdec.h"
#include "smalloc.h"
#include "vec.h"
#include "physics.h"

/* Master switch for the force-field modifications in bondfree.c and
 * nonbonded.c.  Those files pick it up as an extern; it is only assigned
 * inside this module, so without -ionize it stays 0 and the force field is
 * unmodified stock GROMACS 4.5.4. */
/* Defined in gmxlib/bondfree.c, where it is read.  Set to 1 by
 * mcionize_init() when the module is active, i.e. under mdrun -ionize, which
 * is what selects the MolDStruct force-field terms. */
extern int moldstruct_altered_ff;

/* ---------------------------------------------------------------- *
 * Physical constants and model parameters
 * ---------------------------------------------------------------- */

#define MCMD_BOLTZMANN        1.380649e-23  /* J/K                       */
#define MCMD_ELEMENTARY_CHARGE 1.60217663e-19 /* C                       */
#define MCMD_EV_TO_JOULE      1.60218e-19   /* J/eV                      */
#define MCMD_KELVIN_TO_EV     8.61732814974493e-5
#define MCMD_PI               3.1415

/* e^2 / (4 pi eps0), in eV nm.
 *
 * The classical over-the-barrier expression for R_crit below is the
 * atomic-unit form: a binding energy in Hartree gives a radius in bohr.  This
 * module carries binding energies in eV (straight out of energy_levels_X.txt)
 * and positions in nm, so the expression needs the Coulomb constant in those
 * units to be dimensionally a length at all.  Without it R_crit came out a
 * factor 1.44 too small - and the capture volume a factor three too small -
 * for every pair.  H donating to a bare proton, for instance, gave 0.220 nm
 * against the correct 0.317 nm (6 bohr). */
#define MCMD_COULOMB_EV_NM 1.43996454

/* Charge transfer is only considered between atoms closer than this, in nm.
 *
 * This is a pre-filter, not physics: the real criterion is R_crit, and a pair
 * is only ever used if R_crit exceeds its separation.  Measured over a full
 * lysozyme run at 1000 eV, the largest separation of any pair actually used
 * was 0.509 nm (median 0.204, 99.9th percentile 0.478) - but that was with
 * R_crit missing the conversion above, so those numbers all scale by 1.44,
 * putting the observed maximum near 0.73 nm.  The cutoff used to be 6 nm,
 * which is larger than the whole sample, so it rejected nothing; 1 nm was
 * then chosen for headroom over the uncorrected maximum and is too tight
 * once it is corrected.  1.5 nm restores roughly the same factor two of
 * headroom at 3.4x the search volume rather than the 27x a 3 nm cutoff would
 * cost.
 *
 * The tail of the distribution reaches further than that in principle: a
 * donor in a weakly bound excited state next to a heavily charged acceptor
 * can push R_crit past 2 nm.  Those states are rare, and the warning in
 * mcionize_done() reports the largest separation actually used, so raise the
 * cutoff if it fires. */
#define MCMD_CT_CUTOFF_NM 1.5

/* Charge transfer dies out long before the run does: the pulse stops, the
 * Auger cascades that keep manufacturing charge differences finish, and the
 * fragments drift apart.  After that the pass costs a fixed 27 cell lookups
 * per atom per step and finds nothing - measured on a 25000 step lysozyme
 * run, 35 s of the 92 s spent in charge transfer came after the last
 * meaningful transfer.
 *
 * So once nothing has happened for mcmd-charge-transfer-idle steps, stop
 * running the pass every step and only look every
 * mcmd-charge-transfer-recheck steps.  It is not switched off
 * outright: Auger decay continues and can create a fresh donor at any time,
 * and the recheck picks that up.  The default window is deliberately generous
 * because late transfers are sparse rather than absent - at one transfer per
 * fifty steps, a short window would trip constantly.
 *
 * Keyed on observed activity, not on time since the pulse peak.  A pulse
 * based window is the wrong variable: measured on the same run, 7.5% of all
 * transfers happened more than 5 sigma after the peak and 0.7% more than
 * 50 sigma after, when the intensity was already 200 orders of magnitude
 * down. */
/* Smallest energy gain, in eV, that counts as downhill under
 * mcmd-charge-transfer-downhill.  The two sides of that test come from
 * different tables - the acceptor's from total_energies_X.txt, the donor's
 * from energy_levels_X.txt - and for two atoms in the same state they agree
 * only to rounding (C ground state: 11.26030000000003 against 11.2603).  A
 * strict comparison therefore passed exact resonances in both directions,
 * and a C0/C+ pair swapped its electron every step: 99.7% of 3.6 million
 * transfers on a 150 fs lysozyme run, against ~3400 with this margin, with
 * no change in charges, kinetic energies or ion directions. */
#define MCMD_DOWNHILL_TOL_EV 1e-6

#define MCMD_CT_IDLE_STEPS_DEFAULT 2000
#define MCMD_CT_RECHECK_DEFAULT     100

/* Per-step output is written through persistent, generously buffered handles
 * rather than reopening each file every step.  The files are flushed every
 * MCMD_OUTPUT_FLUSH_STEPS so that tailing them during a run still works and a
 * crash loses at most that many steps of analysis output; the simulation
 * state itself is in the checkpoint, not here. */
#define MCMD_OUTPUT_FLUSH_STEPS 1000
#define MCMD_OUTPUT_BUFSIZE     (64 * 1024)

/* FWHM = 2 sqrt(2 ln 2) sigma, for converting the mdp's mcmd-pulse-fwhm into
 * the standard deviation the Gaussian maths uses. */
#define MCMD_FWHM_PER_SIGMA 2.3548200450309493

/* eV per particle -> kJ/mol */
#define MCMD_EV_TO_KJMOL 96.4853321

/* Energy deposited in the plasma by one trapped electron, in eV. */
#define MCMD_TRAPPED_ELECTRON_EV 25.0


/* Energy difference assumed between collisional levels, in eV. */
#define MCMD_COLL_DELTA_E 200.0

/* PB grid constants, settled by the box and speed tests of 2026-10-23..26
 * (environment switches until the 2026-10-27 clean-up):
 * the adaptive box is at least MCMD_PBG_BOXMIN R and shrinks only to where at
 * most MCMD_PBG_SCUT of the electrons are outside (safe shrink), resizing to
 * MCMD_PBG_BOXMARG x the wanted or safe size; the kernel smoothing is
 * s = max(MCMD_PBG_SFAC dx, s_e); outside ions use sub-lattices of stride
 * dx <= distance / MCMD_PBG_OUTLAT_D, except those within MCMD_PBG_MARG cells
 * of the box, which go into one larger FFT. */
#define MCMD_PBG_BOXMIN   1.3
#define MCMD_PBG_SCUT     0.001
#define MCMD_PBG_BOXMARG  1.05
#define MCMD_PBG_SFAC     2.0
#define MCMD_PBG_OUTLAT_D 3.0
#define MCMD_PBG_MARG     16

/* Output directory and files, relative to the run directory. */
#define MCMD_OUTDIR "./simulation_output"

/* Restart file written at the end of every run.  Holds everything the
 * module cannot re-derive from the trajectory: the electronic state, the
 * free-electron plasma, and the sample geometry the plasma volume is
 * measured against. */
#define MCMD_STATE_FILE  "mcionize_state.bin"
#define MCMD_STATE_MAGIC 0x4d434d44 /* "MCMD" */
#define MCMD_STATE_VERSION 3 /* 2: sample radius and R_g about the COM; 3: written at the checkpoints, with the
                              * PB electrons, the step counters, rand() and the output lengths */
#define MCMD_STATE_PREV  "mcionize_state_prev.bin"

/* ---------------------------------------------------------------- *
 * Module state
 * ---------------------------------------------------------------- */

/* One atom as the charge-transfer candidate scan sees it.  Kept as a struct
 * rather than parallel arrays because every field is read together, so a
 * candidate costs one cache line rather than four streams.  q and shell are
 * copies and are refreshed by mcmd_cell_sync() whenever a transfer changes
 * them mid-pass. */
typedef struct
{
    real x, y, z;   /* position                                   */
    real q;         /* net charge                                 */
    int  shell;     /* outermost occupied shell, -1 if stripped   */
    int  idx;       /* global atom index                          */
} t_mcmd_cell_atom;

/* One collisional channel out of a state (mcmd-collisions). */
typedef struct
{
    int    final_state[MCMD_NSHELL];
    int    recomb;     /* 0: ionization (MC type 3), 1: three-body recombination (type 4) */
    int    excite;     /* 0: ionization/recombination; 1: excitation (MC type 5),
                        * 2: de-excitation (type 6), mcmd-coll-excitation */
    double c[4], s0;   /* Sampson coefficients of the ionizing direction            */
    double a0;         /* colex: c[0..3], a0 of the CRETIN 'colex2' fit             */
    double dE;         /* transition energy, eV (total energies)                     */
    double gratio;     /* recombination: g(final, more electrons) / g(initial)       */
    double rate;       /* current rate, 1/s, set once per step                       */
    /* tabulated (collisional_*_table_X.txt): g = ln <sigma v> + dE/T on the grid lnT, the
     * log of the coefficient without its Boltzmann factor (what the fits give); NULL = fit */
    const double *g, *lnT;
    int           nT;
    double        lotz;  /* Lotz ionization: 6.7e-7 a q (cm^3 s^-1 eV^3/2); 0 = not a Lotz channel */
} t_mcmd_cchan;

/* Continuation: a named scalar of the state file, and the file as read */
typedef struct
{
    char   name[40];
    double v;
} t_mcmd_kv;
#define MCMD_KV_MAX 192

typedef struct t_mcmd_rst
{
    int             version;
    gmx_large_int_t step;
    double          geom[7];
    int             counters[4];
    int             nkv;
    t_mcmd_kv      *kv;
    char            rng[128];
    int            *conf, *gs; /* natoms x MCMD_NSHELL */
    double         *q;
    char            path[512];
} t_mcmd_rst;

struct t_mcionize
{
    /* Configuration, cached from the inputrec */
    gmx_bool bChargeTransfer; /* mcmd-charge-transfer                  */
    gmx_bool bAllowHCT;       /* mcmd-allow-h-ct                       */
    gmx_bool ct_h_donor_only; /* mcmd-ct-h-donor-only = 1: only hydrogen
                               * donates, as in the published version  */
    gmx_bool bDownhill;       /* mcmd-charge-transfer-downhill         */
    gmx_bool bInitCharges;    /* mcmd-initial-charges                  */
    gmx_bool bDetailedOutput; /* mcmd-detailed-output                */

    /* Pulse, in SI */
    double t_mean_s;   /* pulse peak time, from mcmd-pulse-peak-time */
    double width_s;    /* Gaussian sigma, from mcmd-pulse-fwhm       */
    double imax;       /* peak intensity, photons/cm^2/s          */
    double photon_energy_eV; /* mcmd-pulse-photon-energy          */

    /* Per-atom electronic state, all of length natoms */
    int natoms;
    int **atom_configurations; /* current K/L/M occupations         */
    int **GS_configurations;   /* ground state, for the net charge  */

    /* Energy of each atom's current electronic state, in eV.  Cached because
     * the charge-transfer pair search would otherwise do one hash lookup per
     * candidate pair, and with the 6 nm cutoff exceeding the size of a
     * typical sample that is essentially every pair, every step. */
    double *energy_of_state;

    /* Per-element atomic data, indexed by the element index (mcmd_atom2idx()) */
    t_mcmd_atomdata *atomdata;

    /* Element index of each atom.  Masses never change during a run, so this
     * is resolved once instead of doing floor(mass + 0.5) and a scan of the
     * element table once per atom per step. */
    int *elem_of_atom;

    /* Element index of hydrogen, resolved once from the element table rather
     * than assumed, so reordering mcmd_elements[] cannot silently break the
     * H-H charge transfer test. */
    int elem_H;

    /* Outermost occupied shell of each atom, or -1 for a fully stripped one.
     * Maintained alongside atom_configurations so the charge-transfer inner
     * loop reads one contiguous array instead of dereferencing a per-atom
     * pointer and walking its three shells for every candidate pair. */
    int *outermost_shell;

    /* Electron-impact ionization and three-body recombination
     * (mcmd-collisions = 1).  The channels out of every state are built
     * once at setup (mcmd_coll_build); their rates depend only on the plasma
     * T_e and n_e, which all atoms share, so they are evaluated once per step
     * for all channels (mcmd_coll_update_rates) and the kernel only reads
     * them.  See the manual, chapter 7. */
    gmx_bool       coll_on;
    gmx_bool       coll_lotz;     /* mcmd-coll-data = lotz: ionization from Lotz's formula, no data file */
    gmx_bool       coll_recomb;          /* include three-body recombination */
    gmx_bool       coll_no_h_recomb;     /* no recombination onto hydrogen   */
    t_mcmd_cchan  *cch[MCMD_ATOM_TABLE_SIZE];       /* channels, by element  */
    double        *cch_g[MCMD_ATOM_TABLE_SIZE][2];  /* g of the tabulated ionization / excitation */
    double         ctab_clamped;                     /* steps with T_e outside a table's grid */
    int           *cch_start[MCMD_ATOM_TABLE_SIZE]; /* [num_rates + 1], by the
                                                     * element's rate state index */
    gmx_bool       coll_active;          /* rates nonzero this step         */
    double         coll_T, coll_ne;      /* plasma state of the current rates */
    double         count_coll_ion, count_coll_rec;
    gmx_bool       coll_excite;      /* mcmd-coll-excitation = 1: collisional (de-)excitation */
    double         count_coll_exc, count_coll_deexc;
    double         coll_ion_shell[MCMD_NSHELL], coll_rec_shell[MCMD_NSHELL]; /* by K, L, M */
    double         coll_energy;          /* net energy taken from the plasma, eV */
    double         coll_pe_energy;       /* cumulative energy that adds to the plasma, eV */
    /* plasma energy budget (diagnostic), cumulative eV: brought in by trapped
     * photo/Auger/field electrons, change from the adiabatic factor (negative
     * = cooling), carried off by evaporation */
    double         eb_dep, eb_adiab, eb_evap;
    FILE          *fcoll;
    gmx_large_int_t coll_step;
    const t_mcmd_cchan **mc_cref;        /* scratch: channel of each MC entry */

    /* Scratch for the Monte Carlo kernel, sized once to the largest number of
     * transitions any loaded state has.  These used to be snew/sfree'd per
     * event, which on a 25000 step run was ~10^8 malloc/free pairs. */
    t_mcmd_transition *mc_transitions;
    double            *mc_dt;
    int                mc_scratch_size;

    /* Visit order for the charge-transfer pass, allocated once. */
    int *ct_order;

    /* Free-electron plasma */
    double num_free_electrons;
    double thermalized_free_electron_energy;
    double escaped_photoelectrons;
    double electron_density;
    double electron_temperature;

    /* Sample geometry */
    double total_mass;
    double radius_of_sample;
    double radius_gyration_equilibrium;
    double sample_volume;
    double current_volume;
    double rg_factor;
    double com[3];        /* centre of mass, nm, updated every step       */
    double radius_now;    /* current equivalent-sphere radius, nm         */
    /* The sample is the uniform ellipsoid of its gyration tensor, semi-axes
     * sqrt(5 lambda_i) along its eigenvectors; radius_now is the
     * volume-equivalent radius (a_1 a_2 a_3)^(1/3). */
    gmx_bool geom_heavy;          /* geometry from non-H atoms only (always) */
    double   ell_axis[3][3];      /* rows: principal directions            */
    double   ell_a[3];            /* semi-axes, nm                          */
    double count_escape_fallback; /* electrons scored without total energies */

    /* The trapped electrons (mcmd-electrons): 0 none, 5 the grid
     * Poisson-Boltzmann cloud, see mcionize_screen_forces(). */
    int     screen_mode;          /* 0 off, 5 pbgrid                          */
    double  screen_E_trapped;     /* kinetic energy of the trapped e-, eV     */
    gmx_bool screen_evap;         /* let trapped e- leave as the sample grows */
    double  screen_lambda_nm;     /* last Debye length, for the log           */
    double  screen_ND;            /* last electrons per Debye sphere          */
    double  screen_evaporated;    /* trapped e- released by expansion         */
    FILE   *fscreen;              /* screening_log.txt                        */

    double  bandwidth_fwhm;       /* relative FWHM of the photon energy        */
    unsigned short erand_state[3];/* own stream: energy draws do not shift the
                                     Monte Carlo's rand() sequence            */

    double   pressure_h_now;      /* electron smoothing width of the last solve, nm (log) */
    double   pressure_work;       /* work of the electron field on the ions, eV */

    /* Self-consistent electrons (mcmd-electrons = pbgrid), see mcmd_pbg_solve() */
    int      pb_nstep;            /* steps between solves                      */
    double  *pb_ne;               /* local electron density w_b rho_b, nm^-3   */
    double   pb_N, pb_h, pb_Dmean;/* N, width and mean depth at the last solve */
    double   pb_W;                /* work on the ions since the last solve, eV */
    double   pb_Eclamp;           /* energy lost to E_kin >= 0, eV             */
    /* The electrons' TOTAL energy
     * E_tot = E_kin + U_es is tracked.  Every event adds exactly what it
     * changes, where it happens: a trapped electron E_emit - D_k, a
     * collisional ionization -I_p - D_k (recombination the reverse), a charge
     * change the ion-electron jump (ledger), evaporation -(kinetic - depth),
     * the electron field's work on the ions -W.  Between solves E_kin
     * (screen_E_trapped) changes only by the kinetic parts; at a solve
     * E_kin = E_tot - U_es of the new cloud, so all relaxation (new electrons
     * falling in, rearrangement, expansion) lands in E_kin.  No mean depth:
     * an ionization and a recombination at the same atom cancel exactly. */
    gmx_bool pb_etot_init;
    double   pb_Etot;
    /* GMX_MCMD_PB_EDIAG=1 (ion-weight, total energy): pb_ediag.txt splits each
     * solve's change of E_kin into the electrostatic steps from the last solve
     * (atoms moving, charges changing, N changing, h changing, rearranging),
     * each minus what was booked for it in between (pbd_*) */
    gmx_bool pb_ediag;
    double   pbd_W, pbd_J, pbd_N, pbd_Ulast;
    real    *pbd_q;               /* charges at the last solve */
    FILE    *fpbd;
    gmx_bool pb_ready;
    gmx_large_int_t pb_step_solve;
    int      pb_nsolve;
    FILE    *fpb, *fpbw;
    gmx_bool pbg_dens_out;        /* GMX_MCMD_PBG_DENSITY=1: write pbg_density.bin (diagnostic) */
    double   tm[10];              /* wall ms per module section on the master, printed at the end (MCMD_TM_*) */
    /* GMX_MCMD_TIMING_DETAIL=1: finer timers (td[], ms) for speed-up estimates */
    gmx_bool tdet;
    double   td[11];
    double   td_rhs2;             /* convolutions spent on a second CG right-hand side (0 since the two-solve CG went, kept for the column) */
    double   td_nsolve, td_nskip, td_newton0, td_newtonex;
    double  *pbg_nprev;           /* density of the solve before the last (extrapolation test) */
    gmx_bool pbg_nprev_ok;
    FILE    *ftd;
    char    *ct_nH;               /* charge transfer, H donors only: neutral hydrogens this pass */
    int     *ct_cellH, ct_cellH_n;/* and their number per cell of the pass's cell list */
    /* and, per cell, the cell_sorted positions of the neutral hydrogens and of
     * the charged atoms (ct_sH/ct_sQ: [ncells+1] offsets into ct_lH/ct_lQ),
     * see mcmd_charge_transfer */
    int     *ct_lH, *ct_lQ, *ct_sH, *ct_sQ, ct_lcap, ct_scap;
    /* energy ledger (energy_ledger.txt): exact energy jumps when charges change,
     * cumulative, eV; see mcmd_ledger_step() */
    float   *led_q;               /* charges at the start of the step          */
    double   led_Jii, led_Jie;    /* ion-ion, ion-electron jumps               */
    double   led_nchg;            /* charge changes counted                    */
    FILE    *fled;
    gmx_bool coll_evlog;          /* GMX_MCMD_COLL_EVLOG=1: coll_event_log.txt, one line per collisional event (pb) */
    FILE    *fcev;
    int      pb_evn;              /* evaporation: sites (pb) or cells (pbgrid)  */
    int     *pb_evlist, pb_evnl;  /* pbgrid: the cells that hold electrons, set at each solve */
    int      pbg_nthreads;        /* FFTW threads for the grid solve            */
    long     park_seq;            /* this rank's count of master-only sections (mcmd_park) */
    double  *pb_evw, *pb_evV;     /* their electrons and potential energy, eV   */
    int      pb_evh_ver;          /* the solve the depth histogram is of */
    double  *pb_evhw, *pb_evhV;   /* depth histogram: electrons and mean potential per bin */

    /* Grid version (mcmd-electrons = pbgrid), see mcmd_pbg_solve() */
    int      pbg_n;               /* cells per side                            */
    int      pbg_nmin, pbg_nmax, pbg_ncap;  /* size per rebuild between nmin and nmax; arrays hold ncap^3 */
    int      pbg_capset;          /* the grid limit is set: 0 = from the sample at the first use, 1 = done */
    int      pbg_nglimit;         /* the most cells per side a sample may get (mcmd-pb-grid-limit) */
    int      pbg_nplan;           /* the size the FFT plans are for            */
    double   pbg_dxt;             /* target spacing at R0, nm: s_e / mcmd-pb-cells-per-smoothing */
    double   pbg_R0;              /* R at the first rebuild: the target grows with R / R0 */
    double   pbg_lo[3], pbg_dx;   /* grid origin and spacing, nm               */
    double   pbg_R, pbg_c[3];     /* R and centre the grid was built for       */
    double   pbg_box;             /* grid half-width / R: the upper limit with a cut, the size without (mcmd-pb-box-max) */
    double   pbg_cut;             /* mcmd-pb-box-cut (0 = fixed box) */
    double   pbg_hbox;            /* the current half-width / R */
    FILE    *fbox;
    double   pbg_se;              /* mcmd-pb-smoothing: smoothing floor s_e, nm: s = max(2 dx, s_e) */
    double   pbg_ntol;            /* mcmd-pb-newton-tol: Newton stops at resid below this */
    int      pbg_outlat;          /* outside ions on sub-lattices (1), checked against the exact sum (2, GMX_MCMD_PBG_OUTLAT_CHECK) */
    FILE    *fol;
    FILE    *fop;
    /* mcmd-pb-long-range: long-range ion-ion Coulomb from the grid (GROMACS does r < rcoulomb) */
    int      pbg_lr, pbg_lr_nstep, lr_nout, lr_headn, lr_nc[3];
    gmx_bool lr_ready;
    gmx_large_int_t lr_step;
    double   lr_rc, lr_W, lr_U, lr_Ush, lr_cs, lr_clo[3];
    double  *pbg_PhiL;            /* master: ion potential on the grid, eV per e */
    float   *pbg_PhiLf;           /* every rank: the same, single precision */
    double  *lr_E;                /* long-range field at the home atoms, eV/nm per e, 3 per atom */
    int     *lr_head, *lr_next;   /* cell list of the charged atoms (cells >= rc) */
    FILE    *flr;
    gmx_bool pbg_have;
    int      pbg_nregrid, pbg_nconv;
    double  *pbg_ne;              /* electrons per cell                        */
    double  *pbg_Phi, *pbg_V, *pbg_g, *pbg_d, *pbg_t, *pbg_Kn, *pbg_nt, *pbg_Vt, *pbg_cg;
    float   *pbg_Ve;              /* electron potential K n, eV (all ranks)     */
    float   *pbg_buf, *pbg_Kh;
    fftwf_complex *pbg_hat;
    fftwf_plan pbg_fwd, pbg_bwd;
    /* ion potential with a margin: ions up to MCMD_PBG_MARG cells outside the box go into one larger FFT,
     * 2 (ng + M) per side, instead of the sub-lattice sums */
    fftwf_plan pbg_fwdE, pbg_bwdE;
    float   *pbg_KhE;
    double   pbg_K0, pbg_ecen[3];
    gmx_large_int_t screen_step;

    /* mcmd-field-ionization: over-the-barrier ionization in the collective
     * field of the PB solution (pbgrid only), see mcmd_pbfi_field */
    gmx_bool pbfi;
    double   pbfi_L;              /* smoothing length, nm (mcmd-field-ionization-smoothing) */
    int      pbfi_kn;             /* grid size, spacing and length the kernel is for */
    double   pbfi_kdx, pbfi_kL;
    float   *pbfi_Kh;             /* Fourier kernel k G(r, L_eff), as pbg_Kh */
    double   pbfi_Gt[13];         /* k G(sqrt(m) dx, L_eff), m = 0..12: own-deposit term */
    double  *pbfi_q, *pbfi_V;     /* source q_ion - n_e, and its potential   */
    double  *pbfi_F;              /* field at every atom, V/nm, last update  */
    gmx_bool pbfi_F_ok;
    double   pbfi_nion, pbfi_Eb;  /* field ionizations, binding energy paid, eV */
    int      pbfi_quiet;          /* updates in a row with no atom near its threshold */
    gmx_bool pbfi_done;           /* stopped for the rest of the run, see mcmd_pbfi_due */
    FILE    *fpbfi_ev;

    /* Current pulse intensity, photons/cm^2/s */
    double intensity;

    /* Largest separation of a pair actually used for a transfer, for the
     * cutoff sanity check at the end of the run. */
    double max_transfer_separation;

    /* Cell list for the charge-transfer pair search, rebuilt every step.
     * Buffers are kept across steps and grown on demand. */
    int    *cell_start;        /* [ncells+1] offsets into cell_sorted   */
    int    *cell_cursor;       /* [ncells] fill positions, build only   */
    int     cell_nalloc;       /* cells the two above are sized for     */
    int    *cell_of_atom;      /* cell each atom was binned into        */
    int    *pos_of_atom;       /* where each atom sits in cell_sorted   */

    /* Atoms gathered into cell order, so the 27-cell scan reads straight
     * through memory instead of chasing a linked list to scattered indices.
     * Holds exactly the fields the candidate rejection tests need. */
    t_mcmd_cell_atom *cell_sorted;

    /* Idle detection for the charge-transfer pass, see
     * MCMD_CT_IDLE_STEPS_DEFAULT.  ct_idle_steps is
     * mcmd-charge-transfer-idle and ct_recheck mcmd-charge-transfer-recheck;
     * ct_idle_steps < 0 means never go idle. */
    gmx_large_int_t ct_idle_steps;
    int             ct_recheck;
    gmx_large_int_t ct_step;
    gmx_large_int_t ct_steps_since_transfer;
    gmx_large_int_t ct_skipped;

    /* Event counters */
    int count_auger;
    int count_fluorescence;
    int count_photoionization;
    int count_chargetransfer;

    /* Output */
    FILE *fcharges;  /* binary, charges_over_time.bin        */
    FILE *fmean;     /* mean_charge_vs_time.txt             */
    FILE *fpulse;    /* pulse_profile.txt                   */
    FILE *felectron; /* electron_data.txt                   */
    FILE *ftrans;    /* electronic_transition_log.txt       */
    char *outbuf[5]; /* stdio buffers, freed with the rest  */
    gmx_large_int_t output_step;
    int             charge_stride;  /* mcmd-charge-output-stride */
    int             text_stride;    /* the per-step text files every nstenergy steps */

    gmx_bool bInitialised; /* set on the first step */
    /* Continuation (-cpi), see mcmd_state_load() */
    t_mcmd_rst *rst;       /* the state file, read at init, applied after   */
    gmx_bool    pbg_cold;  /* the first solve after it starts from the ions  */
    gmx_bool    pbg_rst_grid; /* and on the grid of the checkpoint           */
    int         nstate;    /* state files written by this run                */
    unsigned int seed;     /* the RNG seed (mcmd-seed, or chosen at start)   */
};

/* ---------------------------------------------------------------- *
 * Small helpers
 * ---------------------------------------------------------------- */

/* Seed for the module's RNG: mcmd-seed, or with -1 one decorrelated per job:
 * time(NULL) alone hands identical seeds to every job started within the same
 * second, which silently correlates trajectories across an ensemble, so the
 * pid is mixed in.  A fixed seed makes a run reproducible.  Chosen once in
 * mcionize_init (mc->seed) and written to the log. */
static unsigned int mcmd_rng_seed(const t_inputrec *ir)
{
    if (ir->mcmd_seed >= 0)
    {
        return (unsigned int)ir->mcmd_seed;
    }
    return (unsigned int)(time(NULL) ^ (getpid() << 16));
}

static void mcionize_free(t_mcionize *mc);

/* Open a file under MCMD_OUTDIR, or die trying. */
/* An mdp real (single precision in this build) as the decimal written in the
 * mdp: 0.05 comes back as 0.05, not 0.0500000007. */
static double mcmd_mdp_real(real r)
{
    char buf[32];

    sprintf(buf, "%.7g", (double)r);
    return atof(buf);
}

static FILE *mcmd_open_output(const char *name, const char *mode)
{
    char  path[512];
    FILE *fp;

    sprintf(path, "%s/%s", MCMD_OUTDIR, name);
    fp = fopen(path, mode);
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open '%s' for mode '%s'.", path, mode);
    }

    return fp;
}

/* Open one of the per-step outputs and give it a large stdio buffer, so a
 * step costs an fprintf into memory rather than open/write/close. */
static FILE *mcmd_open_buffered(const char *name, const char *mode, char **buf)
{
    FILE *fp = mcmd_open_output(name, mode);

    snew(*buf, MCMD_OUTPUT_BUFSIZE);
    setvbuf(fp, *buf, _IOFBF, MCMD_OUTPUT_BUFSIZE);

    return fp;
}

/* transitions.bin (replaced electronic_transition_log.txt on 2026-10-27, ~7x
 * smaller): a 16-byte header ("MCMDTRN1", int32 version, int32 record size),
 * then one 20-byte little-endian record per event:
 *   float32 t [fs], int32 atom (0-based), int32 partner (the donor of a
 *   charge transfer, else -1), uint8 type, uint8 K L M before, uint8 K L M
 *   after, uint8 0.
 * type: 0 Auger, 1 fluorescence, 2 photoionization, 3 collisional
 * ionization, 4 recombination, 5 excitation, 6 de-excitation (the rate
 * table's types), 7 charge transfer (atom = acceptor).  Read it with
 * tools/mcmd_output.py. */
#define MCMD_TRN_CT      7
#define MCMD_TRN_VERSION 1
#define MCMD_TRN_RECSIZE 20

static void mcmd_trn_header(FILE *fp)
{
    const int hdr[2] = {MCMD_TRN_VERSION, MCMD_TRN_RECSIZE};

    fwrite("MCMDTRN1", 1, 8, fp);
    fwrite(hdr, sizeof(int), 2, fp);
}

static void mcmd_trn_write(FILE *fp, double t_ps, int atom, int partner, int type, const int *before,
                           const int *after)
{
    unsigned char rec[MCMD_TRN_RECSIZE];
    const float   t = (float)(t_ps * 1e3);
    int           d;

    memcpy(rec, &t, 4);
    memcpy(rec + 4, &atom, 4);
    memcpy(rec + 8, &partner, 4);
    rec[12] = (unsigned char)type;
    for (d = 0; d < 3; d++)
    {
        rec[13 + d] = (unsigned char)before[d];
        rec[16 + d] = (unsigned char)after[d];
    }
    rec[19] = 0;
    fwrite(rec, 1, MCMD_TRN_RECSIZE, fp);
}

/* Resolve an atom's element from its mass, turning the -1 that
 * mcmd_atom2idx() returns for an unlisted mass into a clear error instead of
 * an out-of-bounds read.  Checks that the element is known
 * and its data loaded.  Only used to fill elem_of_atom[] at setup; the hot
 * paths read that array instead. */
static int mcmd_resolve_element(const t_mcionize *mc, const t_mdatoms *mdatoms,
                                int atom)
{
    int idx = mcmd_atom2idx(mdatoms->massT[atom]);

    if (idx < 0)
    {
        gmx_fatal(FARGS,
                  "Atom %d has mass %g, which is not in the element table.  "
                  "Add it to mcmd_elements[] in atomic_data.c, together with "
                  "its ground-state configuration and its Atomic_data files.",
                  atom + 1, mdatoms->massT[atom]);
    }
    if (!mc->atomdata[idx].bPresent)
    {
        gmx_fatal(FARGS, "No atomic data loaded for element index %d.", idx);
    }

    return idx;
}

static int mcmd_atom_element(const t_mcionize *mc, const t_mdatoms *mdatoms,
                             int atom)
{
    return mc->elem_of_atom[atom];
}

/* Net charge implied by an atom's configuration relative to its ground
 * state: positive means electrons have been lost. */
static double mcmd_net_charge(const t_mcionize *mc, int atom)
{
    return (double)((mc->GS_configurations[atom][0] - mc->atom_configurations[atom][0]) +
                    (mc->GS_configurations[atom][1] - mc->atom_configurations[atom][1]) +
                    (mc->GS_configurations[atom][2] - mc->atom_configurations[atom][2]));
}

/* Refresh the cached state energy of one atom.  Fatal if the atom's current
 * configuration has no entry in its energy level file, since the
 * over-the-barrier criterion cannot be evaluated without it. */
static void mcmd_update_state_energy(t_mcionize *mc, const t_mdatoms *mdatoms,
                                     int atom)
{
    int    elem = mcmd_atom_element(mc, mdatoms, atom);
    double E    = mcmd_dict_get(mc->atomdata[elem].energy_levels,
                                mc->atom_configurations[atom]);

    if ((int)E == -1)
    {
        gmx_fatal(FARGS,
                  "No energy level for atom %d (mass %g) in state %d %d %d.",
                  atom + 1, mdatoms->massT[atom],
                  mc->atom_configurations[atom][0],
                  mc->atom_configurations[atom][1],
                  mc->atom_configurations[atom][2]);
    }
    mc->energy_of_state[atom] = E;
}

static int mcmd_outermost_occupied(const int config[MCMD_NSHELL]);

/* Bring both per-atom caches back in step with atom_configurations[atom].
 * Call this everywhere a configuration changes; between such changes the
 * caches stay valid, which is what lets the charge-transfer pass drop the
 * loop that used to re-derive the energy of every atom on every step.
 *
 * The energy cache is only maintained when charge transfer is on: it is the
 * only consumer, and computing it unconditionally would turn a missing energy
 * level into a fatal error in runs that previously completed without one. */
static void mcmd_refresh_atom(t_mcionize *mc, const t_mdatoms *mdatoms,
                              int atom)
{
    mc->outermost_shell[atom] =
        mcmd_outermost_occupied(mc->atom_configurations[atom]);

    if (mc->bChargeTransfer)
    {
        mcmd_update_state_energy(mc, mdatoms, atom);
    }
}

/* Fisher-Yates shuffle.  Draws size-1 values from rand(); deliberately does
 * not reseed, since srand() on every call resets the stream and made repeated
 * calls within the same second return the same permutation. */
static void mcmd_shuffle(int arr[], int size)
{
    int i, j, tmp;

    for (i = size - 1; i > 0; i--)
    {
        j = rand() % (i + 1);

        tmp    = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

/* Index of the outermost occupied shell, or -1 if the atom has no electrons
 * left. */
static int mcmd_outermost_occupied(const int config[MCMD_NSHELL])
{
    int shell;

    for (shell = MCMD_NSHELL - 1; shell >= 0; shell--)
    {
        if (config[shell] > 0)
        {
            return shell;
        }
    }

    return -1;
}

/* ---------------------------------------------------------------- *
 * Charge transfer
 *
 * Classical over-the-barrier model.  For a pair of atoms with
 * different charge the less charged one is the donor, and the
 * transfer is allowed when
 *
 *     R_crit = (Q_D + 1 + 2 sqrt((Q_D + 1) Q_A)) / E_donor
 *
 * exceeds their separation.  Atoms are visited in random order, and
 * each takes the nearest partner that qualifies.
 *
 * Every atom is considered as a partner.  Until 2026-08-20 the partner
 * search ran over hydrogens only, which - combined with the guard that
 * drops fully stripped atoms - meant a neutral hydrogen was always the
 * donor and charge could never move directly between two heavy atoms.
 * ---------------------------------------------------------------- */

/* Bin the atoms into a uniform grid of cells at least MCMD_CT_CUTOFF_NM
 * across, so the pair search only has to look at an atom's own cell and its
 * 26 neighbours instead of every other atom.  Because the cell edge is never
 * smaller than the cutoff, that 3x3x3 block provably contains every atom
 * within the cutoff, so this changes which pairs are *examined* and nothing
 * about which are accepted.
 *
 * The sample expands by an order of magnitude during a run, so the cell count
 * is capped and the cell edge grown to match; that keeps the memory O(N) and
 * costs only a few extra candidates per cell once the system is dilute.
 *
 * Returns the grid origin, dimensions and edge length through its pointers. */
static void mcmd_build_cell_list(t_mcionize *mc, const t_mdatoms *mdatoms,
                                 rvec x[], double origin[3], int ncell[3],
                                 double *edge)
{
    const int max_cells = 20 * mdatoms->nr + 1000;
    double    lo[3], hi[3];
    double    cell = MCMD_CT_CUTOFF_NM;
    double    inv_cell;
    int       i, d, total;

    for (d = 0; d < 3; d++)
    {
        lo[d] = hi[d] = x[0][d];
    }
    for (i = 1; i < mdatoms->nr; i++)
    {
        for (d = 0; d < 3; d++)
        {
            if (x[i][d] < lo[d]) { lo[d] = x[i][d]; }
            if (x[i][d] > hi[d]) { hi[d] = x[i][d]; }
        }
    }

    for (;;)
    {
        double fcells = 1.0;

        for (d = 0; d < 3; d++)
        {
            ncell[d] = (int)((hi[d] - lo[d]) / cell) + 1;
            if (ncell[d] < 1)
            {
                ncell[d] = 1;
            }
            fcells *= (double)ncell[d];
        }
        /* Compare in floating point: the product can overflow int for a
         * badly exploded sample before the cap would catch it. */
        if (fcells <= (double)max_cells)
        {
            break;
        }
        cell *= 2.0;
    }

    total = ncell[0] * ncell[1] * ncell[2];

    if (total + 1 > mc->cell_nalloc)
    {
        srenew(mc->cell_start,  total + 1);
        srenew(mc->cell_cursor, total + 1);
        mc->cell_nalloc = total + 1;
    }

    inv_cell = 1.0 / cell;

    /* Counting sort into cell order.  Pass one bins the atoms and counts
     * them; cell_start is built one slot high so the prefix sum below turns
     * the counts into offsets in place. */
    for (i = 0; i <= total; i++)
    {
        mc->cell_start[i] = 0;
    }

    for (i = 0; i < mdatoms->nr; i++)
    {
        int c[3], idx;

        for (d = 0; d < 3; d++)
        {
            c[d] = (int)((x[i][d] - lo[d]) * inv_cell);
            if (c[d] < 0)          { c[d] = 0; }
            if (c[d] >= ncell[d])  { c[d] = ncell[d] - 1; }
        }
        idx = (c[2] * ncell[1] + c[1]) * ncell[0] + c[0];

        /* Kept so the pair search can recover each atom's cell without
         * repeating the division. */
        mc->cell_of_atom[i] = idx;
        mc->cell_start[idx + 1]++;
    }

    for (i = 0; i < total; i++)
    {
        mc->cell_start[i + 1] += mc->cell_start[i];
        mc->cell_cursor[i]     = mc->cell_start[i];
    }

    /* Pass two places the atoms.  Descending i reproduces the order the old
     * head-insertion linked list gave - within a cell, highest index first -
     * so the candidate scan visits pairs in exactly the same sequence and the
     * results stay bit for bit identical. */
    for (i = mdatoms->nr - 1; i >= 0; i--)
    {
        int p = mc->cell_cursor[mc->cell_of_atom[i]]++;

        mc->pos_of_atom[i]      = p;
        mc->cell_sorted[p].x     = x[i][0];
        mc->cell_sorted[p].y     = x[i][1];
        mc->cell_sorted[p].z     = x[i][2];
        mc->cell_sorted[p].q     = mdatoms->chargeA[i];
        mc->cell_sorted[p].shell = mc->outermost_shell[i];
        mc->cell_sorted[p].idx   = i;
    }

    for (d = 0; d < 3; d++)
    {
        origin[d] = lo[d];
    }
    *edge = cell;
}

/* Push an atom's charge and shell back into its gathered copy.  Positions do
 * not change during a pass, so they never need this. */
static void mcmd_cell_sync(t_mcionize *mc, const t_mdatoms *mdatoms, int atom)
{
    t_mcmd_cell_atom *ca = &mc->cell_sorted[mc->pos_of_atom[atom]];

    ca->q     = mdatoms->chargeA[atom];
    ca->shell = mc->outermost_shell[atom];
}

/* Pick the acceptor orbital the transferred electron enters.
 *
 * Boll et al., Nat. Phys. 18, 423 (2022) put the electron not into the
 * acceptor's outermost occupied shell, nor even into its outermost vacancy,
 * but into whichever vacancy has the binding energy closest to that of the
 * donor orbital the electron came from.
 *
 * The binding energy an electron would have in shell s is a difference of
 * total energies, T(c + e_s) - T(c), which is defined for every shell and
 * every configuration - including the autoionizing ones energy_levels_X.txt
 * can only mark with a placeholder.  That is why total_energies_X.txt is
 * required whenever charge transfer is on: priced from energy_levels_X.txt
 * instead, an acceptor holding an inner-shell hole cannot be scored at all
 * and is refused, which at hard photon energies is most of them.
 *
 * A deeper vacancy - a K hole under an occupied L shell - is skipped either
 * way.  Without total energies it would read back the L binding energy, two
 * orders of magnitude too small, making every K hole look like an ideal match
 * for any valence donor.  With them it could be priced correctly but stays
 * excluded on physics: it is the exclusion the paper applies by hand to the M
 * and N shells of iodine, because refilling the shell that absorbed the
 * photon short-circuits the Auger cascade the rate data already models.
 *
 * On top of the match, mcmd-charge-transfer-downhill optionally requires the
 * transfer to be downhill.  It is off by default: the paper imposes no such
 * condition, and on the lysozyme example turning it on changes the transfer
 * count by an order of magnitude without moving any observable.  Charge
 * transfer there saturates, so the extra firings only reshuffle a
 * distribution the energetics have already fixed.  The switch is kept because
 * that is a property of the regime rather than a theorem.
 *
 * "Downhill" is not E_acceptor > E_donor.  Those are isolated-atom binding
 * energies, and the electron is moving between two ions a fraction of a
 * nanometre apart, so the change in their mutual Coulomb energy is part of
 * the balance: the pair goes from Q_D Q_A to (Q_D + 1)(Q_A - 1), releasing
 * (Q_A - Q_D - 1) e^2 / R.  Comparing the bare binding energies instead makes
 * transfer to a highly charged acceptor look endothermic whenever its outer
 * shell is shallower than the donor's, which is the opposite of what the
 * over-the-barrier picture describes.
 *
 * The role assignment upstream guarantees Q_A > Q_D, so the term is never
 * negative: it can only permit a transfer the bare comparison refused.  It
 * vanishes exactly when Q_A = Q_D + 1, the resonant case, where the acceptor
 * ends up in precisely the state the donor left - so resonant pairs cancel to
 * E_acceptor = E_donor and are excluded whenever the option is on.
 *
 * Returns the chosen shell, or -1 if no vacancy is both scorable and
 * downhill. */
static int mcmd_acceptor_shell(t_mcionize *mc, const t_mdatoms *mdatoms,
                               int acceptor, double E_donor, double dE_coul)
{
    int          elem = mcmd_atom_element(mc, mdatoms, acceptor);
    int         *cfg  = mc->atom_configurations[acceptor];
    t_mcmd_dict *total = mc->atomdata[elem].total_energies;
    int          best = -1;
    double       best_gap = 0;
    double       T_now = 0;
    int          s;

    T_now = mcmd_dict_get(total, cfg);

    if ((int)T_now == -1)
    {
        /* The two tables are generated together from one model, so a state
         * present in energy_levels_X.txt and absent here means the pair does
         * not match.  Silently mispricing the acceptor would be worse. */
        gmx_fatal(FARGS,
                  "No total energy for atom %d (mass %g) in state %d %d %d.  "
                  "total_energies_X.txt and energy_levels_X.txt must be "
                  "generated from the same atomic model.",
                  acceptor + 1, mdatoms->massT[acceptor],
                  cfg[0], cfg[1], cfg[2]);
    }

    for (s = 0; s < MCMD_NSHELL; s++)
    {
        int    trial[MCMD_NSHELL];
        int    k;
        double E_acceptor, gap;

        if (cfg[s] + 1 > mc->GS_configurations[acceptor][s])
        {
            continue;   /* shell is already full: no vacancy here */
        }

        for (k = 0; k < MCMD_NSHELL; k++)
        {
            trial[k] = cfg[k];
        }
        trial[s] += 1;

        if (mcmd_outermost_occupied(trial) != s)
        {
            continue;   /* deeper vacancy; see above */
        }

        {
            double T_new = mcmd_dict_get(total, trial);

            if ((int)T_new == -1)
            {
                continue;
            }

            /* What the electron gains by being bound here, whatever shell
             * this is and whether or not the result autoionizes. */
            E_acceptor = T_new - T_now;
        }


        /* Downhill only, counting the Coulomb release; see the header.  Off
         * by default: the reference method does not impose it, and measured
         * on a 30 fs lysozyme run it changes the transfer count twelvefold
         * while leaving every observable alone.  The margin keeps exact
         * resonances out; see MCMD_DOWNHILL_TOL_EV. */
        if (mc->bDownhill &&
            E_acceptor + dE_coul <= E_donor + MCMD_DOWNHILL_TOL_EV)
        {
            continue;
        }

        /* The same shift applies to the effective depth of every shell on
         * this acceptor, so it belongs in the matched quantity too. */
        gap = fabs(E_acceptor + dE_coul - E_donor);

        if (best < 0 || gap < best_gap)
        {
            best     = s;
            best_gap = gap;
        }
    }

    return best;
}

static double mcmd_pb_wallms(void);
/* GMX_MCMD_TIMING_DETAIL: mc->td[] */
#define MCMD_TD_CONV    0  /* grid convolutions, total                             */
#define MCMD_TD_FFT     1  /* of which the FFTs                                    */
#define MCMD_TD_CG      2  /* grid CG, total (convolutions included)               */
#define MCMD_TD_CGCONV  3  /* convolutions inside the CG                           */
#define MCMD_TD_NEWTON  4  /* grid Newton loop, total                              */
#define MCMD_TD_CTBUILD 5  /* charge transfer: cell-list build                     */
#define MCMD_TD_CTSHUF  6  /* charge transfer: shuffle                             */
#define MCMD_TD_CTLOOP  7  /* charge transfer: pair loop                           */
#define MCMD_TD_EVAP    8  /* evaporation                                          */
#define MCMD_TD_EXTRAP  9  /* the extrapolation test itself (one convolution)      */
#define MCMD_TD_OUTPOT 10  /* potential of the ions outside the grid box           */

static void mcmd_charge_transfer(t_mcionize *mc, t_mdatoms *mdatoms,
                                 t_state *state, double t)
{
    int    *idx_map = mc->ct_order;
    int     i, j, j2;
    int     number_of_charge_transfers = 0;
    int     ct_before[3];               /* the acceptor's K L M before a transfer (transitions.bin) */
    double  origin[3], edge;
    int     ncell[3];

    /* energy_of_state[] and outermost_shell[] are kept current by
     * mcmd_refresh_atom() wherever a configuration changes, so the loop that
     * used to re-derive the energy of every atom here - one hash lookup per
     * atom per step, ~15 million of them on a 25000 step run - is gone. */

    /* With hydrogen donors only (mcmd-ct-h-donor-only), a pair needs a neutral
     * hydrogen: flag them once.  Hydrogen ionizes early and is not
     * recombined (mcmd-coll-h-recombination = 0), so for most of a run there are
     * none, and then the pass is skipped: the shuffle's nr-1 draws are made
     * all the same, so the random stream is that of a full pass. */
    if (mc->ct_h_donor_only)
    {
        int nH = 0;

        if (mc->ct_nH == NULL)
        {
            snew(mc->ct_nH, mdatoms->nr);
        }
        for (j = 0; j < mdatoms->nr; j++)
        {
            mc->ct_nH[j] = (mc->elem_of_atom[j] == mc->elem_H && mc->outermost_shell[j] >= 0);
            nH          += mc->ct_nH[j];
        }
        if (nH == 0)
        {
            for (j = mdatoms->nr - 1; j > 0; j--)
            {
                (void)rand();
            }
            return;
        }
    }

    /* Positions do not change during this pass, so one build serves it all. */
    double tct0 = mc->tdet ? mcmd_pb_wallms() : 0;

    mcmd_build_cell_list(mc, mdatoms, state->x, origin, ncell, &edge);
    if (mc->tdet) { mc->td[MCMD_TD_CTBUILD] += mcmd_pb_wallms() - tct0; tct0 = mcmd_pb_wallms(); }

    /* H donors only: neutral hydrogens per cell, so that an atom with none in
     * its 27 cells (it can only pair with one) is skipped whole, and empty
     * cells are skipped in the search - the same pairs as the full scan */
    if (mc->ct_h_donor_only)
    {
        const int nc = ncell[0] * ncell[1] * ncell[2];

        if (nc > mc->ct_cellH_n)
        {
            srenew(mc->ct_cellH, nc);
            mc->ct_cellH_n = nc;
        }
        memset(mc->ct_cellH, 0, sizeof(int) * nc);
        for (j = 0; j < mdatoms->nr; j++)
        {
            if (mc->ct_nH[j]) { mc->ct_cellH[mc->cell_of_atom[j]]++; }
        }
    }

    /* H donors only, no H-H transfer: a pair is a neutral hydrogen (donor)
     * and an atom of higher charge (acceptor).  Within a pass no atom becomes
     * a neutral hydrogen (hydrogens only donate) and no atom becomes charged
     * that could pair (only hydrogens gain charge, and H-H is off), so the
     * neutral hydrogens and the charged atoms at the start of the pass hold
     * every partner the pass can find.  A neutral hydrogen then scans only
     * the charged atoms of its cells, a charged atom only the neutral
     * hydrogens, and an uncharged heavy atom nothing.  The sub-lists keep the
     * cell_sorted order, so the nearest partner, ties included, is the one
     * the full scan picks: the same pairs, bit for bit. */
    const gmx_bool fast = mc->ct_h_donor_only && !mc->bAllowHCT;

    if (fast)
    {
        const int nc = ncell[0] * ncell[1] * ncell[2];
        int       c, n, kH = 0, kQ = 0;

        if (nc + 1 > mc->ct_scap)
        {
            srenew(mc->ct_sH, nc + 1);
            srenew(mc->ct_sQ, nc + 1);
            mc->ct_scap = nc + 1;
        }
        if (mdatoms->nr > mc->ct_lcap)
        {
            srenew(mc->ct_lH, mdatoms->nr);
            srenew(mc->ct_lQ, mdatoms->nr);
            mc->ct_lcap = mdatoms->nr;
        }
        for (c = 0; c < nc; c++)
        {
            mc->ct_sH[c] = kH;
            mc->ct_sQ[c] = kQ;
            for (n = mc->cell_start[c]; n < mc->cell_start[c + 1]; n++)
            {
                if (mc->ct_nH[mc->cell_sorted[n].idx]) { mc->ct_lH[kH++] = n; }
                if (mc->cell_sorted[n].q > 0)          { mc->ct_lQ[kQ++] = n; }
            }
        }
        mc->ct_sH[nc] = kH;
        mc->ct_sQ[nc] = kQ;
    }

    for (j = 0; j < mdatoms->nr; j++)
    {
        idx_map[j] = j;
    }
    mcmd_shuffle(idx_map, mdatoms->nr);
    if (mc->tdet) { mc->td[MCMD_TD_CTSHUF] += mcmd_pb_wallms() - tct0; tct0 = mcmd_pb_wallms(); }

    for (j2 = 0; j2 < mdatoms->nr; j2++)
    {
        double R_min = 1e7;
        int    R_min_idx = -1; /* acceptor of the closest qualifying pair */
        int    donor_idx  = -1; /* and its donor - the two must stay paired */
        int    R_min_shell = -1; /* and the orbital settled on for it */
        int    INDEX1 = 100, INDEX2 = 100;
        int    shell;

        int  jc[3], dx, dy, dz, cj;
        real xj, yj, zj, qj;
        gmx_bool bJisH;

        j = idx_map[j2];

        /* Fixed for this atom's whole search, so it is resolved once here
         * rather than per candidate. */
        bJisH = (mc->elem_of_atom[j] == mc->elem_H);

        /* Nothing is applied until this atom's search finishes, so its own
         * position and charge are fixed for the duration and can be hoisted
         * out of the candidate loop. */
        xj = state->x[j][0];
        yj = state->x[j][1];
        zj = state->x[j][2];
        qj = mdatoms->chargeA[j];

        /* The cell list build already binned this atom; unpack that rather
         * than repeating the three divisions. */
        cj    = mc->cell_of_atom[j];
        jc[0] = cj % ncell[0];
        jc[1] = (cj / ncell[0]) % ncell[1];
        jc[2] = cj / (ncell[0] * ncell[1]);

        if (fast && !mc->ct_nH[j] && (qj <= 0 || bJisH))
        {
            continue;    /* no partner: see above (an ionized hydrogen only pairs with hydrogens) */
        }
        if (mc->ct_h_donor_only && !mc->ct_nH[j])
        {
            int nH = 0;

            for (dz = -1; dz <= 1; dz++)
            for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++)
            {
                int cx = jc[0] + dx, cy = jc[1] + dy, cz = jc[2] + dz;

                if (cx >= 0 && cx < ncell[0] && cy >= 0 && cy < ncell[1] && cz >= 0 && cz < ncell[2])
                {
                    nH += mc->ct_cellH[(cz * ncell[1] + cy) * ncell[0] + cx];
                }
            }
            if (nH == 0)
            {
                continue;
            }
        }

        /* Own cell plus the 26 around it: the cell edge is at least the
         * cutoff, so nothing within range can lie outside that block. */
        for (dz = -1; dz <= 1; dz++)
        for (dy = -1; dy <= 1; dy++)
        for (dx = -1; dx <= 1; dx++)
        {
            int cx = jc[0] + dx, cy = jc[1] + dy, cz = jc[2] + dz;

            if (cx < 0 || cx >= ncell[0] ||
                cy < 0 || cy >= ncell[1] ||
                cz < 0 || cz >= ncell[2])
            {
                continue;
            }

            if (mc->ct_h_donor_only && !mc->ct_nH[j] && mc->ct_cellH[(cz * ncell[1] + cy) * ncell[0] + cx] == 0)
            {
                continue;
            }
        {
        const int  c       = (cz * ncell[1] + cy) * ncell[0] + cx;
        const int *lst     = !fast ? NULL : (mc->ct_nH[j] ? mc->ct_lQ : mc->ct_lH);
        const int *st      = !fast ? mc->cell_start : (mc->ct_nH[j] ? mc->ct_sQ : mc->ct_sH);
        const int  cell_lo = st[c];
        const int  cell_hi = st[c + 1];
        int        kk;

        for (kk = cell_lo; kk < cell_hi; kk++)
        {
            const int               n  = fast ? lst[kk] : kk;
            const t_mcmd_cell_atom *ca = &mc->cell_sorted[n];
            double Q_D, Q_A, R_cob, R_crit, E_donor, dr2;
            int    cand_donor, cand_acceptor, cand_shell;

            i = ca->idx;

            if (i == j)
            {
                continue;
            }

            /* H donors only: no neutral hydrogen in the pair, no donor (the
             * same pairs the donor test below rejects, before the distance) */
            if (mc->ct_h_donor_only && !mc->ct_nH[i] && !mc->ct_nH[j])
            {
                continue;
            }

            /* A fully stripped atom has nothing to give, but it can still
             * receive: a bare nucleus is all vacancy.  So it is only useless
             * as a candidate when the pairing below would cast it as the
             * donor, which is exactly when it is the less charged of the two.
             * The gathered copy puts shell, charge and position on the same
             * cache line, so a rejected candidate costs one line rather than
             * four scattered loads. */
            if (ca->shell < 0 && ca->q < qj)
            {
                continue;
            }

            /* Equal charges: nothing to transfer.  Charges are integral, so
             * the integer abs() here is exact. */
            if (abs((int)(qj - ca->q)) < 1e-5)
            {
                continue;
            }

            /* Reject on the squared distance so the square root is only
             * taken for the few pairs that survive the cutoff. */
            dr2 = sqr(ca->x - xj) +
                  sqr(ca->y - yj) +
                  sqr(ca->z - zj);

            if (dr2 > MCMD_CT_CUTOFF_NM * MCMD_CT_CUTOFF_NM)
            {
                continue;
            }

            R_cob = sqrt(dr2);

            /* Hydrogen-hydrogen transfer between two ground-state
             * hydrogens moves an electron between wells of identical depth,
             * so it is resonant: the pair swaps roles and the charge
             * distribution is unchanged.  Since transfer is applied whenever
             * the over-the-barrier geometry allows rather than at a physical
             * rate, such a pair fires every step, and it dominates the
             * transfer count for no effect.
             *
             * Off by default on that cost argument.  It is not currently a
             * physics argument: on the lysozyme example enabling it moves no
             * observable, because hydrogen ionizes fully within about a
             * femtosecond and there are then no neutral hydrogens left to
             * shuffle.  Set mcmd-allow-h-ct = 1 to enable. */
            if (!mc->bAllowHCT && bJisH &&
                mc->elem_of_atom[i] == mc->elem_H)
            {
                continue;
            }

            /* The more charged atom is the acceptor.  Which of the pair
             * donates therefore depends on their charges, so the "has an
             * electron to give" test has to be applied after this, not just
             * to i. */
            if (ca->q > qj)
            {
                Q_D            = qj;
                Q_A            = ca->q;
                cand_donor     = j;
                cand_acceptor  = i;
            }
            else
            {
                Q_D            = ca->q;
                Q_A            = qj;
                cand_donor     = i;
                cand_acceptor  = j;
            }

            /* A fully stripped donor has nothing to give, and its state
             * energy is 0, which would make R_crit infinite and let the pair
             * win the nearest-partner contest at any distance only to be
             * rejected later - blocking a transfer another partner could
             * have supplied.
             *
             * A stripped acceptor used to be rejected here too, on the
             * grounds that it had "nowhere to put the electron".  That was
             * wrong: an atom with no occupied shell has a vacancy in every
             * shell, and is the easiest acceptor to place an electron in, not
             * the hardest.  Hydrogen has only a K shell, so the rule made
             * every ionized hydrogen a permanent electron sink - it could
             * donate once and never recombine. */
            if (mc->outermost_shell[cand_donor] < 0)
            {
                continue;
            }

            /* The published version only ever had a neutral hydrogen
             * as the donor (the partner search ran over hydrogens only), so
             * charge never moved between heavy atoms and an ionized hydrogen
             * never took an electron back.  mcmd-ct-h-donor-only = 1 restores
             * that: hydrogen, which barely photoionizes at these energies,
             * still loses its electron to charged neighbours, but bound
             * electrons are not shuffled between heavy atoms. */
            if (mc->ct_h_donor_only && mc->elem_of_atom[cand_donor] != mc->elem_H)
            {
                continue;
            }

            E_donor = mc->energy_of_state[cand_donor];

            R_crit = MCMD_COULOMB_EV_NM *
                     (Q_D + 1 + 2 * sqrt((Q_D + 1) * Q_A)) / E_donor;

            if (R_crit < R_cob)
            {
                continue;
            }

            /* Which orbital the electron would enter - and whether any is
             * both scorable and downhill - is a property of the pair, so it
             * has to be settled here rather than after the nearest partner
             * has won.  Deciding it afterwards let a nearest neighbour that
             * fails the energy test block a farther one that would have
             * passed, and abandoned the atom's transfer for that step
             * entirely. */
            cand_shell = mcmd_acceptor_shell(mc, mdatoms, cand_acceptor,
                                             E_donor,
                                             MCMD_COULOMB_EV_NM *
                                             (Q_A - Q_D - 1) / R_cob);

            if (cand_shell < 0)
            {
                continue;
            }

            if (R_cob < R_min)
            {
                R_min       = R_cob;
                R_min_idx   = cand_acceptor;
                donor_idx   = cand_donor;
                R_min_shell = cand_shell;
            }
        }
        }
        }   /* neighbour cell */

        if (R_min_idx == -1)
        {
            continue;
        }

        /* The electron leaves the donor's outermost occupied shell. */
        shell  = mc->outermost_shell[donor_idx];
        INDEX1 = (shell < 0) ? 100 : shell;

        if (INDEX1 == 100 || mc->atom_configurations[donor_idx][INDEX1] - 1 < 0)
        {
            continue;
        }

        /* R_min_idx is by construction the more charged of the pair and
         * donor_idx the less charged, so the direction is already fixed; the
         * acceptor orbital was settled during the search, since a pair that
         * has none is not a qualifying pair at all.
         *
         * There used to be a fallback here that reverted to the original
         * targeting - outermost occupied shell, if it has room - so that the
         * energy match could never refuse a transfer the old rule allowed.
         * That fallback is what made the match inert: it carried 43.9% of
         * transfers, and every one of them skipped the energy test entirely.
         * Refusing is the correct outcome, and it is a decision rather than a
         * gap: every acceptor target an atom can reach is present in the
         * tables, so a refusal means the electron would either be unbound on
         * the acceptor or less bound than it already is. */
        INDEX2 = R_min_shell;

        {
            int d;

            for (d = 0; d < 3; d++) { ct_before[d] = mc->atom_configurations[R_min_idx][d]; }
        }

        /* Move one electron from the donor's outermost occupied shell to the
         * acceptor's, and apply it immediately.  The previous version
         * accumulated the move into a per-step delta array that was applied
         * in full on every pair but never reset, so an atom taking part in
         * more than one transfer in the same step had its delta applied
         * repeatedly and net charge was not conserved.  That could not happen
         * while partners were restricted to hydrogens, because a hydrogen
         * became fully stripped after donating once and dropped out of the
         * search. */
        mc->atom_configurations[donor_idx][INDEX1]  -= 1;
        mc->atom_configurations[R_min_idx][INDEX2] += 1;

        mdatoms->chargeA[donor_idx] = mcmd_net_charge(mc, donor_idx);
        mdatoms->chargeA[R_min_idx] = mcmd_net_charge(mc, R_min_idx);

        /* Both configurations just changed, so their caches must be refreshed
         * before the next pair is considered. */
        mcmd_refresh_atom(mc, mdatoms, donor_idx);
        mcmd_refresh_atom(mc, mdatoms, R_min_idx);

        /* The gathered copies are what later atoms in this pass will read, so
         * they have to follow the change too.  Missing this would let the
         * rest of the pass see pre-transfer charges. */
        mcmd_cell_sync(mc, mdatoms, donor_idx);
        mcmd_cell_sync(mc, mdatoms, R_min_idx);

        if (R_min > mc->max_transfer_separation)
        {
            mc->max_transfer_separation = R_min;
        }

        number_of_charge_transfers += 1;
        mc->count_chargetransfer   += 1;

        if (mc->bDetailedOutput)
        {
            mcmd_trn_write(mc->ftrans, t, R_min_idx, donor_idx, MCMD_TRN_CT, ct_before,
                           mc->atom_configurations[R_min_idx]);
        }
    }

}

/* ---------------------------------------------------------------- *
 * Kinetic Monte Carlo over electronic transitions
 *
 * First-reaction method: for every transition available from the
 * atom's current state draw a waiting time dt = -ln(u)/rate, take the
 * smallest, and fire it if it lands inside the MD timestep.  Repeat
 * until the accumulated time passes the MD timestep.
 * ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- *
 * Sample geometry and the fate of emitted electrons
 * ---------------------------------------------------------------- */

/* Eigen-decomposition of a symmetric 3x3 matrix by cyclic Jacobi rotations.
 * A is destroyed; lam gets the eigenvalues and the rows of V the matching
 * unit eigenvectors. */
static void mcmd_jacobi3(double A[3][3], double lam[3], double V[3][3])
{
    int i, j, k, sweep;

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++) { V[i][j] = (i == j); }
    }
    for (sweep = 0; sweep < 50; sweep++)
    {
        double off = sqr(A[0][1]) + sqr(A[0][2]) + sqr(A[1][2]);

        if (off < 1e-30 * (sqr(A[0][0]) + sqr(A[1][1]) + sqr(A[2][2]) + 1e-300))
        {
            break;
        }
        for (i = 0; i < 2; i++)
        {
            for (j = i + 1; j < 3; j++)
            {
                double th, t, c, sn;

                if (A[i][j] == 0) { continue; }
                th = (A[j][j] - A[i][i]) / (2.0 * A[i][j]);
                t  = (th >= 0 ? 1.0 : -1.0) / (fabs(th) + sqrt(th * th + 1.0));
                c  = 1.0 / sqrt(t * t + 1.0);
                sn = t * c;
                for (k = 0; k < 3; k++)            /* A <- A G (columns i, j) */
                {
                    double aki = A[k][i], akj = A[k][j];

                    A[k][i] = c * aki - sn * akj;
                    A[k][j] = sn * aki + c * akj;
                }
                for (k = 0; k < 3; k++)            /* A <- G^T A (rows i, j) */
                {
                    double aik = A[i][k], ajk = A[j][k];

                    A[i][k] = c * aik - sn * ajk;
                    A[j][k] = sn * aik + c * ajk;
                }
                for (k = 0; k < 3; k++)            /* eigenvectors as rows */
                {
                    double vik = V[i][k], vjk = V[j][k];

                    V[i][k] = c * vik - sn * vjk;
                    V[j][k] = sn * vik + c * vjk;
                }
            }
        }
    }
    for (i = 0; i < 3; i++) { lam[i] = A[i][i]; }
}

/* Is x inside the sample ellipsoid (the sphere of radius_now before the first
 * geometry pass)?  Gates the collisions until the first PB solve gives n_e. */
static gmx_bool mcmd_in_cloud(const t_mcionize *mc, const rvec x)
{
    double d[3] = {x[XX] - mc->com[0], x[YY] - mc->com[1], x[ZZ] - mc->com[2]};

    if (mc->ell_a[0] > 0)
    {
        double u = 0;
        int    i;

        for (i = 0; i < 3; i++)
        {
            double y = mc->ell_axis[i][0] * d[0] + mc->ell_axis[i][1] * d[1] + mc->ell_axis[i][2] * d[2];

            u += y * y / (mc->ell_a[i] * mc->ell_a[i]);
        }
        return u < 1.0;
    }
    return d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < sqr(mc->radius_now);
}

/* Centre of mass (stored in mc->com) and mass-weighted radius of gyration
 * about it, in nm. */
/* Weight of atom i in the sample geometry: its mass, or 0 for hydrogen
 * (geom_heavy, always on: fast protons leave the electron cloud behind and
 * would otherwise stretch R and the ellipsoid). */
#define MCMD_GEOM_W(i) ((mc->geom_heavy && mdatoms->massT[i] < 1.5) ? 0.0 : (double)mdatoms->massT[i])

static void mcmd_sample_geometry(t_mcionize *mc, const t_mdatoms *mdatoms,
                                 const t_state *state, double *rg)
{
    double mtot = 0, sum = 0;
    int    i, d;

    mc->com[0] = mc->com[1] = mc->com[2] = 0;
    for (i = 0; i < mdatoms->nr; i++)
    {
        for (d = 0; d < 3; d++)
        {
            mc->com[d] += MCMD_GEOM_W(i) * state->x[i][d];
        }
        mtot += MCMD_GEOM_W(i);
    }
    for (d = 0; d < 3; d++)
    {
        mc->com[d] /= mtot;
    }
    for (i = 0; i < mdatoms->nr; i++)
    {
        sum += MCMD_GEOM_W(i) * (sqr(state->x[i][0] - mc->com[0]) +
                                    sqr(state->x[i][1] - mc->com[1]) +
                                    sqr(state->x[i][2] - mc->com[2]));
    }
    *rg = sqrt(sum / mtot);

    {
        double S[3][3] = {{0}}, lam[3];
        int    e;

        for (i = 0; i < mdatoms->nr; i++)
        {
            double dx[3];

            for (d = 0; d < 3; d++) { dx[d] = state->x[i][d] - mc->com[d]; }
            for (d = 0; d < 3; d++)
            {
                for (e = 0; e < 3; e++) { S[d][e] += MCMD_GEOM_W(i) * dx[d] * dx[e]; }
            }
        }
        for (d = 0; d < 3; d++)
        {
            for (e = 0; e < 3; e++) { S[d][e] /= mtot; }
        }
        mcmd_jacobi3(S, lam, mc->ell_axis);
        for (d = 0; d < 3; d++)
        {
            /* a uniform ellipsoid has <x_i^2> = a_i^2 / 5; keep a floor so a
             * flat sample still has a finite thickness */
            mc->ell_a[d] = sqrt(5.0 * (lam[d] > 1e-4 * sqr(*rg) ? lam[d] : 1e-4 * sqr(*rg)));
        }
    }
}

/* Kinetic energy, in eV, of the electron emitted when an atom goes from
 * state `from` to `to` by a transition of the given type.  Taken from the
 * total binding energies, which is what makes it right for every shell:
 *   photoionization  E_ph - (T_from - T_to)   e.g. C 1s at 600 eV: 301 eV
 *   Auger            T_to - T_from            e.g. C KLL:          263 eV
 * Without total_energies_X.txt (allowed when charge transfer is off) the
 * photoelectron falls back to the outermost binding energy, which
 * overestimates core photoelectrons, and an Auger electron is scored as
 * trapped; mc->count_escape_fallback counts those. */
static double mcmd_emitted_electron_energy(t_mcionize *mc,
                                           const t_mcmd_atomdata *ad,
                                           const int from[MCMD_NSHELL],
                                           const int to[MCMD_NSHELL], int type)
{
    double T_from = MCMD_NOT_FOUND, T_to = MCMD_NOT_FOUND;

    if (ad->total_energies != NULL)
    {
        T_from = mcmd_dict_get(ad->total_energies, from);
        T_to   = mcmd_dict_get(ad->total_energies, to);
        if (to[0] == 0 && to[1] == 0 && to[2] == 0)
        {
            T_to = 0.0; /* a bare nucleus has no binding energy */
        }
    }

    if (T_from != MCMD_NOT_FOUND && T_to != MCMD_NOT_FOUND)
    {
        switch (type)
        {
            case 2:
            {
                /* The photon energy jitters with the XFEL bandwidth; a
                 * Box-Muller Gaussian from the module's own erand48 stream. */
                double E_ph = mc->photon_energy_eV;

                if (mc->bandwidth_fwhm > 0)
                {
                    double u1 = erand48(mc->erand_state), u2 = erand48(mc->erand_state);

                    E_ph += mc->bandwidth_fwhm * mc->photon_energy_eV / MCMD_FWHM_PER_SIGMA *
                            sqrt(-2.0 * log(u1 > 0 ? u1 : 1e-300)) * cos(2.0 * M_PI * u2);
                }
                return E_ph - (T_from - T_to);
            }
            case 0: return T_to - T_from;
            default: return 0.0;
        }
    }

    mc->count_escape_fallback += 1;
    if (type == 2)
    {
        double Ip = mcmd_dict_get(ad->energy_levels, from);

        return mc->photon_energy_eV - (Ip == MCMD_NOT_FOUND ? 0.0 : Ip);
    }
    return 0.0;
}

/* Electron temperature of the trapped cloud, eV: its kinetic energy is
 * 3/2 N T (the potential energy is the PB solution's). */
static double mcmd_plasma_Te(const t_mcionize *mc)
{
    if (mc->num_free_electrons <= 0)
    {
        return 0.0;
    }
    return (2.0 / 3.0) * mc->screen_E_trapped / mc->num_free_electrons;
}

#define MCMD_COLL_K    1.0899e-6   /* sqrt(8/pi m_e) pi a0^2 Ry^2, cm^3 s^-1 eV^3/2 */
#define MCMD_LOTZ_K    6.7e-7      /* Lotz (1968) Maxwellian coefficient, cm^3 s^-1 eV^3/2 per a (1e-14 cm^2 eV^2) */
#define MCMD_LOTZ_A    4.5
#define MCMD_SAHA_K    1.656e-22   /* (h^2/2 pi m_e)^3/2 / 2 at 1 eV, cm^3 eV^3/2   */

/* exp(x^2) erfc(x), x >= 0, without overflow. */
static double mcmd_erfcx(double x)
{
    if (x < 10.0)
    {
        return exp(x * x) * erfc(x);
    }
    else
    {
        double y = 1.0 / (x * x);

        return (1.0 - 0.5 * y + 0.75 * y * y - 1.875 * y * y * y + 6.5625 * y * y * y * y) /
               (x * sqrt(M_PI));
    }
}

/* exp(x) E1(x), x > 0: series below 1, continued fraction above
 * (Numerical Recipes expint, n = 1). */
static double mcmd_e1x(double x)
{
    int i;

    if (x <= 0)
    {
        return 0;
    }
    if (x <= 1.0)
    {
        double sum = 0, term = 1;

        for (i = 1; i < 60; i++)
        {
            term *= -x / i;               /* (-x)^i / i! */
            sum  -= term / i;             /* sum_i (-1)^(i+1) x^i / (i i!) */
            if (fabs(term / i) < 1e-17 * fabs(sum)) { break; }
        }
        return exp(x) * (-0.5772156649015329 - log(x) + sum);
    }
    else
    {
        double b = x + 1.0, c = 1e300, d = 1.0 / b, h = d;

        for (i = 1; i < 200; i++)
        {
            double an = -(double)i * i, del;

            b += 2.0;
            d  = 1.0 / (an * d + b);
            c  = b + an / c;
            del = c * d;
            h  *= del;
            if (fabs(del - 1.0) < 1e-15) { break; }
        }
        return h;                         /* = e^x E1(x) */
    }
}

/* q(b) of the Maxwellian rate; 0 if not finite or not positive. */
/* CRETIN 'colex' collision strength W(x) = c0 ln x + c1 + c2/(x+a0) +
 * c3/(x+a0)^2 integrated over a Maxwellian (manual p. 118; its b/(a0+1) term
 * carries c3, not c2 as printed - checked against numerical integration):
 *   q(b) = c0 E1(b)e^b + c1 + c3 b/A + (c2 - c3 b) b E1(Ab) e^(Ab),  A = a0+1 */
#define MCMD_COLEX_K 8.01e-8      /* cm^3 s^-1 eV^1/2 */
static double mcmd_colex_q(const t_mcmd_cchan *ch, double b)
{
    const double A = ch->a0 + 1.0;
    double       q = ch->c[0] * mcmd_e1x(b) + ch->c[1] + ch->c[3] * b / A +
                     (ch->c[2] - ch->c[3] * b) * b * mcmd_e1x(A * b);

    return (isfinite(q) && q > 0) ? q : 0.0;
}

static double mcmd_coll_q(const t_mcmd_cchan *ch, double b)
{
    double q = ch->c[0] + ch->c[1] * sqrt(M_PI * b) * mcmd_erfcx(sqrt(b)) +
               (ch->c[2] * b + ch->c[3]) * mcmd_e1x(b);

    return (isfinite(q) && q > 0) ? q : 0.0;
}

/* Build, once, the collisional channels out of every state of every element
 * present: ionization from collisional_parameters_X.txt, and for each of those
 * the inverse recombination attached to the final state.  Transition energies
 * from the total energies, weights from statistical_weight_X.txt. */
static void mcmd_coll_build(t_mcionize *mc, FILE *fplog)
{
    int el, maxch = 0;

    for (el = 0; el < MCMD_ATOM_TABLE_SIZE; el++)
    {
        t_mcmd_atomdata *ad = &mc->atomdata[el];
        int             *count;
        int              pass, s, t, n, q;
        int              nlotz_skip = 0;

        mc->cch[el] = NULL; mc->cch_start[el] = NULL;
        if (!ad->bPresent)
        {
            continue;
        }
        if (ad->total_energies == NULL)
        {
            gmx_fatal(FARGS, "mcmd-collisions needs total_energies_%s.txt.",
                      mcmd_idx2symbol(ad->elem));
        }
        snew(count, ad->num_rates);
        snew(mc->cch_start[el], ad->num_rates + 1);
        /* pass 0 counts, pass 1 fills */
        for (pass = 0; pass < 2; pass++)
        {
            int *fill = NULL;

            if (pass == 1)
            {
                n = 0;
                for (s = 0; s < ad->num_rates; s++)
                {
                    mc->cch_start[el][s] = n;
                    n                   += count[s];
                    if (count[s] > maxch) { maxch = count[s]; }
                }
                mc->cch_start[el][ad->num_rates] = n;
                snew(mc->cch[el], n > 0 ? n : 1);
                snew(fill, ad->num_rates);
            }
            for (s = 0; s < ad->num_coll; s++)
            {
                const t_mcmd_coll *co = &ad->coll[s];
                int                si = mcmd_state_index_get(&ad->rate_index, co->initial_state);
                double             Ti = mcmd_dict_get(ad->total_energies, co->initial_state);
                double             gi = mcmd_weight_for_state(&ad->weights, co->initial_state);

                if (si < 0)
                {
                    gmx_fatal(FARGS, "Collisional state %d %d %d of %s has no rate data.",
                              co->initial_state[0], co->initial_state[1], co->initial_state[2],
                              mcmd_idx2symbol(ad->elem));
                }
                for (t = 0; t < co->num_transitions; t++)
                {
                    const int *fs = co->final_states[t];
                    int        sf = mcmd_state_index_get(&ad->rate_index, fs);
                    double     Tf = mcmd_dict_get(ad->total_energies, fs);
                    double     gf = mcmd_weight_for_state(&ad->weights, fs);
                    double     dE = Ti - Tf;
                    int        k;

                    if (sf < 0 || Ti == MCMD_NOT_FOUND || Tf == MCMD_NOT_FOUND || dE <= 0 ||
                        gi <= 0 || gf <= 0)
                    {
                        gmx_fatal(FARGS, "Collisional transition %d %d %d -> %d %d %d of %s: "
                                  "missing rate data, total energy or weight (dE %g, g %g %g).",
                                  co->initial_state[0], co->initial_state[1], co->initial_state[2],
                                  fs[0], fs[1], fs[2], mcmd_idx2symbol(ad->elem), dE, gi, gf);
                    }
                    if (pass == 0)
                    {
                        count[si]++;
                        if (mc->coll_recomb && !(mc->coll_no_h_recomb && el == mc->elem_H)) { count[sf]++; }
                        continue;
                    }
                    for (k = 0; k < 2; k++)
                    {
                        t_mcmd_cchan *ch;
                        int           from = (k == 0) ? si : sf;

                        if (k == 1 && !mc->coll_recomb) { break; }
                        /* mcmd-coll-h-recombination = 0: a proton stays a proton.
                         * Charge transfer fires as soon as the geometry allows,
                         * so a hydrogen that recombines hands the electron on
                         * at once and shuttles plasma electrons into heavy ions;
                         * at these T_e and n_e Saha keeps hydrogen ionized
                         * anyway, and the published model never let an ionized
                         * hydrogen take an electron back. */
                        if (k == 1 && mc->coll_no_h_recomb && el == mc->elem_H) { break; }
                        ch = &mc->cch[el][mc->cch_start[el][from] + fill[from]++];
                        ch->recomb = k;
                        ch->excite = 0;
                        ch->final_state[0] = (k == 0) ? fs[0] : co->initial_state[0];
                        ch->final_state[1] = (k == 0) ? fs[1] : co->initial_state[1];
                        ch->final_state[2] = (k == 0) ? fs[2] : co->initial_state[2];
                        ch->c[0] = co->coll_rates[t][0]; ch->c[1] = co->coll_rates[t][1];
                        ch->c[2] = co->coll_rates[t][2]; ch->c[3] = co->coll_rates[t][3];
                        ch->s0   = co->coll_rates[t][4];
                        ch->dE   = dE;
                        ch->gratio = gi / gf;
                        ch->rate = 0;
                    }
                }
            }
            /* mcmd-coll-data = lotz (default since 2026-10-28): electron-impact ionization from Lotz's
             * formula (1967/68, b = 0, a = 4.5e-14 cm^2 eV^2) instead of a data file.  Every state with
             * transition data loses one electron from each occupied shell, wherever that state exists:
             * q = the shell's occupation, P = dE from the total energies.  Gives the same channels as the
             * CRETIN fit files (H C N O S) and coefficients 13-25 % above them; the explosion moves less
             * (ion KE -2.5 % at 600 eV 1e7, the rest within the seed noise).  Recombination by detailed
             * balance as for the fits. */
            for (s = 0; mc->coll_lotz && s < ad->num_rates; s++)
            {
                const int *is = ad->rates[s].initial_state;
                int        sh;

                if (mcmd_state_index_get(&ad->rate_index, is) != s)
                {
                    continue;                                   /* a repeated initial state: the first wins */
                }
                for (sh = 0; sh < MCMD_NSHELL; sh++)
                {
                    int    fs[MCMD_NSHELL], sf, kk;
                    double Ti, Tf, gi, gf, dE;

                    if (is[sh] <= 0)
                    {
                        continue;
                    }
                    fs[0] = is[0]; fs[1] = is[1]; fs[2] = is[2];
                    fs[sh]--;
                    sf = mcmd_state_index_get(&ad->rate_index, fs);
                    Ti = mcmd_dict_get(ad->total_energies, is);
                    Tf = mcmd_dict_get(ad->total_energies, fs);
                    gi = mcmd_weight_for_state(&ad->weights, is);
                    gf = mcmd_weight_for_state(&ad->weights, fs);
                    dE = Ti - Tf;
                    if (sf < 0 || Ti == MCMD_NOT_FOUND || Tf == MCMD_NOT_FOUND || dE <= 0 || gi <= 0 || gf <= 0)
                    {
                        if (pass == 0) { nlotz_skip++; }
                        continue;
                    }
                    if (pass == 0)
                    {
                        count[s]++;
                        if (mc->coll_recomb && !(mc->coll_no_h_recomb && el == mc->elem_H)) { count[sf]++; }
                        continue;
                    }
                    for (kk = 0; kk < 2; kk++)
                    {
                        t_mcmd_cchan *ch;
                        int           from = (kk == 0) ? s : sf;

                        if (kk == 1 && (!mc->coll_recomb || (mc->coll_no_h_recomb && el == mc->elem_H)))
                        {
                            break;
                        }
                        ch = &mc->cch[el][mc->cch_start[el][from] + fill[from]++];
                        ch->recomb = kk;
                        ch->excite = 0;
                        ch->final_state[0] = (kk == 0) ? fs[0] : is[0];
                        ch->final_state[1] = (kk == 0) ? fs[1] : is[1];
                        ch->final_state[2] = (kk == 0) ? fs[2] : is[2];
                        ch->c[0] = ch->c[1] = ch->c[2] = ch->c[3] = 0;
                        ch->s0     = ch->a0 = 0;
                        ch->dE     = dE;
                        ch->gratio = gi / gf;
                        ch->g      = NULL;
                        ch->lotz   = MCMD_LOTZ_K * MCMD_LOTZ_A * is[sh];
                        ch->rate   = 0;
                    }
                }
            }
            /* Tabulated coefficients (collisional_ionization_table_X.txt and
             * collisional_excitation_table_X.txt, 2026-10-28), the generator-
             * independent alternative to the fits: the same channels, with
             * g = ln <sigma v> + dE/T on the table's grid in place of the fit
             * (see mcmd_coll_update_rates).  q = 0: ionization, 1: excitation. */
            for (q = 0; q < 2; q++)
            {
                const t_mcmd_ctab *tb = (q == 0) ? ad->coll_tab : (mc->coll_excite ? ad->colex_tab : NULL);

                if (tb == NULL)
                {
                    continue;
                }
                if (pass == 1 && mc->cch_g[el][q] == NULL)
                {
                    snew(mc->cch_g[el][q], (size_t)tb->ntr * tb->nT);
                }
                for (t = 0; t < tb->ntr; t++)
                {
                    const int *is = tb->ini[t], *fs = tb->fin[t];
                    int        si = mcmd_state_index_get(&ad->rate_index, is);
                    int        sf = mcmd_state_index_get(&ad->rate_index, fs);
                    double     Ti = mcmd_dict_get(ad->total_energies, is);
                    double     Tf = mcmd_dict_get(ad->total_energies, fs);
                    double     gi = mcmd_weight_for_state(&ad->weights, is);
                    double     gf = mcmd_weight_for_state(&ad->weights, fs);
                    double     dE = Ti - Tf;
                    int        dn = (fs[0] + fs[1] + fs[2]) - (is[0] + is[1] + is[2]);
                    int        kk, j;
                    double    *g;

                    if (si < 0 || sf < 0 || Ti == MCMD_NOT_FOUND || Tf == MCMD_NOT_FOUND || dE <= 0 ||
                        gi <= 0 || gf <= 0 || dn != (q == 0 ? -1 : 0))
                    {
                        gmx_fatal(FARGS, "Tabulated collisional %s %d %d %d -> %d %d %d of %s: missing rate data, "
                                  "total energy or weight, or a wrong electron count (dE %g, g %g %g, change %d).",
                                  q == 0 ? "ionization" : "excitation", is[0], is[1], is[2], fs[0], fs[1], fs[2],
                                  mcmd_idx2symbol(ad->elem), dE, gi, gf, dn);
                    }
                    if (pass == 0)
                    {
                        count[si]++;
                        if (q == 1 || (mc->coll_recomb && !(mc->coll_no_h_recomb && el == mc->elem_H)))
                        {
                            count[sf]++;
                        }
                        continue;
                    }
                    g = mc->cch_g[el][q] + (size_t)t * tb->nT;
                    {
                        const double *lk = tb->lnk + (size_t)t * tb->nT;
                        int           j0 = -1, j1 = -1;

                        for (j = 0; j < tb->nT; j++)
                        {
                            g[j] = lk[j] + dE / exp(tb->lnT[j]);
                            if (isfinite(lk[j]))
                            {
                                if (j0 < 0) { j0 = j; } else if (j1 < 0) { j1 = j; }
                            }
                        }
                        /* a 0 (too small to write) is extrapolated linearly in ln T from the first two
                         * values given (g is smooth there); a transition that is 0 throughout never fires */
                        for (j = 0; j < tb->nT; j++)
                        {
                            if (isfinite(lk[j]))
                            {
                                continue;
                            }
                            if (j0 < 0)
                            {
                                g[j] = -1e30;
                            }
                            else if (j1 < 0)
                            {
                                g[j] = g[j0];
                            }
                            else
                            {
                                g[j] = g[j0] + (g[j1] - g[j0]) / (tb->lnT[j1] - tb->lnT[j0]) * (tb->lnT[j] - tb->lnT[j0]);
                            }
                        }
                    }
                    for (kk = 0; kk < 2; kk++)
                    {
                        t_mcmd_cchan *ch;
                        int           from = (kk == 0) ? si : sf;

                        if (q == 0 && kk == 1 && (!mc->coll_recomb || (mc->coll_no_h_recomb && el == mc->elem_H)))
                        {
                            break;
                        }
                        ch = &mc->cch[el][mc->cch_start[el][from] + fill[from]++];
                        ch->recomb = (q == 0) ? kk : 0;
                        ch->excite = (q == 0) ? 0 : kk + 1;
                        ch->final_state[0] = (kk == 0) ? fs[0] : is[0];
                        ch->final_state[1] = (kk == 0) ? fs[1] : is[1];
                        ch->final_state[2] = (kk == 0) ? fs[2] : is[2];
                        ch->c[0] = ch->c[1] = ch->c[2] = ch->c[3] = 0;
                        ch->s0   = ch->a0 = 0;
                        ch->dE   = dE;
                        ch->gratio = gi / gf;
                        ch->g    = g;
                        ch->lnT  = tb->lnT;
                        ch->nT   = tb->nT;
                        ch->rate = 0;
                    }
                }
            }
            /* mcmd-coll-excitation: excitation lower -> upper from the
             * colex2 fit, and its inverse de-excitation attached to the upper
             * state.  dE = T(lower) - T(upper) > 0 from the total energies. */
            for (s = 0; mc->coll_excite && s < ad->num_colex; s++)
            {
                const t_mcmd_coll *co = &ad->colex[s];
                int                si = mcmd_state_index_get(&ad->rate_index, co->initial_state);
                double             Ti = mcmd_dict_get(ad->total_energies, co->initial_state);
                double             gi = mcmd_weight_for_state(&ad->weights, co->initial_state);

                for (t = 0; t < co->num_transitions; t++)
                {
                    const int *fs = co->final_states[t];
                    int        sf = mcmd_state_index_get(&ad->rate_index, fs);
                    double     Tf = mcmd_dict_get(ad->total_energies, fs);
                    double     gf = mcmd_weight_for_state(&ad->weights, fs);
                    double     dE = Ti - Tf;
                    int        k;

                    if (si < 0 || sf < 0 || Ti == MCMD_NOT_FOUND || Tf == MCMD_NOT_FOUND ||
                        dE <= 0 || gi <= 0 || gf <= 0)
                    {
                        gmx_fatal(FARGS, "Collisional excitation %d %d %d -> %d %d %d of %s: "
                                  "missing rate data, total energy or weight (dE %g, g %g %g).",
                                  co->initial_state[0], co->initial_state[1], co->initial_state[2],
                                  fs[0], fs[1], fs[2], mcmd_idx2symbol(ad->elem), dE, gi, gf);
                    }
                    if (pass == 0)
                    {
                        count[si]++;
                        count[sf]++;
                        continue;
                    }
                    for (k = 0; k < 2; k++)
                    {
                        int           from = (k == 0) ? si : sf;
                        t_mcmd_cchan *ch   = &mc->cch[el][mc->cch_start[el][from] + fill[from]++];

                        ch->recomb = 0;
                        ch->excite = k + 1;
                        ch->final_state[0] = (k == 0) ? fs[0] : co->initial_state[0];
                        ch->final_state[1] = (k == 0) ? fs[1] : co->initial_state[1];
                        ch->final_state[2] = (k == 0) ? fs[2] : co->initial_state[2];
                        ch->c[0] = co->coll_rates[t][0]; ch->c[1] = co->coll_rates[t][1];
                        ch->c[2] = co->coll_rates[t][2]; ch->c[3] = co->coll_rates[t][3];
                        ch->a0   = co->coll_rates[t][4];
                        ch->s0   = 0;
                        ch->dE   = dE;
                        ch->gratio = gi / gf;   /* de-excitation: g(lower) / g(upper) */
                        ch->rate = 0;
                    }
                }
            }
            sfree(fill);
        }
        sfree(count);
        if (ad->max_transitions + maxch > mc->mc_scratch_size)
        {
            mc->mc_scratch_size = ad->max_transitions + maxch;
        }
        if (fplog != NULL)
        {
            fprintf(fplog, "mcionize: collisional channels for %s: %d (ionization from %s%s)\n",
                    mcmd_idx2symbol(ad->elem), mc->cch_start[el][ad->num_rates],
                    mc->coll_lotz ? "Lotz's formula" : (ad->coll_tab ? "the table" : "the fit file"),
                    nlotz_skip ? "; some shells without a final state, total energy or weight skipped" : "");
        }
    }
    srenew(mc->mc_transitions, mc->mc_scratch_size);
    srenew(mc->mc_dt, mc->mc_scratch_size);
    snew(mc->mc_cref, mc->mc_scratch_size);
}

/* Evaluate every channel's rate for the current plasma, once per step.  The
 * electrons fill the sample sphere, n_e = N_trapped / (4/3 pi R^3). */
static void mcmd_coll_update_rates(t_mcionize *mc)
{
    const double N = mc->num_free_electrons, R = mc->radius_now;
    const double T = mcmd_plasma_Te(mc);
    double       ne;
    int          el, i;
    gmx_bool     clamped = FALSE;

    mc->coll_active = (N >= 1.0 && T > 0 && R > 0);
    mc->coll_T      = T;
    mc->coll_ne     = ne = mc->coll_active ? N / (4.0 / 3.0 * M_PI * R * R * R * 1e-21) : 0;
    if (!mc->coll_active)
    {
        return;
    }
    for (el = 0; el < MCMD_ATOM_TABLE_SIZE; el++)
    {
        int n;

        if (mc->cch[el] == NULL) { continue; }
        n = mc->cch_start[el][mc->atomdata[el].num_rates];
        for (i = 0; i < n; i++)
        {
            t_mcmd_cchan *ch = &mc->cch[el][i];
            double        b  = ch->dE / T;
            double        k;

            if (ch->lotz > 0)
            {
                /* Lotz, without its Boltzmann factor: 6.7e-7 a q T^-3/2 e^b E1(b) / b */
                k        = ch->lotz * pow(T, -1.5) * mcmd_e1x(b) / b;
                ch->rate = ch->recomb ? ne * ne * k * ch->gratio * MCMD_SAHA_K * pow(T, -1.5) : ne * k * exp(-b);
                continue;
            }
            if (ch->g != NULL)
            {
                /* tabulated: g = ln <sigma v> + dE/T is smooth in ln T (the Boltzmann factor taken out),
                 * so it is interpolated linearly there; outside the grid it is held at the edge */
                const double *x = ch->lnT;
                double        lt = log(T), w;
                int           lo = 0, hi = ch->nT - 1;

                if (lt <= x[0] || lt >= x[hi])
                {
                    clamped = TRUE;
                    lo      = (lt <= x[0]) ? 0 : hi - 1;
                    w       = (lt <= x[0]) ? 0.0 : 1.0;
                }
                else
                {
                    while (hi - lo > 1)
                    {
                        int mid = (lo + hi) / 2;

                        if (x[mid] <= lt) { lo = mid; } else { hi = mid; }
                    }
                    w = (lt - x[lo]) / (x[lo + 1] - x[lo]);
                }
                k = exp((1.0 - w) * ch->g[lo] + w * ch->g[lo + 1]);   /* coefficient without e^-b */
                if (ch->excite)
                {
                    ch->rate = (ch->excite == 1) ? ne * k * exp(-b) : ne * k * ch->gratio;
                }
                else
                {
                    ch->rate = ch->recomb ? ne * ne * k * ch->gratio * MCMD_SAHA_K * pow(T, -1.5)
                                          : ne * k * exp(-b);
                }
                continue;
            }
            if (ch->excite)
            {
                /* CRETIN colex: k_up = 8.01e-8 q(b) e^-b / sqrt(T) cm^3/s per ion
                 * in the lower level; de-excitation by detailed balance,
                 * k_down = k_up (g_lower/g_upper) e^b (the e^-b cancels). */
                k        = MCMD_COLEX_K * mcmd_colex_q(ch, b) / sqrt(T);
                ch->rate = (ch->excite == 1) ? ne * k * exp(-b) : ne * k * ch->gratio;
                continue;
            }
            k = MCMD_COLL_K * ch->s0 * mcmd_coll_q(ch, b) / (ch->dE * sqrt(T));

            ch->rate = ch->recomb ? ne * ne * k * ch->gratio * MCMD_SAHA_K * pow(T, -1.5)
                                  : ne * k * exp(-b);
        }
    }
    if (clamped)
    {
        if (mc->ctab_clamped == 0)
        {
            fprintf(stderr, "\nmcionize: T_e = %g eV is outside a collisional table's grid; its coefficients are "
                    "held at the grid edge (counted in the log at the end)\n", T);
        }
        mc->ctab_clamped += 1;
    }
}

/* ---------------------------------------------------------------- *
 * Self-consistent trapped electrons (mcmd-electrons = pbgrid): the shared
 * pieces.  The grid solve itself is described above mcmd_pbg_conv().
 *
 * The trapped electrons are a classical Boltzmann gas at T_e in the field of
 * the ions, the density that minimizes the free energy
 *     F[n] = -int n Phi + 1/2 int int n k G n + T int n (ln n - 1),   int n = N,
 * Phi = k G * q_ion, with every charge smoothed: G(r, s) = erf(r/s) / r
 * (2 / (s sqrt(pi)) at r = 0).  F is convex; its minimum is
 * n ~ exp(-V / T), V = -Phi + k G * n: Poisson-Boltzmann.  Cold, the
 * electrons fill the core until it is neutral and the net charge sits on the
 * outside; hot, they spread out.  Solved every mcmd-pb-interval (default
 * 0.2 fs since 2026-10-26), on the master.  With the current model 0.2 fs
 * (with MD 2 as and Newton 1e-2, 4 seeds) gives the 0.1 fs results within
 * the seed noise in ion energy, charge, escape and ion directions; the 3 %
 * drop seen on 2026-10-23 was with the older model.  From 0.5 fs the cloud
 * stays too hot and the ion energy drops by up to 20 % (2026-10-22 scan).
 *
 * The trapped electrons' total energy (kinetic 3/2 N T plus their
 * electrostatic energy) changes only by explicit exchanges: the work of their
 * field on the moving ions, electrons born or captured (with their kinetic
 * energy and the potential energy at their atom), binding energy paid in
 * collisions, the jump in the ion-electron energy when a charge changes, and
 * evaporation.  T_e is the temperature at which the equilibrium density
 * reproduces that energy.
 *
 * The escape test uses the same potential: an electron leaving atom k
 * escapes if E > D_k, the depth at k with k's own ion left out.
 * Evaporation: the Maxwell tail above each cell's depth, at the e-e
 * collision rate.
 * ---------------------------------------------------------------- */
#define MCMD_PB_SQRTPI 1.7724538509055159

/* Wall-clock timers of the module's sections, on the master (mc->tm[], ms):
 * where the "Rest" of GROMACS's cycle accounting goes.  The force sections are
 * timed on the master rank only; under -pd the others do their share at the
 * same time. */
#define MCMD_TM_CT     0   /* charge-transfer pass                                */
#define MCMD_TM_PLASMA 1   /* geometry, plasma update, evaporation, rates         */
#define MCMD_TM_MC     2   /* Monte Carlo (photo, Auger, collisional, escape)     */
#define MCMD_TM_OUT    3   /* per-step output                                     */
#define MCMD_TM_STEP   4   /* the whole master part of mcionize_step              */
#define MCMD_TM_PBF    5   /* PB forces, solves included                          */
#define MCMD_TM_SCREEN 8   /* the whole of mcionize_screen_forces                 */
#define MCMD_TM_SOLVE  9   /* PB solves (master)                                  */

static double mcmd_pb_G(double r, double s)
{
    return (r < 1e-6 * s) ? 2.0 / (s * MCMD_PB_SQRTPI) : erf(r / s) / r;
}

/* dG/dr */
static double mcmd_pb_dG(double r, double s)
{
    double u = r / s;

    if (u > 6.0)
    {
        return -1.0 / (r * r);
    }
    return (2.0 / MCMD_PB_SQRTPI * u * exp(-u * u) - erf(u)) / (r * r);
}

static double mcmd_pb_d2(const rvec a, const rvec b)
{
    return sqr(a[XX] - b[XX]) + sqr(a[YY] - b[YY]) + sqr(a[ZZ] - b[ZZ]);
}

static double mcmd_pbg_depth(const t_mcionize *mc, const t_mdatoms *md, rvec *x, int k);

/* Depth, eV (> 0 bound), of an electron at atom k: ions other than k and the
 * electron cloud, see mcmd_pbg_depth(). */
static double mcmd_pb_depth(const t_mcionize *mc, const t_mdatoms *md, rvec *x, int k)
{
    return mcmd_pbg_depth(mc, md, x, k);
}

/* Master-only sections (the PB solve, the long-range ion potential): the
 * other thread-MPI ranks would spin in the broadcast that follows (sched_yield,
 * all their cores busy) and take the cores from the FFTW threads of the
 * master.  They sleep on a condition variable instead until the master is
 * done.  Each rank counts the sections it has passed (park_seq); the master
 * publishes its count when done, so a rank that arrives late does not wait.
 * Thread-MPI ranks are threads of one process; with other MPI (several
 * processes) this is a no-op. */
#if defined(GMX_THREADS)
static pthread_mutex_t mcmd_park_mx   = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  mcmd_park_cv   = PTHREAD_COND_INITIALIZER;
static long            mcmd_park_done = 0;
#endif

static void mcmd_park(t_mcionize *mc, const t_commrec *cr)
{
#if defined(GMX_THREADS)
    if (!PAR(cr))
    {
        return;
    }
    mc->park_seq++;
    pthread_mutex_lock(&mcmd_park_mx);
    if (MASTER(cr))
    {
        mcmd_park_done = mc->park_seq;
        pthread_cond_broadcast(&mcmd_park_cv);
    }
    else
    {
        while (mcmd_park_done < mc->park_seq)
        {
            pthread_cond_wait(&mcmd_park_cv, &mcmd_park_mx);
        }
    }
    pthread_mutex_unlock(&mcmd_park_mx);
#endif
}

static double mcmd_pb_wallms(void)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return 1e3 * tv.tv_sec + 1e-3 * tv.tv_usec;
}

/* ---------------------------------------------------------------- *
 * Self-consistent trapped electrons on a 3D grid (mcmd-electrons = pbgrid):
 * the full Poisson-Boltzmann solve.
 *
 * The electrons are a density on a cubic grid (cells per side chosen at each
 * rebuild from mcmd-pb-cells-per-smoothing, at least 24, at most
 * mcmd-pb-grid-limit) over at most +-mcmd-pb-box-max R about the heavy-atom
 * centre, cell counts n_i
 * (sum N).  The ions are deposited on the same grid (cloud-in-cell; ions
 * outside the box are added point by point).  Every cell is a Gaussian of
 * width sigma = the grid spacing, so two cells interact with
 *     K(r) = k erf(r / 2 sigma) / r,
 * an open-boundary convolution (Hockney: zero-padded (2n)^3 grid, FFTW).
 * The electrons minimize the same free energy as on the atoms,
 *     F(n) = -sum_i n_i Phi_i + 1/2 n.K n + T sum_i n_i [ln(n_i / dV) - 1],
 * Phi = K * q the ions' attraction: n_i ~ exp(-V_i / T), V = -Phi + K n.
 * Newton with Jacobi-preconditioned conjugate gradients (one convolution per
 * CG step), warm-started from the last solution; the grid is rebuilt (and
 * the last solution interpolated onto it) when R has changed by 10 % or the
 * centre has moved by a cell.
 *
 * Unlike the atom version the electrons do not move with the atoms between
 * solves: the electron potential K n stays on the grid (scaled by
 * N / N_solve) and the ions move through it.  The force on an ion is
 * q grad(K n) at its position, from the whole cloud: there is no cutoff, so
 * run with a Coulomb cutoff that covers the sample (the ion-ion forces must
 * be long-range too).  Outside the grid the cloud acts as a point charge at
 * its centroid.  Energy bookkeeping, escape test, deposits, collisions at
 * the local density and evaporation (per cell) are as in the ion-weight version (removed 2026-10-22).
 * ---------------------------------------------------------------- */
#define MCMD_PBG_IDX(i, j, k) ((((size_t)(i)) * ng + (j)) * ng + (k))

/* out = K * a on the grid (open boundaries) */
static void mcmd_pbg_conv(t_mcionize *mc, const double *a, double *out)
{
    const int    ng = mc->pbg_n, m = 2 * ng;
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    size_t       c;
    int          i, j, k;
    const double tc0 = mc->tdet ? mcmd_pb_wallms() : 0;
    double       tf;

#pragma omp parallel for
    for (i = 0; i < m; i++)
    {
        memset(mc->pbg_buf + (size_t)i * m * m, 0, sizeof(float) * (size_t)m * m);
    }
#pragma omp parallel for private(j, k)
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                mc->pbg_buf[((size_t)i * m + j) * m + k] = (float)a[MCMD_PBG_IDX(i, j, k)];
            }
        }
    }
    tf = mc->tdet ? mcmd_pb_wallms() : 0;
    fftwf_execute(mc->pbg_fwd);
    if (mc->tdet) { mc->td[MCMD_TD_FFT] += mcmd_pb_wallms() - tf; }
#pragma omp parallel for
    for (c = 0; c < nh; c++)
    {
        mc->pbg_hat[c][0] *= mc->pbg_Kh[c];
        mc->pbg_hat[c][1] *= mc->pbg_Kh[c];
    }
    tf = mc->tdet ? mcmd_pb_wallms() : 0;
    fftwf_execute(mc->pbg_bwd);
    if (mc->tdet) { mc->td[MCMD_TD_FFT] += mcmd_pb_wallms() - tf; }
#pragma omp parallel for private(j, k)
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                out[MCMD_PBG_IDX(i, j, k)] = mc->pbg_buf[((size_t)i * m + j) * m + k];
            }
        }
    }
    mc->pbg_nconv++;
    if (mc->tdet) { mc->td[MCMD_TD_CONV] += mcmd_pb_wallms() - tc0; }
}

/* Smoothing length of the grid kernel k erf(r/s)/r: 2 dx, or the physical
 * floor s_e (mcmd-pb-smoothing) if larger. */
static double mcmd_pbg_s(const t_mcionize *mc)
{
    const double s = MCMD_PBG_SFAC * mc->pbg_dx;

    return (mc->pbg_se > s) ? mc->pbg_se : s;
}

/* Kernel of the current spacing, Fourier transformed, normalization included. */
static void mcmd_pbg_kernel(t_mcionize *mc)
{
    const int    ng = mc->pbg_n, m = 2 * ng;
    const double s  = mcmd_pbg_s(mc);
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    size_t       c;
    int          i, j, k;

#pragma omp parallel for private(j, k)
    for (i = 0; i < m; i++)
    {
        double xi = ((i < m - i) ? i : m - i) * mc->pbg_dx;

        for (j = 0; j < m; j++)
        {
            double yj = ((j < m - j) ? j : m - j) * mc->pbg_dx;

            for (k = 0; k < m; k++)
            {
                double zk = ((k < m - k) ? k : m - k) * mc->pbg_dx;

                mc->pbg_buf[((size_t)i * m + j) * m + k] =
                    (float)(MCMD_COULOMB_EV_NM * mcmd_pb_G(sqrt(xi * xi + yj * yj + zk * zk), s));
            }
        }
    }
    mc->pbg_K0 = MCMD_COULOMB_EV_NM * mcmd_pb_G(0, s);
    fftwf_execute(mc->pbg_fwd);
#pragma omp parallel for
    for (c = 0; c < nh; c++)
    {
        mc->pbg_Kh[c] = mc->pbg_hat[c][0] / ((double)m * m * m);
    }
}

/* The same kernel on the FFT of 2 (ng + M) per side (ion potential with the
 * margin, mcmd_pbg_conv_margin). */
static void mcmd_pbg_kernel_margin(t_mcionize *mc)
{
    const int    ng = mc->pbg_n, m = 2 * (ng + MCMD_PBG_MARG);
    const double s  = mcmd_pbg_s(mc);
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    size_t       c;
    int          i, j, k;

    for (i = 0; i < m; i++)
    {
        double xi = ((i < m - i) ? i : m - i) * mc->pbg_dx;

        for (j = 0; j < m; j++)
        {
            double yj = ((j < m - j) ? j : m - j) * mc->pbg_dx;

            for (k = 0; k < m; k++)
            {
                double zk = ((k < m - k) ? k : m - k) * mc->pbg_dx;

                mc->pbg_buf[((size_t)i * m + j) * m + k] =
                    (float)(MCMD_COULOMB_EV_NM * mcmd_pb_G(sqrt(xi * xi + yj * yj + zk * zk), s));
            }
        }
    }
    fftwf_execute(mc->pbg_fwdE);
    for (c = 0; c < nh; c++)
    {
        mc->pbg_KhE[c] = mc->pbg_hat[c][0] / ((double)m * m * m);
    }
}

/* Cloud-in-cell weights of point x: corner cell i0 and fractions f; FALSE if
 * outside the grid. */
static gmx_bool mcmd_pbg_cic(const t_mcionize *mc, const real *x, int i0[3], double f[3])
{
    int d;

    for (d = 0; d < 3; d++)
    {
        double u = (x[d] - mc->pbg_lo[d]) / mc->pbg_dx;

        i0[d] = (int)floor(u);
        f[d]  = u - i0[d];
        if (i0[d] < 0 || i0[d] >= mc->pbg_n - 1)
        {
            return FALSE;
        }
    }
    return TRUE;
}

/* Ion attraction Phi (eV) on the grid for charges q at positions x; returns
 * the number of charged ions outside the box. */
static int mcmd_pbg_ionpot_q(t_mcionize *mc, const t_mdatoms *md, rvec *x, const real *qq, double *Phi);

static int mcmd_pbg_ionpot(t_mcionize *mc, const t_mdatoms *md, rvec *x, double *Phi)
{
    return mcmd_pbg_ionpot_q(mc, md, x, md->chargeA, Phi);
}

/* Potential (eV per e) of the charged ions that CIC cannot deposit (outside
 * the box or in its last cell), added to Phi, with the grid's smoothed kernel
 * G(r, s), s = 2 dx (bare 1/r beyond 4 s, where they agree).  stride 0: exact on every grid
 * point; otherwise each ion on the sub-lattice of stride S (S dx <= d / D, d
 * its distance from the box, S a power of two <= 16), interpolated. */
static void mcmd_pbg_outpot(t_mcionize *mc, const t_mdatoms *md, rvec *x, const real *qq, double *Phi, int lat,
                            const gmx_bool *skip)
{
    const int    ng = mc->pbg_n;
    const double dx = mc->pbg_dx, s = mcmd_pbg_s(mc);
    double      *c = NULL, *fr;
    int         *li, *node, *lev, L, S, nc, a, i, j, k, d, i0[3];
    double       f[3];
    int          nlev[5] = {0, 0, 0, 0, 0};
    double       tsum[5] = {0, 0, 0, 0, 0}, tint[5] = {0, 0, 0, 0, 0}, tl;

    snew(lev, md->nr);
    for (a = 0; a < md->nr; a++)
    {
        double dd = 0;

        lev[a] = -1;
        if (qq[a] == 0 || (skip != NULL && skip[a]) || mcmd_pbg_cic(mc, x[a], i0, f))
        {
            continue;
        }
        for (d = 0; d < 3; d++)
        {
            double u = mc->pbg_lo[d] - x[a][d], v = x[a][d] - (mc->pbg_lo[d] + (ng - 1) * dx);

            u   = (u > v) ? u : v;
            dd += (u > 0) ? u * u : 0;
        }
        dd     = sqrt(dd);
        lev[a] = 0;
        while (lat && lev[a] < 4 && (2 << lev[a]) * dx <= dd / MCMD_PBG_OUTLAT_D) { lev[a]++; }
        nlev[lev[a]]++;
    }
    snew(li, ng); snew(fr, ng); snew(node, ng + 1);
    for (L = 0; L <= 4; L++)
    {
        S  = 1 << L;
        nc = (ng - 1 + S - 1) / S + 1;
        for (a = 0; a < md->nr && lev[a] != L; a++) { }
        if (a == md->nr)
        {
            continue;
        }
        tl = mcmd_pb_wallms();
        for (i = 0; i < nc; i++) { node[i] = (i * S < ng - 1) ? i * S : ng - 1; }
        snew(c, (size_t)nc * nc * nc);
        for (; a < md->nr; a++)
        {
            if (lev[a] != L)
            {
                continue;
            }
            for (i = 0; i < nc; i++)
            {
                double ex = mc->pbg_lo[0] + node[i] * dx - x[a][0];

                for (j = 0; j < nc; j++)
                {
                    double ey = mc->pbg_lo[1] + node[j] * dx - x[a][1];

                    for (k = 0; k < nc; k++)
                    {
                        double ez = mc->pbg_lo[2] + node[k] * dx - x[a][2];

                        double r = sqrt(ex * ex + ey * ey + ez * ez);

                        c[((size_t)i * nc + j) * nc + k] += MCMD_COULOMB_EV_NM * qq[a] * (r < 4 * s ? mcmd_pb_G(r, s) : 1.0 / r);
                    }
                }
            }
        }
        tsum[L] = mcmd_pb_wallms() - tl;
        tl      = mcmd_pb_wallms();
        if (S == 1)
        {
            for (i = 0; i < ng * ng * ng; i++) { Phi[i] += c[i]; }
        }
        else
        {
            for (i = 0; i < ng; i++)
            {
                li[i] = (i / S < nc - 1) ? i / S : nc - 2;
                fr[i] = (double)(i - node[li[i]]) / (node[li[i] + 1] - node[li[i]]);
            }
            for (i = 0; i < ng; i++)
            {
                for (j = 0; j < ng; j++)
                {
                    for (k = 0; k < ng; k++)
                    {
                        double v = 0;
                        int    u, w, z;

                        for (u = 0; u < 2; u++)
                        {
                            for (w = 0; w < 2; w++)
                            {
                                for (z = 0; z < 2; z++)
                                {
                                    v += (u ? fr[i] : 1 - fr[i]) * (w ? fr[j] : 1 - fr[j]) * (z ? fr[k] : 1 - fr[k]) *
                                         c[((size_t)(li[i] + u) * nc + li[j] + w) * nc + li[k] + z];
                                }
                            }
                        }
                        Phi[MCMD_PBG_IDX(i, j, k)] += v;
                    }
                }
            }
        }
        tint[L] = mcmd_pb_wallms() - tl;
        sfree(c);
    }
    sfree(li); sfree(fr); sfree(node); sfree(lev);
    if (mc->pbg_outlat == 2 && lat)
    {
        if (mc->fop == NULL)
        {
            mc->fop = mcmd_open_output("pbg_outpot.txt", "w");
            fprintf(mc->fop, "# per call: step ng dx[nm] ions at stride 1 2 4 8 16 | sum time [ms] at 1 2 4 8 16 | "
                    "interpolation/add time [ms] at 1 2 4 8 16  (GMX_MCMD_PBG_OUTLAT_CHECK)\n");
        }
        fprintf(mc->fop, "%d %d %.4f  %d %d %d %d %d  %.1f %.1f %.1f %.1f %.1f  %.1f %.1f %.1f %.1f %.1f\n",
                (int)mc->screen_step, ng, dx, nlev[0], nlev[1], nlev[2], nlev[3], nlev[4], tsum[0], tsum[1], tsum[2],
                tsum[3], tsum[4], tint[0], tint[1], tint[2], tint[3], tint[4]);
    }
}

/* Phi = potential on the box of the ions inside it and within MCMD_PBG_MARG = M
 * cells outside: deposited (CIC) on the grid extended by M cells per side and
 * convolved with the kernel on an FFT of 2 (ng + M) per side, which keeps the
 * open boundary for every pair (offsets up to ng - 1 + M).  dep[a] marks the
 * deposited ions; *nout counts the charged ions outside the box, as without
 * the margin.  Returns the number of charged ions beyond the margin. */
static int mcmd_pbg_conv_margin(t_mcionize *mc, const t_mdatoms *md, rvec *x, const real *qq, double *Phi,
                                gmx_bool *dep, int *nout)
{
    const int    ng = mc->pbg_n, M = MCMD_PBG_MARG, m = 2 * (ng + M);
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    const double tc0 = mc->tdet ? mcmd_pb_wallms() : 0;
    int          a, d, i, j, k, i0[3], nfar = 0;
    double       f[3], tf;
    size_t       c;

    for (i = 0; i < m; i++)
    {
        memset(mc->pbg_buf + (size_t)i * m * m, 0, sizeof(float) * (size_t)m * m);
    }
    for (a = 0; a < md->nr; a++)
    {
        gmx_bool ok = TRUE;
        int      u, v, w;

        dep[a] = FALSE;
        if (qq[a] == 0)
        {
            continue;
        }
        if (!mcmd_pbg_cic(mc, x[a], i0, f))
        {
            (*nout)++;
        }
        for (d = 0; d < 3; d++)
        {
            double r = (x[a][d] - mc->pbg_lo[d]) / mc->pbg_dx;

            i0[d] = (int)floor(r);
            f[d]  = r - i0[d];
            ok    = ok && i0[d] >= -M && i0[d] <= ng - 2 + M;
        }
        if (!ok)
        {
            nfar++;
            continue;
        }
        dep[a] = TRUE;
        for (u = 0; u < 2; u++)
        {
            for (v = 0; v < 2; v++)
            {
                for (w = 0; w < 2; w++)
                {
                    mc->pbg_buf[((size_t)(i0[0] + M + u) * m + i0[1] + M + v) * m + i0[2] + M + w] +=
                        (float)(qq[a] * (u ? f[0] : 1 - f[0]) * (v ? f[1] : 1 - f[1]) * (w ? f[2] : 1 - f[2]));
                }
            }
        }
    }
    tf = mc->tdet ? mcmd_pb_wallms() : 0;
    fftwf_execute(mc->pbg_fwdE);
    if (mc->tdet) { mc->td[MCMD_TD_FFT] += mcmd_pb_wallms() - tf; }
    for (c = 0; c < nh; c++)
    {
        mc->pbg_hat[c][0] *= mc->pbg_KhE[c];
        mc->pbg_hat[c][1] *= mc->pbg_KhE[c];
    }
    tf = mc->tdet ? mcmd_pb_wallms() : 0;
    fftwf_execute(mc->pbg_bwdE);
    if (mc->tdet) { mc->td[MCMD_TD_FFT] += mcmd_pb_wallms() - tf; }
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                Phi[MCMD_PBG_IDX(i, j, k)] = mc->pbg_buf[((size_t)(i + M) * m + j + M) * m + k + M];
            }
        }
    }
    mc->pbg_nconv++;
    if (mc->tdet) { mc->td[MCMD_TD_CONV] += mcmd_pb_wallms() - tc0; }
    return nfar;
}

static int mcmd_pbg_ionpot_q(t_mcionize *mc, const t_mdatoms *md, rvec *x, const real *qq, double *Phi)
{
    const int    ng = mc->pbg_n;
    const size_t G  = (size_t)ng * ng * ng;
    int          nout = 0, nfar;
    gmx_bool    *dep = NULL;

    snew(dep, md->nr);
    nfar = mcmd_pbg_conv_margin(mc, md, x, qq, Phi, dep, &nout);
    if (nfar > 0)
    {
        const double to0 = mc->tdet ? mcmd_pb_wallms() : 0;

        if (mc->pbg_outlat == 2)
        {
            /* check: sub-lattices against the exact sum, both on zero */
            double *ex, *ap, emax = 0, e2 = 0, pmax = 0;
            size_t  c;

            snew(ex, G); snew(ap, G);
            mcmd_pbg_outpot(mc, md, x, qq, ex, 0, dep);
            mcmd_pbg_outpot(mc, md, x, qq, ap, 1, dep);
            for (c = 0; c < G; c++)
            {
                double e = fabs(ap[c] - ex[c]);

                emax  = (e > emax) ? e : emax;
                e2   += e * e;
                pmax  = (fabs(ex[c]) > pmax) ? fabs(ex[c]) : pmax;
                Phi[c] += ap[c];
            }
            if (mc->fol == NULL)
            {
                mc->fol = mcmd_open_output("pbg_outlat.txt", "w");
                fprintf(mc->fol, "# per ion-potential call: ions_outside box[nm] dx[nm] max|Phi_out|[eV] "
                        "max_err[eV] rms_err[eV]  (stride dx <= distance / %g)\n", MCMD_PBG_OUTLAT_D);
            }
            fprintf(mc->fol, "%d %.4f %.5f %.4e %.3e %.3e\n", nout, 0.5 * (ng - 1) * mc->pbg_dx, mc->pbg_dx, pmax, emax,
                    sqrt(e2 / G));
            fflush(mc->fol);
            sfree(ex); sfree(ap);
        }
        else
        {
            mcmd_pbg_outpot(mc, md, x, qq, Phi, mc->pbg_outlat, dep);
        }
        if (mc->tdet) { mc->td[MCMD_TD_OUTPOT] += mcmd_pb_wallms() - to0; }
    }
    sfree(dep);
    return nout;
}

/* F(n) and V = -Phi + K n */
static double mcmd_pbg_free(t_mcionize *mc, const double *Phi, double T, const double *n, double *V)
{
    const int    ng = mc->pbg_n;
    const size_t G  = (size_t)ng * ng * ng;
    const double dV = mc->pbg_dx * mc->pbg_dx * mc->pbg_dx;
    double       F  = 0;
    size_t       c;

    mcmd_pbg_conv(mc, n, mc->pbg_Kn);
#pragma omp parallel for reduction(+:F)
    for (c = 0; c < G; c++)
    {
        V[c] = -Phi[c] + mc->pbg_Kn[c];
        F   += n[c] * (-Phi[c] + 0.5 * mc->pbg_Kn[c] + T * (log(n[c] / dV) - 1.0));
    }
    return F;
}

/* Electrostatic energy of the electrons n with the ions (Phi) and themselves. */
static double mcmd_pbg_energy(t_mcionize *mc, const double *Phi, const double *n, double sc)
{
    const int    ng = mc->pbg_n;
    const size_t G  = (size_t)ng * ng * ng;
    double       U  = 0;
    size_t       c;

    mcmd_pbg_conv(mc, n, mc->pbg_Kn);
#pragma omp parallel for reduction(+:U)
    for (c = 0; c < G; c++)
    {
        U += sc * n[c] * (-Phi[c] + 0.5 * sc * mc->pbg_Kn[c]);
    }
    return U;
}

/* The constrained Newton step in one solve: H d + mu 1 = -g, sum(d) = 0,
 * H = K + T diag(1/n), by projected Jacobi-preconditioned CG.  The
 * preconditioned residual z = M^-1 (r - lambda 1) with lambda such that
 * sum(z) = 0, so all directions (and d) keep the electron count.  Stops when
 * the projected residual r - lambda 1 is down by 1e-4 (squared 1e-8, as the
 * two-solve CG). */
static int mcmd_pbg_pcg(t_mcionize *mc, double T, const double *n, const double *g, double *dsol)
{
    const int    ng = mc->pbg_n;
    const size_t G  = (size_t)ng * ng * ng;
    double      *p = mc->pbg_cg, *q = p + G, *r = q + G, *z = r + G;
    double       sm = 0, rz = 0, rr0 = 0, lam, sr;
    size_t       c;
    int          it;

#pragma omp parallel for reduction(+:sm)
    for (c = 0; c < G; c++) { sm += 1.0 / (mc->pbg_K0 + T / n[c]); }
#pragma omp parallel for
    for (c = 0; c < G; c++) { dsol[c] = 0; r[c] = -g[c]; }

    /* z = P M^-1 r; returns r.z, and the squared projected residual in *rrp */
#define MCMD_PCG_PRECOND(rrp)                                                      \
    {                                                                              \
        double sz = 0, rrl = 0, srl = 0;                                           \
        _Pragma("omp parallel for reduction(+:sz)")                                \
        for (c = 0; c < G; c++) { sz += r[c] / (mc->pbg_K0 + T / n[c]); }          \
        lam = sz / sm;                                                             \
        _Pragma("omp parallel for reduction(+:srl, rrl)")                          \
        for (c = 0; c < G; c++)                                                    \
        {                                                                          \
            const double rp = r[c] - lam;                                          \
                                                                                   \
            z[c] = rp / (mc->pbg_K0 + T / n[c]);                                   \
            srl += r[c] * z[c];                                                    \
            rrl += rp * rp;                                                        \
        }                                                                          \
        sr     = srl;                                                              \
        *(rrp) = rrl;                                                              \
    }

    MCMD_PCG_PRECOND(&rr0);
    rz = sr;
#pragma omp parallel for
    for (c = 0; c < G; c++) { p[c] = z[c]; }
    for (it = 0; it < 500 && rr0 > 0; it++)
    {
        double pq = 0, a, rzn, rr;

        mcmd_pbg_conv(mc, p, q);
#pragma omp parallel for reduction(+:pq)
        for (c = 0; c < G; c++)
        {
            q[c] += T * p[c] / n[c];
            pq   += p[c] * q[c];
        }
        a = rz / pq;
#pragma omp parallel for
        for (c = 0; c < G; c++)
        {
            dsol[c] += a * p[c];
            r[c]    -= a * q[c];
        }
        MCMD_PCG_PRECOND(&rr);
        rzn = sr;
        if (rr < 1e-8 * rr0)
        {
            it++;
            break;
        }
#pragma omp parallel for
        for (c = 0; c < G; c++) { p[c] = z[c] + rzn / rz * p[c]; }
        rz = rzn;
    }
#undef MCMD_PCG_PRECOND
    return it;
}

/* The grid limit of this sample (master; once, before the arrays are made):
 * the cells the target spacing needs at the largest box, at R0.  The arrays
 * are made just before the first rebuild, which takes R0 = the radius now, so
 * R0 is fixed here; the spacing grows with R / R0 and the box is at most
 * pbg_box R, so the rebuilds never need more cells than this. */
static void mcmd_pbg_setcap(t_mcionize *mc)
{
    int ng;

    if (mc->pbg_capset)
    {
        return;
    }
    if (mc->pbg_R0 <= 0) { mc->pbg_R0 = mc->radius_now; }
    ng = (int)ceil(2.0 * mc->pbg_box * mc->pbg_R0 / mc->pbg_dxt) + 1;
    ng = 8 * ((ng + 7) / 8);
    ng = (ng < mc->pbg_nmin) ? mc->pbg_nmin : ng;
    if (ng > mc->pbg_nglimit)
    {
        fprintf(stderr, "\nmcionize: WARNING the PB grid would need %d cells per side at %g R; limited to %d "
                "(mcmd-pb-grid-limit), so the spacing is coarser than s_e / %g late in the run\n",
                ng, mc->pbg_box, mc->pbg_nglimit, mc->pbg_se / mc->pbg_dxt);
        ng = mc->pbg_nglimit;
    }
    mc->pbg_nmax   = ng;
    mc->pbg_ncap   = ng;
    mc->pbg_n      = ng;
    mc->pbg_capset = 1;
    fprintf(stderr, "\nmcionize: PB grid up to %d cells per side (R0 %.4g nm, spacing %.3g nm at R0, box up to %g R)\n",
            ng, mc->pbg_R0, mc->pbg_dxt, mc->pbg_box);
}

/* FFT plans for the grid size and the kernels for its spacing (master) */
static void mcmd_pbg_plans(t_mcionize *mc)
{
    if (mc->pbg_nplan != mc->pbg_n)
    {
        const int m = 2 * mc->pbg_n;

        if (mc->pbg_nplan > 0)
        {
            fftwf_destroy_plan(mc->pbg_fwd);
            fftwf_destroy_plan(mc->pbg_bwd);
            fftwf_destroy_plan(mc->pbg_fwdE);
            fftwf_destroy_plan(mc->pbg_bwdE);
        }
        mc->pbg_fwd   = fftwf_plan_dft_r2c_3d(m, m, m, mc->pbg_buf, mc->pbg_hat, FFTW_ESTIMATE);
        mc->pbg_bwd   = fftwf_plan_dft_c2r_3d(m, m, m, mc->pbg_hat, mc->pbg_buf, FFTW_ESTIMATE);
        {
            const int mE = 2 * (mc->pbg_n + MCMD_PBG_MARG);

            mc->pbg_fwdE = fftwf_plan_dft_r2c_3d(mE, mE, mE, mc->pbg_buf, mc->pbg_hat, FFTW_ESTIMATE);
            mc->pbg_bwdE = fftwf_plan_dft_c2r_3d(mE, mE, mE, mc->pbg_hat, mc->pbg_buf, FFTW_ESTIMATE);
        }
        mc->pbg_nplan = mc->pbg_n;
    }
    mcmd_pbg_kernel(mc);
    mcmd_pbg_kernel_margin(mc);
}

/* The grid arrays and FFT plans (master; once). */
static void mcmd_pbg_alloc(t_mcionize *mc, int n)
{
    int    ng, m;
    size_t G;

    if (mc->pbg_ne != NULL)
    {
        return;
    }
    mcmd_pbg_setcap(mc);
    ng = mc->pbg_ncap;
    m  = 2 * ng;
    G  = (size_t)ng * ng * ng;
    snew(mc->pbg_ne, G); snew(mc->pbg_Phi, G); snew(mc->pbg_V, G); snew(mc->pbg_g, G); snew(mc->pbg_d, G);
    snew(mc->pbg_t, G); snew(mc->pbg_Kn, G);
    snew(mc->pbg_nt, G); snew(mc->pbg_Vt, G); snew(mc->pbg_cg, 6 * G);
    {
        /* the FFT buffers also hold the ion potential's FFT with the margin */
        const int mE = 2 * (ng + MCMD_PBG_MARG);

        mc->pbg_buf = fftwf_malloc(sizeof(float) * (size_t)mE * mE * mE);
        mc->pbg_hat = fftwf_malloc(sizeof(fftwf_complex) * (size_t)mE * mE * (mE / 2 + 1));
        snew(mc->pbg_KhE, (size_t)mE * mE * (mE / 2 + 1));
    }
    snew(mc->pbg_Kh, (size_t)m * m * (m / 2 + 1));
#ifdef MCMD_FFTW_THREADS
    /* the other ranks wait at the broadcast while the master solves, so
     * the FFTs can use their cores (GMX_MCMD_PBG_FFT_THREADS, default the
     * number of ranks) */
    fftwf_init_threads();
    fftwf_plan_with_nthreads(mc->pbg_nthreads > 0 ? mc->pbg_nthreads : 1);
    fprintf(stderr, "\nmcionize: PB grid FFTs on %d threads\n", mc->pbg_nthreads > 0 ? mc->pbg_nthreads : 1);
#endif
#ifdef _OPENMP
    /* the grid loops (padding, CG, Newton) on as many OpenMP threads, in the
     * master's parallel regions */
    omp_set_num_threads(mc->pbg_nthreads > 0 ? mc->pbg_nthreads : 1);
#endif
    mc->pbg_nplan = 0;            /* plans made by the first rebuild, for its size */
    if (mc->pb_ne == NULL)    { snew(mc->pb_ne, n); }
    if (mc->pbg_rst_grid)
    {
        /* a continuation: the grid of the checkpoint (size, spacing, origin,
         * restored by mcmd_state_apply), so that the first solve's energy is
         * on the same smoothing as before */
        mcmd_pbg_plans(mc);
        mc->pbg_have     = TRUE;
        mc->pbg_rst_grid = FALSE;
    }
}

/* (Re)build the grid about the current centre and radius; the last solution
 * is carried over by trilinear interpolation. */
static void mcmd_pbg_grid(t_mcionize *mc)
{
    const int    ng_old = mc->pbg_n;
    const double half   = mc->pbg_hbox * mc->radius_now;
    double       lo_old[3], dx_old = mc->pbg_dx, dV, dV_old = dx_old * dx_old * dx_old;
    int          ng, i, j, k, d;
    size_t       G;

    {
        /* adaptive size: spacing near the target, which grows with the
         * sample (R / R0), a multiple of 8 */
        double dxt;

        if (mc->pbg_R0 <= 0) { mc->pbg_R0 = mc->radius_now; }
        dxt = mc->pbg_dxt * ((mc->radius_now > mc->pbg_R0) ? mc->radius_now / mc->pbg_R0 : 1.0);
        ng  = (int)ceil(2.0 * half / dxt) + 1;
        ng = 8 * ((ng + 7) / 8);
        ng = (ng < mc->pbg_nmin) ? mc->pbg_nmin : ((ng > mc->pbg_nmax) ? mc->pbg_nmax : ng);
    }
    G = (size_t)ng * ng * ng;
    for (d = 0; d < 3; d++) { lo_old[d] = mc->pbg_lo[d]; }
    mc->pbg_dx = 2.0 * half / (ng - 1);
    dV         = mc->pbg_dx * mc->pbg_dx * mc->pbg_dx;
    for (d = 0; d < 3; d++) { mc->pbg_lo[d] = mc->com[d] - half; }
    mc->pbg_R = mc->radius_now;
    for (d = 0; d < 3; d++) { mc->pbg_c[d] = mc->com[d]; }
    if (mc->pbg_have)
    {
        /* density of the old grid at the new cell centres */
        for (i = 0; i < ng; i++)
        {
            for (j = 0; j < ng; j++)
            {
                for (k = 0; k < ng; k++)
                {
                    double u[3] = {(mc->pbg_lo[0] + i * mc->pbg_dx - lo_old[0]) / dx_old,
                                   (mc->pbg_lo[1] + j * mc->pbg_dx - lo_old[1]) / dx_old,
                                   (mc->pbg_lo[2] + k * mc->pbg_dx - lo_old[2]) / dx_old};
                    int    i0[3], ok = 1, a, b, c;
                    double f[3], s = 0;

                    for (d = 0; d < 3; d++)
                    {
                        i0[d] = (int)floor(u[d]);
                        f[d]  = u[d] - i0[d];
                        if (i0[d] < 0 || i0[d] >= ng_old - 1) { ok = 0; }
                    }
                    if (ok)
                    {
                        for (a = 0; a < 2; a++)
                        {
                            for (b = 0; b < 2; b++)
                            {
                                for (c = 0; c < 2; c++)
                                {
                                    /* the old density in the old layout */
                                    s += (a ? f[0] : 1 - f[0]) * (b ? f[1] : 1 - f[1]) * (c ? f[2] : 1 - f[2]) *
                                         mc->pbg_ne[(((size_t)(i0[0] + a)) * ng_old + (i0[1] + b)) * ng_old + (i0[2] + c)];
                                }
                            }
                        }
                    }
                    mc->pbg_t[MCMD_PBG_IDX(i, j, k)] = s / dV_old * dV;
                }
            }
        }
        memcpy(mc->pbg_ne, mc->pbg_t, sizeof(double) * G);
    }
    mc->pbg_n = ng;
    mcmd_pbg_plans(mc);
    mc->pbg_have = TRUE;
    mc->pbg_nregrid++;
}

/* Adaptive box (mcmd-pb-box-cut): the half-width / R that leaves a
 * fraction pbg_cut of the electrons outside.  Inside the current box the
 * density is the grid's; beyond it, up to the upper limit pbg_box, the
 * Boltzmann density in the field of the net charge Q at the centre,
 * A exp(k Q / (T r)), A fitted to the outermost layer of cells.  *tail: the
 * electrons so estimated outside the current box; *hsafe: the half-width / R
 * that leaves only MCMD_PBG_SCUT outside (the size a shrink may go to). */
static double mcmd_pbg_boxwant(t_mcionize *mc, double T, double *tail, double *hsafe)
{
    const int    ng = mc->pbg_n;
    const double R = mc->radius_now, w = 0.5 * mc->pbg_dx, dV = mc->pbg_dx * mc->pbg_dx * mc->pbg_dx;
    const double hmax = mc->pbg_box * R, half = 0.5 * (ng - 1) * mc->pbg_dx;
    const double kQT = MCMD_COULOMB_EV_NM * mc->escaped_photoelectrons / T;
    const int    nb = (int)(hmax / w) + ng + 4;
    double      *h, tot = 0, sn = 0, sw = 0, A, cum, want;
    int          i, j, k, b, bs;

    snew(h, nb);
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                const double v = mc->pbg_ne[MCMD_PBG_IDX(i, j, k)];
                const int    ix[3] = {i, j, k};
                double       dm = 0, r2 = 0;
                int          e = ng, d;

                for (d = 0; d < 3; d++)
                {
                    double u = mc->pbg_lo[d] + ix[d] * mc->pbg_dx - mc->com[d];

                    r2 += u * u;
                    u   = fabs(u);
                    dm  = (u > dm) ? u : dm;
                    e   = (ix[d] < e) ? ix[d] : e;
                    e   = (ng - 1 - ix[d] < e) ? ng - 1 - ix[d] : e;
                }
                b = (int)(dm / w);
                h[b < nb ? b : nb - 1] += v;
                tot += v;
                if (e == 0)
                {
                    /* exponent relative to the box half-width: <= 0 here and outside */
                    sn += v;
                    sw += dV * exp(kQT * (1.0 / sqrt(r2) - 1.0 / half));
                }
            }
        }
    }
    *tail  = 0;
    *hsafe = mc->pbg_hbox;
    if (tot <= 0)
    {
        sfree(h);
        return mc->pbg_hbox;
    }
    A = (sw > 0) ? sn / sw : 0;
    if (A > 0 && hmax > 1.001 * half)
    {
        const double hs = 2.0 * mc->pbg_dx, vs = hs * hs * hs;
        const int    M  = (int)ceil(hmax / hs);

        for (i = -M; i <= M; i++)
        {
            for (j = -M; j <= M; j++)
            {
                for (k = -M; k <= M; k++)
                {
                    const double u[3] = {i * hs, j * hs, k * hs};
                    double       dm = 0, v;
                    int          d;

                    for (d = 0; d < 3; d++) { dm = (fabs(u[d]) > dm) ? fabs(u[d]) : dm; }
                    if (dm <= half || dm > hmax)
                    {
                        continue;
                    }
                    v = A * vs * exp(kQT * (1.0 / sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) - 1.0 / half));
                    b = (int)(dm / w);
                    h[b < nb ? b : nb - 1] += v;
                    *tail += v;
                }
            }
        }
    }
    /* the cuts count the electrons outside against all of them (counting
     * them against the net charge Q was tested 2026-10-25 and not kept) */
    {
        const double ref = tot + *tail;

        for (i = 0, b = -1, bs = -1, cum = 0; i < nb && bs < 0; i++)
        {
            cum += h[i];
            if (b < 0 && (tot + *tail - cum) / ref <= mc->pbg_cut) { b = i; }
            if ((tot + *tail - cum) / ref <= MCMD_PBG_SCUT) { bs = i; }
        }
    }
    b      = (b < 0) ? nb - 1 : b;
    want   = ((b + 1) * w + mc->pbg_dx) / R;
    *hsafe = (bs < 0) ? mc->pbg_box : ((bs + 1) * w + mc->pbg_dx) / R;
    *hsafe = (*hsafe > mc->pbg_box) ? mc->pbg_box : ((*hsafe < MCMD_PBG_BOXMIN) ? MCMD_PBG_BOXMIN : *hsafe);
    *tail /= tot;
    sfree(h);
    want = (want > mc->pbg_box) ? mc->pbg_box : want;
    return (want < MCMD_PBG_BOXMIN) ? MCMD_PBG_BOXMIN : want;
}

/* Master: new density for N trapped electrons at T_e, and the energy
 * bookkeeping of the interval since the last solve. */
static void mcmd_pbg_solve(t_mcionize *mc, const t_inputrec *ir, const t_mdatoms *md, rvec *x, double t)
{
    const int    n  = md->nr;
    int          ng = mc->pbg_n;
    size_t       G  = (size_t)ng * ng * ng;
    const double N  = mc->num_free_electrons;
    const double Tmin = 1.0;
    double       T  = (mcmd_plasma_Te(mc) > Tmin) ? mcmd_plasma_Te(mc) : Tmin;
    int          nT = 1;
    double       t0 = mcmd_pb_wallms(), Unew, dU = 0, F, resid = 0, sum, dV, Dmean = 0, Dmax = -1e30;
    double       ecore = 0, qcore = 0;
    int          it, cgit = 0, nout, i, j, k, a, d;
    size_t       c;

    mcmd_pbg_alloc(mc, n);
    mc->pbg_nconv = 0;

    /* total energy: the electrostatic steps from the last solve, on the grid
     * the last density lives on (energy split, GMX_MCMD_PB_EDIAG), and U3 =
     * the energy now (atoms moved, charges changed, N scaled) before a rebuild */
    const gmx_bool etot_on = mc->pb_etot_init && mc->pbg_have && mc->pb_nsolve > 0 && !mc->pbg_cold;
    double         U1 = 0, U2 = 0, U3 = 0, U4 = 0;
    gmx_bool       rebuilt = FALSE;

    if (etot_on)
    {
        if (mc->pb_ediag && mc->pbd_q != NULL)
        {
            mcmd_pbg_ionpot_q(mc, md, x, mc->pbd_q, mc->pbg_Phi);
            U1 = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, 1.0);
            mcmd_pbg_ionpot(mc, md, x, mc->pbg_Phi);
            U2 = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, 1.0);
        }
        else
        {
            mcmd_pbg_ionpot(mc, md, x, mc->pbg_Phi);
        }
        U3 = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, N / mc->pb_N);
        U4 = U3;
    }

    /* adaptive box: the half-width / R the last density asks for */
    double   hwant = mc->pbg_hbox, hsafe = mc->pbg_hbox, btail = 0;
    gmx_bool bgrow, bshrink;

    if (mc->pbg_cut > 0 && mc->pbg_have && mc->pb_nsolve > 0 && N > 0 && !mc->pbg_cold)
    {
        hwant = mcmd_pbg_boxwant(mc, T, &btail, &hsafe);
    }
    /* grow at once, 5 % beyond what is wanted (no lag behind a spreading
     * cloud).  A shrink folds the outer electrons into the box and the
     * rebuild's energy shift is booked as numerical: a box that grew and
     * shrank on every solve (600 eV 1e7, up to 36 times in a run) drained
     * ~0.6 MeV from the cloud.  So shrink only to where at most MCMD_PBG_SCUT
     * (<< pbg_cut) is outside, and by at least 10 %: the box then never
     * flips back at the next solve */
    bgrow   = (hwant > mc->pbg_hbox);
    bshrink = (!bgrow && MCMD_PBG_BOXMARG * hsafe < 0.9 * mc->pbg_hbox);

    /* rebuild the grid first (the last density interpolated onto it, same
     * count), so that the energy before and after are both on this grid: a
     * rebuild changes the smoothing, which is not the electrons' energy */
    if (!mc->pbg_have || (!mc->pbg_cold && (mc->radius_now > 1.1 * mc->pbg_R || mc->radius_now < 0.9 * mc->pbg_R ||
        sqrt(sqr(mc->com[0] - mc->pbg_c[0]) + sqr(mc->com[1] - mc->pbg_c[1]) + sqr(mc->com[2] - mc->pbg_c[2])) >
        mc->pbg_dx || bgrow || bshrink)))
    {
        gmx_bool had = mc->pbg_have;

        hwant = bshrink ? MCMD_PBG_BOXMARG * hsafe : (!bgrow ? mc->pbg_hbox : MCMD_PBG_BOXMARG * hwant);
        hwant = (hwant < mc->pbg_box) ? hwant : mc->pbg_box;
        hwant = (hwant > MCMD_PBG_BOXMIN) ? hwant : MCMD_PBG_BOXMIN;
        mc->pbg_hbox = hwant;
        mcmd_pbg_grid(mc);
        ng = mc->pbg_n;
        G  = (size_t)ng * ng * ng;
        if (had && mc->pb_nsolve > 0)
        {
            sum = 0;
            for (c = 0; c < G; c++) { sum += mc->pbg_ne[c]; }
            for (c = 0; c < G; c++) { mc->pbg_ne[c] = (sum > 0) ? mc->pbg_ne[c] * mc->pb_N / sum : mc->pb_N / G; }
        }
        rebuilt = had;
    }
    if (mc->pbg_cut > 0)
    {
        if (mc->fbox == NULL)
        {
            mc->fbox = mcmd_open_output("pbg_box.txt", "w");
            fprintf(mc->fbox, "# t[ps] R[nm] box[R] dx[nm] outside_est(before rebuild) rebuilt N  (mcmd-pb-box-cut %g, "
                    "box %g..%g R, shrink cut %g, margin %g)\n",
                    mc->pbg_cut, MCMD_PBG_BOXMIN, mc->pbg_box, MCMD_PBG_SCUT, MCMD_PBG_BOXMARG);
        }
        fprintf(mc->fbox, "%.5f %.4f %.3f %.5f %.3e %d %d\n", t, mc->radius_now, mc->pbg_hbox, mc->pbg_dx,
                btail / (1.0 + btail), (int)rebuilt, ng);
        fflush(mc->fbox);
    }
    if (etot_on && rebuilt)
    {
        /* the rebuild's change of the energy (smoothing, interpolation) is
         * numerical: E_tot moves with it, E_kin does not */
        mcmd_pbg_ionpot(mc, md, x, mc->pbg_Phi);
        U4           = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, N / mc->pb_N);
        mc->pb_Etot += U4 - U3;
    }
    dV = mc->pbg_dx * mc->pbg_dx * mc->pbg_dx;
    nout = mcmd_pbg_ionpot(mc, md, x, mc->pbg_Phi);

    /* starting density: the last one (scaled), or on the ion charge (also
     * the first solve after a restart, which has no density) */
    if (mc->pb_nsolve == 0 || mc->pbg_cold)
    {
        mc->pbg_cold = FALSE;
        memset(mc->pbg_ne, 0, sizeof(double) * G);
        for (a = 0; a < n; a++)
        {
            int    i0[3], u, v, w;
            double f[3];

            if (md->chargeA[a] <= 0 || !mcmd_pbg_cic(mc, x[a], i0, f)) { continue; }
            for (u = 0; u < 2; u++)
            {
                for (v = 0; v < 2; v++)
                {
                    for (w = 0; w < 2; w++)
                    {
                        mc->pbg_ne[MCMD_PBG_IDX(i0[0] + u, i0[1] + v, i0[2] + w)] +=
                            md->chargeA[a] * (u ? f[0] : 1 - f[0]) * (v ? f[1] : 1 - f[1]) * (w ? f[2] : 1 - f[2]);
                    }
                }
            }
        }
    }
    sum = 0;
    for (c = 0; c < G; c++) { sum += mc->pbg_ne[c]; }
    for (c = 0; c < G; c++)
    {
        mc->pbg_ne[c] = mc->pbg_ne[c] / sum * N * (1.0 - 1e-3) + 1e-3 * N / G;
    }

    /* timing detail: would a guess extrapolated from the last two solves start
     * closer?  (its residual, one extra convolution; behaviour unchanged) */
    double resid_ex = -1, tn0 = mc->tdet ? mcmd_pb_wallms() : 0;

    if (mc->tdet && mc->pbg_nprev_ok && !rebuilt && mc->pb_nsolve > 1)
    {
        const double tx0 = mcmd_pb_wallms();
        double       st = 0, gb = 0, r = 0;

        for (c = 0; c < G; c++)
        {
            mc->pbg_nt[c] = 2.0 * mc->pbg_ne[c] - mc->pbg_nprev[c] * (N / mc->pb_N);
            if (mc->pbg_nt[c] < 1e-300) { mc->pbg_nt[c] = 1e-300; }
            st += mc->pbg_nt[c];
        }
        for (c = 0; c < G; c++) { mc->pbg_nt[c] *= N / st; }
        mcmd_pbg_free(mc, mc->pbg_Phi, T, mc->pbg_nt, mc->pbg_Vt);
        for (c = 0; c < G; c++)
        {
            mc->pbg_g[c] = mc->pbg_Vt[c] + T * log(mc->pbg_nt[c] / dV);
            gb          += mc->pbg_nt[c] * mc->pbg_g[c];
        }
        gb /= N;
        for (c = 0; c < G; c++) { r += mc->pbg_nt[c] * sqr(mc->pbg_g[c] - gb); }
        resid_ex = sqrt(r / N) / T;
        mc->td[MCMD_TD_EXTRAP] += mcmd_pb_wallms() - tx0;
    }

    /* the density at T (Newton, inner); T is then moved until it is the
     * electrons' own temperature: the root of g(T) = E_tot - U(n_T) - 3/2 N T,
     * which decreases with T (outer; bracketed secant, geometric bisection) */
    const gmx_bool tsc = mc->pb_etot_init && mc->pb_nsolve > 0 && N > 0;
    double         resid0 = -1, Tlo = 0, glo = 0, Thi = 0, ghi = 0;
    gmx_bool       have_lo = FALSE, have_hi = FALSE;

    for (nT = 1;; nT++)
    {
        /* Newton */
        F = mcmd_pbg_free(mc, mc->pbg_Phi, T, mc->pbg_ne, mc->pbg_V);
        for (it = 0; it < 60; it++)
        {
            double gbar = 0, alpha, Ft = 0;
            int    ls;
            gmx_bool accept;

#pragma omp parallel for reduction(+:gbar)
            for (c = 0; c < G; c++)
            {
                mc->pbg_g[c] = mc->pbg_V[c] + T * log(mc->pbg_ne[c] / dV);
                gbar        += mc->pbg_ne[c] * mc->pbg_g[c];
            }
            gbar /= N;
            resid = 0;
#pragma omp parallel for reduction(+:resid)
            for (c = 0; c < G; c++) { resid += mc->pbg_ne[c] * sqr(mc->pbg_g[c] - gbar); }
            resid = sqrt(resid / N) / T;
            if (it == 0) { resid0 = resid; }
            if (resid < mc->pbg_ntol)
            {
                break;
            }
#pragma omp parallel for
            for (c = 0; c < G; c++) { mc->pbg_g[c] -= gbar; }
            {
                const double tg0 = mc->tdet ? mcmd_pb_wallms() : 0, tcv0 = mc->td[MCMD_TD_CONV];

                cgit += mcmd_pbg_pcg(mc, T, mc->pbg_ne, mc->pbg_g, mc->pbg_d);
                if (mc->tdet)
                {
                    mc->td[MCMD_TD_CG]     += mcmd_pb_wallms() - tg0;
                    mc->td[MCMD_TD_CGCONV] += mc->td[MCMD_TD_CONV] - tcv0;
                }
            }
            /* line search: take the step when F or the residual falls (no
             * convolution: V_t is known), at most 8 halvings.  Until 2026-10-25
             * only a fall of F counted, halving up to 30 times; near the target
             * the fall of F (~ N T resid^2) is below the single-precision noise
             * of F, so a tenth of the solves stalled after ~400 convolutions. */
            accept = FALSE;
            for (alpha = 1.0, ls = 0; ls < 8; ls++, alpha *= 0.5)
            {
                double st = 0;

#pragma omp parallel for reduction(+:st)
                for (c = 0; c < G; c++)
                {
                    double u = alpha * mc->pbg_d[c] / mc->pbg_ne[c];

                    if (u > 3.0) { u = 3.0; }
                    if (u < -3.0) { u = -3.0; }
                    mc->pbg_nt[c] = mc->pbg_ne[c] * exp(u);
                    if (mc->pbg_nt[c] < 1e-300) { mc->pbg_nt[c] = 1e-300; }
                    st += mc->pbg_nt[c];
                }
#pragma omp parallel for
                for (c = 0; c < G; c++) { mc->pbg_nt[c] *= N / st; }
                Ft = mcmd_pbg_free(mc, mc->pbg_Phi, T, mc->pbg_nt, mc->pbg_Vt);
                if (Ft <= F)
                {
                    accept = TRUE;
                    break;
                }
                {
                    double gt = 0, rt = 0;

#pragma omp parallel for reduction(+:gt)
                    for (c = 0; c < G; c++)
                    {
                        mc->pbg_g[c] = mc->pbg_Vt[c] + T * log(mc->pbg_nt[c] / dV);
                        gt          += mc->pbg_nt[c] * mc->pbg_g[c];
                    }
                    gt /= N;
#pragma omp parallel for reduction(+:rt)
                    for (c = 0; c < G; c++) { rt += mc->pbg_nt[c] * sqr(mc->pbg_g[c] - gt); }
                    if (sqrt(rt / N) / T < resid)
                    {
                        accept = TRUE;
                        break;
                    }
                }
            }
            if (!accept)
            {
                break;      /* no descent left at single-precision convolutions */
            }
            F = Ft;
            memcpy(mc->pbg_ne, mc->pbg_nt, sizeof(double) * G);
            memcpy(mc->pbg_V, mc->pbg_Vt, sizeof(double) * G);
        }

        if (!tsc || nT >= 12)
        {
            break;
        }
        {
            const double U  = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, 1.0);
            const double g  = mc->pb_Etot - U - 1.5 * N * T;
            const double Tt = (mc->pb_Etot - U) / (1.5 * N);
            double       Tn;

            if (fabs(g) <= 1e-2 * 1.5 * N * T || (T <= Tmin && g < 0))
            {
                break;
            }
            if (g > 0) { Tlo = T; glo = g; have_lo = TRUE; }
            else       { Thi = T; ghi = g; have_hi = TRUE; }
            if (have_lo && have_hi)
            {
                Tn = Tlo + glo * (Thi - Tlo) / (glo - ghi);
                if (!(Tn > Tlo + 0.05 * (Thi - Tlo) && Tn < Thi - 0.05 * (Thi - Tlo)))
                {
                    Tn = sqrt(Tlo * Thi);
                }
            }
            else if (g > 0)
            {
                Tn = (Tt < 4.0 * T) ? Tt : 4.0 * T;
            }
            else
            {
                Tn = (Tt > 0.25 * T) ? Tt : 0.25 * T;
                if (Tn < Tmin) { Tn = Tmin; }
            }
            T = Tn;
        }
    }

    if (mc->tdet)
    {
        mc->td[MCMD_TD_NEWTON] += mcmd_pb_wallms() - tn0;
        mc->td_nsolve          += 1;
        mc->td_nskip           += (resid0 >= 0 && resid0 < mc->pbg_ntol);
        if (resid0 >= 0 && resid_ex >= 0) { mc->td_newton0 += resid0; mc->td_newtonex += resid_ex; }
        if (mc->ftd == NULL)
        {
            mc->ftd = mcmd_open_output("pbg_timing.txt", "w");
            fprintf(mc->ftd, "# t[ps] newton cg convolutions conv_rhs2_total t_fft_total[ms] t_conv_total[ms] "
                    "t_cg_total[ms] t_cgconv_total[ms] t_newton_total[ms] resid0 resid_extrapolated rebuilt\n");
        }
        fprintf(mc->ftd, "%.5f %d %d %d %.0f %.0f %.0f %.0f %.0f %.0f %.3e %.3e %d\n", t, it, cgit, mc->pbg_nconv,
                mc->td_rhs2, mc->td[MCMD_TD_FFT], mc->td[MCMD_TD_CONV], mc->td[MCMD_TD_CG], mc->td[MCMD_TD_CGCONV],
                mc->td[MCMD_TD_NEWTON], resid0, resid_ex, rebuilt);
        fflush(mc->ftd);
        if (mc->pbg_nprev == NULL) { snew(mc->pbg_nprev, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
        memcpy(mc->pbg_nprev, mc->pbg_ne, sizeof(double) * G);
        mc->pbg_nprev_ok = TRUE;
    }

    /* the electron potential on the grid (for the forces and the escape test),
     * per-cell weights and depths (evaporation), centroid, local n_e at the
     * atoms (collisions), diagnostics */
    mcmd_pbg_conv(mc, mc->pbg_ne, mc->pbg_Kn);
    if (mc->pbg_Ve == NULL) { snew(mc->pbg_Ve, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
    for (d = 0; d < 3; d++) { mc->pbg_ecen[d] = 0; }
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                double xc[3] = {mc->pbg_lo[0] + i * mc->pbg_dx, mc->pbg_lo[1] + j * mc->pbg_dx,
                                mc->pbg_lo[2] + k * mc->pbg_dx};

                c               = MCMD_PBG_IDX(i, j, k);
                mc->pbg_Ve[c]   = (float)mc->pbg_Kn[c];
                Dmean          -= mc->pbg_ne[c] * mc->pbg_V[c];
                if (-mc->pbg_V[c] > Dmax && mc->pbg_ne[c] > 1e-6 * N / G) { Dmax = -mc->pbg_V[c]; }
                for (d = 0; d < 3; d++) { mc->pbg_ecen[d] += mc->pbg_ne[c] * xc[d]; }
                if (sqr(xc[0] - mc->com[0]) + sqr(xc[1] - mc->com[1]) + sqr(xc[2] - mc->com[2]) <
                    0.25 * sqr(mc->radius_now))
                {
                    ecore += mc->pbg_ne[c];
                }
            }
        }
    }
    Dmean /= N;
    for (d = 0; d < 3; d++) { mc->pbg_ecen[d] /= N; }
    mc->pb_evn = (int)G;
    mc->pb_evw = mc->pbg_ne;
    mc->pb_evV = mc->pbg_V;
    /* evaporation runs every step: only over the cells that hold electrons
     * (the rest together hold < 3e-5 N) */
    if (mc->pb_evlist == NULL) { snew(mc->pb_evlist, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
    mc->pb_evnl = 0;
    for (c = 0; c < G; c++)
    {
        if (mc->pbg_ne[c] > 1e-10 * N) { mc->pb_evlist[mc->pb_evnl++] = (int)c; }
    }
    for (a = 0; a < n; a++)
    {
        int    i0[3];
        double f[3];

        mc->pb_ne[a] = 0;
        if (mcmd_pbg_cic(mc, x[a], i0, f))
        {
            int u, v, w;

            for (u = 0; u < 2; u++)
            {
                for (v = 0; v < 2; v++)
                {
                    for (w = 0; w < 2; w++)
                    {
                        mc->pb_ne[a] += (u ? f[0] : 1 - f[0]) * (v ? f[1] : 1 - f[1]) * (w ? f[2] : 1 - f[2]) *
                                        mc->pbg_ne[MCMD_PBG_IDX(i0[0] + u, i0[1] + v, i0[2] + w)] / dV;
                    }
                }
            }
        }
        if (sqr(x[a][0] - mc->com[0]) + sqr(x[a][1] - mc->com[1]) + sqr(x[a][2] - mc->com[2]) <
            0.25 * sqr(mc->radius_now))
        {
            qcore += md->chargeA[a];
        }
    }

    /* energy bookkeeping */
    {
        const double K0 = mc->screen_E_trapped;

        Unew = mcmd_pbg_energy(mc, mc->pbg_Phi, mc->pbg_ne, 1.0);
        if (etot_on && mc->pb_ediag && mc->pbd_q != NULL)
        {
            if (mc->fpbd == NULL)
            {
                mc->fpbd = mcmd_open_output("pb_ediag.txt", "w");
                fprintf(mc->fpbd, "# t[ps] N dK_solve  m1_motion m2_charges m3_N m4_regrid(not charged) m5_rearrange  "
                        "booked_W booked_J booked_N  [eV; dK_solve = -(m1+m2+m3+m5)]\n");
            }
            fprintf(mc->fpbd, "%.5f %.0f %.1f  %.1f %.1f %.1f %.1f %.1f  %.1f %.1f %.1f\n", t, N,
                    (mc->pb_Etot - Unew) - K0,
                    (U1 - mc->pbd_Ulast) - mc->pbd_W, (U2 - U1) - mc->pbd_J, (U3 - U2) - mc->pbd_N,
                    U4 - U3, Unew - U4, mc->pbd_W, mc->pbd_J, mc->pbd_N);
            fflush(mc->fpbd);
        }
        if (mc->pb_ediag)
        {
            if (mc->pbd_q == NULL) { snew(mc->pbd_q, n); }
            for (a = 0; a < n; a++) { mc->pbd_q[a] = md->chargeA[a]; }
            mc->pbd_Ulast = Unew;
            mc->pbd_W = mc->pbd_J = mc->pbd_N = 0;
        }
        if (!mc->pb_etot_init)
        {
            mc->pb_Etot      = K0 + Unew;
            mc->pb_etot_init = TRUE;
        }
        mc->screen_E_trapped = mc->pb_Etot - Unew;
        if (mc->screen_E_trapped < 0)
        {
            mc->pb_Eclamp        -= mc->screen_E_trapped;
            mc->screen_E_trapped  = 0;
            mc->pb_Etot           = Unew;
        }
        dU = mc->screen_E_trapped - K0;
    }
    mc->pb_h     = mc->pbg_dx;
    mc->pb_N     = N;
    mc->pb_Dmean = Dmean;

    if (mc->fpb == NULL)
    {
        mc->fpb = mcmd_open_output("pb_log.txt", "w");
        fprintf(mc->fpb, "# t[ps] N Q T_e[eV] dx[nm] ions_outside newton cg resid wall[ms] D_mean[eV] D_max[eV] "
                "W_interval[eV] dU_interval[eV] E_kin[eV] E_clamped[eV] core_e/q(r<R/2) convolutions regrids T_iterations\n");
        if (mc->pbg_dens_out)
        {
            mc->fpbw = mcmd_open_output("pbg_density.bin", "wb");
        }
    }
    fprintf(mc->fpb, "%.5f %.0f %.0f %.2f %.4f %d %d %d %.2e %.0f %.1f %.1f %.1f %.1f %.1f %.1f %.4f %d %d %d\n", t, N,
            mc->escaped_photoelectrons, T, mc->pbg_dx, nout, it, cgit, resid, mcmd_pb_wallms() - t0, Dmean, Dmax,
            mc->pb_W, dU, mc->screen_E_trapped, mc->pb_Eclamp, qcore > 0 ? ecore / qcore : 0.0, mc->pbg_nconv,
            mc->pbg_nregrid, nT);
    fflush(mc->fpb);
    /* GMX_MCMD_PBG_DENSITY=1 (diagnostic, off since 2026-10-27: ~80 % of a run's output): the density every
     * 10th solve: header t[ps] N T dx lo[3] ng, then ng^3 floats n_e [nm^-3], x slowest */
    if (mc->fpbw != NULL && mc->pb_nsolve % 10 == 0)
    {
        float hdr[8] = {(float)t, (float)N, (float)T, (float)mc->pbg_dx, (float)mc->pbg_lo[0], (float)mc->pbg_lo[1],
                        (float)mc->pbg_lo[2], (float)ng};

        fwrite(hdr, sizeof(float), 8, mc->fpbw);
        for (c = 0; c < G; c++)
        {
            float v = (float)(mc->pbg_ne[c] / dV);

            fwrite(&v, sizeof(float), 1, mc->fpbw);
        }
        fflush(mc->fpbw);
    }
    mc->pb_nsolve++;
    mc->tm[MCMD_TM_SOLVE] += mcmd_pb_wallms() - t0;
    mc->pb_W = 0;
}

/* electron potential (eV, repulsive for an electron) at x and its gradient;
 * outside the grid the cloud is a point charge at its centroid */
static double mcmd_pbg_Ve(const t_mcionize *mc, const real *x, double grad[3])
{
    const int ng = mc->pbg_n;
    int       i0[3], d;
    double    f[3];

    if (mcmd_pbg_cic(mc, x, i0, f) && i0[0] >= 1 && i0[1] >= 1 && i0[2] >= 1 &&
        i0[0] < ng - 2 && i0[1] < ng - 2 && i0[2] < ng - 2)
    {
        const float *V = mc->pbg_Ve;
        double       v = 0;
        int          a, b, c;

        for (d = 0; d < 3; d++) { grad[d] = 0; }
        for (a = 0; a < 2; a++)
        {
            for (b = 0; b < 2; b++)
            {
                for (c = 0; c < 2; c++)
                {
                    const int    i = i0[0] + a, j = i0[1] + b, k = i0[2] + c;
                    const double wt = (a ? f[0] : 1 - f[0]) * (b ? f[1] : 1 - f[1]) * (c ? f[2] : 1 - f[2]);

                    v       += wt * V[MCMD_PBG_IDX(i, j, k)];
                    grad[0] += wt * (V[MCMD_PBG_IDX(i + 1, j, k)] - V[MCMD_PBG_IDX(i - 1, j, k)]);
                    grad[1] += wt * (V[MCMD_PBG_IDX(i, j + 1, k)] - V[MCMD_PBG_IDX(i, j - 1, k)]);
                    grad[2] += wt * (V[MCMD_PBG_IDX(i, j, k + 1)] - V[MCMD_PBG_IDX(i, j, k - 1)]);
                }
            }
        }
        for (d = 0; d < 3; d++) { grad[d] /= 2.0 * mc->pbg_dx; }
        return v;
    }
    {
        double r[3], r2 = 0, rr;

        for (d = 0; d < 3; d++) { r[d] = x[d] - mc->pbg_ecen[d]; r2 += r[d] * r[d]; }
        rr = sqrt(r2);
        for (d = 0; d < 3; d++) { grad[d] = -MCMD_COULOMB_EV_NM * mc->pb_N * r[d] / (r2 * rr); }
        return MCMD_COULOMB_EV_NM * mc->pb_N / rr;
    }
}

/* Depth, eV, of an electron at atom k: ions other than k (direct, smoothed
 * over the grid kernel's length, mcmd_pbg_s), the electron cloud from the grid. */
static double mcmd_pbg_depth(const t_mcionize *mc, const t_mdatoms *md, rvec *x, int k)
{
    const double s  = mc->pbg_have ? mcmd_pbg_s(mc) : 0.1;
    const double sc = (mc->pb_ready && mc->pb_N > 0) ? mc->num_free_electrons / mc->pb_N : 0.0;
    double       U  = 0, grad[3];
    int          j;

    for (j = 0; j < md->nr; j++)
    {
        if (j != k && md->chargeA[j] != 0)
        {
            U += MCMD_COULOMB_EV_NM * md->chargeA[j] * mcmd_pb_G(sqrt(mcmd_pb_d2(x[k], x[j])), s);
        }
    }
    if (sc > 0 && mc->pbg_Ve != NULL)
    {
        U -= sc * mcmd_pbg_Ve(mc, x[k], grad);
    }
    return U;
}


/* ---------------------------------------------------------------- *
 * mcmd-pb-long-range: the long-range ion-ion Coulomb force from the grid.
 * GROMACS does the bare force within rcoulomb (rc).  The grid has every
 * pair, smoothed (kernel k G(r, s), s = 2 dx, after cloud-in-cell deposit),
 * and the smoothed pairs within rc are taken out again:
 *   E_a = -grad Phi(x_a) + sum_{b != a, r < rc} k q_b G'(r) (x_a - x_b) / r,
 * so pairs beyond rc come from the grid alone (error ~ erfc(rc / s)).  Every
 * pair within rc is taken out, excluded ones too: GROMACS gives those nothing
 * (or its 1-4 term), as without a cutoff.  The field is computed every
 * PB solve (mcmd-pb-interval) for every home atom and held, times the atom's
 * current charge, in between.  Atoms outside the grid stencil: the bare
 * pairs beyond rc directly.  GMX_MCMD_PBG_LR_CHECK=1 compares, at each refresh,
 * with the exact all-pairs field (master, O(N^2)).
 * ---------------------------------------------------------------- */

/* Value and gradient of the float grid g at x (cloud-in-cell of central
 * differences); FALSE when the stencil leaves the grid. */
static gmx_bool mcmd_pbg_gradf(const t_mcionize *mc, const float *g, const real *x, double *val, double grad[3])
{
    const int ng = mc->pbg_n;
    int       i0[3], a, b, c, d;
    double    f[3];

    if (!(mcmd_pbg_cic(mc, x, i0, f) && i0[0] >= 1 && i0[1] >= 1 && i0[2] >= 1 &&
          i0[0] < ng - 2 && i0[1] < ng - 2 && i0[2] < ng - 2))
    {
        return FALSE;
    }
    *val = 0;
    for (d = 0; d < 3; d++) { grad[d] = 0; }
    for (a = 0; a < 2; a++)
    {
        for (b = 0; b < 2; b++)
        {
            for (c = 0; c < 2; c++)
            {
                const int    i = i0[0] + a, j = i0[1] + b, k = i0[2] + c;
                const double wt = (a ? f[0] : 1 - f[0]) * (b ? f[1] : 1 - f[1]) * (c ? f[2] : 1 - f[2]);

                *val    += wt * g[MCMD_PBG_IDX(i, j, k)];
                grad[0] += wt * (g[MCMD_PBG_IDX(i + 1, j, k)] - g[MCMD_PBG_IDX(i - 1, j, k)]);
                grad[1] += wt * (g[MCMD_PBG_IDX(i, j + 1, k)] - g[MCMD_PBG_IDX(i, j - 1, k)]);
                grad[2] += wt * (g[MCMD_PBG_IDX(i, j, k + 1)] - g[MCMD_PBG_IDX(i, j, k - 1)]);
            }
        }
    }
    for (d = 0; d < 3; d++) { grad[d] /= 2.0 * mc->pbg_dx; }
    return TRUE;
}

/* The grid potential at x of a unit charge deposited at x itself (its own
 * term in Phi, eV per e). */
static double mcmd_pbg_selfpot(const t_mcionize *mc, const real *x)
{
    const double s = mcmd_pbg_s(mc);
    int          i0[3], u, v, d;
    double       f[3], w[8], p[8][3], U = 0;

    if (!mcmd_pbg_cic(mc, x, i0, f))
    {
        return 0;
    }
    for (u = 0; u < 8; u++)
    {
        const int o[3] = {(u >> 2) & 1, (u >> 1) & 1, u & 1};

        w[u] = 1;
        for (d = 0; d < 3; d++)
        {
            w[u]   *= o[d] ? f[d] : 1 - f[d];
            p[u][d] = o[d] * mc->pbg_dx;
        }
    }
    for (u = 0; u < 8; u++)
    {
        for (v = 0; v < 8; v++)
        {
            U += w[u] * w[v] * mcmd_pb_G(sqrt(sqr(p[u][0] - p[v][0]) + sqr(p[u][1] - p[v][1]) + sqr(p[u][2] - p[v][2])), s);
        }
    }
    return MCMD_COULOMB_EV_NM * U;
}

/* Cell list of the charged atoms, cells of at least rc (at most 64 a side). */
static void mcmd_lr_cells(t_mcionize *mc, const t_mdatoms *md, rvec *x)
{
    double hi[3];
    int    a, d, nct;

    if (mc->lr_next == NULL) { snew(mc->lr_next, md->nr); }
    for (d = 0; d < 3; d++) { mc->lr_clo[d] = 1e30; hi[d] = -1e30; }
    for (a = 0; a < md->nr; a++)
    {
        if (md->chargeA[a] == 0) { continue; }
        for (d = 0; d < 3; d++)
        {
            if (x[a][d] < mc->lr_clo[d]) { mc->lr_clo[d] = x[a][d]; }
            if (x[a][d] > hi[d])         { hi[d] = x[a][d]; }
        }
    }
    mc->lr_cs = mc->lr_rc;
    for (d = 0; d < 3; d++)
    {
        if (hi[d] < mc->lr_clo[d]) { mc->lr_clo[d] = hi[d] = 0; }
        if ((hi[d] - mc->lr_clo[d]) / 64.0 > mc->lr_cs) { mc->lr_cs = (hi[d] - mc->lr_clo[d]) / 64.0; }
    }
    for (nct = 1, d = 0; d < 3; d++)
    {
        mc->lr_nc[d] = (int)((hi[d] - mc->lr_clo[d]) / mc->lr_cs) + 1;
        nct         *= mc->lr_nc[d];
    }
    if (nct > mc->lr_headn)
    {
        srenew(mc->lr_head, nct);
        mc->lr_headn = nct;
    }
    for (a = 0; a < nct; a++) { mc->lr_head[a] = -1; }
    for (a = 0; a < md->nr; a++)
    {
        int ci[3], c;

        if (md->chargeA[a] == 0) { continue; }
        for (d = 0; d < 3; d++)
        {
            ci[d] = (int)((x[a][d] - mc->lr_clo[d]) / mc->lr_cs);
            if (ci[d] >= mc->lr_nc[d]) { ci[d] = mc->lr_nc[d] - 1; }
        }
        c              = (ci[0] * mc->lr_nc[1] + ci[1]) * mc->lr_nc[2] + ci[2];
        mc->lr_next[a] = mc->lr_head[c];
        mc->lr_head[c] = a;
    }
}

/* Long-range field E (eV/nm per e) and potential phi (eV per e) at atom a
 * from the other ions; qnear = the charge of the other ions within rc. */
static void mcmd_lr_field(const t_mcionize *mc, const t_mdatoms *md, rvec *x, int a, double E[3], double *phi,
                          double *qnear)
{
    const double rc2 = mc->lr_rc * mc->lr_rc, s = mcmd_pbg_s(mc), k = MCMD_COULOMB_EV_NM;
    double       grad[3], val;
    int          b, d;

    if (mcmd_pbg_gradf(mc, mc->pbg_PhiLf, x[a], &val, grad))
    {
        int ci[3], u, v, w;

        for (d = 0; d < 3; d++) { E[d] = -grad[d]; }
        *phi   = val - md->chargeA[a] * mcmd_pbg_selfpot(mc, x[a]);
        *qnear = 0;
        for (d = 0; d < 3; d++)
        {
            ci[d] = (int)floor((x[a][d] - mc->lr_clo[d]) / mc->lr_cs);
            if (ci[d] < 0)                { ci[d] = 0; }
            if (ci[d] >= mc->lr_nc[d])    { ci[d] = mc->lr_nc[d] - 1; }
        }
        for (u = ci[0] - 1; u <= ci[0] + 1; u++)
        {
            for (v = ci[1] - 1; v <= ci[1] + 1; v++)
            {
                for (w = ci[2] - 1; w <= ci[2] + 1; w++)
                {
                    if (u < 0 || v < 0 || w < 0 || u >= mc->lr_nc[0] || v >= mc->lr_nc[1] || w >= mc->lr_nc[2])
                    {
                        continue;
                    }
                    for (b = mc->lr_head[(u * mc->lr_nc[1] + v) * mc->lr_nc[2] + w]; b >= 0; b = mc->lr_next[b])
                    {
                        double dr[3], d2, r, g;

                        if (b == a) { continue; }
                        for (d = 0; d < 3; d++) { dr[d] = x[a][d] - x[b][d]; }
                        d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                        if (d2 >= rc2) { continue; }
                        r = sqrt(d2);
                        g = mcmd_pb_dG(r, s);
                        for (d = 0; d < 3; d++) { E[d] += k * md->chargeA[b] * g * dr[d] / r; }
                        *phi   -= k * md->chargeA[b] * mcmd_pb_G(r, s);
                        *qnear += md->chargeA[b];
                    }
                }
            }
        }
    }
    else
    {
        /* outside the grid: the bare pairs beyond rc */
        for (d = 0; d < 3; d++) { E[d] = 0; }
        *phi   = 0;
        *qnear = 0;
        for (b = 0; b < md->nr; b++)
        {
            double dr[3], d2, r;

            if (b == a || md->chargeA[b] == 0) { continue; }
            for (d = 0; d < 3; d++) { dr[d] = x[a][d] - x[b][d]; }
            d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (d2 < rc2) { *qnear += md->chargeA[b]; continue; }
            r = sqrt(d2);
            for (d = 0; d < 3; d++) { E[d] += k * md->chargeA[b] * dr[d] / (d2 * r); }
            *phi += k * md->chargeA[b] / r;
        }
    }
}

/* Every rank, every step (grid mode with mcmd-pb-long-range): refresh the field
 * when due, then add q E to the home atoms and book its work. */
static void mcmd_pbg_lr(t_mcionize *mc, const t_inputrec *ir, const t_mdatoms *md, rvec *x, rvec *v, rvec *f,
                        const t_commrec *cr)
{
    size_t       G;
    double       work = 0;
    int          a, d;

    if (!mc->lr_ready || mc->screen_step - mc->lr_step >= mc->pbg_lr_nstep)
    {
        const double t  = mc->screen_step * ir->delta_t + ir->init_t;
        double       hb[6], U[2] = {0, 0}, cr_rms = -1, cr_max = -1;
        size_t       c;

        if (mc->pbg_nthreads == 0)
        {
            const char *e = getenv("GMX_MCMD_PBG_FFT_THREADS");

            mc->pbg_nthreads = (e != NULL && atoi(e) > 0) ? atoi(e) : (PAR(cr) ? cr->nnodes : 1);
        }
        if (MASTER(cr))
        {
            mcmd_pbg_alloc(mc, md->nr);
            if (mc->pbg_PhiLf == NULL) { snew(mc->pbg_PhiLf, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
            if (mc->pbg_PhiL == NULL) { snew(mc->pbg_PhiL, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
            /* the grid follows the sample here only while the PB solve does
             * not (no solve yet, or no electrons): the solve owns the rebuilds
             * and their energy bookkeeping */
            if (!mc->pbg_have ||
                ((mc->pb_nsolve == 0 || mc->num_free_electrons <= 0) &&
                 (mc->radius_now > 1.1 * mc->pbg_R || mc->radius_now < 0.9 * mc->pbg_R ||
                  sqrt(sqr(mc->com[0] - mc->pbg_c[0]) + sqr(mc->com[1] - mc->pbg_c[1]) +
                       sqr(mc->com[2] - mc->pbg_c[2])) > mc->pbg_dx)))
            {
                mcmd_pbg_grid(mc);
            }
            G           = (size_t)mc->pbg_n * mc->pbg_n * mc->pbg_n;
            mc->lr_nout = mcmd_pbg_ionpot(mc, md, x, mc->pbg_PhiL);
            for (c = 0; c < G; c++) { mc->pbg_PhiLf[c] = (float)mc->pbg_PhiL[c]; }
        }
        mcmd_park(mc, cr);
        hb[0] = mc->pbg_dx;
        for (d = 0; d < 3; d++) { hb[1 + d] = mc->pbg_lo[d]; }
        hb[4] = mc->pbg_have; hb[5] = mc->lr_nout;
        if (PAR(cr))
        {
            /* the other ranks learn the grid limit (set on the master from the
             * sample) here, before their copy of the potential is made */
            gmx_bcast(sizeof(int), &mc->pbg_n, cr);
            gmx_bcast(sizeof(int), &mc->pbg_ncap, cr);
            if (mc->pbg_PhiLf == NULL) { snew(mc->pbg_PhiLf, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
        }
        G = (size_t)mc->pbg_n * mc->pbg_n * mc->pbg_n;
        if (PAR(cr))
        {
            gmx_bcast(G * sizeof(float), mc->pbg_PhiLf, cr);
            gmx_bcast(sizeof(hb), hb, cr);
        }
        mc->pbg_dx = hb[0];
        for (d = 0; d < 3; d++) { mc->pbg_lo[d] = hb[1 + d]; }
        mc->pbg_have = (gmx_bool)hb[4]; mc->lr_nout = (int)hb[5];

        mcmd_lr_cells(mc, md, x);
        if (mc->lr_E == NULL) { snew(mc->lr_E, 3 * md->nr); }
        for (a = md->start; a < md->start + md->homenr; a++)
        {
            double phi, qn;

            mcmd_lr_field(mc, md, x, a, mc->lr_E + 3 * a, &phi, &qn);
            U[0] += 0.5 * md->chargeA[a] * phi;
            U[1] += 0.5 * md->chargeA[a] * MCMD_COULOMB_EV_NM * qn / mc->lr_rc;
        }
        if (PAR(cr))
        {
            gmx_sumd(2, U, cr);
        }
        mc->lr_U   = U[0];
        mc->lr_Ush = U[1];

        if (mc->pbg_lr == 2 && MASTER(cr))
        {
            /* exact all-pairs bare field against long-range + the bare pairs
             * within rc (what GROMACS adds, exclusions aside), charged atoms */
            double num = 0, den = 0;

            for (a = 0; a < md->nr; a++)
            {
                double El[3], Ex[3] = {0, 0, 0}, En[3] = {0, 0, 0}, phi, qn, e2 = 0, x2 = 0;
                int    b;

                if (md->chargeA[a] == 0) { continue; }
                mcmd_lr_field(mc, md, x, a, El, &phi, &qn);
                for (b = 0; b < md->nr; b++)
                {
                    double dr[3], d2, r3;

                    if (b == a || md->chargeA[b] == 0) { continue; }
                    for (d = 0; d < 3; d++) { dr[d] = x[a][d] - x[b][d]; }
                    d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                    r3 = d2 * sqrt(d2);
                    for (d = 0; d < 3; d++)
                    {
                        Ex[d] += MCMD_COULOMB_EV_NM * md->chargeA[b] * dr[d] / r3;
                        if (d2 < mc->lr_rc * mc->lr_rc)
                        {
                            En[d] += MCMD_COULOMB_EV_NM * md->chargeA[b] * dr[d] / r3;
                        }
                    }
                }
                for (d = 0; d < 3; d++)
                {
                    e2 += sqr(El[d] + En[d] - Ex[d]);
                    x2 += sqr(Ex[d]);
                }
                num += sqr(md->chargeA[a]) * e2;
                den += sqr(md->chargeA[a]) * x2;
                if (x2 > 0 && sqrt(e2 / x2) > cr_max) { cr_max = sqrt(e2 / x2); }
            }
            cr_rms = (den > 0) ? sqrt(num / den) : 0;
        }
        if (MASTER(cr))
        {
            if (mc->flr == NULL)
            {
                mc->flr = mcmd_open_output("pbg_lr.txt", "w");
                fprintf(mc->flr, "# t[ps] rc[nm] s=2dx[nm] erfc(rc/s) ions_outside U_lr[eV] W_lr[eV]  "
                        "check: rms_rel max_rel (|q|-weighted force error against all pairs; -1 off)  "
                        "U_shift[eV] (k q q / rc over pairs within rc: GROMACS's cutoff shift)\n");
            }
            fprintf(mc->flr, "%.5f %.3f %.4f %.2e %d %.1f %.1f  %.3e %.3e  %.1f\n", t, mc->lr_rc, mcmd_pbg_s(mc),
                    erfc(mc->lr_rc / mcmd_pbg_s(mc)), mc->lr_nout, mc->lr_U, mc->lr_W, cr_rms, cr_max,
                    mc->lr_Ush);
            fflush(mc->flr);
        }
        mc->lr_ready = TRUE;
        mc->lr_step  = mc->screen_step;
    }
    for (a = md->start; a < md->start + md->homenr; a++)
    {
        const double q = md->chargeA[a];

        if (q == 0) { continue; }
        for (d = 0; d < 3; d++)
        {
            const double F = q * mc->lr_E[3 * a + d] * MCMD_EV_TO_KJMOL;     /* kJ/mol/nm */

            f[a][d] += F;
            work    += F * v[a][d] * ir->delta_t;
        }
    }
    if (PAR(cr))
    {
        gmx_sumd(1, &work, cr);
    }
    mc->lr_W += work / MCMD_EV_TO_KJMOL;
}

/* ---------------------------------------------------------------- *
 * mcmd-field-ionization = 1: over-the-barrier (Bethe) ionization in the
 * collective field of the PB solution.  An atom is stripped while
 *     F >= Ip^2 / (4 k Z),
 * Z the charge after the step, Ip from energy_levels.
 *
 * The collective potential of an electron at y is -P_L(y),
 *     P_L = k G_L * (q_ion - n_e),   G_L = erf(r / L) / r,
 * every charge spread over L (L_eff = max(L, s) of the grid kernel,
 * mcmd-field-ionization-smoothing, default 0.5 nm), the atom's own charge left out.
 * Neighbour ions closer than ~L are smeared: their barrier lowering is the
 * charge-transfer channel's.  F = |grad P_L| at the atom, V/nm.
 *
 * On the grid (trapped electrons present): ions deposited by cloud-in-cell
 * and the cloud, convolved with the Fourier kernel of L_eff, gradient by the
 * same stencil as the forces; the atom's own deposit is subtracted exactly
 * (its stencil response to the kernel, from a table of the 13 distances a
 * corner pair can have).  Ions outside the grid by their bare field.  Atoms
 * outside the grid, or no grid: direct sum over the ions with G_L, and the
 * cloud as a point charge at its centroid.
 *
 * The field is computed on the master after every solve (before trapping,
 * every pb_nstep steps) and applied every step in the Monte Carlo part, see
 * mcmd_pbfi_apply.  Every ionization is logged to pbfield_events.txt.
 *
 * Field ionization happens while the sample charges up; afterwards every
 * threshold is far above the field.  So the field is only computed within
 * [peak - 6 sigma, peak + 4 sigma] of the pulse (all lysozyme events at
 * 600 eV, 1e5 and 1e7: -3.8 to +2.6 sigma), and stopped for good once no
 * atom has been within a factor 2 of its threshold for MCMD_PBFI_QUIET
 * updates in a row, but never before peak + sigma.
 * ---------------------------------------------------------------- */

#define MCMD_PBFI_QUIET 5

/* Whether the field is due at t (ps), see above */
static gmx_bool mcmd_pbfi_due(const t_mcionize *mc, double t)
{
    const double u = (t * 1e-12 - mc->t_mean_s) / mc->width_s;

    return mc->pbfi && !mc->pbfi_done && u >= -6.0 && u <= 4.0;
}

/* Bethe threshold (V/nm) of atom i's current state; -1 without data */
static double mcmd_pbfi_thr(const t_mcionize *mc, const t_mdatoms *md, int i)
{
    const t_mcmd_atomdata *ad = &mc->atomdata[mc->elem_of_atom[i]];
    double                 Ip;

    if (mcmd_outermost_occupied(mc->atom_configurations[i]) < 0)
    {
        return -1;
    }
    Ip = mcmd_dict_get(ad->energy_levels, mc->atom_configurations[i]);
    if (Ip == MCMD_NOT_FOUND || Ip <= 0)
    {
        return -1;
    }
    return Ip * Ip / (4.0 * MCMD_COULOMB_EV_NM * (md->chargeA[i] + 1.0));
}

static void mcmd_pbfi_open(t_mcionize *mc)
{
    if (mc->fpbfi_ev == NULL)
    {
        mc->fpbfi_ev = mcmd_open_output("pbfield_events.txt", "w");
        fprintf(mc->fpbfi_ev, "# field ionization in the PB field, L = %g nm\n"
                "# t[ps] atom mass q_after r[nm] F[V/nm] Ip[eV] D_k[eV]\n", mc->pbfi_L);
    }
}

/* out = Kh * a on the grid (Kh a Fourier kernel as pbg_Kh) */
static void mcmd_pbfi_conv(t_mcionize *mc, const float *Kh, const double *a, double *out)
{
    const int    ng = mc->pbg_n, m = 2 * ng;
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    size_t       c;
    int          i, j, k;

    memset(mc->pbg_buf, 0, sizeof(float) * (size_t)m * m * m);
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                mc->pbg_buf[((size_t)i * m + j) * m + k] = (float)a[MCMD_PBG_IDX(i, j, k)];
            }
        }
    }
    fftwf_execute(mc->pbg_fwd);
    for (c = 0; c < nh; c++)
    {
        mc->pbg_hat[c][0] *= Kh[c];
        mc->pbg_hat[c][1] *= Kh[c];
    }
    fftwf_execute(mc->pbg_bwd);
    for (i = 0; i < ng; i++)
    {
        for (j = 0; j < ng; j++)
        {
            for (k = 0; k < ng; k++)
            {
                out[MCMD_PBG_IDX(i, j, k)] = mc->pbg_buf[((size_t)i * m + j) * m + k];
            }
        }
    }
}

/* Fourier kernel k G(r, L) of the current grid into Kh */
static void mcmd_pbfi_kernel(t_mcionize *mc, double L, float *Kh)
{
    const int    ng = mc->pbg_n, m = 2 * ng;
    const size_t nh = (size_t)m * m * (m / 2 + 1);
    size_t       c;
    int          i, j, k;

    for (i = 0; i < m; i++)
    {
        double xi = ((i < m - i) ? i : m - i) * mc->pbg_dx;

        for (j = 0; j < m; j++)
        {
            double yj = ((j < m - j) ? j : m - j) * mc->pbg_dx;

            for (k = 0; k < m; k++)
            {
                double zk = ((k < m - k) ? k : m - k) * mc->pbg_dx;

                mc->pbg_buf[((size_t)i * m + j) * m + k] =
                    (float)(MCMD_COULOMB_EV_NM * mcmd_pb_G(sqrt(xi * xi + yj * yj + zk * zk), L));
            }
        }
    }
    fftwf_execute(mc->pbg_fwd);
    for (c = 0; c < nh; c++)
    {
        Kh[c] = mc->pbg_hat[c][0] / ((double)m * m * m);
    }
}

/* grad P_L at atom a by direct sum: ions with G_L (a itself left out; bare
 * beyond 4 L, where erf is 1 to 1e-8), the cloud of N electrons as a point
 * charge at its centroid. */
static void mcmd_pbfi_direct(const t_mcionize *mc, const t_mdatoms *md, rvec *x, int a, double L, double N,
                             double g[3])
{
    int b, d;

    for (d = 0; d < 3; d++) { g[d] = 0; }
    for (b = 0; b < md->nr; b++)
    {
        double r2, r, dG;

        if (b == a || md->chargeA[b] == 0)
        {
            continue;
        }
        r2 = mcmd_pb_d2(x[a], x[b]);
        if (r2 <= 0)
        {
            continue;
        }
        r  = sqrt(r2);
        dG = (r > 4.0 * L) ? -1.0 / r2 : mcmd_pb_dG(r, L);
        for (d = 0; d < 3; d++)
        {
            g[d] += MCMD_COULOMB_EV_NM * md->chargeA[b] * dG * (x[a][d] - x[b][d]) / r;
        }
    }
    if (N > 0)
    {
        double r2 = 0, r;

        for (d = 0; d < 3; d++) { r2 += sqr(x[a][d] - mc->pbg_ecen[d]); }
        r = sqrt(r2);
        if (r > 0)
        {
            for (d = 0; d < 3; d++)
            {
                g[d] += MCMD_COULOMB_EV_NM * N * (x[a][d] - mc->pbg_ecen[d]) / (r2 * r);
            }
        }
    }
}

/* The field F (V/nm) at every atom with electrons left, into pbfi_F, at t
 * (ps).  Master only; grid: a solution is current (after mcmd_pbg_solve). */
static void mcmd_pbfi_field(t_mcionize *mc, const t_mdatoms *md, rvec *x, gmx_bool grid, double t)
{
    const int n = md->nr;
    int       a, d;
    gmx_bool  near = FALSE;

    if (mc->pbfi_F == NULL)
    {
        snew(mc->pbfi_F, n);
    }
    if (grid)
    {
        const int     ng = mc->pbg_n;
        const size_t  G  = (size_t)ng * ng * ng;
        const double  s  = mcmd_pbg_s(mc), L = (mc->pbfi_L > s) ? mc->pbfi_L : s;
        const double  N  = mc->num_free_electrons;
        const double *P;
        int          *out, nout = 0, i0[3], m;
        double        f[3];
        size_t        c;

        if (mc->pbfi_q == NULL)
        {
            const size_t Gc = (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap;
            const int    m2 = 2 * mc->pbg_ncap;

            snew(mc->pbfi_q, Gc);
            snew(mc->pbfi_V, Gc);
            snew(mc->pbfi_Kh, (size_t)m2 * m2 * (m2 / 2 + 1));
        }
        if (mc->pbfi_kn != ng || mc->pbfi_kdx != mc->pbg_dx || mc->pbfi_kL != L)
        {
            mcmd_pbfi_kernel(mc, L, mc->pbfi_Kh);
            for (m = 0; m <= 12; m++)
            {
                mc->pbfi_Gt[m] = MCMD_COULOMB_EV_NM * mcmd_pb_G(sqrt((double)m) * mc->pbg_dx, L);
            }
            mc->pbfi_kn = ng; mc->pbfi_kdx = mc->pbg_dx; mc->pbfi_kL = L;
        }
        /* source: ions by cloud-in-cell minus the cloud */
        snew(out, n);
        for (c = 0; c < G; c++) { mc->pbfi_q[c] = -mc->pbg_ne[c]; }
        for (a = 0; a < n; a++)
        {
            int u, v, w;

            if (md->chargeA[a] == 0)
            {
                continue;
            }
            if (!mcmd_pbg_cic(mc, x[a], i0, f))
            {
                out[nout++] = a;
                continue;
            }
            for (u = 0; u < 2; u++)
            {
                for (v = 0; v < 2; v++)
                {
                    for (w = 0; w < 2; w++)
                    {
                        mc->pbfi_q[MCMD_PBG_IDX(i0[0] + u, i0[1] + v, i0[2] + w)] +=
                            md->chargeA[a] * (u ? f[0] : 1 - f[0]) * (v ? f[1] : 1 - f[1]) * (w ? f[2] : 1 - f[2]);
                    }
                }
            }
        }
        mcmd_pbfi_conv(mc, mc->pbfi_Kh, mc->pbfi_q, mc->pbfi_V);
        P = mc->pbfi_V;
        for (a = 0; a < n; a++)
        {
            double g[3] = {0, 0, 0};
            int    b, u, v;

            mc->pbfi_F[a] = 0;
            if (mcmd_outermost_occupied(mc->atom_configurations[a]) < 0)
            {
                continue;    /* bare: nothing to strip */
            }
            if (mcmd_pbg_cic(mc, x[a], i0, f) && i0[0] >= 1 && i0[1] >= 1 && i0[2] >= 1 &&
                i0[0] < ng - 2 && i0[1] < ng - 2 && i0[2] < ng - 2)
            {
                double wt[8];
                int    o[8][3];

                for (u = 0; u < 8; u++)
                {
                    o[u][0] = (u >> 2) & 1; o[u][1] = (u >> 1) & 1; o[u][2] = u & 1;
                    wt[u]   = (o[u][0] ? f[0] : 1 - f[0]) * (o[u][1] ? f[1] : 1 - f[1]) *
                              (o[u][2] ? f[2] : 1 - f[2]);
                }
                for (u = 0; u < 8; u++)
                {
                    const int i = i0[0] + o[u][0], j = i0[1] + o[u][1], k = i0[2] + o[u][2];

                    g[0] += wt[u] * (P[MCMD_PBG_IDX(i + 1, j, k)] - P[MCMD_PBG_IDX(i - 1, j, k)]);
                    g[1] += wt[u] * (P[MCMD_PBG_IDX(i, j + 1, k)] - P[MCMD_PBG_IDX(i, j - 1, k)]);
                    g[2] += wt[u] * (P[MCMD_PBG_IDX(i, j, k + 1)] - P[MCMD_PBG_IDX(i, j, k - 1)]);
                }
                /* own deposit: evaluation corner u, deposit corner v */
                if (md->chargeA[a] != 0)
                {
                    for (u = 0; u < 8; u++)
                    {
                        for (v = 0; v < 8; v++)
                        {
                            for (d = 0; d < 3; d++)
                            {
                                int mp = 0, mm = 0, e;

                                for (e = 0; e < 3; e++)
                                {
                                    const int off = o[u][e] - o[v][e];

                                    mp += (off + (e == d)) * (off + (e == d));
                                    mm += (off - (e == d)) * (off - (e == d));
                                }
                                g[d] -= wt[u] * wt[v] * md->chargeA[a] * (mc->pbfi_Gt[mp] - mc->pbfi_Gt[mm]);
                            }
                        }
                    }
                }
                for (d = 0; d < 3; d++) { g[d] /= 2.0 * mc->pbg_dx; }
                /* ions outside the grid, bare */
                for (b = 0; b < nout; b++)
                {
                    const int    ob = out[b];
                    const double r  = sqrt(mcmd_pb_d2(x[a], x[ob]));

                    for (d = 0; d < 3; d++)
                    {
                        g[d] -= MCMD_COULOMB_EV_NM * md->chargeA[ob] * (x[a][d] - x[ob][d]) / (r * r * r);
                    }
                }
            }
            else
            {
                mcmd_pbfi_direct(mc, md, x, a, L, N, g);
            }
            mc->pbfi_F[a] = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        }
        sfree(out);
    }
    else
    {
        for (a = 0; a < n; a++)
        {
            double g[3];

            mc->pbfi_F[a] = 0;
            if (mcmd_outermost_occupied(mc->atom_configurations[a]) < 0)
            {
                continue;
            }
            mcmd_pbfi_direct(mc, md, x, a, mc->pbfi_L, 0, g);
            mc->pbfi_F[a] = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        }
    }
    mc->pbfi_F_ok = TRUE;

    for (a = 0; a < n && !near; a++)
    {
        const double thr = mcmd_pbfi_thr(mc, md, a);

        near = (thr > 0 && mc->pbfi_F[a] >= 0.5 * thr);
    }
    mc->pbfi_quiet = near ? 0 : mc->pbfi_quiet + 1;
    if (mc->pbfi_quiet >= MCMD_PBFI_QUIET && t * 1e-12 >= mc->t_mean_s + mc->width_s)
    {
        mc->pbfi_done = TRUE;
        mcmd_pbfi_open(mc);
        fprintf(mc->fpbfi_ev, "# %.5f stopped: no atom within a factor 2 of its threshold for %d updates\n", t,
                MCMD_PBFI_QUIET);
        fflush(mc->fpbfi_ev);
    }
}

/* Strip every atom while the field of the last update is over the Bethe
 * threshold of its current state.  Called every step from the Monte Carlo
 * part (master), so the charges reach the other ranks with the rest of the
 * step's changes.
 *
 * Energy: the bound electron (energy -Ip relative to its ion) leaves over the
 * saddle and falls down the collective potential; it joins the cloud with
 * total energy -D_k - Ip, D_k the depth at its atom (own ion left out, as
 * for every birth).  So the cloud pays the binding energy, as in an impact
 * ionization with no kinetic energy brought in; the change in the cloud-ion
 * energy from the new charge is the ledger's (mcmd_ledger). */
static void mcmd_pbfi_apply(t_mcionize *mc, t_mdatoms *md, t_state *state, double t)
{
    int i;

    if (!mc->pbfi_F_ok || !mcmd_pbfi_due(mc, t))
    {
        return;
    }
    mcmd_pbfi_open(mc);
    for (i = 0; i < md->nr; i++)
    {
        const t_mcmd_atomdata *ad = &mc->atomdata[mc->elem_of_atom[i]];

        for (;;)
        {
            int    shell = mcmd_outermost_occupied(mc->atom_configurations[i]);
            int    next[MCMD_NSHELL], d;
            double Ip, Z, D_k;

            if (shell < 0)
            {
                break;
            }
            Ip = mcmd_dict_get(ad->energy_levels, mc->atom_configurations[i]);
            if (Ip == MCMD_NOT_FOUND || Ip <= 0)
            {
                break;
            }
            Z = md->chargeA[i] + 1.0;
            if (mc->pbfi_F[i] < Ip * Ip / (4.0 * MCMD_COULOMB_EV_NM * Z))
            {
                break;
            }
            for (d = 0; d < MCMD_NSHELL; d++) { next[d] = mc->atom_configurations[i][d]; }
            next[shell]--;
            if (!(next[0] == 0 && next[1] == 0 && next[2] == 0) &&
                mcmd_state_index_get(&ad->rate_index, next) == -1)
            {
                break;
            }
            D_k = mcmd_pb_depth(mc, md, state->x, i);
            mc->atom_configurations[i][shell]--;
            md->chargeA[i] += 1;
            mcmd_refresh_atom(mc, md, i);
            mc->num_free_electrons               += 1;
            mc->thermalized_free_electron_energy += MCMD_TRAPPED_ELECTRON_EV;
            mc->pbfi_nion                        += 1;
            mc->pbfi_Eb                          += Ip;
            mc->pb_Etot -= D_k + Ip;
            mc->pbd_N   -= D_k;
            fprintf(mc->fpbfi_ev, "%.5f %d %.0f %.0f %.4f %.1f %.2f %.2f\n", t, i + 1, md->massT[i], md->chargeA[i],
                    sqrt(sqr(state->x[i][XX] - mc->com[0]) + sqr(state->x[i][YY] - mc->com[1]) +
                         sqr(state->x[i][ZZ] - mc->com[2])), mc->pbfi_F[i], Ip, D_k);
        }
    }
    fflush(mc->fpbfi_ev);
}

/* Every rank: solve when due (master, broadcast the electron potential), then
 * q grad(K n) on the home ions; the work comes out of the electrons' energy. */
static void mcmd_pbg_forces(t_mcionize *mc, const t_inputrec *ir, const t_mdatoms *md, rvec *x, rvec *v,
                            rvec *f, const t_commrec *cr)
{
    size_t       G;
    const double N  = mc->num_free_electrons;
    double       work = 0, sc;
    int          a, d;

    if (N <= 0)
    {
        const double t = mc->screen_step * ir->delta_t + ir->init_t;

        if (MASTER(cr) && mc->screen_step % mc->pb_nstep == 0 && mcmd_pbfi_due(mc, t))
        {
            mcmd_pbfi_field(mc, md, x, FALSE, t);
        }
        return;
    }
    if (mc->pbg_nthreads == 0)
    {
        const char *e = getenv("GMX_MCMD_PBG_FFT_THREADS");

        mc->pbg_nthreads = (e != NULL && atoi(e) > 0) ? atoi(e) : (PAR(cr) ? cr->nnodes : 1);
    }
    if (!mc->pb_ready || mc->screen_step - mc->pb_step_solve >= mc->pb_nstep || N > 1.2 * mc->pb_N + 10)
    {
        double hb[14];

        if (MASTER(cr))
        {
            mcmd_pbg_solve(mc, ir, md, x, mc->screen_step * ir->delta_t + ir->init_t);
            const double t = mc->screen_step * ir->delta_t + ir->init_t;

            if (mcmd_pbfi_due(mc, t))
            {
                mcmd_pbfi_field(mc, md, x, TRUE, t);
            }
        }
        mcmd_park(mc, cr);
        if (PAR(cr))
        {
            /* the grid limit is set on the master at its first solve */
            gmx_bcast(sizeof(int), &mc->pbg_n, cr);
            gmx_bcast(sizeof(int), &mc->pbg_ncap, cr);
        }
        if (mc->pbg_Ve == NULL) { snew(mc->pbg_Ve, (size_t)mc->pbg_ncap * mc->pbg_ncap * mc->pbg_ncap); }
        G = (size_t)mc->pbg_n * mc->pbg_n * mc->pbg_n;
        hb[0] = mc->pb_N; hb[1] = mc->pb_Dmean; hb[2] = mc->screen_E_trapped; hb[3] = mc->pbg_dx;
        for (d = 0; d < 3; d++) { hb[4 + d] = mc->pbg_lo[d]; hb[7 + d] = mc->pbg_ecen[d]; }
        hb[10] = mc->pbg_have; hb[11] = mc->pb_h; hb[12] = mc->pb_Etot; hb[13] = mc->pb_etot_init;
        if (PAR(cr))
        {
            gmx_bcast(G * sizeof(float), mc->pbg_Ve, cr);
            gmx_bcast(sizeof(hb), hb, cr);
        }
        mc->pb_N = hb[0]; mc->pb_Dmean = hb[1]; mc->screen_E_trapped = hb[2]; mc->pbg_dx = hb[3];
        for (d = 0; d < 3; d++) { mc->pbg_lo[d] = hb[4 + d]; mc->pbg_ecen[d] = hb[7 + d]; }
        mc->pbg_have = (gmx_bool)hb[10]; mc->pb_h = hb[11]; mc->pb_Etot = hb[12]; mc->pb_etot_init = (gmx_bool)hb[13];
        mc->pb_ready      = TRUE;
        mc->pb_step_solve = mc->screen_step;
    }
    sc = N / mc->pb_N;
    for (a = md->start; a < md->start + md->homenr; a++)
    {
        double grad[3];

        if (md->chargeA[a] == 0)
        {
            continue;
        }
        mcmd_pbg_Ve(mc, x[a], grad);
        for (d = 0; d < 3; d++)
        {
            double F = md->chargeA[a] * sc * grad[d] * MCMD_EV_TO_KJMOL;    /* kJ/mol/nm */

            f[a][d] += F;
            work    += F * v[a][d] * ir->delta_t;
        }
    }
    if (PAR(cr))
    {
        gmx_sumd(1, &work, cr);
    }
    work /= MCMD_EV_TO_KJMOL;
    mc->pressure_work    += work;
    mc->pb_W             += work;
    mc->pb_Etot -= work;
    mc->pbd_W   -= work;
}

static void mcmd_monte_carlo(t_mcionize *mc, const t_inputrec *ir,
                             t_mdatoms *mdatoms, t_state *state, double t, int k0, int k1)
{
    const double dt_md_s = ir->delta_t * 1e-12;
    int          k;

    for (k = k0; k < k1; k++)
    {
        const t_mcmd_atomdata *ad = &mc->atomdata[mcmd_atom_element(mc, mdatoms, k)];
        double                 DT_current = 0.0;
        int                    match[MCMD_NSHELL];

        while (DT_current < dt_md_s)
        {
            t_mcmd_transition *possible_transitions;
            double            *dt_processes;
            double             min_value;
            int                min_value_index = -1;
            int                total_transitions;
            int                match_index;
            int                offset;
            double             net_charge;
            int                i;

            match[0] = mc->atom_configurations[k][0];
            match[1] = mc->atom_configurations[k][1];
            match[2] = mc->atom_configurations[k][2];

            /* Array read into the per-element table built at setup, in place
             * of a scan over every state of the element - 324 of them for
             * sulphur, and one scan per atom per event per step. */
            match_index = mcmd_state_index_get(&ad->rate_index, match);

            if (match_index == -1)
            {
                if (match[0] == 0 && match[1] == 0 && match[2] == 0)
                {
                    /* Fully stripped.  Nothing can decay, but recombination
                     * could still bring an electron back, so the state is
                     * kept rather than treated as an error. */
                    break;
                }
                gmx_fatal(FARGS,
                          "No transition data for atom %d (mass %g) in state "
                          "%d %d %d.",
                          k + 1, mdatoms->massT[k],
                          match[0], match[1], match[2]);
            }

            offset            = ad->rates[match_index].num_transitions;
            total_transitions = offset;

            if (total_transitions > mc->mc_scratch_size)
            {
                /* Cannot happen: the buffers are sized to the largest
                 * transition count of any loaded state. */
                gmx_fatal(FARGS,
                          "Monte Carlo scratch too small (%d needed, %d "
                          "allocated).",
                          total_transitions, mc->mc_scratch_size);
            }
            possible_transitions = mc->mc_transitions;
            dt_processes         = mc->mc_dt;

            /* Radiative and Auger channels, straight from the data file. */
            for (i = 0; i < offset; i++)
            {
                possible_transitions[i].final_state[0] = ad->rates[match_index].final_states[i][0];
                possible_transitions[i].final_state[1] = ad->rates[match_index].final_states[i][1];
                possible_transitions[i].final_state[2] = ad->rates[match_index].final_states[i][2];
                possible_transitions[i].type           = ad->rates[match_index].types[i];

                /* Only photoionization is driven by the pulse; Auger and
                 * fluorescence rates are intrinsic. */
                if (possible_transitions[i].type == 2)
                {
                    possible_transitions[i].rate =
                        mc->intensity * ad->rates[match_index].rates[i];
                }
                else
                {
                    possible_transitions[i].rate = ad->rates[match_index].rates[i];
                }
            }

            mc->electron_density = mc->num_free_electrons / mc->current_volume;

            /* Collisional channels, for atoms inside the electron cloud.
             * pb: at the electron density at the atom instead of the mean
             * N/V (T_e is one temperature): ionization and excitation scale
             * with n_e, three-body recombination with n_e^2; every atom the
             * electrons reach takes part. */
            double pb_f = 1.0;

            if (mc->coll_on && mc->coll_active && mc->pb_ready && mc->pb_ne != NULL)
            {
                pb_f = mc->pb_ne[k] * (mc->num_free_electrons / mc->pb_N) * 1e21 / mc->coll_ne;
            }
            if (mc->coll_on && mc->coll_active && pb_f > 0 &&
                (mc->pb_ready || mcmd_in_cloud(mc, state->x[k])))
            {
                const int     el = mcmd_atom_element(mc, mdatoms, k);
                const int     c0 = mc->cch_start[el][match_index];
                const int     c1 = mc->cch_start[el][match_index + 1];
                int           c;

                for (c = c0; c < c1; c++)
                {
                    const t_mcmd_cchan *ch = &mc->cch[el][c];

                    /* the rates are frozen for the step; never recombine an
                     * electron the plasma no longer has */
                    if (ch->recomb && mc->num_free_electrons < 1.0) { continue; }
                    possible_transitions[total_transitions].final_state[0] = ch->final_state[0];
                    possible_transitions[total_transitions].final_state[1] = ch->final_state[1];
                    possible_transitions[total_transitions].final_state[2] = ch->final_state[2];
                    possible_transitions[total_transitions].rate = ch->rate * (ch->recomb ? pb_f * pb_f : pb_f);
                    possible_transitions[total_transitions].type =
                        ch->excite ? 4 + ch->excite : (ch->recomb ? 4 : 3);
                    mc->mc_cref[total_transitions] = ch;
                    total_transitions++;
                }
            }

            /* First-reaction step: one uniform deviate per channel. */
            min_value = 1.1 * dt_md_s; /* a bit beyond the MD timestep */

            for (i = 0; i < total_transitions; i++)
            {
                double u = (double)rand() / (RAND_MAX + 1.0);

                /* rand() can return 0, so u is drawn from [0,1) inclusive of
                 * zero, and -log(0) is +inf.  A zero deviate formally means an
                 * infinite waiting time -- the same "never fires" outcome as a
                 * zero rate -- so it is folded into that branch rather than
                 * clamped to an epsilon.  Rare per draw (4.7e-10 at
                 * RAND_MAX = 2^31-1) but near-certain over a large system: a
                 * 50000-atom run draws ~2e5 deviates per step and expects
                 * ~3 zeros in 30000 steps.  The isinf() check below stays as
                 * an invariant on the remaining arithmetic. */
                if (possible_transitions[i].rate < 1e-100 || u <= 0.0)
                {
                    /* Rate of zero, or a zero deviate: never fires. */
                    dt_processes[i] = 1e6;
                }
                else
                {
                    dt_processes[i] = -log(u) / possible_transitions[i].rate;
                }

                if (dt_processes[i] < min_value)
                {
                    min_value       = dt_processes[i];
                    min_value_index = i;
                }

                if (isinf(dt_processes[i]))
                {
                    gmx_fatal(FARGS,
                              "Infinite Monte Carlo waiting time for atom %d "
                              "(mass %g) in state %d %d %d, rate %g, "
                              "intensity %g.",
                              k + 1, mdatoms->massT[k],
                              match[0], match[1], match[2],
                              possible_transitions[i].rate, mc->intensity);
                }
            }

            DT_current += min_value;

            if (DT_current > dt_md_s)
            {
                /* The next event falls outside this MD step. */
                break;
            }

            net_charge =
                (double)(match[0] - possible_transitions[min_value_index].final_state[0] +
                         match[1] - possible_transitions[min_value_index].final_state[1] +
                         match[2] - possible_transitions[min_value_index].final_state[2]);

            /* An emitted electron either escapes the sample's space charge or
             * is trapped into the free-electron plasma.  Fluorescence emits
             * none (net_charge 0); a negative net charge is a capture, which
             * takes an electron out of the plasma. */
            if (possible_transitions[min_value_index].type >= 3)
            {
                /* Collisional event.  Ionization: a plasma electron gives
                 * up dE and frees one more, which stays in the plasma.
                 * Three-body recombination: one plasma electron is bound and
                 * the binding energy dE goes to the third electron. */
                const t_mcmd_cchan *ch = mc->mc_cref[min_value_index];
                /* the plasma pays dE for ionization and excitation, and gets
                 * it back from recombination and de-excitation */
                double              dE = (net_charge > 0 || ch->excite == 1) ? -ch->dE : ch->dE;

                /* The electron is freed (or captured) at its atom, so the
                 * cloud's total energy gains (loses) its potential energy
                 * there, -D_k. */
                if (net_charge != 0)
                {
                    double D_k = mcmd_pb_depth(mc, mdatoms, state->x, k);

                    mc->pb_Etot -= net_charge * D_k;
                    mc->pbd_N   -= net_charge * D_k;
                    if (mc->coll_evlog && mc->pb_ready)
                    {
                        /* the solve's own electron potential at the atom (all ions,
                         * own included, smoothed as for D_mean), interpolated
                         * from the cells */
                        double D_site = 0, sc = mc->num_free_electrons / mc->pb_N, rr;
                        int    i0[3];
                        double f[3];

                        if (mcmd_pbg_cic(mc, state->x[k], i0, f))
                        {
                            const int ng = mc->pbg_n;
                            int       u, v, w;

                            for (u = 0; u < 2; u++)
                            {
                                for (v = 0; v < 2; v++)
                                {
                                    for (w = 0; w < 2; w++)
                                    {
                                        D_site -= (u ? f[0] : 1 - f[0]) * (v ? f[1] : 1 - f[1]) *
                                                  (w ? f[2] : 1 - f[2]) *
                                                  mc->pbg_V[MCMD_PBG_IDX(i0[0] + u, i0[1] + v, i0[2] + w)];
                                    }
                                }
                            }
                        }
                        if (mc->fcev == NULL)
                        {
                            mc->fcev = mcmd_open_output("coll_event_log.txt", "w");
                            fprintf(mc->fcev, "# t[ps] atom net_charge(+1 ionization, -1 recombination) q_atom_after "
                                    "D_k[eV] D_site[eV] D_mean[eV] n_e_local[nm^-3] r/R\n");
                        }
                        rr = sqrt(sqr(state->x[k][0] - mc->com[0]) + sqr(state->x[k][1] - mc->com[1]) +
                                  sqr(state->x[k][2] - mc->com[2])) / mc->radius_now;
                        fprintf(mc->fcev, "%.5f %d %+.0f %.0f %.2f %.2f %.2f %.4g %.3f\n", t, k + 1, net_charge,
                                (double)mdatoms->chargeA[k] + net_charge, D_k, D_site, mc->pb_Dmean,
                                (mc->pb_ne != NULL) ? mc->pb_ne[k] * sc : 0.0, rr);
                    }
                }
                mc->num_free_electrons               += net_charge;
                mc->thermalized_free_electron_energy += net_charge * MCMD_TRAPPED_ELECTRON_EV;
                mc->screen_E_trapped                 += dE;
                mc->coll_energy                      -= dE;
                mc->pb_Etot                          += dE;
                if (mc->screen_E_trapped < 0) { mc->screen_E_trapped = 0; }
            }
            else if (net_charge > 0)
            {
                double E = mcmd_emitted_electron_energy(
                        mc, ad, match,
                        possible_transitions[min_value_index].final_state,
                        possible_transitions[min_value_index].type);

                /* none: every emitted electron leaves.  pbgrid: it escapes if
                 * it can climb out of the self-consistent potential of ions and
                 * electrons at its atom, else it joins the cloud with its
                 * kinetic energy and its potential energy there, -D_k (the
                 * kinetic part alone until the next solve). */
                if (mc->screen_mode == 0)
                {
                    mc->escaped_photoelectrons += 1;
                }
                else
                {
                    const double D_k = mcmd_pb_depth(mc, mdatoms, state->x, k);

                    if (E > D_k)
                    {
                        mc->escaped_photoelectrons += 1;
                    }
                    else
                    {
                        const double dep = (E > 0 ? E : 0);

                        mc->num_free_electrons               += 1;
                        mc->thermalized_free_electron_energy += MCMD_TRAPPED_ELECTRON_EV;
                        mc->pb_Etot                          += dep - D_k;
                        mc->pbd_N                            -= D_k;
                        mc->screen_E_trapped                 += dep;
                        mc->eb_dep                           += dep;
                    }
                }
            }
            else if (net_charge < 0)
            {
                mc->num_free_electrons               += net_charge;
                mc->thermalized_free_electron_energy +=
                    net_charge * MCMD_TRAPPED_ELECTRON_EV;
            }

            mc->electron_density = mc->num_free_electrons / mc->current_volume;
            mc->electron_temperature =
                MCMD_KELVIN_TO_EV *
                ((2.0 / (3.0 * MCMD_BOLTZMANN)) *
                 ((MCMD_EV_TO_JOULE * mc->thermalized_free_electron_energy) /
                  mc->num_free_electrons));

            if (mc->bDetailedOutput)
            {
                mcmd_trn_write(mc->ftrans, t, k, -1, possible_transitions[min_value_index].type, match,
                               possible_transitions[min_value_index].final_state);
            }

            mc->atom_configurations[k][0] = possible_transitions[min_value_index].final_state[0];
            mc->atom_configurations[k][1] = possible_transitions[min_value_index].final_state[1];
            mc->atom_configurations[k][2] = possible_transitions[min_value_index].final_state[2];

            mdatoms->chargeA[k] += net_charge;

            switch (possible_transitions[min_value_index].type)
            {
                case 0: mc->count_auger += 1; break;
                case 1: mc->count_fluorescence += 1; break;
                case 2: mc->count_photoionization += 1; break;
                case 3:
                case 4:
                {
                    /* diagnostic: which shell lost (ionization) or gained
                     * (recombination) the electron */
                    int sh;

                    for (sh = 0; sh < MCMD_NSHELL; sh++)
                    {
                        int d = possible_transitions[min_value_index].final_state[sh] - match[sh];

                        if (d < 0) { mc->coll_ion_shell[sh] += 1; }
                        if (d > 0) { mc->coll_rec_shell[sh] += 1; }
                    }
                    if (possible_transitions[min_value_index].type == 3) { mc->count_coll_ion += 1; }
                    else { mc->count_coll_rec += 1; }
                    break;
                }
                case 5: mc->count_coll_exc += 1; break;
                case 6: mc->count_coll_deexc += 1; break;
                default: break;
            }

            /* The configuration changed, so the caches the charge-transfer
             * pass reads must follow it. */
            mcmd_refresh_atom(mc, mdatoms, k);
        }
    }
}

/* ---------------------------------------------------------------- *
 * Per-step bookkeeping and output
 * ---------------------------------------------------------------- */

static void mcmd_update_volume(t_mcionize *mc, t_mdatoms *mdatoms,
                               t_state *state)
{
    double radius_of_gyration;

    mcmd_sample_geometry(mc, mdatoms, state, &radius_of_gyration);

    mc->rg_factor  = radius_of_gyration / mc->radius_gyration_equilibrium;
    mc->radius_now = mc->radius_of_sample * mc->rg_factor;
    /* Sphere volume, so the linear expansion factor enters cubed. */
    mc->current_volume = mc->sample_volume *
                         mc->rg_factor * mc->rg_factor * mc->rg_factor;
    {
        /* Volume-equivalent radius of the gyration ellipsoid (equal to
         * sqrt(5/3) R_g for a sphere); the volume in cm^3 as above. */
        mc->radius_now     = cbrt(mc->ell_a[0] * mc->ell_a[1] * mc->ell_a[2]);
        mc->current_volume = 4.0 / 3.0 * M_PI * 1e-21 * mc->ell_a[0] * mc->ell_a[1] * mc->ell_a[2];
    }
}


static void mcmd_write_step_output(t_mcionize *mc, t_mdatoms *mdatoms, double t)
{
    double mean_charge = 0;
    int    i;

    for (i = 0; i < mdatoms->nr; i++)
    {
        mean_charge += mdatoms->chargeA[i];
    }
    mean_charge /= (double)mdatoms->nr;

    if (!mc->bDetailedOutput)
    {
        return;
    }

    if ((mc->output_step % mc->charge_stride) == 0)
    {
        fwrite(mdatoms->chargeA, sizeof(mdatoms->chargeA[0]), mdatoms->nr,
               mc->fcharges);
    }

    if ((mc->output_step % mc->text_stride) == 0)
    {
        fprintf(mc->fmean,  "%lf %lf\n", t, mean_charge);
        fprintf(mc->fpulse, "%lf %lf\n", t, mc->intensity);

        /* Plain numeric columns so the file loads directly with np.loadtxt:
         *   time [ps]   e- density [cm^-3]   e- temperature [eV]   R_g / R_g(0)
         * (The Debye length the manual used to promise was never computed.) */
        fprintf(mc->felectron, "%g %g %g %g\n",
                t, mc->electron_density, mc->electron_temperature, mc->rg_factor);
    }

    /* Periodic flush: keeps the files usable while a run is in flight without
     * paying a write() every step. */
    if ((++mc->output_step % MCMD_OUTPUT_FLUSH_STEPS) == 0)
    {
        fflush(mc->fcharges);
        fflush(mc->fmean);
        fflush(mc->fpulse);
        fflush(mc->felectron);
        fflush(mc->ftrans);
    }
}

/* ---------------------------------------------------------------- *
 * Restart state
 *
 * mcionize_checkpoint() writes ./simulation_output/mcionize_state.bin at the
 * start of every step GROMACS writes a checkpoint on, before the module
 * advances, so it is the electronic state that goes with the .cpt's x and v:
 * mdrun -cpi resumes at that step and runs it in full.  The file before is
 * kept as mcionize_state_prev.bin (as GROMACS keeps state_prev.cpt); at a
 * restart the one whose step is the .cpt's is used.
 *
 * Version 3 adds, after the version 2 blocks (geometry, event counters), the
 * module's scalar state by name (the PB electrons' total and kinetic energy,
 * the box and grid sizes, the step counters, the cumulative diagnostics, the
 * energy draws' stream), rand()'s state, and the length of every output file
 * at the checkpoint: a restart cuts each file back to that length and
 * appends, so the steps a killed run did after its last checkpoint are not
 * doubled.  The PB density itself is not saved (up to 384^3 cells): the first
 * step after a restart solves afresh (from the ion charges) on the grid of the
 * checkpoint, rebuilt from its saved size, spacing and origin, and since the electrons' total energy is restored their kinetic energy
 * follows from it, E_kin = E_tot - U, as at any solve.  That one extra solve
 * (and long-range update) is the only difference from an uninterrupted run
 * with the same seed.
 * ---------------------------------------------------------------- */

/* rand()'s state lives in our own buffer (initstate / setstate, glibc), so it
 * can be saved.  128 bytes is glibc's default generator (TYPE_3), so a seed
 * gives the same sequence as srand() did.  Two buffers, because setstate()
 * first stores the position into the buffer it leaves. */
static char mcmd_rngbuf[2][128];
static int  mcmd_rngcur = -1;

static void mcmd_srand(unsigned int seed)
{
    if (mcmd_rngcur < 0)
    {
        mcmd_rngcur = 0;
        initstate(seed, mcmd_rngbuf[0], sizeof(mcmd_rngbuf[0]));
    }
    else
    {
        srand(seed);
    }
}

static void mcmd_rng_save(char out[128])
{
    setstate(mcmd_rngbuf[mcmd_rngcur]);
    memcpy(out, mcmd_rngbuf[mcmd_rngcur], 128);
}

static void mcmd_rng_load(const char in[128])
{
    const int nx = 1 - mcmd_rngcur;

    memcpy(mcmd_rngbuf[nx], in, 128);
    setstate(mcmd_rngbuf[nx]);
    mcmd_rngcur = nx;
}

/* The output files, for their lengths at a checkpoint.  The first
 * MCMD_NLOG_INIT are opened by mcionize_init (mcmd-detailed-output), the rest
 * on their first write. */
static const struct
{
    const char *name;
    size_t      off;
} mcmd_logs[] = {
    {"charges_over_time.bin", offsetof(t_mcionize, fcharges)},
    {"mean_charge_vs_time.txt", offsetof(t_mcionize, fmean)},
    {"pulse_profile.txt", offsetof(t_mcionize, fpulse)},
    {"electron_data.txt", offsetof(t_mcionize, felectron)},
    {"transitions.bin", offsetof(t_mcionize, ftrans)},
    {"pb_log.txt", offsetof(t_mcionize, fpb)},
    {"pbg_density.bin", offsetof(t_mcionize, fpbw)},
    {"pbg_box.txt", offsetof(t_mcionize, fbox)},
    {"pb_ediag.txt", offsetof(t_mcionize, fpbd)},
    {"pbg_lr.txt", offsetof(t_mcionize, flr)},
    {"pbg_timing.txt", offsetof(t_mcionize, ftd)},
    {"pbg_outlat.txt", offsetof(t_mcionize, fol)},
    {"pbg_outpot.txt", offsetof(t_mcionize, fop)},
    {"pbfield_events.txt", offsetof(t_mcionize, fpbfi_ev)},
    {"coll_event_log.txt", offsetof(t_mcionize, fcev)},
    {"collisional_log.txt", offsetof(t_mcionize, fcoll)},
    {"energy_ledger.txt", offsetof(t_mcionize, fled)},
    {"screening_log.txt", offsetof(t_mcionize, fscreen)}
};
#define MCMD_NLOG         ((int)(sizeof(mcmd_logs) / sizeof(mcmd_logs[0])))
#define MCMD_NLOG_INIT    5
#define MCMD_LOGFP(mc, k) ((FILE **)((char *)(mc) + mcmd_logs[k].off))

static gmx_bool mcmd_kv_get(const t_mcmd_kv *kv, int nkv, const char *name, double *v)
{
    int i;

    for (i = 0; i < nkv; i++)
    {
        if (strcmp(kv[i].name, name) == 0)
        {
            *v = kv[i].v;
            return TRUE;
        }
    }
    return FALSE;
}

/* The scalar state by name: bWrite copies it into kv and returns the count;
 * otherwise sets it from kv, counting the names kv lacks in *nmiss. */
static int mcmd_state_kv(t_mcionize *mc, t_mcmd_kv *kv, int nkv, gmx_bool bWrite, int *nmiss)
{
    int    n = 0;
    double v;

#define MCMD_KV(f) \
    do \
    { \
        if (bWrite) \
        { \
            strncpy(kv[n].name, #f, sizeof(kv[n].name) - 1); \
            kv[n].name[sizeof(kv[n].name) - 1] = '\0'; \
            kv[n].v = (double)(mc->f); \
        } \
        else if (mcmd_kv_get(kv, nkv, #f, &v)) { mc->f = v; } \
        else { (*nmiss)++; } \
        n++; \
    } while (0)
    /* step counters */
    MCMD_KV(screen_step); MCMD_KV(ct_step); MCMD_KV(ct_steps_since_transfer); MCMD_KV(ct_skipped);
    MCMD_KV(coll_step); MCMD_KV(output_step);
    /* the plasma and the PB electrons */
    MCMD_KV(electron_temperature); MCMD_KV(electron_density); MCMD_KV(rg_factor); MCMD_KV(current_volume);
    MCMD_KV(radius_now);
    MCMD_KV(screen_E_trapped); MCMD_KV(pb_Etot); MCMD_KV(pb_etot_init); MCMD_KV(pb_Eclamp); MCMD_KV(pb_W);
    MCMD_KV(pressure_work); MCMD_KV(screen_evaporated);
    MCMD_KV(eb_dep); MCMD_KV(eb_adiab); MCMD_KV(eb_evap);
    MCMD_KV(pb_N); MCMD_KV(pb_h); MCMD_KV(pb_Dmean); MCMD_KV(pb_nsolve);
    /* grid and box */
    MCMD_KV(pbg_R0); MCMD_KV(pbg_hbox); MCMD_KV(pbg_n); MCMD_KV(pbg_ncap); MCMD_KV(pbg_nmax); MCMD_KV(pbg_capset);
    MCMD_KV(pbg_nregrid); MCMD_KV(pbg_have); MCMD_KV(pbg_dx); MCMD_KV(pbg_R);
    MCMD_KV(pbg_lo[0]); MCMD_KV(pbg_lo[1]); MCMD_KV(pbg_lo[2]);
    MCMD_KV(pbg_c[0]); MCMD_KV(pbg_c[1]); MCMD_KV(pbg_c[2]);
    /* collisional processes, energy ledger, field ionization: cumulative */
    MCMD_KV(count_coll_ion); MCMD_KV(count_coll_rec); MCMD_KV(count_coll_exc); MCMD_KV(count_coll_deexc);
    MCMD_KV(coll_ion_shell[0]); MCMD_KV(coll_ion_shell[1]); MCMD_KV(coll_ion_shell[2]);
    MCMD_KV(coll_rec_shell[0]); MCMD_KV(coll_rec_shell[1]); MCMD_KV(coll_rec_shell[2]);
    MCMD_KV(coll_energy); MCMD_KV(coll_pe_energy);
    MCMD_KV(led_Jii); MCMD_KV(led_Jie); MCMD_KV(led_nchg);
    MCMD_KV(pbfi_nion); MCMD_KV(pbfi_Eb); MCMD_KV(pbfi_quiet); MCMD_KV(pbfi_done);
    MCMD_KV(max_transfer_separation); MCMD_KV(count_escape_fallback); MCMD_KV(ctab_clamped);
    /* the energy draws' own stream */
    MCMD_KV(erand_state[0]); MCMD_KV(erand_state[1]); MCMD_KV(erand_state[2]);
#undef MCMD_KV
    return n;
}

static void mcmd_state_write(t_mcionize *mc, const t_mdatoms *mdatoms)
{
    char            tmp[512], cur[512], prev[512], rng[128];
    FILE           *fp;
    int             header[3], counters[4], nkv, k, i;
    gmx_large_int_t step = mc->screen_step;
    double          geom[7];
    t_mcmd_kv      *kv;

    snew(kv, MCMD_KV_MAX);
    nkv = mcmd_state_kv(mc, kv, 0, TRUE, NULL);
    if (nkv + MCMD_NLOG > MCMD_KV_MAX)
    {
        gmx_fatal(FARGS, "mcionize: %d state variables, raise MCMD_KV_MAX", nkv + MCMD_NLOG);
    }
    for (k = 0; k < MCMD_NLOG; k++)
    {
        FILE *f = *MCMD_LOGFP(mc, k);

        sprintf(kv[nkv].name, "log:%s", mcmd_logs[k].name);
        if (f != NULL)
        {
            fflush(f);
        }
        kv[nkv++].v = (f != NULL) ? (double)ftello(f) : -1.0;
    }
    mcmd_rng_save(rng);

    sprintf(tmp, "%s/%s.tmp", MCMD_OUTDIR, MCMD_STATE_FILE);
    sprintf(cur, "%s/%s", MCMD_OUTDIR, MCMD_STATE_FILE);
    sprintf(prev, "%s/%s", MCMD_OUTDIR, MCMD_STATE_PREV);
    fp = fopen(tmp, "wb");
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open '%s' for writing.", tmp);
    }
    header[0] = MCMD_STATE_MAGIC;
    header[1] = MCMD_STATE_VERSION;
    header[2] = mc->natoms;
    fwrite(header, sizeof(int), 3, fp);
    fwrite(&step, sizeof(step), 1, fp);

    geom[0] = mc->total_mass;
    geom[1] = mc->radius_of_sample;
    geom[2] = mc->radius_gyration_equilibrium;
    geom[3] = mc->sample_volume;
    geom[4] = mc->num_free_electrons;
    geom[5] = mc->thermalized_free_electron_energy;
    geom[6] = mc->escaped_photoelectrons;
    fwrite(geom, sizeof(double), 7, fp);

    counters[0] = mc->count_auger;
    counters[1] = mc->count_fluorescence;
    counters[2] = mc->count_photoionization;
    counters[3] = mc->count_chargetransfer;
    fwrite(counters, sizeof(int), 4, fp);

    fwrite(&nkv, sizeof(int), 1, fp);
    fwrite(kv, sizeof(t_mcmd_kv), nkv, fp);
    fwrite(rng, 1, sizeof(rng), fp);

    for (i = 0; i < mc->natoms; i++)
    {
        double charge = (double)mdatoms->chargeA[i];

        fwrite(mc->atom_configurations[i], sizeof(int), MCMD_NSHELL, fp);
        fwrite(mc->GS_configurations[i], sizeof(int), MCMD_NSHELL, fp);
        fwrite(&charge, sizeof(double), 1, fp);
    }
    /* complete or not at all: a run killed while writing keeps the last one */
    if (ferror(fp) || fclose(fp) != 0)
    {
        gmx_fatal(FARGS, "Writing '%s' failed (disk full?).", tmp);
    }
    rename(cur, prev);
    if (rename(tmp, cur) != 0)
    {
        gmx_fatal(FARGS, "Could not rename '%s' to '%s'.", tmp, cur);
    }
    sfree(kv);
    mc->nstate++;
}

/* Read the state file that goes with the checkpoint being resumed (its step
 * is the .cpt's, ir->init_step), and cut the output files back to their
 * length at that checkpoint.  mcmd_state_apply() sets the state once the
 * module is set up.  Master only, before the output files are opened. */
static void mcmd_state_load(t_mcionize *mc, const t_inputrec *ir)
{
    const char *names[2] = {MCMD_STATE_FILE, MCMD_STATE_PREV};
    char        seen[2][64];
    t_mcmd_rst *r;
    int         f, i, k, header[3];

    snew(r, 1);
    for (f = 0; f < 2; f++)
    {
        FILE *fp;

        sprintf(r->path, "%s/%s", MCMD_OUTDIR, names[f]);
        fp = fopen(r->path, "rb");
        if (fp == NULL)
        {
            sprintf(seen[f], "missing");
            continue;
        }
        if (fread(header, sizeof(int), 3, fp) != 3 || header[0] != MCMD_STATE_MAGIC)
        {
            gmx_fatal(FARGS, "'%s' is not an mcionize state file.", r->path);
        }
        if (header[2] != mc->natoms)
        {
            gmx_fatal(FARGS, "'%s' holds %d atoms but this system has %d.", r->path, header[2], mc->natoms);
        }
        r->version = header[1];
        if (r->version == 2)
        {
            /* the end-of-run file of older builds: no step, no PB electrons */
            if (mc->screen_mode != 0)
            {
                gmx_fatal(FARGS, "'%s' is a version 2 state file (an older build), without the state of the "
                          "PB electrons: a run with mcmd-electrons = pbgrid cannot be continued from it.", r->path);
            }
            r->step = ir->init_step;
        }
        else if (r->version == MCMD_STATE_VERSION)
        {
            if (fread(&r->step, sizeof(r->step), 1, fp) != 1)
            {
                gmx_fatal(FARGS, "'%s' is truncated.", r->path);
            }
        }
        else
        {
            gmx_fatal(FARGS, "'%s' has version %d, this build reads versions 2 and %d.", r->path, r->version,
                      MCMD_STATE_VERSION);
        }
        if (r->step != ir->init_step)
        {
            sprintf(seen[f], "from step %lld", (long long)r->step);
            fclose(fp);
            continue;
        }
        if (fread(r->geom, sizeof(double), 7, fp) != 7 || fread(r->counters, sizeof(int), 4, fp) != 4)
        {
            gmx_fatal(FARGS, "'%s' is truncated.", r->path);
        }
        if (r->version >= 3)
        {
            if (fread(&r->nkv, sizeof(int), 1, fp) != 1 || r->nkv < 0 || r->nkv > MCMD_KV_MAX)
            {
                gmx_fatal(FARGS, "'%s' is damaged.", r->path);
            }
            snew(r->kv, r->nkv);
            if (fread(r->kv, sizeof(t_mcmd_kv), r->nkv, fp) != (size_t)r->nkv ||
                fread(r->rng, 1, sizeof(r->rng), fp) != sizeof(r->rng))
            {
                gmx_fatal(FARGS, "'%s' is truncated.", r->path);
            }
        }
        snew(r->conf, (size_t)mc->natoms * MCMD_NSHELL);
        snew(r->gs, (size_t)mc->natoms * MCMD_NSHELL);
        snew(r->q, mc->natoms);
        for (i = 0; i < mc->natoms; i++)
        {
            if (fread(r->conf + (size_t)i * MCMD_NSHELL, sizeof(int), MCMD_NSHELL, fp) != MCMD_NSHELL ||
                fread(r->gs + (size_t)i * MCMD_NSHELL, sizeof(int), MCMD_NSHELL, fp) != MCMD_NSHELL ||
                fread(&r->q[i], sizeof(double), 1, fp) != 1)
            {
                gmx_fatal(FARGS, "'%s' is truncated at atom %d.", r->path, i + 1);
            }
        }
        fclose(fp);
        break;
    }
    if (f == 2)
    {
        gmx_fatal(FARGS, "Continuing from step %lld needs the mcionize state written at that checkpoint, but "
                  "%s/%s is %s and %s is %s.  (Checkpoints of builds before state version 3 have no state file "
                  "of their own.)", (long long)ir->init_step, MCMD_OUTDIR, MCMD_STATE_FILE, seen[0],
                  MCMD_STATE_PREV, seen[1]);
    }
    mc->rst = r;

    /* the output files back to the checkpoint; one not yet open then was
     * started after it */
    if (r->version >= 3)
    {
        for (k = 0; k < MCMD_NLOG; k++)
        {
            char        path[512], key[48];
            double      off = -1;
            struct stat st;

            sprintf(path, "%s/%s", MCMD_OUTDIR, mcmd_logs[k].name);
            sprintf(key, "log:%s", mcmd_logs[k].name);
            if (stat(path, &st) != 0)
            {
                if (mcmd_kv_get(r->kv, r->nkv, key, &off) && off > 0)
                {
                    gmx_fatal(FARGS, "'%s' was %.0f bytes at the checkpoint but is gone.", path, off);
                }
                continue;
            }
            if (!mcmd_kv_get(r->kv, r->nkv, key, &off) || off < 0)
            {
                off = 0;
            }
            if ((double)st.st_size < off)
            {
                gmx_fatal(FARGS, "'%s' is %.0f bytes, shorter than at the checkpoint (%.0f).", path,
                          (double)st.st_size, off);
            }
            if (truncate(path, (off_t)off) != 0)
            {
                gmx_fatal(FARGS, "Could not cut '%s' back to %.0f bytes.", path, off);
            }
        }
    }
}

/* Set the state read by mcmd_state_load() (master). */
static void mcmd_state_apply(t_mcionize *mc, t_mdatoms *mdatoms)
{
    t_mcmd_rst *r = mc->rst;
    int         i, j, k, nmiss = 0;

    mc->total_mass                       = r->geom[0];
    mc->radius_of_sample                 = r->geom[1];
    mc->radius_gyration_equilibrium      = r->geom[2];
    mc->sample_volume                    = r->geom[3];
    mc->num_free_electrons               = r->geom[4];
    mc->thermalized_free_electron_energy = r->geom[5];
    mc->escaped_photoelectrons           = r->geom[6];

    mc->count_auger           = r->counters[0];
    mc->count_fluorescence    = r->counters[1];
    mc->count_photoionization = r->counters[2];
    mc->count_chargetransfer  = r->counters[3];

    for (i = 0; i < mc->natoms; i++)
    {
        for (j = 0; j < MCMD_NSHELL; j++)
        {
            mc->atom_configurations[i][j] = r->conf[(size_t)i * MCMD_NSHELL + j];
            mc->GS_configurations[i][j]   = r->gs[(size_t)i * MCMD_NSHELL + j];
        }
        mdatoms->chargeA[i] = (real)r->q[i];
    }

    if (r->version >= 3)
    {
        mcmd_state_kv(mc, r->kv, r->nkv, FALSE, &nmiss);
        if (nmiss > 0)
        {
            gmx_fatal(FARGS, "'%s' lacks %d of the module's state variables.", r->path, nmiss);
        }
        mcmd_rng_load(r->rng);
        /* the first solve starts from the ion charges, on the grid of the
         * checkpoint, rebuilt when the arrays are allocated (mcmd_pbg_alloc) */
        mc->pbg_rst_grid = mc->pbg_have;
        mc->pbg_have     = FALSE;
        mc->pbg_cold     = TRUE;
        /* the files that were open at the checkpoint: append, the headers are
         * there already */
        for (k = MCMD_NLOG_INIT; k < MCMD_NLOG; k++)
        {
            char   key[48];
            double off;

            sprintf(key, "log:%s", mcmd_logs[k].name);
            if (mcmd_kv_get(r->kv, r->nkv, key, &off) && off > 0)
            {
                *MCMD_LOGFP(mc, k) = mcmd_open_output(mcmd_logs[k].name,
                                                      strstr(mcmd_logs[k].name, ".bin") ? "ab" : "a");
            }
        }
    }

    sfree(r->kv);
    sfree(r->conf);
    sfree(r->gs);
    sfree(r->q);
    sfree(r);
    mc->rst = NULL;

    /* Restored configurations, so the caches have to be rebuilt from them. */
    for (i = 0; i < mc->natoms; i++)
    {
        mcmd_refresh_atom(mc, mdatoms, i);
    }
}

/* ---------------------------------------------------------------- *
 * First-step initialisation
 * ---------------------------------------------------------------- */

static void mcmd_first_step(t_mcionize *mc, t_mdatoms *mdatoms, t_state *state)
{
    double radius_of_gyration = 0;
    int    i, j;

    mc->electron_temperature               = 0.0;
    mc->num_free_electrons                 = 0.0;
    mc->thermalized_free_electron_energy   = 0.0;
    mc->escaped_photoelectrons             = 0.0;

    for (i = 0; i < mdatoms->nr; i++)
    {
        mc->total_mass += mdatoms->massT[i];
        /* The model tracks integer net charge, so the force field's partial
         * charges are discarded. */
        mdatoms->chargeA[i] = 0.0;
    }

    /* Sample size about the centre of mass.  This used to be measured from
     * the box origin - the radius as the distance to the corner of the
     * largest coordinates, R_g about (0,0,0) - which for a protein sitting
     * in the middle of its box overestimated the radius several times over
     * and made R_g follow the box placement rather than the expansion. */
    mcmd_sample_geometry(mc, mdatoms, state, &radius_of_gyration);
    mc->radius_gyration_equilibrium = radius_of_gyration;

    /* Radius of the uniform sphere with the same R_g. */
    mc->radius_of_sample = sqrt(5.0 / 3.0) * radius_of_gyration;
    mc->radius_now       = mc->radius_of_sample;

    mcmd_srand(mc->seed);

    mc->num_free_electrons               = 0;
    mc->thermalized_free_electron_energy = 0;
    /* Sphere volume in cm^3 (1 nm^3 = 1e-21 cm^3), so the electron density
     * below comes out per cm^3. */
    mc->sample_volume = 4.0 / 3.0 * M_PI * 1e-21 *
                        mc->radius_of_sample * mc->radius_of_sample *
                        mc->radius_of_sample;

    /* Ground states */
    for (i = 0; i < mdatoms->nr; i++)
    {
        const int *gs = mcmd_idx_ground_state(mcmd_atom2idx(mdatoms->massT[i]));

        if (gs != NULL)
        {
            for (j = 0; j < MCMD_NSHELL; j++)
            {
                mc->atom_configurations[i][j] = gs[j];
                mc->GS_configurations[i][j]   = gs[j];
            }
        }
    }

    if (mc->bInitCharges)
    {
        char  line[100];
        FILE *fp = fopen("charges.txt", "r");

        if (fp == NULL)
        {
            gmx_fatal(FARGS, "Could not open 'charges.txt' for reading.");
        }

        while (fgets(line, sizeof(line), fp))
        {
            int charge_index, temp_charge;

            if (sscanf(line, "%d %d", &charge_index, &temp_charge) != 2)
            {
                gmx_fatal(FARGS,
                          "Malformed line in 'charges.txt': expected "
                          "'<index> <charge>', got '%s'.", line);
            }
            if (charge_index < 1 || charge_index > mdatoms->nr)
            {
                gmx_fatal(FARGS,
                          "Atom index %d in 'charges.txt' is outside "
                          "[1, %d].", charge_index, mdatoms->nr);
            }
            mdatoms->chargeA[charge_index - 1] = (real)temp_charge;
        }
        fclose(fp);

        /* Strip the corresponding electrons, outermost shell first. */
        for (i = 0; i < mdatoms->nr; i++)
        {
            for (j = 0; j < (int)mdatoms->chargeA[i]; j++)
            {
                int shell = mcmd_outermost_occupied(mc->atom_configurations[i]);

                if (shell < 0)
                {
                    gmx_fatal(FARGS,
                              "Atom %d was given charge %g in 'charges.txt', "
                              "but it does not have that many electrons.",
                              i + 1, (double)mdatoms->chargeA[i]);
                }
                mc->atom_configurations[i][shell] -= 1;
            }
        }
    }

    /* Seed the per-atom caches now that the configurations are final. */
    for (i = 0; i < mdatoms->nr; i++)
    {
        mcmd_refresh_atom(mc, mdatoms, i);
    }
}

/* ---------------------------------------------------------------- *
 * Public entry points
 *
 * mdp options read by this module.  All of them only take effect under
 * mdrun -ionize, which is also what enables the altered force field.
 *
 *   mcmd-charge-transfer          enable charge transfer
 *   mcmd-charge-transfer-idle     quiet steps before the pass goes idle
 *                                 (0 = built-in default, negative = never)
 *   mcmd-charge-transfer-recheck  step interval to recheck on once idle
 *   mcmd-autostop                 stop once E_kin/E_tot exceeds the threshold
 *   mcmd-autostop-threshold       that threshold
 *   mcmd-initial-charges          read initial charges from charges.txt
 *   mcmd-detailed-output          log electronic dynamics (heavy I/O)
 *   mcmd-collisional-ionization   NOT IMPLEMENTED
 *   mcmd-pulse-peak-time          [fs]
 *   mcmd-pulse-fwhm               [fs], full width at half maximum
 *   mcmd-pulse-photons            total photons in the pulse
 *   mcmd-pulse-focal-diameter     [nm]
 *   mcmd-pulse-photon-energy      [eV]
 * ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- *
 * Run summary (2026-10-28): what the module will do, written to md.log
 * and stderr when a run starts, and warnings for settings outside the
 * validated model.  Warnings only: the run goes on.
 * ---------------------------------------------------------------- */

/* The converged production settings the warnings compare against. */
#define MCMD_REF_DT_AS        2.0     /* MD step, as       */
#define MCMD_REF_PB_INTERVAL  0.2     /* fs between solves */
#define MCMD_REF_NEWTON_TOL   1e-2
#define MCMD_REF_CELLS        2.5     /* cells per s_e     */
#define MCMD_REF_SMOOTHING    0.25    /* s_e, nm           */
#define MCMD_LR_AUTO_ATOMS    10000   /* auto long range above this */

static void mcmd_print_summary(FILE *fp, const t_mcionize *mc, const t_inputrec *ir,
                               const t_mdatoms *md, const t_commrec *cr, gmx_bool bRestart,
                               int nHarmonicBonds)
{
    const char *bar  = "==================================================================================";
    const double dt_as = ir->delta_t * 1e6;
    int          count[MCMD_ATOM_TABLE_SIZE] = { 0 };
    double       qion = 0;
    char         warn[16][320];
    int          nw = 0, i, n = 0;
    gmx_bool     te_all = TRUE, sw_all = TRUE;

    for (i = 0; i < md->nr; i++)
    {
        count[mc->elem_of_atom[i]]++;
        qion += md->chargeA[i];
    }
    fprintf(fp, "\n%s\n MolDStruct MC/MD%s\n", bar, bRestart ? " (continued from a checkpoint)" : "");
    fprintf(fp, " System        %d atoms (", md->nr);
    for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
    {
        if (count[i] > 0)
        {
            fprintf(fp, "%s%s %d", n++ ? ", " : "", mcmd_idx2symbol(i), count[i]);
            te_all = te_all && (mc->atomdata[i].total_energies != NULL);
            sw_all = sw_all && mc->atomdata[i].bWeightsFile;
        }
    }
    /* At a fresh start chargeA still holds the topology's partial charges, which the module replaces by the
     * ionization charges (all 0) at the first step; after a restart it holds the ionization charges. */
    if (bRestart)
    {
        fprintf(fp, "), ionization charge %+.0f\n", qion);
    }
    else
    {
        fprintf(fp, "), all neutral at the start\n");   /* partial charges replaced at step 0 */
    }
    {
        const double d_um  = ir->mcmd_pulse_focal_diameter * 1e-3;
        const double fluen = ir->mcmd_pulse_photons / (M_PI * 0.25 * d_um * d_um);

        fprintf(fp, " Pulse         %g eV, Gaussian, FWHM %g fs, peak at %g fs, %g photons into %g nm\n"
                "               (%.3g photons/um^2); bandwidth %g\n",
                ir->mcmd_pulse_photon_energy, ir->mcmd_pulse_fwhm, ir->mcmd_pulse_peak_time,
                ir->mcmd_pulse_photons, ir->mcmd_pulse_focal_diameter, fluen, ir->mcmd_pulse_bandwidth);
    }
    fprintf(fp, " Run           %s steps x %g as = %g fs; seed %u (%s); %d rank%s\n",
            gmx_step_str(ir->nsteps, warn[0]), dt_as, ir->nsteps * dt_as * 1e-3, mc->seed,
            ir->mcmd_seed >= 0 ? "mcmd-seed" : "chosen at start, mcmd-seed = -1", cr->nnodes,
            cr->nnodes > 1 ? "s" : "");
    if (ir->mcmd_atomic_data == 1)
    {
        fprintf(fp, " Atomic data   built in: screened-hydrogenic shells (K, L, M) at %g eV, mcmd-atomic-data = hydrogenic\n"
                "               (written to %s/Atomic_data)\n", ir->mcmd_pulse_photon_energy, MCMD_OUTDIR);
    }
    else
    {
        fprintf(fp, " Atomic data   ./Atomic_data (mcmd-atomic-data = files); total energies %s, statistical weights %s\n",
                te_all ? "yes" : "no", !mc->coll_on ? "not used" : sw_all ? "from the files" : "computed");
    }

    fprintf(fp, "\n Processes\n");
    fprintf(fp, "   %-40s on (always)\n", "photoionization, Auger, fluorescence");
    if (mc->bChargeTransfer)
    {
        fprintf(fp, "   %-40s on   %s, H-H %s, downhill %s\n", "charge transfer",
                mc->ct_h_donor_only ? "hydrogen donors only" : "any donor",
                mc->bAllowHCT ? "on" : "off", mc->bDownhill ? "on" : "off");
    }
    else
    {
        fprintf(fp, "   %-40s off\n", "charge transfer");
    }
    if (mc->screen_mode == 5)
    {
        fprintf(fp, "   %-40s on   mcmd-electrons = pbgrid\n", "trapped electrons (Poisson-Boltzmann)");
        fprintf(fp, "   %-40s s_e %g nm, %g cells/s_e, box <= %g R, <= %d cells,\n"
                "   %-40s solve every %g fs, Newton tol %g\n", "    grid",
                mc->pbg_se, mc->pbg_se / mc->pbg_dxt, mc->pbg_box, mc->pbg_nglimit, "",
                mc->pb_nstep * dt_as * 1e-3, mc->pbg_ntol);
        fprintf(fp, "   %-40s %s\n", "    evaporation", mc->screen_evap ? "on" : "off");
        if (mc->pbg_lr)
        {
            fprintf(fp, "   %-40s %g nm cutoff + grid long range\n", "ion-ion forces", mc->lr_rc);
        }
        else
        {
            fprintf(fp, "   %-40s all pairs\n", "ion-ion forces");
        }
        if (mc->coll_on)
        {
            fprintf(fp, "   %-40s on   %s\n", "electron-impact ionization",
                    mc->coll_lotz ? "Lotz, built in (mcmd-coll-data = lotz)"
                    : "from files (mcmd-coll-data = files)");
            fprintf(fp, "   %-40s %s\n", "three-body recombination",
                    mc->coll_recomb ? (mc->coll_no_h_recomb ? "on (not onto H)" : "on") : "off");
            fprintf(fp, "   %-40s %s\n", "collisional excitation", mc->coll_excite ? "on" : "off");
        }
        else
        {
            fprintf(fp, "   %-40s off\n", "electron-impact ionization");
        }
        if (mc->pbfi)
        {
            fprintf(fp, "   %-40s on   L = %g nm\n", "field ionization", mc->pbfi_L);
        }
        else
        {
            fprintf(fp, "   %-40s off\n", "field ionization");
        }
    }
    else
    {
        fprintf(fp, "   %-40s none (mcmd-electrons = none: every electron leaves)\n", "trapped electrons");
    }
    fprintf(fp, "\n Output        simulation_output/: %s, charges_over_time.bin every %d steps;\n"
            "               autostop %s\n",
            mc->bDetailedOutput ? "transitions.bin (mcmd-detailed-output = 1)" : "summary files only",
            mc->charge_stride, ir->mcmd_autostop ? "on" : "off");

    /* warnings */
#define MCMD_WARN(...) do { if (nw < 16) { snprintf(warn[nw++], sizeof(warn[0]), __VA_ARGS__); } } while (0)
    if (dt_as > MCMD_REF_DT_AS * (1 + 1e-6))
    {
        MCMD_WARN("dt = %g as: the model is converged at <= %g as", dt_as, MCMD_REF_DT_AS);
    }
    if (nHarmonicBonds > 0)
    {
        MCMD_WARN("%d harmonic bonds: they neither weaken nor break under -ionize; set morse = yes in the mdp\n"
                  "     so that grompp turns them into Morse bonds", nHarmonicBonds);
    }
    if (ir->mcmd_pulse_photons <= 0 || ir->mcmd_pulse_photon_energy <= 0)
    {
        MCMD_WARN("mcmd-pulse-photons or mcmd-pulse-photon-energy is 0: nothing will be ionized");
    }
    if (mc->bAllowHCT)
    {
        MCMD_WARN("mcmd-allow-h-ct = 1: resonant H-H charge transfer, not fully validated (see README)");
    }
    if (mc->bChargeTransfer && !mc->ct_h_donor_only)
    {
        MCMD_WARN("mcmd-ct-h-donor-only = 0: not fully validated");
    }
    if (mc->screen_mode == 5)
    {
        if (ir->mcmd_pb_interval > MCMD_REF_PB_INTERVAL * (1 + 1e-6) ||
            ir->mcmd_pb_newton_tol > MCMD_REF_NEWTON_TOL * (1 + 1e-6) ||
            ir->mcmd_pb_cells_per_smoothing < MCMD_REF_CELLS * (1 - 1e-6))
        {
            MCMD_WARN("PB solver coarser than the converged settings (mcmd-pb-interval <= %g fs, "
                      "mcmd-pb-newton-tol <= %g, mcmd-pb-cells-per-smoothing >= %g)",
                      MCMD_REF_PB_INTERVAL, MCMD_REF_NEWTON_TOL, MCMD_REF_CELLS);
        }
        if (fabs(ir->mcmd_pb_smoothing - MCMD_REF_SMOOTHING) > 1e-6)
        {
            MCMD_WARN("mcmd-pb-smoothing = %g nm (default %g): the ion kinetic energy depends on it\n"
                      "     (charge, escape and directions do not)",
                      ir->mcmd_pb_smoothing, MCMD_REF_SMOOTHING);
        }
        if (!mc->screen_evap)
        {
            MCMD_WARN("mcmd-evaporation = 0: not fully validated");
        }
        if (mc->coll_on && mc->coll_recomb)
        {
            MCMD_WARN("mcmd-coll-recombination = 1: three-body recombination is not fully validated");
        }
        if (mc->coll_on && mc->coll_excite)
        {
            MCMD_WARN("mcmd-coll-excitation = 1: collisional excitation is not fully validated");
        }
        if (!mc->pbg_lr && md->nr > MCMD_LR_AUTO_ATOMS)
        {
            MCMD_WARN("all pairs with %d atoms: O(N^2), about 2x slower than mcmd-pb-long-range = auto "
                      "(GroEL, 55k atoms: 1.9x)", md->nr);
        }
    }
#undef MCMD_WARN
    if (nw > 0)
    {
        fprintf(fp, "\n WARNINGS (%d)\n", nw);
        for (i = 0; i < nw; i++)
        {
            fprintf(fp, "   ! %s\n", warn[i]);
        }
    }
    else
    {
        fprintf(fp, " Settings      all at the validated production values\n");
    }
    fprintf(fp, "%s\n\n", bar);
    fflush(fp);
}

t_mcionize *mcionize_init(FILE *fplog, const t_inputrec *ir,
                          t_mdatoms *mdatoms, const t_commrec *cr,
                          gmx_bool bStartFromCpt, int nHarmonicBonds)
{
    t_mcionize *mc;
    struct stat st;
    double      nphot, rho;
    int         i;

    if (DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS,
                  "The MC/MD ionization module does not support domain "
                  "decomposition: it needs a global pair search over all "
                  "coordinates, and DD renumbers atoms into rank-local "
                  "indices.  Run with -pd for particle decomposition, or "
                  "with -nt 1.");
    }

    snew(mc, 1);

    mc->bChargeTransfer = (ir->mcmd_charge_transfer != 0);
    mc->bAllowHCT       = (ir->mcmd_allow_H_CT != 0);
    mc->bDownhill       = (ir->mcmd_charge_transfer_downhill != 0);
    mc->elem_H          = mcmd_mass2idx(1);
    mc->bInitCharges    = (ir->mcmd_initial_charges != 0);
    mc->bDetailedOutput        = (ir->mcmd_detailed_output != 0);

    {
        /* Moved to the mdp (tpr versions 77 and 78) or removed in the
         * 2026-10-22 and 2026-10-27 clean-ups: refuse them, so an old launcher
         * cannot run a different model without noticing. */
        static const char *gone[][2] = {
            {"GMX_MCMD_SEED", "mcmd-seed"}, {"GMX_MCMD_PBG_NEWTON_TOL", "mcmd-pb-newton-tol"},
            {"GMX_MCMD_PBG_DX_FACTOR", "mcmd-pb-cells-per-smoothing"}, {"GMX_MCMD_PBG_BOXMAX", "mcmd-pb-box-max"},
            {"GMX_MCMD_PBG_NGLIMIT", "mcmd-pb-grid-limit"},
            {"GMX_MCMD_PBG_OUTLAT_LOG", "GMX_MCMD_PBG_OUTLAT_CHECK=1 (it also writes pbg_outpot.txt)"},
            {"GMX_MCMD_PARK_STEP", NULL}, {"GMX_MCMD_PBG_LS", NULL}, {"GMX_MCMD_PBG_CUT_Q", NULL},
            {"GMX_MCMD_PBG_GROW_MARGIN", NULL}, {"GMX_MCMD_PBG_SHRINK_MARGIN", NULL},
            {"GMX_MCMD_PBG_NOSHRINK_EXPAND", NULL}, {"GMX_MCMD_PBG_SHRINK_CUT", NULL}, {"GMX_MCMD_PBG_BOXMIN", NULL},
            {"GMX_MCMD_PBG_OUT_MARGIN", NULL}, {"GMX_MCMD_PBG_NGMAX", NULL}, {"GMX_MCMD_PBG_SFAC", NULL},
            {"GMX_MCMD_SCREEN", "mcmd-electrons"}, {"GMX_MCMD_SCREEN_EVAP", "mcmd-evaporation"},
            {"GMX_MCMD_PB_NSTEP", "mcmd-pb-interval (in fs)"}, {"GMX_MCMD_PBG_SE", "mcmd-pb-smoothing"},
            {"GMX_MCMD_PBG_BOX_CUT", "mcmd-pb-box-cut"}, {"GMX_MCMD_PBG_LR", "mcmd-pb-long-range"},
            {"GMX_MCMD_BANDWIDTH", "mcmd-pulse-bandwidth"}, {"GMX_MCMD_PB_FIELDION", "mcmd-field-ionization"},
            {"GMX_MCMD_PB_FIELDION_L", "mcmd-field-ionization-smoothing"},
            {"GMX_MCMD_CT_H_DONOR", "mcmd-ct-h-donor-only"}, {"GMX_MCMD_COLLISIONAL", "mcmd-collisions"},
            {"GMX_MCMD_COLL_EXCITE", "mcmd-coll-excitation"}, {"GMX_MCMD_COLL_RECOMB", "mcmd-coll-recombination"},
            {"GMX_MCMD_COLL_NO_H_RECOMB", "mcmd-coll-h-recombination (the opposite sense)"},
            {"GMX_MCMD_PBG_OUTLAT", "GMX_MCMD_PBG_OUTLAT_CHECK=1 for the check"},
            {"GMX_MCMD_PBG_N", NULL}, {"GMX_MCMD_PBG_NAUTO", NULL}, {"GMX_MCMD_PBG_NMIN", NULL},
            {"GMX_MCMD_PBG_NMAX", NULL}, {"GMX_MCMD_PBG_DXT", NULL}, {"GMX_MCMD_PBG_BOX", NULL},
            {"GMX_MCMD_PBG_BOX_MIN", NULL}, {"GMX_MCMD_PBG_OUTLAT_D", NULL}, {"GMX_MCMD_PBG_LR_NSTEP", NULL},
            {"GMX_MCMD_PB_ESC_S", NULL}, {"GMX_MCMD_PBG_PCG", NULL}, {"GMX_MCMD_PBG_TSC", NULL},
            {"GMX_MCMD_PB_ETOT", NULL}, {"GMX_MCMD_EVAP_HIST", NULL}, {"GMX_MCMD_COLL_PE", NULL},
            {"GMX_MCMD_GEOM_HEAVY", NULL}, {"GMX_MCMD_ELLIPSOID", NULL}, {"GMX_MCMD_PLASMA", NULL},
            {"GMX_MCMD_PRESSURE", NULL}, {"GMX_MCMD_PRESSURE_H", NULL}, {"GMX_MCMD_PRESSURE_CAP", NULL},
            {"GMX_MCMD_PRESSURE_PDV", NULL}, {"GMX_MCMD_PRESSURE_SYM", NULL}, {"GMX_MCMD_ADIABATIC", NULL},
            {"GMX_MCMD_SCREEN_TE", NULL}, {"GMX_MCMD_SCREEN_GATE", NULL}, {"GMX_MCMD_SCREEN_WS", NULL},
            {"GMX_MCMD_SCREEN_NOFORCE", NULL}, {"GMX_MCMD_WELL_FLAT", NULL}, {"GMX_MCMD_FIELD_ION", NULL},
            {"GMX_MCMD_FIELD_SCALE", NULL}, {"GMX_MCMD_FIELD_STRIDE", NULL}, {"GMX_MCMD_FIELD_SCREENED", NULL},
            {"GMX_MCMD_MC_PAR", NULL}, {"GMX_MCMD_COLL_AVG_R", NULL}
        };
        const char *e;
        int         g;

        for (g = 0; g < (int)(sizeof(gone) / sizeof(gone[0])); g++)
        {
            if (getenv(gone[g][0]) != NULL)
            {
                if (gone[g][1] != NULL)
                {
                    gmx_fatal(FARGS, "%s is no longer read: set %s in the mdp instead (and unset %s).",
                              gone[g][0], gone[g][1], gone[g][0]);
                }
                gmx_fatal(FARGS, "%s was removed (a 2026-10-22 or 2026-10-27 clean-up); unset it.", gone[g][0]);
            }
        }

        mc->seed = mcmd_rng_seed(ir);
        if (fplog != NULL)
        {
            fprintf(fplog, "mcionize: RNG seed %u%s\n", mc->seed, ir->mcmd_seed >= 0 ? " (mcmd-seed)" : " (chosen, mcmd-seed = -1)");
        }

        /* Free electrons: none (plain MolDStruct: every emitted electron
         * leaves) or the self-consistent grid Poisson-Boltzmann cloud. */
        mc->screen_mode = (ir->mcmd_electrons == 1) ? 5 : 0;
        mc->screen_evap = (ir->mcmd_evaporation != 0);
        mc->pb_evh_ver  = -1;
        mc->pb_nstep    = (int)floor(mcmd_mdp_real(ir->mcmd_pb_interval) / (ir->delta_t * 1e3) + 0.5);
        if (mc->pb_nstep < 1)
        {
            mc->pb_nstep = 1;
        }
        /* diagnostics (output only; the physics is the same with or without) */
        e = getenv("GMX_MCMD_TIMING_DETAIL");
        mc->tdet = (e != NULL && atoi(e) != 0);
        e = getenv("GMX_MCMD_PB_EDIAG");
        mc->pb_ediag = (e != NULL && atoi(e) != 0);
        e = getenv("GMX_MCMD_COLL_EVLOG");
        mc->coll_evlog = (e != NULL && atoi(e) != 0);
        e = getenv("GMX_MCMD_PBG_DENSITY");
        mc->pbg_dens_out = (e != NULL && atoi(e) != 0);
        /* Grid: adaptive size (spacing s_e / mcmd-pb-cells-per-smoothing at R0,
         * at least 24 cells per side, at most what the sample needs at the
         * largest box: set from its radius at the first solve, mcmd_pbg_setcap,
         * capped at mcmd-pb-grid-limit) and box (grows when more than pbg_cut
         * of the electrons are outside, shrinks only to where at most
         * MCMD_PBG_SCUT are, MCMD_PBG_BOXMIN..mcmd-pb-box-max R). */
        mc->pbg_nmin  = 24;
        mc->pbg_nmax  = 96;             /* placeholders until the limit is set */
        mc->pbg_ncap  = mc->pbg_nmax;
        mc->pbg_n     = mc->pbg_nmax;   /* until the first rebuild chooses */
        mc->pbg_capset = 0;
        mc->pbg_nglimit = 8 * ((ir->mcmd_pb_grid_limit + 7) / 8);
        mc->pbg_box   = mcmd_mdp_real(ir->mcmd_pb_box_max);
        mc->pbg_cut   = mcmd_mdp_real(ir->mcmd_pb_box_cut);
        /* the adaptive box starts at the minimum; a fixed box (cut 0) is the
         * maximum */
        mc->pbg_hbox  = (mc->pbg_cut > 0) ? MCMD_PBG_BOXMIN : mc->pbg_box;
        mc->pbg_se    = mcmd_mdp_real(ir->mcmd_pb_smoothing);
        mc->pbg_dxt   = mc->pbg_se / mcmd_mdp_real(ir->mcmd_pb_cells_per_smoothing);
        mc->pbg_ntol  = mcmd_mdp_real(ir->mcmd_pb_newton_tol);
        /* the ions outside the box on sub-lattices; GMX_MCMD_PBG_OUTLAT_CHECK=1
         * (diagnostic) checks them against the exact sum (pbg_outlat.txt) and
         * logs the ions and time per stride level (pbg_outpot.txt) */
        e = getenv("GMX_MCMD_PBG_OUTLAT_CHECK");
        mc->pbg_outlat   = (e != NULL && atoi(e) != 0) ? 2 : 1;
        /* GMX_MCMD_PBG_LR_CHECK=1 (diagnostic): check against all pairs */
        e = getenv("GMX_MCMD_PBG_LR_CHECK");
        mc->pbg_lr = (ir->mcmd_pb_long_range == 0) ? 0 : (e != NULL && atoi(e) != 0) ? 2 : 1;
        mc->pbg_lr_nstep = mc->pb_nstep;
        mc->lr_rc = ir->rcoulomb;
        mc->bandwidth_fwhm = mcmd_mdp_real(ir->mcmd_pulse_bandwidth);
        mc->pbfi           = (ir->mcmd_field_ionization != 0);
        if (mc->pbfi)
        {
            mc->pbfi_L = mcmd_mdp_real(ir->mcmd_field_ionization_smoothing);
        }
        mc->ct_h_donor_only = (ir->mcmd_ct_h_donor_only != 0);
        mc->coll_on          = (ir->mcmd_collisions != 0);
        mc->coll_lotz        = (ir->mcmd_coll_data == 1);
        mc->coll_recomb      = (ir->mcmd_coll_recombination != 0);
        mc->coll_no_h_recomb = (ir->mcmd_coll_h_recombination == 0);
        mc->coll_excite      = (ir->mcmd_coll_excitation != 0);
        /* sample geometry (centre, R_g, ellipsoid from the gyration tensor)
         * from the non-hydrogen atoms: fast protons leave the cloud behind */
        mc->geom_heavy = TRUE;
        {
            unsigned int sd = mc->seed ^ 0x9E3779B9u;

            mc->erand_state[0] = 0x330E;
            mc->erand_state[1] = (unsigned short)(sd & 0xFFFF);
            mc->erand_state[2] = (unsigned short)(sd >> 16);
        }

        if (mc->screen_mode != 0 && DOMAINDECOMP(cr))
        {
            gmx_fatal(FARGS, "The electron cloud needs -pd or -nt 1.");
        }
    }

    /* charges_over_time.bin is one real per atom per step, which dwarfs
     * everything else the module writes - 196 MB for a 25000 step run of a
     * 1960 atom protein, and it grows with both.  Writing only every Nth step
     * cuts that by N.  Clamped to at least 1.  The three per-step text files
     * (mean charge, pulse, electron data) follow nstenergy since 2026-10-27
     * (0.1 fs in production; every step they were ~4.5 MB a lysozyme run). */
    mc->charge_stride = (ir->mcmd_charge_output_stride > 0 ?
                         ir->mcmd_charge_output_stride : 1);
    mc->text_stride   = (ir->nstenergy > 0 ? ir->nstenergy : 1);

    /* The module is active, so the altered force field is too.  This used to
     * be a separate mdp switch (userint1); it is now implied by -ionize,
     * which is the only way to reach this function.  Set here rather than in
     * mcionize_step() so it is in place before the first do_force(). */
    moldstruct_altered_ff = 1;

    /* Idle skipping for the charge-transfer pass.  0 means "use the built-in
     * default", which is also what a tpr generated before these two knobs
     * existed carries, so old inputs keep the tuned behaviour instead of
     * silently getting a zero-length window and a modulo by zero.  A negative
     * mcmd-charge-transfer-idle disables the skip and runs the pass on
     * every step. */
    mc->ct_idle_steps = (ir->mcmd_charge_transfer_idle != 0 ?
                         (gmx_large_int_t)ir->mcmd_charge_transfer_idle :
                         (gmx_large_int_t)MCMD_CT_IDLE_STEPS_DEFAULT);
    mc->ct_recheck    = (ir->mcmd_charge_transfer_recheck > 0 ?
                         ir->mcmd_charge_transfer_recheck :
                         MCMD_CT_RECHECK_DEFAULT);

    mc->natoms = mdatoms->nr;

    /* Pulse, converted to SI once here rather than every step. */
    mc->t_mean_s         = ir->mcmd_pulse_peak_time * 1e-15;   /* fs -> s */
    nphot                = (double)ir->mcmd_pulse_photons;
    /* The mdp gives the FWHM; everything below works in the standard
     * deviation, so convert once here. */
    mc->width_s          = (ir->mcmd_pulse_fwhm * 1e-15) / MCMD_FWHM_PER_SIGMA;
    rho                  = ((double)ir->mcmd_pulse_focal_diameter)
                           * 1e-9 * 100;                        /* nm -> cm */
    mc->photon_energy_eV = (double)ir->mcmd_pulse_photon_energy;

    if (nphot > 0)
    {
        mc->imax = (nphot / (MCMD_PI * sqr(rho / 2.0))) /
                   (mc->width_s * sqrt(2.0 * MCMD_PI));
    }
    else
    {
        mc->imax = 0;
    }

    snew(mc->atom_configurations, mc->natoms);
    snew(mc->GS_configurations, mc->natoms);
    snew(mc->energy_of_state, mc->natoms);
    snew(mc->outermost_shell, mc->natoms);
    snew(mc->elem_of_atom, mc->natoms);
    snew(mc->ct_order, mc->natoms);
    snew(mc->cell_of_atom, mc->natoms);
    snew(mc->pos_of_atom, mc->natoms);
    snew(mc->cell_sorted, mc->natoms);
    mc->cell_start   = NULL;
    mc->cell_cursor  = NULL;
    mc->cell_nalloc  = 0;
    for (i = 0; i < mc->natoms; i++)
    {
        snew(mc->atom_configurations[i], MCMD_NSHELL);
        snew(mc->GS_configurations[i], MCMD_NSHELL);
    }

    /* Only the master writes output: every rank runs this function, but the
     * physics and all the I/O happen on the master alone. */
    if (MASTER(cr))
    {
        /* continuation: the state of the checkpoint, and the output files cut
         * back to it, before they are opened for appending */
        if (bStartFromCpt)
        {
            mcmd_state_load(mc, ir);
        }
        if (stat(MCMD_OUTDIR, &st) == -1)
        {
            if (mkdir(MCMD_OUTDIR, 0700) != 0)
            {
                gmx_fatal(FARGS, "Could not create the '%s' directory.",
                          MCMD_OUTDIR);
            }
        }

        /* Only the detailed outputs are opened here; without them the module
         * writes nothing per step.  Opening with "w" also truncates, so a
         * fresh run does not append to the output of the last one. */
        if (mc->bDetailedOutput)
        {
            const char *tmode = bStartFromCpt ? "a" : "w";

            mc->fcharges  = mcmd_open_buffered("charges_over_time.bin",
                                               bStartFromCpt ? "ab" : "wb",
                                               &mc->outbuf[0]);
            mc->fmean     = mcmd_open_buffered("mean_charge_vs_time.txt",
                                               tmode, &mc->outbuf[1]);
            mc->fpulse    = mcmd_open_buffered("pulse_profile.txt",
                                               tmode, &mc->outbuf[2]);
            mc->felectron = mcmd_open_buffered("electron_data.txt",
                                               tmode, &mc->outbuf[3]);
            mc->ftrans    = mcmd_open_buffered("transitions.bin",
                                               bStartFromCpt ? "ab" : "wb", &mc->outbuf[4]);
            if (!bStartFromCpt)
            {
                mcmd_trn_header(mc->ftrans);
            }
        }

        {
            FILE *fmasses = mcmd_open_output("masses.bin", "wb");

            fwrite(mdatoms->massT, sizeof(mdatoms->massT[0]), mdatoms->nr,
                   fmasses);
            fclose(fmasses);
        }
    }


    /* Load atomic data only for the elements actually present.  Loading all
     * MCMD_ATOM_TABLE_SIZE slots was what produced a page of
     * "rate_transitions__.txt" errors on every run. */
    snew(mc->atomdata, MCMD_ATOM_TABLE_SIZE);
    for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
    {
        mc->atomdata[i].bPresent = FALSE;
    }
    for (i = 0; i < mdatoms->nr; i++)
    {
        int idx = mcmd_atom2idx(mdatoms->massT[i]);

        if (idx < 0)
        {
            gmx_fatal(FARGS,
                      "Atom %d has mass %g, which is not in the element "
                      "table.  Add it to mcmd_elements[] in atomic_data.c.",
                      i + 1, mdatoms->massT[i]);
        }
        if (!mc->atomdata[idx].bPresent)
        {
            mcmd_atomdata_read(&mc->atomdata[idx], idx, !mc->coll_on ? 0 : (mc->coll_lotz ? 2 : 1),
                               ir->mcmd_atomic_data, ir->mcmd_pulse_photon_energy);
            if (fplog != NULL)
            {
                fprintf(fplog, "mcionize: %s atomic data for %s\n",
                        ir->mcmd_atomic_data == 1 ? "computed" : "loaded", mcmd_idx2symbol(idx));
            }
        }
    }

    /* The built-in data as files, so that a run can be inspected, and repeated with mcmd-atomic-data = files */
    if (ir->mcmd_atomic_data == 1 && MASTER(cr))
    {
        const char *dir = MCMD_OUTDIR "/Atomic_data";

        if (stat(dir, &st) == -1 && mkdir(dir, 0700) != 0)
        {
            gmx_fatal(FARGS, "Could not create the '%s' directory.", dir);
        }
        for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
        {
            if (mc->atomdata[i].bPresent)
            {
                mcmd_atomdata_write(&mc->atomdata[i], dir);
            }
        }
    }

    /* Charge transfer needs to price an electron entering any shell of the
     * acceptor, which only total energies can do.  Priced from
     * energy_levels_X.txt instead, an acceptor holding an inner-shell hole
     * cannot be scored and is refused - at hard photon energies that is most
     * of them, and transfer is suppressed by more than an order of magnitude
     * while every number still looks reasonable.  A silent result that wrong
     * is worse than not running, so this is fatal rather than a warning.
     *
     * Only when transfer is actually on: a pure ionization run never reads
     * the table. */
    if (mc->bChargeTransfer)
    {
        char missing[MCMD_NUM_ELEMENTS * 4];
        int  nmissing = 0;

        missing[0] = '\0';

        for (i = 0; i < MCMD_NUM_ELEMENTS; i++)
        {
            if (mc->atomdata[i].bPresent &&
                mc->atomdata[i].total_energies == NULL)
            {
                if (nmissing > 0)
                {
                    strcat(missing, ", ");
                }
                strcat(missing, mcmd_idx2symbol(i));
                nmissing++;
            }
        }

        if (nmissing > 0)
        {
            gmx_fatal(FARGS,
                      "Charge transfer is on but Atomic_data has no "
                      "total_energies_X.txt for %s.  Regenerate the atomic "
                      "data (see ATOMIC_DATA.md) or set "
                      "mcmd-charge-transfer = 0.",
                      missing);
        }
    }

    /* Masses are fixed for the run, so resolve every atom's element once here
     * and let the hot paths index elem_of_atom[] instead. */
    for (i = 0; i < mdatoms->nr; i++)
    {
        mc->elem_of_atom[i] = mcmd_resolve_element(mc, mdatoms, i);
    }

    /* Size the Monte Carlo scratch to the worst case over the elements
     * actually present, once, so the kernel never allocates. */
    mc->mc_scratch_size = 1;
    for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
    {
        if (mc->atomdata[i].bPresent &&
            mc->atomdata[i].max_transitions > mc->mc_scratch_size)
        {
            mc->mc_scratch_size = mc->atomdata[i].max_transitions;
        }
    }
    snew(mc->mc_transitions, mc->mc_scratch_size);
    snew(mc->mc_dt, mc->mc_scratch_size);
    if (mc->coll_on)
    {
        mcmd_coll_build(mc, fplog);
    }

    if (fplog != NULL && mc->bDetailedOutput)
    {
        fprintf(fplog,
                "mcionize: detailed output on; charges_over_time.bin written "
                "every %d steps, the text files every step\n",
                mc->charge_stride);
    }

    if (fplog != NULL && mc->bChargeTransfer)
    {
        if (mc->ct_idle_steps < 0)
        {
            fprintf(fplog, "mcionize: charge transfer runs on every step "
                    "(idle skipping disabled by mcmd-charge-transfer-idle)\n");
        }
        else
        {
            fprintf(fplog, "mcionize: charge transfer drops to every %d steps "
                    "after %d steps without a transfer\n",
                    mc->ct_recheck, (int)mc->ct_idle_steps);
        }
    }

    mcmd_srand(mc->seed);

    if (bStartFromCpt)
    {
        /* Continuation: recover the electronic state rather than resetting to
         * ground state.  The reference radius of gyration in particular must
         * come from the file - recomputing it from an already-expanded sample
         * would reset rg_factor to 1 and shrink the plasma volume back to its
         * t=0 value. */
        if (MASTER(cr))
        {
            mcmd_state_apply(mc, mdatoms);
        }
        if (PAR(cr))
        {
            gmx_bcast(mdatoms->nr * sizeof(mdatoms->chargeA[0]),
                      mdatoms->chargeA, cr);
            /* every rank schedules the PB solves by it */
            gmx_bcast(sizeof(mc->screen_step), &mc->screen_step, cr);
        }
        mc->bInitialised = TRUE;

        if (fplog != NULL)
        {
            fprintf(fplog,
                    "mcionize: restarted at step %lld (%.0f trapped electrons, "
                    "%.0f escaped)\n", (long long)mc->screen_step,
                    mc->num_free_electrons, mc->escaped_photoelectrons);
        }
    }
    else
    {
        mc->bInitialised = FALSE;
    }

    if (MASTER(cr))
    {
        if (fplog != NULL)
        {
            mcmd_print_summary(fplog, mc, ir, mdatoms, cr, bStartFromCpt, nHarmonicBonds);
        }
        mcmd_print_summary(stderr, mc, ir, mdatoms, cr, bStartFromCpt, nHarmonicBonds);
    }

    return mc;
}

/* Energy ledger: at the end of the master's part of the step, the exact energy
 * jumps of the charge changes made in it (photo, Auger, collisional, field
 * ionization, charge transfer), at this step's positions, all pairs:
 *   ion-ion       dq_k k sum_{j != k} q_j / r_kj     (atoms taken one by one,
 *                                                     so pairs that both changed
 *                                                     are counted exactly)
 *   ion-electron  -dq_k * (electron potential at k): pb the clouds
 *                 k sum_c w_c G(r, sqrt2 h), pbgrid sc K n at x_k.
 * GROMACS's total energy then changes by W (the electron field's work) + J_ii
 * up to integration error and charge-dependent bonded terms. */
static void mcmd_ledger_start(t_mcionize *mc, const t_mdatoms *md)
{
    int a;

    if (mc->led_q == NULL)
    {
        snew(mc->led_q, md->nr);
    }
    for (a = 0; a < md->nr; a++) { mc->led_q[a] = md->chargeA[a]; }
}

static void mcmd_ledger_step(t_mcionize *mc, const t_mdatoms *md, rvec *x)
{
    const double sc = (mc->pb_ready && mc->pb_N > 0) ? mc->num_free_electrons / mc->pb_N : 0.0;
    int          k, j;

    if (mc->led_q == NULL)
    {
        return;
    }
    for (k = 0; k < md->nr; k++)
    {
        const double dq = md->chargeA[k] - mc->led_q[k];
        double       u  = 0;

        if (dq == 0)
        {
            continue;
        }
        for (j = 0; j < md->nr; j++)
        {
            if (j != k && mc->led_q[j] != 0)
            {
                u += mc->led_q[j] / sqrt(mcmd_pb_d2(x[k], x[j]));
            }
        }
        mc->led_Jii += dq * MCMD_COULOMB_EV_NM * u;
        if (sc > 0 && mc->pbg_Ve != NULL)
        {
            double grad[3];

            {
                const double dj = -dq * sc * mcmd_pbg_Ve(mc, x[k], grad);

                mc->led_Jie += dj;
                mc->pb_Etot += dj;
                mc->pbd_J   += dj;
            }
        }
        mc->led_q[k]  = md->chargeA[k];
        mc->led_nchg += 1;
    }
}

static void mcmd_bcast_plasma(t_mcionize *mc, const t_commrec *cr)
{
        double buf[21];
        int    i, j;

        buf[0] = mc->num_free_electrons; buf[1] = mc->escaped_photoelectrons;
        buf[2] = mc->screen_E_trapped;   buf[3] = mc->radius_now;
        buf[4] = mc->com[0]; buf[5] = mc->com[1]; buf[6] = mc->com[2];
        /* the ellipsoid too: the PB solve, which runs on every rank, picks its
         * sites by ellipsoidal distance */
        for (i = 0; i < 3; i++)
        {
            buf[7 + i] = mc->ell_a[i];
            for (j = 0; j < 3; j++) { buf[10 + 3 * i + j] = mc->ell_axis[i][j]; }
        }
        buf[19] = mc->pb_Etot;
        buf[20] = mc->pb_etot_init;
        gmx_bcast(sizeof(buf), buf, cr);
        mc->pb_Etot      = buf[19];
        mc->pb_etot_init = (gmx_bool)buf[20];
        mc->num_free_electrons = buf[0]; mc->escaped_photoelectrons = buf[1];
        mc->screen_E_trapped   = buf[2]; mc->radius_now = buf[3];
        mc->com[0] = buf[4]; mc->com[1] = buf[5]; mc->com[2] = buf[6];
        for (i = 0; i < 3; i++)
        {
            mc->ell_a[i] = buf[7 + i];
            for (j = 0; j < 3; j++) { mc->ell_axis[i][j] = buf[10 + 3 * i + j]; }
        }
}

/* Evaporation tail at T over sites w (electrons), V (potential energy, eV),
 * the n entries of list (or 0..n-1 without one): the Maxwell tail above each
 * site's depth, Q(3/2, D/T) of the electrons, Q(5/2, D/T) of their energy,
 * and their depth. */
static void mcmd_evap_tail(const double *w, const double *V, const int *list, int n, double T,
                           double *dN, double *dE, double *dND)
{
    int ib;

    *dN = *dE = *dND = 0;
    for (ib = 0; ib < n; ib++)
    {
        const int b = (list != NULL) ? list[ib] : ib;
        double    xb, ex3, q32;

        if (w[b] <= 0) { continue; }
        xb  = -V[b] / T;
        if (xb < 0) { xb = 0; }
        ex3 = exp(-xb);
        q32 = erfc(sqrt(xb)) + 2.0 * sqrt(xb / M_PI) * ex3;
        *dN  += w[b] * q32;
        *dND += w[b] * q32 * (-V[b]);
        *dE  += w[b] * 1.5 * T * (q32 + 4.0 * xb * sqrt(xb) * ex3 / (3.0 * sqrt(M_PI)));
    }
}

#define MCMD_EVH_NB 4000
/* The sites binned by depth: bin 0 holds depth <= 0 (potential set to 0, as
 * the tail clamps it), bins 1.. MCMD_EVH_NB the depths up to the deepest. */
static void mcmd_evap_hist(t_mcionize *mc)
{
    const int n = (mc->pb_evlist != NULL) ? mc->pb_evnl : mc->pb_evn;
    double    Dmax = 0;
    int       ib, k;

    if (mc->pb_evhw == NULL) { snew(mc->pb_evhw, MCMD_EVH_NB + 1); snew(mc->pb_evhV, MCMD_EVH_NB + 1); }
    for (k = 0; k <= MCMD_EVH_NB; k++) { mc->pb_evhw[k] = mc->pb_evhV[k] = 0; }
    for (ib = 0; ib < n; ib++)
    {
        const int b = (mc->pb_evlist != NULL) ? mc->pb_evlist[ib] : ib;

        if (mc->pb_evw[b] > 0 && -mc->pb_evV[b] > Dmax) { Dmax = -mc->pb_evV[b]; }
    }
    for (ib = 0; ib < n; ib++)
    {
        const int    b = (mc->pb_evlist != NULL) ? mc->pb_evlist[ib] : ib;
        const double D = -mc->pb_evV[b];

        if (mc->pb_evw[b] <= 0) { continue; }
        if (D <= 0)
        {
            mc->pb_evhw[0] += mc->pb_evw[b];
            continue;
        }
        k = 1 + (int)(D / Dmax * MCMD_EVH_NB);
        if (k > MCMD_EVH_NB) { k = MCMD_EVH_NB; }
        mc->pb_evhw[k] += mc->pb_evw[b];
        mc->pb_evhV[k] += mc->pb_evw[b] * mc->pb_evV[b];
    }
    for (k = 1; k <= MCMD_EVH_NB; k++)
    {
        if (mc->pb_evhw[k] > 0) { mc->pb_evhV[k] /= mc->pb_evhw[k]; }
    }
    mc->pb_evh_ver = mc->pb_nsolve;
}

void mcionize_step(t_mcionize *mc, FILE *fplog, const t_inputrec *ir,
                   t_mdatoms *mdatoms, t_state *state,
                   double t, const t_commrec *cr, t_nrnb *nrnb)
{
    if (mc == NULL)
    {
        return;
    }

    mc->intensity = mc->imax *
                    exp(-0.5 * sqr((t * 1e-12 - mc->t_mean_s) / mc->width_s));

    /* Under particle decomposition the integrator has only advanced this
     * rank's home atoms, so the rest of state->x is a step stale.  The charge
     * transfer pair search, the PB solve and the geometry need all of them, so
     * gather first, over the whole ring.  move_x() (what do_force() uses) only
     * runs pd_shift pulses, enough for the half-shell pair list: for 1B7U on 6
     * ranks that is 4, and the master's copy of the last block was then only
     * refreshed when the trajectory was written (nstxout).  GMX_MCMD_PD_CHECK
     * checks it. */
    if (PAR(cr))
    {
        move_rvecs(cr, FALSE, FALSE, GMX_LEFT, GMX_RIGHT, state->x, NULL, cr->nnodes - 1, nrnb);
    }

    /* GMX_MCMD_PD_CHECK=1 (diagnostic): every 50 steps, each rank's sum of
     * its home coordinates against the master's copy of that block */
    if (PAR(cr) && getenv("GMX_MCMD_PD_CHECK") != NULL && mc->screen_step % 50 == 0)
    {
        double s[4 * 64] = {0}, own[3];
        int    nn = cr->nnodes, r, a, d;
        int   *ix = pd_index(cr);

        for (a = mdatoms->start; a < mdatoms->start + mdatoms->homenr; a++)
        {
            for (d = 0; d < 3; d++) { s[4 * cr->nodeid + d] += state->x[a][d]; }
        }
        s[4 * cr->nodeid + 3] = 1;
        gmx_sumd(4 * nn, s, cr);
        if (MASTER(cr))
        {
            fprintf(stderr, "PD_CHECK step %d shift %d bshift %d:", (int)mc->screen_step, pd_shift(cr), pd_bshift(cr));
            for (r = 0; r < nn; r++)
            {
                own[0] = own[1] = own[2] = 0;
                for (a = ix[r]; a < ix[r + 1]; a++)
                {
                    for (d = 0; d < 3; d++) { own[d] += state->x[a][d]; }
                }
                fprintf(stderr, "  [%d: %d..%d] %.3e", r, ix[r], ix[r + 1],
                        fabs(own[0] - s[4 * r]) + fabs(own[1] - s[4 * r + 1]) + fabs(own[2] - s[4 * r + 2]));
            }
            fprintf(stderr, "\n");
        }
    }

    /* The module itself is serial, on the master.  It is ~3% of the step, and
     * the electron escape test reads escaped_photoelectrons back inside the
     * same atom loop that updates it, so splitting the loop across ranks
     * would change which electrons are trapped.  Keeping it on one rank
     * preserves the serial physics exactly and still lets the force loop -
     * which is the actual cost - use every rank. */
    if (MASTER(cr))
    {
        double tm_step0 = mcmd_pb_wallms(), tm_x;

        if (mc->bInitialised)
        {
            mcmd_ledger_start(mc, mdatoms);
        }

        if (!mc->bInitialised)
        {
            /* Keyed off the flag rather than t == 0: on a continuation tinit
             * is non-zero, and the old test meant the module was never
             * initialised at all. */
            mcmd_first_step(mc, mdatoms, state);
            mc->bInitialised = TRUE;
        }
        else if (mc->bChargeTransfer)
        {
            mc->ct_step++;

            if (mc->ct_idle_steps < 0 ||
                mc->ct_steps_since_transfer < mc->ct_idle_steps ||
                (mc->ct_step % mc->ct_recheck) == 0)
            {
                int before = mc->count_chargetransfer;

                tm_x = mcmd_pb_wallms();
                mcmd_charge_transfer(mc, mdatoms, state, t);
                mc->tm[MCMD_TM_CT] += mcmd_pb_wallms() - tm_x;

                if (mc->count_chargetransfer > before)
                {
                    mc->ct_steps_since_transfer = 0;
                }
                else
                {
                    mc->ct_steps_since_transfer++;
                }
            }
            else
            {
                int r;

                /* Draw and discard exactly what the skipped pass would have
                 * consumed (mcmd_shuffle takes nr-1 values), so the random
                 * stream is the same as in a run that never skipped.  Without
                 * this the skip would shift every later Monte Carlo draw and
                 * behave like a different seed, which would also make it
                 * impossible to tell a genuinely missed transfer from stream
                 * drift when checking the optimisation. */
                for (r = mdatoms->nr - 1; r > 0; r--)
                {
                    (void)rand();
                }

                mc->ct_steps_since_transfer++;
                mc->ct_skipped++;
            }
        }

        tm_x = mcmd_pb_wallms();
        mcmd_update_volume(mc, mdatoms, state);

        /* Plasma update: evaporation of the trapped electrons. */
        if (mc->screen_mode != 0 && mc->num_free_electrons > 0)
        {
            const double Q = mc->escaped_photoelectrons, R = mc->radius_now;

            if (Q > 0)
            {
                double T, n_cm3, lnL, nu, frac;

                /* At each cell the Maxwell tail above its depth leaves; the
                 * tail refills at the e-e collision rate (NRL) of the mean
                 * density, so per step only a fraction nu*dt of it goes. */
                const double tev0 = mc->tdet ? mcmd_pb_wallms() : 0;

                if (mc->screen_evap && mc->pb_ready && mc->pb_evw != NULL)
                {
                    /* the tail above each depth, Q(3/2, D/T) of the electrons,
                     * Q(5/2, D/T) of their energy, from the histogram of the
                     * cells' depths (rebuilt after every solve) */
                    double dN = 0, dE = 0, dND = 0, sc = mc->num_free_electrons / mc->pb_N;

                    T = mcmd_plasma_Te(mc);
                    if (T > 0)
                    {
                        if (mc->pb_evh_ver != mc->pb_nsolve)
                        {
                            mcmd_evap_hist(mc);
                        }
                        mcmd_evap_tail(mc->pb_evhw, mc->pb_evhV, NULL, MCMD_EVH_NB + 1, T, &dN, &dE, &dND);
                        n_cm3 = mc->num_free_electrons / (4.0 / 3.0 * M_PI * R * R * R * 1e-21);
                        lnL   = 24.0 - log(sqrt(n_cm3) / T);
                        if (lnL < 2.0) { lnL = 2.0; }
                        nu    = 2.91e-6 * n_cm3 * lnL * pow(T, -1.5);
                        frac  = nu * ir->delta_t * 1e-12;
                        if (frac > 1.0) { frac = 1.0; }
                        dN *= sc * frac; dE *= sc * frac; dND *= sc * frac;
                        if (dN > mc->num_free_electrons) { dN = mc->num_free_electrons; }
                        if (dE > mc->screen_E_trapped)   { dE = mc->screen_E_trapped; }
                        mc->pb_Etot                -= dE - dND;    /* kinetic energy minus depth of the leavers */
                        mc->pbd_N                  += dND;
                        mc->num_free_electrons     -= dN;
                        mc->screen_E_trapped       -= dE;
                        mc->eb_evap                += dE;
                        mc->escaped_photoelectrons += dN;
                        mc->screen_evaporated      += dN;
                    }
                    if (mc->tdet) { mc->td[MCMD_TD_EVAP] += mcmd_pb_wallms() - tev0; }
                }
            }
        }

        if (mc->pbfi)
        {
            mcmd_pbfi_apply(mc, mdatoms, state, t);
        }

        if (mc->coll_on)
        {
            mcmd_coll_update_rates(mc);
            if (mc->fcoll == NULL)
            {
                mc->fcoll = mcmd_open_output("collisional_log.txt", "w");
                fprintf(mc->fcoll, "# t[ps] T_e[eV] n_e[cm^-3] N_trapped E_plasma[eV] "
                        "n_coll_ion n_coll_rec energy_to_ionization[eV] E_in_trapped[eV] "
                        "E_adiabatic[eV] E_evaporated[eV] R[nm] Q_escaped E_coll_PE[eV]\n");
            }
            if (mc->coll_step++ % 100 == 0)
            {
                fprintf(mc->fcoll, "%g %.3f %.4g %.0f %.1f %.0f %.0f %.1f %.1f %.1f %.1f %.4f %.1f %.1f\n", t,
                        mc->coll_T, mc->coll_ne, mc->num_free_electrons, mc->screen_E_trapped,
                        mc->count_coll_ion, mc->count_coll_rec, mc->coll_energy, mc->eb_dep,
                        mc->eb_adiab, mc->eb_evap, mc->radius_now, mc->escaped_photoelectrons,
                        mc->coll_pe_energy);
            }
        }
        mc->tm[MCMD_TM_PLASMA] += mcmd_pb_wallms() - tm_x;
        tm_x                    = mcmd_pb_wallms();
        mcmd_monte_carlo(mc, ir, mdatoms, state, t, 0, mdatoms->nr);
        mc->tm[MCMD_TM_MC]     += mcmd_pb_wallms() - tm_x;
        tm_x                    = mcmd_pb_wallms();
        mcmd_write_step_output(mc, mdatoms, t);
        mc->tm[MCMD_TM_OUT]    += mcmd_pb_wallms() - tm_x;
        mcmd_ledger_step(mc, mdatoms, state->x);
        mc->tm[MCMD_TM_STEP]   += mcmd_pb_wallms() - tm_step0;
    }

    /* Screening: the trapped-electron state lives on the master, but
     * under -pd every rank corrects the forces on its own home atoms. */
    if (PAR(cr) && mc->screen_mode != 0)
    {
        mcmd_bcast_plasma(mc, cr);
    }

    /* Publish the new charges before do_force() runs.  Without this every
     * pair spanning two ranks would be evaluated against a stale charge, in
     * both the forces and the energies. */
    if (PAR(cr))
    {
        gmx_bcast(mdatoms->nr * sizeof(mdatoms->chargeA[0]),
                  mdatoms->chargeA, cr);
    }
}

/* ---------------------------------------------------------------- *
 * The forces of the trapped electrons on the ions.
 *
 * The MD sees bare ion charges, as if every emitted electron had left.  Here
 * the forces do_force() just produced get the field of the trapped electron
 * cloud of the last PB solve (mcmd_pbg_forces), and with mcmd-pb-long-range the
 * long-range ion-ion force beyond rcoulomb (mcmd_pbg_lr).  Also writes the
 * energy ledger and screening_log.txt.
 * ---------------------------------------------------------------- */
void mcionize_screen_forces(t_mcionize *mc, const t_inputrec *ir,
                            t_mdatoms *mdatoms, rvec *x, rvec *v, rvec *f,
                            const t_commrec *cr)
{
    const double N_tr    = (mc != NULL) ? mc->num_free_electrons : 0;
    double       f_trap  = 0;
    int          applied = 0;
    double       tm_sf, tm_x;

    if (mc == NULL || mc->screen_mode == 0)
    {
        return;
    }
    mc->screen_step++;
    tm_sf = mcmd_pb_wallms();

    if (N_tr > 0)
    {
        /* lambda_D of the mean density, for the log only */
        double R  = mc->radius_now, Te = mcmd_plasma_Te(mc);
        double ne = N_tr / (4.0 / 3.0 * M_PI * R * R * R);
        double lambda;

        f_trap               = N_tr / (N_tr + mc->escaped_photoelectrons);
        lambda               = (Te > 0) ? 0.23508 * sqrt(Te / ne) : 0;
        mc->screen_ND        = ne * 4.0 / 3.0 * M_PI * lambda * lambda * lambda;
        mc->screen_lambda_nm = lambda;
        applied              = 1;
    }
    tm_x = mcmd_pb_wallms();
    mcmd_pbg_forces(mc, ir, mdatoms, x, v, f, cr);
    if (mc->pbg_lr)
    {
        mcmd_pbg_lr(mc, ir, mdatoms, x, v, f, cr);
    }
    mc->pressure_h_now      = mc->pb_h;
    mc->tm[MCMD_TM_PBF]    += mcmd_pb_wallms() - tm_x;
    mc->tm[MCMD_TM_SCREEN] += mcmd_pb_wallms() - tm_sf;

    /* energy ledger, every 100 steps, on GROMACS's energy steps (nstenergy) */
    if (MASTER(cr) && ((mc->screen_step - 1) % 100) == 0)
    {
        if (mc->fled == NULL)
        {
            mc->fled = mcmd_open_output("energy_ledger.txt", "w");
            fprintf(mc->fled, "# step t[ps] J_ii[eV] J_ie[eV] W[eV] E_kin_e[eV] N_trapped charge_changes "
                    "(cumulative: ion-ion and ion-electron energy jumps of charge changes, electron-field work on "
                    "the ions)\n");
        }
        fprintf(mc->fled, "%ld %.5f %.3f %.3f %.3f %.3f %.1f %.0f\n", (long)(mc->screen_step - 1),
                (mc->screen_step - 1) * ir->delta_t + ir->init_t, mc->led_Jii, mc->led_Jie, mc->pressure_work,
                mc->screen_E_trapped, mc->num_free_electrons, mc->led_nchg);
        fflush(mc->fled);
    }

    if (MASTER(cr) && (mc->screen_step % 100) == 1)
    {
        if (mc->fscreen == NULL)
        {
            mc->fscreen = mcmd_open_output("screening_log.txt", "w");
            fprintf(mc->fscreen, "# step escaped trapped trapped_frac mean_E_trapped[eV] "
                    "R[nm] lambda_D[nm] N_Debye applied evaporated T_e[eV] pressure_work[eV] "
                    "pressure_h[nm]\n");
        }
        fprintf(mc->fscreen, "%ld %.0f %.0f %.4f %.2f %.4f %.5f %.4g %d %.0f %.2f %.1f %.4f\n",
                (long)mc->screen_step, mc->escaped_photoelectrons, N_tr, f_trap,
                N_tr > 0 ? mc->screen_E_trapped / N_tr : 0.0, mc->radius_now,
                mc->screen_lambda_nm, mc->screen_ND, applied, mc->screen_evaporated,
                mcmd_plasma_Te(mc), mc->pressure_work, mc->pressure_h_now);
    }
}

gmx_bool mcionize_autostop(const t_mcionize *mc, const t_inputrec *ir,
                           double E_kin, double E_tot, gmx_large_int_t step,
                           gmx_bool bEnergiesGlobal, const t_commrec *cr)
{
    int stop = 0;

    /* mc == NULL already means -ionize was not given, which is what gates the
     * whole module.  All three of these tests are rank-independent, so every
     * rank returns here together and none is left behind at the broadcast. */
    if (mc == NULL || ir->mcmd_autostop == 0)
    {
        return FALSE;
    }

    /* enerd->term[] is only summed across ranks on global-communication
     * steps.  On any other step it carries this rank's own partial energies,
     * and the ratio computed from them is not the system's - under particle
     * decomposition that made the test fire on one rank and not the others. */
    if (!bEnergiesGlobal)
    {
        return FALSE;
    }

    if (MASTER(cr))
    {
        double ratio = E_kin / E_tot;

        if (ratio > ir->mcmd_autostop_threshold && step > 10000)
        {
            fprintf(stderr,
                    "\nSimulation terminated: ratio %f is over threshold %f\n",
                    ratio, (double)ir->mcmd_autostop_threshold);
            stop = 1;
        }
    }

    /* Decide once and tell everyone.  The threshold is crossed by a hair
     * (0.990001 against 0.990000 in the shipped example), so even identical
     * inputs cannot be relied on to give every rank the same answer. */
    if (PAR(cr))
    {
        gmx_bcast(sizeof(stop), &stop, cr);
    }

    return (stop != 0);
}

void mcionize_checkpoint(t_mcionize *mc, const t_mdatoms *mdatoms, gmx_large_int_t step, const t_commrec *cr)
{
    if (mc == NULL || !MASTER(cr) || !mc->bInitialised)
    {
        return;
    }
    if (mc->screen_mode != 0 && mc->screen_step != step)
    {
        gmx_fatal(FARGS, "mcionize: the module is at step %lld at the checkpoint of step %lld (a run that did not "
                  "start at step 0?); its state could not be resumed.", (long long)mc->screen_step, (long long)step);
    }
    mc->screen_step = step;
    mcmd_state_write(mc, mdatoms);
}

void mcionize_done(t_mcionize *mc, t_mdatoms *mdatoms, const t_commrec *cr)
{
    FILE *fp;
    int   i;

    if (mc == NULL)
    {
        return;
    }

    if (MASTER(cr))
    {
        static const char *nm[10] = {"charge transfer", "geometry + plasma update + rates", "Monte Carlo",
                                     "per-step output", "mcionize_step (master, total)", "PB forces (solves incl.)",
                                     NULL, NULL, "mcionize_screen_forces (total)", "PB solves"};
        FILE *ft = mcmd_open_output("timing_mcionize.txt", "w");
        int   k;

        fprintf(stderr, "\nmcionize: wall time per module section (master rank):\n");
        for (k = 0; k < 10; k++)
        {
            if (nm[k] == NULL)
            {
                continue;    /* 6, 7: the removed Debye and pressure forces */
            }
            fprintf(stderr, "  %-36s %10.1f s\n", nm[k], 1e-3 * mc->tm[k]);
            if (ft) { fprintf(ft, "%-36s %10.1f s\n", nm[k], 1e-3 * mc->tm[k]); }
        }
        if (mc->tdet && ft)
        {
            fprintf(ft, "-- detail (GMX_MCMD_TIMING_DETAIL)\n");
            fprintf(ft, "%-36s %10.1f s\n", "grid convolutions, total", 1e-3 * mc->td[MCMD_TD_CONV]);
            fprintf(ft, "%-36s %10.1f s\n", "  of which FFTs", 1e-3 * mc->td[MCMD_TD_FFT]);
            fprintf(ft, "%-36s %10.1f s\n", "grid CG, total", 1e-3 * mc->td[MCMD_TD_CG]);
            fprintf(ft, "%-36s %10.1f s\n", "  of which convolutions", 1e-3 * mc->td[MCMD_TD_CGCONV]);
            fprintf(ft, "%-36s %10.0f\n", "  convolutions for the 2nd RHS", mc->td_rhs2);
            fprintf(ft, "%-36s %10.1f s\n", "grid Newton loop, total", 1e-3 * mc->td[MCMD_TD_NEWTON]);
            fprintf(ft, "%-36s %10.0f of %.0f\n", "solves already converged at start", mc->td_nskip, mc->td_nsolve);
            fprintf(ft, "%-36s %10.3e vs %.3e\n", "sum of start residuals: now vs extrap", mc->td_newton0,
                    mc->td_newtonex);
            fprintf(ft, "%-36s %10.1f s\n", "extrapolation test (extra cost)", 1e-3 * mc->td[MCMD_TD_EXTRAP]);
            fprintf(ft, "%-36s %10.1f s\n", "charge transfer: cell-list build", 1e-3 * mc->td[MCMD_TD_CTBUILD]);
            fprintf(ft, "%-36s %10.1f s\n", "charge transfer: shuffle", 1e-3 * mc->td[MCMD_TD_CTSHUF]);
            fprintf(ft, "%-36s %10.1f s\n", "evaporation", 1e-3 * mc->td[MCMD_TD_EVAP]);
            fprintf(ft, "%-36s %10.1f s\n", "potential of ions outside the box", 1e-3 * mc->td[MCMD_TD_OUTPOT]);
        }
        if (ft) { fclose(ft); }
    }

    if (!MASTER(cr))
    {
        mcionize_free(mc);
        return;
    }

    /* The state for a continuation is written at the checkpoints
     * (mcionize_checkpoint), not here: the last one is at the start of the
     * last step, which is where mdrun -cpi resumes. */
    if (mc->nstate == 0)
    {
        fprintf(stderr, "\nmcionize: no checkpoint was written, so this run cannot be continued\n");
    }

    if (mc->ct_skipped > 0)
    {
        fprintf(stderr,
                "\nmcionize: charge transfer went idle; the pass ran on %.0f%% "
                "of steps (sampled every %d once idle).\n",
                100.0 * (double)(mc->ct_step - mc->ct_skipped) /
                (double)mc->ct_step, mc->ct_recheck);
    }

    fprintf(stderr,
            "\nmcionize: %.0f electrons escaped, %.0f trapped in the sample "
            "(final radius %.2f nm).\n",
            mc->escaped_photoelectrons, mc->num_free_electrons, mc->radius_now);
    if (mc->ctab_clamped > 0)
    {
        fprintf(stderr, "mcionize: %.0f steps with T_e outside a collisional table's grid (coefficients held at "
                "the edge)\n", mc->ctab_clamped);
    }
    if (mc->count_escape_fallback > 0)
    {
        fprintf(stderr,
                "mcionize: %.0f electron energies had no total_energies_X.txt "
                "entry and used the outer-shell fallback.\n",
                mc->count_escape_fallback);
    }

    if (mc->max_transfer_separation > 0.8 * MCMD_CT_CUTOFF_NM)
    {
        fprintf(stderr,
                "\nWARNING: a charge transfer used a pair %.3f nm apart, "
                "against a search cutoff of %.3f nm.  Transfers may be being "
                "lost to the cutoff; raise MCMD_CT_CUTOFF_NM in mcionize.c.\n",
                mc->max_transfer_separation, (double)MCMD_CT_CUTOFF_NM);
    }

    /* "w", not "a".  Both files describe the whole run rather than a segment
     * of it - the event counters are saved into mcionize_state.bin and
     * restored by -cpi, and the charges are the final ones either way - so a
     * second run in the same directory has nothing to add, only a stale copy
     * to leave behind.  Appending also ran the counts of one run into the
     * header of the next, since the last line carried no newline. */
    fp = mcmd_open_output("procces_statistics.txt", "w");
    fprintf(fp, "Processes during the simulation \n| Auger decay | "
                "Fluorescence | Photo-ionization | charge transfer |\n");
    fprintf(fp, "%d %d %d %d\n", mc->count_auger, mc->count_fluorescence,
            mc->count_photoionization, mc->count_chargetransfer);
    fclose(fp);

    if (mc->fscreen != NULL)
    {
        fclose(mc->fscreen);
        mc->fscreen = NULL;
    }
    if (mc->coll_on)
    {
        fprintf(stderr, "mcionize: collisional: %.0f ionizations, %.0f three-body "
                "recombinations, %.0f eV net to ionization\n",
                mc->count_coll_ion, mc->count_coll_rec, mc->coll_energy);
        fprintf(stderr, "mcionize: collisional by shell (K L M): ionizations %.0f %.0f %.0f, "
                "recombinations %.0f %.0f %.0f\n",
                mc->coll_ion_shell[0], mc->coll_ion_shell[1], mc->coll_ion_shell[2],
                mc->coll_rec_shell[0], mc->coll_rec_shell[1], mc->coll_rec_shell[2]);
        if (mc->coll_excite)
        {
            fprintf(stderr, "mcionize: collisional: %.0f excitations, %.0f de-excitations\n",
                    mc->count_coll_exc, mc->count_coll_deexc);
        }
        if (mc->fcoll != NULL)
        {
            fclose(mc->fcoll);
            mc->fcoll = NULL;
        }
    }

    fp = mcmd_open_output("charges.txt", "w");
    for (i = 1; i < mc->natoms + 1; i++)
    {
        fprintf(fp, "%d %d\n", i, (int)mdatoms->chargeA[i - 1]);
    }
    fclose(fp);

    mcionize_free(mc);
}

static void mcionize_free(t_mcionize *mc)
{
    int i;

    if (mc->fcharges != NULL)  { fclose(mc->fcharges); }
    if (mc->fmean != NULL)     { fclose(mc->fmean); }
    if (mc->fpulse != NULL)    { fclose(mc->fpulse); }
    if (mc->felectron != NULL) { fclose(mc->felectron); }
    if (mc->ftrans != NULL)    { fclose(mc->ftrans); }
    for (i = 0; i < 5; i++)
    {
        sfree(mc->outbuf[i]);
    }

    for (i = 0; i < mc->natoms; i++)
    {
        sfree(mc->atom_configurations[i]);
        sfree(mc->GS_configurations[i]);
    }
    sfree(mc->atom_configurations);
    sfree(mc->GS_configurations);
    sfree(mc->energy_of_state);
    sfree(mc->outermost_shell);
    sfree(mc->elem_of_atom);
    sfree(mc->ct_order);
    sfree(mc->mc_transitions);
    sfree(mc->mc_dt);
    sfree(mc->mc_cref);
    {
        int el;

        for (el = 0; el < MCMD_ATOM_TABLE_SIZE; el++)
        {
            sfree(mc->cch[el]);
            sfree(mc->cch_start[el]);
            sfree(mc->cch_g[el][0]);
            sfree(mc->cch_g[el][1]);
        }
    }
    sfree(mc->cell_of_atom);
    sfree(mc->pos_of_atom);
    sfree(mc->cell_sorted);
    sfree(mc->cell_start);
    sfree(mc->cell_cursor);

    for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
    {
        mcmd_atomdata_done(&mc->atomdata[i], mc->coll_on);
    }
    sfree(mc->atomdata);

    sfree(mc);
}
