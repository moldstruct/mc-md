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
#include <sys/stat.h>
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

/* Energy deposited in the plasma by one trapped electron, in eV. */
#define MCMD_TRAPPED_ELECTRON_EV 25.0

/* Ionization energy used in the electron escape test, in eV. */
#define MCMD_IONIZATION_ENERGY_EV 871.0

/* Energy difference assumed between collisional levels, in eV. */
#define MCMD_COLL_DELTA_E 200.0

/* Output directory and files, relative to the run directory. */
#define MCMD_OUTDIR "./simulation_output"

/* Restart file written at the end of every run.  Holds everything the
 * module cannot re-derive from the trajectory: the electronic state, the
 * free-electron plasma, and the sample geometry the plasma volume is
 * measured against. */
#define MCMD_STATE_FILE  "mcionize_state.bin"
#define MCMD_STATE_MAGIC 0x4d434d44 /* "MCMD" */
#define MCMD_STATE_VERSION 1

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

struct t_mcionize
{
    /* Configuration, cached from the inputrec */
    gmx_bool bChargeTransfer; /* mcmd-charge-transfer                  */
    gmx_bool bAllowHCT;       /* mcmd-allow-h-ct                       */
    gmx_bool bDownhill;       /* mcmd-charge-transfer-downhill         */
    gmx_bool bCollisional;    /* mcmd-collisional-ionization, see below */
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

    /* Per-element atomic data, indexed by mcmd_mass2idx() */
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

    gmx_bool bInitialised; /* set on the first step */
};

/* ---------------------------------------------------------------- *
 * Small helpers
 * ---------------------------------------------------------------- */

/* Seed for the module's RNG.  Normally decorrelated per job: time(NULL) alone
 * hands identical seeds to every job started within the same second, which
 * silently correlates trajectories across an ensemble, so the pid is mixed
 * in.  Setting GMX_MCMD_SEED pins the stream instead, which is what makes a
 * run reproducible for debugging and regression testing. */
static unsigned int mcmd_rng_seed(void)
{
    const char *env = getenv("GMX_MCMD_SEED");

    if (env != NULL)
    {
        return (unsigned int)strtoul(env, NULL, 10);
    }
    return (unsigned int)(time(NULL) ^ (getpid() << 16));
}

static void mcionize_free(t_mcionize *mc);

/* Open a file under MCMD_OUTDIR, or die trying. */
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

/* Truncate a text output file so a fresh run does not append to the last. */
static void mcmd_truncate_output(const char *name)
{
    fclose(mcmd_open_output(name, "w"));
}

/* Resolve an atom's element from its mass, turning the -1 that
 * mcmd_mass2idx() returns for an unlisted mass into a clear error instead of
 * an out-of-bounds read.  Checks that the element is known
 * and its data loaded.  Only used to fill elem_of_atom[] at setup; the hot
 * paths read that array instead. */
static int mcmd_resolve_element(const t_mcionize *mc, const t_mdatoms *mdatoms,
                                int atom)
{
    int mass = (int)floor(mdatoms->massT[atom] + 0.5);
    int idx  = mcmd_mass2idx(mass);

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
         * while leaving every observable alone. */
        if (mc->bDownhill && E_acceptor + dE_coul <= E_donor)
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

static void mcmd_charge_transfer(t_mcionize *mc, t_mdatoms *mdatoms,
                                 t_state *state, double t)
{
    int    *idx_map = mc->ct_order;
    int     i, j, j2;
    int     number_of_charge_transfers = 0;
    FILE   *fp = NULL;
    double  origin[3], edge;
    int     ncell[3];

    /* energy_of_state[] and outermost_shell[] are kept current by
     * mcmd_refresh_atom() wherever a configuration changes, so the loop that
     * used to re-derive the energy of every atom here - one hash lookup per
     * atom per step, ~15 million of them on a 25000 step run - is gone. */

    /* Positions do not change during this pass, so one build serves it all. */
    mcmd_build_cell_list(mc, mdatoms, state->x, origin, ncell, &edge);

    for (j = 0; j < mdatoms->nr; j++)
    {
        idx_map[j] = j;
    }
    mcmd_shuffle(idx_map, mdatoms->nr);

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

        {
        const int cell_lo = mc->cell_start[(cz * ncell[1] + cy) * ncell[0] + cx];
        const int cell_hi = mc->cell_start[(cz * ncell[1] + cy) * ncell[0] + cx + 1];
        int       n;

        for (n = cell_lo; n < cell_hi; n++)
        {
            const t_mcmd_cell_atom *ca = &mc->cell_sorted[n];
            double Q_D, Q_A, R_cob, R_crit, E_donor, dr2;
            int    cand_donor, cand_acceptor, cand_shell;

            i = ca->idx;

            if (i == j)
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

        if (mc->bDetailedOutput)
        {
            fprintf(mc->ftrans,
                    "Charge transfer at t: %lf, Previous state: "
                    "[%d, %d, %d]. Charge init: %lf ",
                    t,
                    mc->atom_configurations[R_min_idx][0],
                    mc->atom_configurations[R_min_idx][1],
                    mc->atom_configurations[R_min_idx][2],
                    (double)mdatoms->chargeA[R_min_idx]);
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
            fprintf(mc->ftrans,
                    " New state: [%d, %d, %d]. With %i charge transfers "
                    "occuring. Charge after: %lf \n",
                    mc->atom_configurations[R_min_idx][0],
                    mc->atom_configurations[R_min_idx][1],
                    mc->atom_configurations[R_min_idx][2],
                    number_of_charge_transfers,
                    (double)mdatoms->chargeA[R_min_idx]);
        }
    }

}

/* ---------------------------------------------------------------- *
 * Collisional rates
 *
 * NOT IMPLEMENTED.  The data plumbing (collisional_parameters_X.txt,
 * statistical_weight_X.txt, the exponential integral) is kept so the
 * channel can be finished later, but the cross sections below are
 * untested and the clamps are arbitrary.  mcmd-collisional-ionization is
 * documented as
 * unimplemented; treat any run with it enabled as experimental.
 * ---------------------------------------------------------------- */

static void mcmd_collisional_rates(t_mcionize *mc,
                                   const t_mcmd_coll *coll, int coll_index,
                                   const t_mcmd_weights *weights,
                                   const int match[MCMD_NSHELL],
                                   int itrans, int offset,
                                   t_mcmd_transition *out)
{
    const double *p = coll[coll_index].coll_rates[itrans - offset];
    const int    *final_state = coll[coll_index].final_states[itrans - offset];
    double        tev, B, qb, gi, gf;
    double        sigma_collisional, sigma_recombination;

    out->final_state[0] = final_state[0];
    out->final_state[1] = final_state[1];
    out->final_state[2] = final_state[2];

    tev = mc->thermalized_free_electron_energy;
    B   = MCMD_COLL_DELTA_E / tev;

    qb = p[0] +
         p[1] * sqrt(MCMD_PI * B) * erfc(B) * exp(B * B) +
         p[2] * B +
         p[3] * mcmd_exp_integral_Ei(B) * exp(B) * p[4];

    gi = mcmd_weight_for_state(weights, match);
    gf = mcmd_weight_for_state(weights, final_state);

    if (fabs(mc->num_free_electrons) < 1e-5 || isnan(qb) || qb < 0.0)
    {
        sigma_collisional   = 0.00001;
        sigma_recombination = 0.00001;
    }
    else
    {
        sigma_collisional = mc->electron_density *
                            (1.09 * 1e-6 * qb * exp(-B)) /
                            (MCMD_COLL_DELTA_E * sqrt(tev));
        sigma_recombination = sigma_collisional * 1.66 * 1e-22 * (gi / gf) *
                              (mc->electron_density / pow(tev, 1.5)) *
                              exp(-MCMD_COLL_DELTA_E / (MCMD_BOLTZMANN * tev));

        if (sigma_collisional < 0 || fabs(sigma_collisional) < 1e-5)
        {
            sigma_collisional = 0.00001;
        }
        if (sigma_recombination < 0 || fabs(sigma_recombination) < 1e-5)
        {
            sigma_recombination = 0.0001;
        }
    }

    /* Losing electrons is a collision, gaining one is recombination. */
    if ((match[0] - final_state[0] +
         match[1] - final_state[1] +
         match[2] - final_state[2]) > 0)
    {
        out->rate = sigma_collisional;
    }
    else
    {
        out->rate = sigma_recombination;
    }

    if (isinf(out->rate) || isnan(out->rate))
    {
        out->rate = 0.0;
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

static void mcmd_monte_carlo(t_mcionize *mc, const t_inputrec *ir,
                             t_mdatoms *mdatoms, t_state *state, double t)
{
    const double dt_md_s = ir->delta_t * 1e-12;
    int          k;
    FILE        *fp;

    for (k = 0; k < mdatoms->nr; k++)
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
            int                match_index, match_index_coll = -1;
            int                offset;
            double             net_charge, critical_potential;
            int                i;

            match[0] = mc->atom_configurations[k][0];
            match[1] = mc->atom_configurations[k][1];
            match[2] = mc->atom_configurations[k][2];

            /* Array read into the per-element table built at setup, in place
             * of a scan over every state of the element - 324 of them for
             * sulphur, and one scan per atom per event per step. */
            match_index = mcmd_state_index_get(&ad->rate_index, match);

            if (mc->bCollisional)
            {
                match_index_coll = mcmd_state_index_get(&ad->coll_index, match);
            }

            if (match_index == -1 ||
                (mc->bCollisional && match_index_coll == -1))
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
            if (mc->bCollisional)
            {
                total_transitions += ad->coll[match_index_coll].num_transitions;
            }

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

            for (i = offset; i < total_transitions; i++)
            {
                mcmd_collisional_rates(mc, ad->coll, match_index_coll,
                                       &ad->weights, match, i, offset,
                                       &possible_transitions[i]);
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

            if (mc->bDetailedOutput)
            {
                fprintf(mc->ftrans,
                        "Making a transition, the time is %lf and current "
                        "Monte-Carlo DT is: %e.\n", t, DT_current);
            }

            net_charge =
                (double)(match[0] - possible_transitions[min_value_index].final_state[0] +
                         match[1] - possible_transitions[min_value_index].final_state[1] +
                         match[2] - possible_transitions[min_value_index].final_state[2]);

            /* Does the ejected electron escape the sample's space charge, or
             * is it trapped into the free-electron plasma? */
            critical_potential =
                100 * (MCMD_ELEMENTARY_CHARGE * MCMD_ELEMENTARY_CHARGE * 4 *
                       MCMD_PI * 9e9 * 1e9 *
                       (mc->escaped_photoelectrons / mc->current_volume) *
                       state->box[XX][XX] * state->box[XX][XX]) / 3;

            if ((mc->photon_energy_eV - MCMD_IONIZATION_ENERGY_EV) *
                MCMD_EV_TO_JOULE > critical_potential)
            {
                mc->escaped_photoelectrons += 1;
            }
            else
            {
                mc->num_free_electrons += net_charge;
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
                fprintf(mc->ftrans, "Transition of type %d \n",
                        possible_transitions[min_value_index].type);
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
    double radius_of_gyration = 0;
    int    i;

    for (i = 0; i < mdatoms->nr; i++)
    {
        radius_of_gyration += mdatoms->massT[i] * (sqr(state->x[i][0]) +
                                                   sqr(state->x[i][1]) +
                                                   sqr(state->x[i][2]));
    }
    radius_of_gyration = sqrt(radius_of_gyration / mc->total_mass);

    mc->rg_factor = radius_of_gyration / mc->radius_gyration_equilibrium;
    /* Sphere volume, so the linear expansion factor enters cubed. */
    mc->current_volume = mc->sample_volume *
                         mc->rg_factor * mc->rg_factor * mc->rg_factor;
}

static void mcmd_write_step_output(t_mcionize *mc, t_mdatoms *mdatoms, double t)
{
    double mean_charge = 0;
    FILE  *fp;
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

    fprintf(mc->fmean,  "%lf %lf\n", t, mean_charge);
    fprintf(mc->fpulse, "%lf %lf\n", t, mc->intensity);

    /* Plain numeric columns so the file loads directly with np.loadtxt:
     *   time [ps]   e- density [cm^-3]   e- temperature [eV]   R_g / R_g(0)
     * (The Debye length the manual used to promise was never computed.) */
    fprintf(mc->felectron, "%g %g %g %g\n",
            t, mc->electron_density, mc->electron_temperature, mc->rg_factor);

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
 * ---------------------------------------------------------------- */

static void mcmd_state_write(const t_mcionize *mc, const t_mdatoms *mdatoms)
{
    FILE *fp = mcmd_open_output(MCMD_STATE_FILE, "wb");
    int   header[3];
    int   counters[4];
    double geom[7];
    int    i;

    header[0] = MCMD_STATE_MAGIC;
    header[1] = MCMD_STATE_VERSION;
    header[2] = mc->natoms;
    fwrite(header, sizeof(int), 3, fp);

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

    for (i = 0; i < mc->natoms; i++)
    {
        double charge = (double)mdatoms->chargeA[i];

        fwrite(mc->atom_configurations[i], sizeof(int), MCMD_NSHELL, fp);
        fwrite(mc->GS_configurations[i], sizeof(int), MCMD_NSHELL, fp);
        fwrite(&charge, sizeof(double), 1, fp);
    }

    fclose(fp);
}

static void mcmd_state_read(t_mcionize *mc, t_mdatoms *mdatoms)
{
    char   path[512];
    FILE  *fp;
    int    header[3];
    int    counters[4];
    double geom[7];
    int    i;

    sprintf(path, "%s/%s", MCMD_OUTDIR, MCMD_STATE_FILE);
    fp = fopen(path, "rb");
    if (fp == NULL)
    {
        gmx_fatal(FARGS,
                  "Restarting an ionizing run needs '%s', which is written at "
                  "the end of every run, but it could not be opened.  Without "
                  "it the electronic state and the free-electron plasma cannot "
                  "be recovered.", path);
    }

    if (fread(header, sizeof(int), 3, fp) != 3 ||
        header[0] != MCMD_STATE_MAGIC)
    {
        gmx_fatal(FARGS, "'%s' is not an mcionize state file.", path);
    }
    if (header[1] != MCMD_STATE_VERSION)
    {
        gmx_fatal(FARGS, "'%s' has version %d, this build writes version %d.",
                  path, header[1], MCMD_STATE_VERSION);
    }
    if (header[2] != mc->natoms)
    {
        gmx_fatal(FARGS,
                  "'%s' holds %d atoms but this system has %d.",
                  path, header[2], mc->natoms);
    }

    if (fread(geom, sizeof(double), 7, fp) != 7 ||
        fread(counters, sizeof(int), 4, fp) != 4)
    {
        gmx_fatal(FARGS, "'%s' is truncated.", path);
    }

    mc->total_mass                       = geom[0];
    mc->radius_of_sample                 = geom[1];
    mc->radius_gyration_equilibrium      = geom[2];
    mc->sample_volume                    = geom[3];
    mc->num_free_electrons               = geom[4];
    mc->thermalized_free_electron_energy = geom[5];
    mc->escaped_photoelectrons           = geom[6];

    mc->count_auger           = counters[0];
    mc->count_fluorescence    = counters[1];
    mc->count_photoionization = counters[2];
    mc->count_chargetransfer  = counters[3];

    for (i = 0; i < mc->natoms; i++)
    {
        double charge;

        if (fread(mc->atom_configurations[i], sizeof(int), MCMD_NSHELL, fp) != MCMD_NSHELL ||
            fread(mc->GS_configurations[i], sizeof(int), MCMD_NSHELL, fp) != MCMD_NSHELL ||
            fread(&charge, sizeof(double), 1, fp) != 1)
        {
            gmx_fatal(FARGS, "'%s' is truncated at atom %d.", path, i + 1);
        }
        mdatoms->chargeA[i] = (real)charge;
    }

    fclose(fp);

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
    double max_x = 0, max_y = 0, max_z = 0;
    double radius_of_gyration = 0;
    int    i, j;

    mc->electron_temperature               = 0.0;
    mc->num_free_electrons                 = 0.0;
    mc->thermalized_free_electron_energy   = 0.0;
    mc->escaped_photoelectrons             = 0.0;

    /* Radius of the sample, estimated from the largest coordinates. */
    for (i = 0; i < mdatoms->nr; i++)
    {
        if (state->x[i][0] > max_x) { max_x = state->x[i][0]; }
        if (state->x[i][1] > max_y) { max_y = state->x[i][1]; }
        if (state->x[i][2] > max_z) { max_z = state->x[i][2]; }
    }
    mc->radius_of_sample = sqrt(max_x * max_x + max_y * max_y + max_z * max_z);

    for (i = 0; i < mdatoms->nr; i++)
    {
        radius_of_gyration += mdatoms->massT[i] * (sqr(state->x[i][0]) +
                                                   sqr(state->x[i][1]) +
                                                   sqr(state->x[i][2]));
        mc->total_mass += mdatoms->massT[i];
        /* The model tracks integer net charge, so the force field's partial
         * charges are discarded. */
        mdatoms->chargeA[i] = 0.0;
    }
    mc->radius_gyration_equilibrium =
        sqrt(radius_of_gyration / mc->total_mass);

    srand(mcmd_rng_seed());

    mc->num_free_electrons               = 0;
    mc->thermalized_free_electron_energy = 0;
    /* Sphere volume in cm^3: the 100^3 and (1e-9)^3 factors together turn
     * a radius in nm into a volume in cm^3, so the electron density
     * below comes out per cm^3. */
    mc->sample_volume = 100 * 100 * 100 * 1e-9 * 1e-9 * 1e-9 *
                        mc->radius_of_sample * mc->radius_of_sample *
                        mc->radius_of_sample;

    /* Ground states */
    for (i = 0; i < mdatoms->nr; i++)
    {
        const int *gs = mcmd_ground_state((int)floor(mdatoms->massT[i] + 0.5));

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

t_mcionize *mcionize_init(FILE *fplog, const t_inputrec *ir,
                          t_mdatoms *mdatoms, const t_commrec *cr,
                          gmx_bool bStartFromCpt)
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
    mc->bCollisional    = (ir->mcmd_collisional_ionization != 0);
    mc->bInitCharges    = (ir->mcmd_initial_charges != 0);
    mc->bDetailedOutput        = (ir->mcmd_detailed_output != 0);

    /* charges_over_time.bin is one real per atom per step, which dwarfs
     * everything else the module writes - 196 MB for a 25000 step run of a
     * 1960 atom protein, and it grows with both.  Writing only every Nth step
     * cuts that by N.  Clamped to at least 1; the three text files stay
     * per-step, since together they are under 100 bytes a step. */
    mc->charge_stride = (ir->mcmd_charge_output_stride > 0 ?
                         ir->mcmd_charge_output_stride : 1);

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
        if (stat(MCMD_OUTDIR, &st) == -1)
        {
            if (mkdir(MCMD_OUTDIR, 0700) != 0)
            {
                gmx_fatal(FARGS, "Could not create the '%s' directory.",
                          MCMD_OUTDIR);
            }
        }

        /* Only the detailed outputs are opened here; without them the module
         * writes nothing per step.  Opening with "w" also truncates, which is
         * what mcmd_truncate_output() used to do at the first step. */
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
            mc->ftrans    = mcmd_open_buffered("electronic_transition_log.txt",
                                               tmode, &mc->outbuf[4]);
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
        int mass = (int)floor(mdatoms->massT[i] + 0.5);
        int idx  = mcmd_mass2idx(mass);

        if (idx < 0)
        {
            gmx_fatal(FARGS,
                      "Atom %d has mass %g, which is not in the element "
                      "table.  Add it to mcmd_elements[] in atomic_data.c.",
                      i + 1, mdatoms->massT[i]);
        }
        if (!mc->atomdata[idx].bPresent)
        {
            mcmd_atomdata_read(&mc->atomdata[idx], idx, mc->bCollisional);
            if (fplog != NULL)
            {
                fprintf(fplog, "mcionize: loaded atomic data for %s\n",
                        mcmd_mass2symbol(mass));
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
                strcat(missing, mcmd_mass2symbol(mcmd_idx2mass(i)));
                nmissing++;
            }
        }

        if (nmissing > 0)
        {
            gmx_fatal(FARGS,
                      "Charge transfer is on but Atomic_data has no "
                      "total_energies_X.txt for %s.  Regenerate the atomic "
                      "data - explode_tools.generate_atomic_data() writes it "
                      "alongside energy_levels_X.txt - or set "
                      "mcmd-charge-transfer = 0.  See \"Supplying atomic "
                      "data\" in the README.",
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

    srand(mcmd_rng_seed());

    if (bStartFromCpt)
    {
        /* Continuation: recover the electronic state rather than resetting to
         * ground state.  The reference radius of gyration in particular must
         * come from the file - recomputing it from an already-expanded sample
         * would reset rg_factor to 1 and shrink the plasma volume back to its
         * t=0 value. */
        if (MASTER(cr))
        {
            mcmd_state_read(mc, mdatoms);
        }
        if (PAR(cr))
        {
            gmx_bcast(mdatoms->nr * sizeof(mdatoms->chargeA[0]),
                      mdatoms->chargeA, cr);
        }
        mc->bInitialised = TRUE;

        if (fplog != NULL)
        {
            fprintf(fplog,
                    "mcionize: restarted from %s/%s (%.0f trapped electrons, "
                    "%.0f escaped)\n",
                    MCMD_OUTDIR, MCMD_STATE_FILE,
                    mc->num_free_electrons, mc->escaped_photoelectrons);
        }
    }
    else
    {
        mc->bInitialised = FALSE;
    }

    return mc;
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
     * transfer pair search needs all of them, so gather first - the same ring
     * exchange do_force() performs a little later. */
    if (PAR(cr))
    {
        move_x(fplog, cr, GMX_LEFT, GMX_RIGHT, state->x, nrnb);
    }

    /* The module itself is serial, on the master.  It is ~3% of the step, and
     * the electron escape test reads escaped_photoelectrons back inside the
     * same atom loop that updates it, so splitting the loop across ranks
     * would change which electrons are trapped.  Keeping it on one rank
     * preserves the serial physics exactly and still lets the force loop -
     * which is the actual cost - use every rank. */
    if (MASTER(cr))
    {
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

                mcmd_charge_transfer(mc, mdatoms, state, t);

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

        mcmd_update_volume(mc, mdatoms, state);
        mcmd_monte_carlo(mc, ir, mdatoms, state, t);
        mcmd_write_step_output(mc, mdatoms, t);
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

void mcionize_done(t_mcionize *mc, t_mdatoms *mdatoms, const t_commrec *cr)
{
    FILE *fp;
    int   i;

    if (mc == NULL)
    {
        return;
    }

    if (!MASTER(cr))
    {
        mcionize_free(mc);
        return;
    }

    mcmd_state_write(mc, mdatoms);

    if (mc->ct_skipped > 0)
    {
        fprintf(stderr,
                "\nmcionize: charge transfer went idle; the pass ran on %.0f%% "
                "of steps (sampled every %d once idle).\n",
                100.0 * (double)(mc->ct_step - mc->ct_skipped) /
                (double)mc->ct_step, mc->ct_recheck);
    }

    if (mc->max_transfer_separation > 0.8 * MCMD_CT_CUTOFF_NM)
    {
        fprintf(stderr,
                "\nWARNING: a charge transfer used a pair %.3f nm apart, "
                "against a search cutoff of %.3f nm.  Transfers may be being "
                "lost to the cutoff; raise MCMD_CT_CUTOFF_NM in mcionize.c.\n",
                mc->max_transfer_separation, (double)MCMD_CT_CUTOFF_NM);
    }

    fp = mcmd_open_output("procces_statistics.txt", "a");
    fprintf(fp, "Processes during the simulation \n| Auger decay | "
                "Fluorescence | Photo-ionization | charge transfer |\n");
    fprintf(fp, "%d %d %d %d", mc->count_auger, mc->count_fluorescence,
            mc->count_photoionization, mc->count_chargetransfer);
    fclose(fp);

    fp = mcmd_open_output("charges.txt", "a");
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
    sfree(mc->cell_of_atom);
    sfree(mc->pos_of_atom);
    sfree(mc->cell_sorted);
    sfree(mc->cell_start);
    sfree(mc->cell_cursor);

    for (i = 0; i < MCMD_ATOM_TABLE_SIZE; i++)
    {
        mcmd_atomdata_done(&mc->atomdata[i], mc->bCollisional);
    }
    sfree(mc->atomdata);

    sfree(mc);
}
