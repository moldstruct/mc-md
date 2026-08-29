/*
 * atomic_data.h
 *
 * Atomic data for the MolDStruct MC/MD ionization module: the table of
 * supported elements, and the per-element data read out of ./Atomic_data/
 * (energy levels, transition rates, and the collisional parameters and
 * statistical weights used by the not-yet-implemented collisional channel).
 *
 * An electronic state is three integers: the K, L and M shell occupations.
 *
 * This file is part of the MolDStruct modifications to GROMACS 4.5.4 and is
 * distributed under the same licence as GROMACS.
 */

#ifndef _atomic_data_h
#define _atomic_data_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Number of shells tracked per atom: K, L, M. */
#define MCMD_NSHELL 3

/* Number of elements the code knows the mass and ground state of. */
#define MCMD_NUM_ELEMENTS 15

/* Upper bound on element indices, used to size the per-element data table.
 * Must be >= MCMD_NUM_ELEMENTS. */
#define MCMD_ATOM_TABLE_SIZE 20

/* Value returned by mcmd_dict_get() and mcmd_weight_for_state() when the
 * requested state is not in the table.  All real energies and weights are
 * positive, so a negative sentinel is unambiguous. */
#define MCMD_NOT_FOUND (-1.0)

/* ---------------------------------------------------------------- *
 * Element table
 * ---------------------------------------------------------------- */

/* Element index for a rounded atomic mass, or -1 if the mass is not in the
 * table.  Callers must check for -1: it used to be passed straight to an
 * array subscript. */
int mcmd_mass2idx(int mass);

/* Rounded atomic mass for an element index, or -1 if out of range. */
int mcmd_idx2mass(int idx);

/* Periodic-table symbol in upper case ("H", "FE", ...) for a rounded mass,
 * or "_" if the mass is not in the table. */
const char *mcmd_mass2symbol(int mass);

/* Ground-state K/L/M occupations for a rounded mass, or NULL if unknown. */
const int *mcmd_ground_state(int mass);

/* ---------------------------------------------------------------- *
 * Transition data
 * ---------------------------------------------------------------- */

/* One transition out of a state.  type is 0 for Auger-Meitner decay, 1 for
 * fluorescence and 2 for photoionization; only type 2 has its rate scaled by
 * the pulse profile. */
typedef struct
{
    int    final_state[MCMD_NSHELL];
    double rate;
    int    type;
} t_mcmd_transition;

/* All transitions out of one initial state, as read from
 * rate_transitions_X.txt. */
typedef struct
{
    int     initial_state[MCMD_NSHELL];
    int     num_transitions;
    int   **final_states;  /* [num_transitions][MCMD_NSHELL] */
    double *rates;         /* [num_transitions] */
    int    *types;         /* [num_transitions] */
} t_mcmd_rate;

/* Collisional counterpart of t_mcmd_rate: same shape, but five fitting
 * parameters per transition instead of a single rate.  Read from
 * collisional_parameters_X.txt.  Only used when userint6 is set, which is
 * not a supported configuration yet - see mcionize.c. */
#define MCMD_NCOLL_PARAM 5

typedef struct
{
    int      initial_state[MCMD_NSHELL];
    int      num_transitions;
    int    **final_states; /* [num_transitions][MCMD_NSHELL] */
    double **coll_rates;   /* [num_transitions][MCMD_NCOLL_PARAM] */
} t_mcmd_coll;

/* Statistical weights per state, from statistical_weight_X.txt.  Also only
 * needed by the collisional channel. */
typedef struct
{
    int     num_states;
    int   **states; /* [num_states][MCMD_NSHELL] */
    double *weight; /* [num_states] */
} t_mcmd_weights;

/* Energy levels, keyed on the K/L/M occupations. */
typedef struct t_mcmd_dict t_mcmd_dict;

/* Direct-index lookup from a K/L/M state to a position in an array of states,
 * replacing a linear scan.  Shell occupations are small and bounded (K <= 2,
 * L <= 8, M <= 18 for the elements in the table), so a dense array indexed by
 * (K * dimL + L) * dimM + M costs a couple of kB per element and answers in
 * one load.  Entries are -1 where no such state exists. */
typedef struct
{
    int  dimK, dimL, dimM;
    int *idx;   /* [dimK * dimL * dimM], -1 where absent */
} t_mcmd_state_index;

/* Look up a state, or -1.  Out-of-range occupations return -1 rather than
 * indexing past the table. */
int mcmd_state_index_get(const t_mcmd_state_index *si,
                         const int state[MCMD_NSHELL]);

/* Everything known about one element. */
typedef struct
{
    int             mass;
    gmx_bool        bPresent; /* FALSE if this element is not in the system */
    t_mcmd_dict    *energy_levels;

    /* Total electronic binding energy per configuration, from
     * total_energies_X.txt.  Optional: NULL when the file is absent, which is
     * every Atomic_data directory generated before it existed.  Where it is
     * present the charge-transfer criterion can score an electron entering
     * any shell, as T(c + e_s) - T(c), instead of only the outermost one. */
    t_mcmd_dict    *total_energies;
    t_mcmd_rate    *rates;
    int             num_rates;
    t_mcmd_coll    *coll;
    int             num_coll;
    t_mcmd_weights  weights;

    /* O(1) replacements for the linear scans over rates[] and coll[].  Built
     * by mcmd_atomdata_read() from the same arrays, so they always agree with
     * mcmd_rate_state_index() / mcmd_coll_state_index(). */
    t_mcmd_state_index rate_index;
    t_mcmd_state_index coll_index;

    /* Largest num_transitions over all states of this element, so callers can
     * size a reusable scratch buffer once instead of allocating per event. */
    int             max_transitions;
} t_mcmd_atomdata;

/* Read the data files for one element into ad.  Reads the collisional
 * parameters and statistical weights as well when bCollisional is TRUE.
 * Fatal if a required file is missing or malformed. */
void mcmd_atomdata_read(t_mcmd_atomdata *ad, int elem_idx,
                        gmx_bool bCollisional);

/* Release everything mcmd_atomdata_read() allocated. */
void mcmd_atomdata_done(t_mcmd_atomdata *ad, gmx_bool bCollisional);

/* Index of the entry whose initial state matches, or -1. */
int mcmd_rate_state_index(const t_mcmd_rate *rates, int num_rates,
                          const int state[MCMD_NSHELL]);
int mcmd_coll_state_index(const t_mcmd_coll *coll, int num_coll,
                          const int state[MCMD_NSHELL]);

/* Energy of a state in eV, or MCMD_NOT_FOUND. */
double mcmd_dict_get(const t_mcmd_dict *dict, const int state[MCMD_NSHELL]);

/* Statistical weight of a state, or MCMD_NOT_FOUND. */
double mcmd_weight_for_state(const t_mcmd_weights *w,
                             const int state[MCMD_NSHELL]);

/* ---------------------------------------------------------------- *
 * Exponential integral Ei(x), used by the collisional cross sections
 * ---------------------------------------------------------------- */

double mcmd_exp_integral_Ei(double x);

#ifdef __cplusplus
}
#endif

#endif /* _atomic_data_h */
