/*
 * atomic_data.c
 *
 * Element table and ./Atomic_data/ file parsing for the MolDStruct MC/MD
 * ionization module.  See atomic_data.h.
 *
 * This file is part of the MolDStruct modifications to GROMACS 4.5.4 and is
 * distributed under the same licence as GROMACS.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic_data.h"
#include "gmx_fatal.h"
#include "smalloc.h"

/* Directories are relative to the directory mdrun is started from. */
#define MCMD_PATH_ENERGY "./Atomic_data/energy_levels_"
#define MCMD_PATH_RATES  "./Atomic_data/rate_transitions_"
#define MCMD_PATH_COLL   "./Atomic_data/collisional_parameters_"
#define MCMD_PATH_WEIGHT "./Atomic_data/statistical_weight_"

/* Longest line accepted in a data file.  A single initial state can have many
 * final states on one line, so this needs headroom. */
#define MCMD_LINEBUF 1024

/* ================================================================ *
 * Element table
 *
 * Adding an element means adding a row to both tables below, and
 * supplying its energy_levels_X.txt and rate_transitions_X.txt.
 * ================================================================ */

typedef struct
{
    int  mass; /* rounded atomic mass, used as the species identifier */
    char symbol[3];
    int  ground_state[MCMD_NSHELL]; /* K, L, M occupations */
} t_mcmd_element;

static const t_mcmd_element mcmd_elements[MCMD_NUM_ELEMENTS] = {
    {  1, "H",  {1, 0, 0}  },
    { 12, "C",  {2, 4, 0}  },
    { 14, "N",  {2, 5, 0}  },
    { 16, "O",  {2, 6, 0}  },
    { 19, "F",  {2, 7, 0}  },
    { 24, "MG", {2, 8, 2}  },
    { 31, "P",  {2, 8, 5}  },
    { 32, "S",  {2, 8, 6}  },
    { 35, "CL", {2, 8, 7}  },
    { 40, "CA", {2, 8, 10} },
    { 56, "FE", {2, 8, 16} },
    { 59, "NI", {2, 8, 18} },
    { 28, "SI", {2, 8, 4}  },
    { 23, "NA", {2, 8, 1}  },
    /* Iodine's configuration is a placeholder copied from nickel: its bonded
     * behaviour is right but its ionization physics is not. */
    { 127, "I", {2, 8, 18} }
};

int mcmd_mass2idx(int mass)
{
    int i;

    for (i = 0; i < MCMD_NUM_ELEMENTS; i++)
    {
        if (mcmd_elements[i].mass == mass)
        {
            return i;
        }
    }
    return -1;
}

int mcmd_idx2mass(int idx)
{
    if (idx >= 0 && idx < MCMD_NUM_ELEMENTS)
    {
        return mcmd_elements[idx].mass;
    }
    return -1;
}

const char *mcmd_mass2symbol(int mass)
{
    int idx = mcmd_mass2idx(mass);

    return (idx >= 0) ? mcmd_elements[idx].symbol : "_";
}

const int *mcmd_ground_state(int mass)
{
    int idx = mcmd_mass2idx(mass);

    return (idx >= 0) ? mcmd_elements[idx].ground_state : NULL;
}

/* Build "<prefix><SYMBOL>.txt".  Caller frees. */
static char *mcmd_data_path(const char *prefix, int mass)
{
    const char *symbol = mcmd_mass2symbol(mass);
    char       *path;
    int         len;

    len = (int)(strlen(prefix) + strlen(symbol) + strlen(".txt") + 1);
    snew(path, len);
    sprintf(path, "%s%s.txt", prefix, symbol);

    return path;
}

/* Number of lines in a file, or -1 if it cannot be opened. */
static int mcmd_count_lines(const char *path)
{
    FILE *fp;
    int   count = 0;
    int   ch;
    int   last_was_newline = 1;

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        return -1;
    }

    while ((ch = fgetc(fp)) != EOF)
    {
        last_was_newline = (ch == '\n');
        if (last_was_newline)
        {
            count++;
        }
    }
    /* A final line without a trailing newline still counts. */
    if (!last_was_newline)
    {
        count++;
    }

    fclose(fp);

    return count;
}

static gmx_bool mcmd_state_equal(const int a[MCMD_NSHELL], const int b[MCMD_NSHELL])
{
    return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}

/* ================================================================ *
 * Energy levels
 *
 * A small open-chaining hash table keyed on the three shell
 * occupations.  Lookups are exact; a missing state returns
 * MCMD_NOT_FOUND.
 * ================================================================ */

#define MCMD_DICT_SIZE 211 /* prime, comfortably above the largest table */

typedef struct t_mcmd_dict_entry
{
    int                       state[MCMD_NSHELL];
    double                    energy;
    struct t_mcmd_dict_entry *next;
} t_mcmd_dict_entry;

struct t_mcmd_dict
{
    t_mcmd_dict_entry *bucket[MCMD_DICT_SIZE];
};

static int mcmd_dict_hash(const int state[MCMD_NSHELL])
{
    unsigned int h = (unsigned int)(state[0] * 73856093) ^
                     (unsigned int)(state[1] * 19349663) ^
                     (unsigned int)(state[2] * 83492791);

    return (int)(h % MCMD_DICT_SIZE);
}

static void mcmd_dict_insert(t_mcmd_dict *dict, const int state[MCMD_NSHELL],
                             double energy)
{
    int                h = mcmd_dict_hash(state);
    t_mcmd_dict_entry *e;

    snew(e, 1);
    e->state[0]   = state[0];
    e->state[1]   = state[1];
    e->state[2]   = state[2];
    e->energy     = energy;
    e->next       = dict->bucket[h];
    dict->bucket[h] = e;
}

double mcmd_dict_get(const t_mcmd_dict *dict, const int state[MCMD_NSHELL])
{
    const t_mcmd_dict_entry *e;

    if (dict == NULL)
    {
        return MCMD_NOT_FOUND;
    }

    for (e = dict->bucket[mcmd_dict_hash(state)]; e != NULL; e = e->next)
    {
        if (mcmd_state_equal(e->state, state))
        {
            return e->energy;
        }
    }

    return MCMD_NOT_FOUND;
}

static t_mcmd_dict *mcmd_dict_read(const char *path)
{
    t_mcmd_dict *dict;
    FILE        *fp;
    int          state[MCMD_NSHELL];
    double       energy;
    int          i;

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open the energy level file '%s'.", path);
    }

    snew(dict, 1);
    for (i = 0; i < MCMD_DICT_SIZE; i++)
    {
        dict->bucket[i] = NULL;
    }

    while (fscanf(fp, "%d %d %d %lf",
                  &state[0], &state[1], &state[2], &energy) == 4)
    {
        mcmd_dict_insert(dict, state, energy);
    }

    fclose(fp);

    return dict;
}

static void mcmd_dict_done(t_mcmd_dict *dict)
{
    int i;

    if (dict == NULL)
    {
        return;
    }

    for (i = 0; i < MCMD_DICT_SIZE; i++)
    {
        t_mcmd_dict_entry *e = dict->bucket[i];

        while (e != NULL)
        {
            t_mcmd_dict_entry *next = e->next;

            sfree(e);
            e = next;
        }
    }
    sfree(dict);
}

/* ================================================================ *
 * Transition rates
 *
 * One initial state per line:
 *   a b c ; a' b' c' rate type ; a'' b'' c'' rate type ; ...
 * Every final state carries its own rate and its own type; a token
 * with fewer than five fields is fatal, since data predating the type
 * column cannot be interpreted (the type is not recoverable from the
 * state change).
 * ================================================================ */

static void mcmd_rate_read_line(FILE *fp, const char *path, t_mcmd_rate *rate)
{
    char  line[MCMD_LINEBUF];
    char *token;
    char *saveptr;
    int   nsemi = 0;
    int   i;

    if (fgets(line, sizeof(line), fp) == NULL)
    {
        gmx_fatal(FARGS, "Unexpected end of file while reading '%s'.", path);
    }

    if (sscanf(line, "%d %d %d;",
               &rate->initial_state[0],
               &rate->initial_state[1],
               &rate->initial_state[2]) != 3)
    {
        gmx_fatal(FARGS, "Could not read an initial state from '%s'.", path);
    }

    /* The initial state is followed by one ';' per final state. */
    for (i = 0; line[i] != '\0'; i++)
    {
        if (line[i] == ';')
        {
            nsemi++;
        }
    }
    rate->num_transitions = nsemi - 1;

    if (rate->num_transitions < 0)
    {
        gmx_fatal(FARGS, "Malformed transition line in '%s'.", path);
    }

    snew(rate->final_states, rate->num_transitions);
    snew(rate->rates, rate->num_transitions);
    snew(rate->types, rate->num_transitions);
    for (i = 0; i < rate->num_transitions; i++)
    {
        snew(rate->final_states[i], MCMD_NSHELL);
    }

    /* Skip the initial-state token, then read one final state per token. */
    token = strtok_r(line, ";", &saveptr);
    token = strtok_r(NULL, ";", &saveptr);

    for (i = 0; i < rate->num_transitions; i++)
    {
        if (token == NULL)
        {
            gmx_fatal(FARGS, "Too few transitions on a line of '%s'.", path);
        }
        if (sscanf(token, "%d %d %d %lf %d",
                   &rate->final_states[i][0],
                   &rate->final_states[i][1],
                   &rate->final_states[i][2],
                   &rate->rates[i],
                   &rate->types[i]) != 5)
        {
            gmx_fatal(FARGS,
                      "Could not read transition data from '%s'.  Each final "
                      "state needs five fields: 'a b c rate type'.", path);
        }
        token = strtok_r(NULL, ";", &saveptr);
    }
}

static t_mcmd_rate *mcmd_rates_read(const char *path, int *num_rates)
{
    t_mcmd_rate *rates;
    FILE        *fp;
    int          i;

    *num_rates = mcmd_count_lines(path);
    if (*num_rates < 0)
    {
        gmx_fatal(FARGS, "Could not open the transition rate file '%s'.", path);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open the transition rate file '%s'.", path);
    }

    snew(rates, *num_rates);
    for (i = 0; i < *num_rates; i++)
    {
        mcmd_rate_read_line(fp, path, &rates[i]);
    }

    fclose(fp);

    return rates;
}

static void mcmd_rates_done(t_mcmd_rate *rates, int num_rates)
{
    int i, j;

    if (rates == NULL)
    {
        return;
    }

    for (i = 0; i < num_rates; i++)
    {
        for (j = 0; j < rates[i].num_transitions; j++)
        {
            sfree(rates[i].final_states[j]);
        }
        sfree(rates[i].final_states);
        sfree(rates[i].rates);
        sfree(rates[i].types);
    }
    sfree(rates);
}

int mcmd_rate_state_index(const t_mcmd_rate *rates, int num_rates,
                          const int state[MCMD_NSHELL])
{
    int i;

    for (i = 0; i < num_rates; i++)
    {
        if (mcmd_state_equal(rates[i].initial_state, state))
        {
            return i;
        }
    }

    return -1;
}

/* ================================================================ *
 * Collisional parameters and statistical weights
 *
 * Only read when userint6 is set.  The collisional channel is not
 * implemented, so this data is parsed but the physics using it is
 * untested - see the comment above mcionize_collisional_rates().
 * ================================================================ */

static void mcmd_coll_read_line(FILE *fp, const char *path, t_mcmd_coll *coll)
{
    char  line[MCMD_LINEBUF];
    char *token;
    char *saveptr;
    int   nsemi = 0;
    int   i;

    if (fgets(line, sizeof(line), fp) == NULL)
    {
        gmx_fatal(FARGS, "Unexpected end of file while reading '%s'.", path);
    }

    if (sscanf(line, "%d %d %d;",
               &coll->initial_state[0],
               &coll->initial_state[1],
               &coll->initial_state[2]) != 3)
    {
        gmx_fatal(FARGS, "Could not read an initial state from '%s'.", path);
    }

    for (i = 0; line[i] != '\0'; i++)
    {
        if (line[i] == ';')
        {
            nsemi++;
        }
    }
    coll->num_transitions = nsemi - 1;

    if (coll->num_transitions < 0)
    {
        gmx_fatal(FARGS, "Malformed collisional line in '%s'.", path);
    }

    snew(coll->final_states, coll->num_transitions);
    snew(coll->coll_rates, coll->num_transitions);
    for (i = 0; i < coll->num_transitions; i++)
    {
        snew(coll->final_states[i], MCMD_NSHELL);
        snew(coll->coll_rates[i], MCMD_NCOLL_PARAM);
    }

    token = strtok_r(line, ";", &saveptr);
    token = strtok_r(NULL, ";", &saveptr);

    for (i = 0; i < coll->num_transitions; i++)
    {
        if (token == NULL)
        {
            gmx_fatal(FARGS, "Too few transitions on a line of '%s'.", path);
        }
        if (sscanf(token, "%d %d %d %lf %lf %lf %lf %lf",
                   &coll->final_states[i][0],
                   &coll->final_states[i][1],
                   &coll->final_states[i][2],
                   &coll->coll_rates[i][0],
                   &coll->coll_rates[i][1],
                   &coll->coll_rates[i][2],
                   &coll->coll_rates[i][3],
                   &coll->coll_rates[i][4]) != 3 + MCMD_NCOLL_PARAM)
        {
            gmx_fatal(FARGS,
                      "Could not read collisional data from '%s'.  Each final "
                      "state needs three occupations and %d parameters.",
                      path, MCMD_NCOLL_PARAM);
        }
        token = strtok_r(NULL, ";", &saveptr);
    }
}

static t_mcmd_coll *mcmd_coll_read(const char *path, int *num_coll)
{
    t_mcmd_coll *coll;
    FILE        *fp;
    int          i;

    *num_coll = mcmd_count_lines(path);
    if (*num_coll < 0)
    {
        gmx_fatal(FARGS, "Could not open the collisional parameter file '%s'.",
                  path);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open the collisional parameter file '%s'.",
                  path);
    }

    snew(coll, *num_coll);
    for (i = 0; i < *num_coll; i++)
    {
        mcmd_coll_read_line(fp, path, &coll[i]);
    }

    fclose(fp);

    return coll;
}

static void mcmd_coll_done(t_mcmd_coll *coll, int num_coll)
{
    int i, j;

    if (coll == NULL)
    {
        return;
    }

    for (i = 0; i < num_coll; i++)
    {
        for (j = 0; j < coll[i].num_transitions; j++)
        {
            sfree(coll[i].final_states[j]);
            sfree(coll[i].coll_rates[j]);
        }
        sfree(coll[i].final_states);
        sfree(coll[i].coll_rates);
    }
    sfree(coll);
}

int mcmd_coll_state_index(const t_mcmd_coll *coll, int num_coll,
                          const int state[MCMD_NSHELL])
{
    int i;

    for (i = 0; i < num_coll; i++)
    {
        if (mcmd_state_equal(coll[i].initial_state, state))
        {
            return i;
        }
    }

    return -1;
}

static void mcmd_weights_read(t_mcmd_weights *w, const char *path)
{
    FILE *fp;
    int   i;

    w->num_states = mcmd_count_lines(path);
    if (w->num_states < 0)
    {
        gmx_fatal(FARGS, "Could not open the statistical weight file '%s'.",
                  path);
    }

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        gmx_fatal(FARGS, "Could not open the statistical weight file '%s'.",
                  path);
    }

    snew(w->states, w->num_states);
    snew(w->weight, w->num_states);

    for (i = 0; i < w->num_states; i++)
    {
        snew(w->states[i], MCMD_NSHELL);
        if (fscanf(fp, "%d %d %d %lf",
                   &w->states[i][0], &w->states[i][1], &w->states[i][2],
                   &w->weight[i]) != 4)
        {
            gmx_fatal(FARGS, "Could not read statistical weights from '%s'.",
                      path);
        }
    }

    fclose(fp);
}

static void mcmd_weights_done(t_mcmd_weights *w)
{
    int i;

    if (w->states == NULL)
    {
        return;
    }

    for (i = 0; i < w->num_states; i++)
    {
        sfree(w->states[i]);
    }
    sfree(w->states);
    sfree(w->weight);

    w->states     = NULL;
    w->weight     = NULL;
    w->num_states = 0;
}

double mcmd_weight_for_state(const t_mcmd_weights *w,
                             const int state[MCMD_NSHELL])
{
    int i;

    for (i = 0; i < w->num_states; i++)
    {
        if (mcmd_state_equal(w->states[i], state))
        {
            return w->weight[i];
        }
    }

    return MCMD_NOT_FOUND;
}

/* ================================================================ *
 * Public entry points
 * ================================================================ */

/* ---------------------------------------------------------------- *
 * Direct-index state lookup
 * ---------------------------------------------------------------- */

int mcmd_state_index_get(const t_mcmd_state_index *si,
                         const int state[MCMD_NSHELL])
{
    if (si->idx == NULL ||
        state[0] < 0 || state[0] >= si->dimK ||
        state[1] < 0 || state[1] >= si->dimL ||
        state[2] < 0 || state[2] >= si->dimM)
    {
        return -1;
    }

    return si->idx[(state[0] * si->dimL + state[1]) * si->dimM + state[2]];
}

static void mcmd_state_index_clear(t_mcmd_state_index *si)
{
    si->dimK = si->dimL = si->dimM = 0;
    si->idx  = NULL;
}

/* Build the table from num_states initial states, laid out with the given
 * stride in ints.  Where a state appears more than once the first occurrence
 * wins, matching the linear scans this replaces. */
static void mcmd_state_index_build(t_mcmd_state_index *si,
                                   const int *first_state, size_t stride,
                                   int num_states)
{
    int    i, total;
    size_t byte_stride = stride;

    mcmd_state_index_clear(si);

    if (num_states <= 0)
    {
        return;
    }

    for (i = 0; i < num_states; i++)
    {
        const int *st = (const int *)((const char *)first_state + i * byte_stride);

        if (st[0] + 1 > si->dimK) { si->dimK = st[0] + 1; }
        if (st[1] + 1 > si->dimL) { si->dimL = st[1] + 1; }
        if (st[2] + 1 > si->dimM) { si->dimM = st[2] + 1; }
    }

    total = si->dimK * si->dimL * si->dimM;
    snew(si->idx, total);
    for (i = 0; i < total; i++)
    {
        si->idx[i] = -1;
    }

    for (i = 0; i < num_states; i++)
    {
        const int *st = (const int *)((const char *)first_state + i * byte_stride);
        int        at = (st[0] * si->dimL + st[1]) * si->dimM + st[2];

        if (si->idx[at] == -1)
        {
            si->idx[at] = i;
        }
    }
}

void mcmd_atomdata_read(t_mcmd_atomdata *ad, int elem_idx,
                        gmx_bool bCollisional)
{
    char *path;
    int   i;

    ad->mass          = mcmd_idx2mass(elem_idx);
    ad->bPresent      = TRUE;
    ad->energy_levels = NULL;
    ad->rates         = NULL;
    ad->num_rates     = 0;
    ad->coll          = NULL;
    ad->num_coll      = 0;
    ad->weights.states     = NULL;
    ad->weights.weight     = NULL;
    ad->weights.num_states = 0;
    ad->max_transitions    = 0;
    mcmd_state_index_clear(&ad->rate_index);
    mcmd_state_index_clear(&ad->coll_index);

    if (ad->mass < 0)
    {
        gmx_fatal(FARGS, "Invalid element index %d.", elem_idx);
    }

    path = mcmd_data_path(MCMD_PATH_ENERGY, ad->mass);
    ad->energy_levels = mcmd_dict_read(path);
    sfree(path);

    path = mcmd_data_path(MCMD_PATH_RATES, ad->mass);
    ad->rates = mcmd_rates_read(path, &ad->num_rates);
    sfree(path);

    if (bCollisional)
    {
        path = mcmd_data_path(MCMD_PATH_COLL, ad->mass);
        ad->coll = mcmd_coll_read(path, &ad->num_coll);
        sfree(path);

        path = mcmd_data_path(MCMD_PATH_WEIGHT, ad->mass);
        mcmd_weights_read(&ad->weights, path);
        sfree(path);
    }

    /* Index the states once, so the per-event lookups are array reads rather
     * than scans over every state of the element (324 of them for sulphur). */
    if (ad->num_rates > 0)
    {
        mcmd_state_index_build(&ad->rate_index, ad->rates[0].initial_state,
                               sizeof(ad->rates[0]), ad->num_rates);
    }
    if (bCollisional && ad->num_coll > 0)
    {
        mcmd_state_index_build(&ad->coll_index, ad->coll[0].initial_state,
                               sizeof(ad->coll[0]), ad->num_coll);
    }

    /* A single event considers the radiative/Auger transitions of a state
     * plus, when enabled, its collisional ones, so the bound a caller needs
     * is the sum of the two maxima rather than the larger of them. */
    {
        int maxr = 0, maxc = 0;

        for (i = 0; i < ad->num_rates; i++)
        {
            if (ad->rates[i].num_transitions > maxr)
            {
                maxr = ad->rates[i].num_transitions;
            }
        }
        for (i = 0; i < ad->num_coll; i++)
        {
            if (ad->coll[i].num_transitions > maxc)
            {
                maxc = ad->coll[i].num_transitions;
            }
        }
        ad->max_transitions = maxr + maxc;
    }
}

void mcmd_atomdata_done(t_mcmd_atomdata *ad, gmx_bool bCollisional)
{
    if (!ad->bPresent)
    {
        return;
    }

    mcmd_dict_done(ad->energy_levels);
    mcmd_rates_done(ad->rates, ad->num_rates);

    if (bCollisional)
    {
        mcmd_coll_done(ad->coll, ad->num_coll);
        mcmd_weights_done(&ad->weights);
    }

    sfree(ad->rate_index.idx);
    sfree(ad->coll_index.idx);
    mcmd_state_index_clear(&ad->rate_index);
    mcmd_state_index_clear(&ad->coll_index);

    ad->bPresent      = FALSE;
    ad->energy_levels = NULL;
    ad->rates         = NULL;
    ad->coll          = NULL;
}

/* ================================================================ *
 * Exponential integral Ei(x)
 *
 * Taken unchanged from the public-domain implementation at
 * mymathlib.com (exponential_integral_Ei.c); only the wrapper name
 * differs.  Used by the collisional cross sections.
 * ================================================================ */

static long double xExponential_Integral_Ei(long double x);
static long double Continued_Fraction_Ei(long double x);
static long double Power_Series_Ei(long double x);
static long double Argument_Addition_Series_Ei(long double x);

static const long double epsilon = 10.0 * LDBL_EPSILON;

double mcmd_exp_integral_Ei(double x)
{
    return (double)xExponential_Integral_Ei((long double)x);
}

static long double xExponential_Integral_Ei(long double x)
{
    if (x < -5.0L)
        return Continued_Fraction_Ei(x);
    if (x == 0.0L)
        return -DBL_MAX;
    if (x < 6.8L)
        return Power_Series_Ei(x);
    if (x < 50.0L)
        return Argument_Addition_Series_Ei(x);
    return Continued_Fraction_Ei(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei(long double x)
{
    long double Am1 = 1.0L;
    long double A0 = 0.0L;
    long double Bm1 = 0.0L;
    long double B0 = 1.0L;
    long double a = expl(x);
    long double b = -x + 1.0L;
    long double Ap1 = b * A0 + a * Am1;
    long double Bp1 = b * B0 + a * Bm1;
    int j = 1;

    a = 1.0L;
    while (fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1))
    {
        if (fabsl(Bp1) > 1.0L)
        {
            Am1 = A0 / Bp1;
            A0 = Ap1 / Bp1;
            Bm1 = B0 / Bp1;
            B0 = 1.0L;
        }
        else
        {
            Am1 = A0;
            A0 = Ap1;
            Bm1 = B0;
            B0 = Bp1;
        }
        a = -j * j;
        b += 2.0L;
        Ap1 = b * A0 + a * Am1;
        Bp1 = b * B0 + a * Bm1;
        j += 1;
    }
    return (-Ap1 / Bp1);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei(long double x)
{
    long double xn = -x;
    long double Sn = -x;
    long double Sm1 = 0.0L;
    long double hsum = 1.0L;
    long double g = 0.5772156649015328606065121L;
    long double y = 1.0L;
    long double factorial = 1.0L;

    if (x == 0.0L)
        return (long double)-DBL_MAX;

    while (fabsl(Sn - Sm1) > epsilon * fabsl(Sm1))
    {
        Sm1 = Sn;
        y += 1.0L;
        xn *= (-x);
        factorial *= y;
        hsum += (1.0 / y);
        Sn += hsum * xn / factorial;
    }
    return (g + logl(fabsl(x)) - expl(x) * Sn);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
    static long double ei[] = {
        1.915047433355013959531e2L, 4.403798995348382689974e2L,
        1.037878290717089587658e3L, 2.492228976241877759138e3L,
        6.071406374098611507965e3L, 1.495953266639752885229e4L,
        3.719768849068903560439e4L, 9.319251363396537129882e4L,
        2.349558524907683035782e5L, 5.955609986708370018502e5L,
        1.516637894042516884433e6L, 3.877904330597443502996e6L,
        9.950907251046844760026e6L, 2.561565266405658882048e7L,
        6.612718635548492136250e7L, 1.711446713003636684975e8L,
        4.439663698302712208698e8L, 1.154115391849182948287e9L,
        3.005950906525548689841e9L, 7.842940991898186370453e9L,
        2.049649711988081236484e10L, 5.364511859231469415605e10L,
        1.405991957584069047340e11L, 3.689732094072741970640e11L,
        9.694555759683939661662e11L, 2.550043566357786926147e12L,
        6.714640184076497558707e12L, 1.769803724411626854310e13L,
        4.669055014466159544500e13L, 1.232852079912097685431e14L,
        3.257988998672263996790e14L, 8.616388199965786544948e14L,
        2.280446200301902595341e15L, 6.039718263611241578359e15L,
        1.600664914324504111070e16L, 4.244796092136850759368e16L,
        1.126348290166966760275e17L, 2.990444718632336675058e17L,
        7.943916035704453771510e17L, 2.111342388647824195000e18L,
        5.614329680810343111535e18L, 1.493630213112993142255e19L,
        3.975442747903744836007e19L, 1.058563689713169096306e20L};
    int k = (int)(x + 0.5);
    int j = 0;
    long double xx = (long double)k;
    long double dx = x - xx;
    long double xxj = xx;
    long double edx = expl(dx);
    long double Sm = 1.0L;
    long double Sn = (edx - 1.0L) / xxj;
    long double term = DBL_MAX;
    long double factorial = 1.0L;
    long double dxj = 1.0L;

    while (fabsl(term) > epsilon * fabsl(Sn))
    {
        j++;
        factorial *= (long double)j;
        xxj *= xx;
        dxj *= (-dx);
        Sm += (dxj / factorial);
        term = (factorial * (edx * Sm - 1.0L)) / xxj;
        Sn += term;
    }

    return ei[k - 7] + Sn * expl(xx);
}
