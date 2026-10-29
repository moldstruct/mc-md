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
#define MCMD_PATH_TOTAL  "./Atomic_data/total_energies_"
#define MCMD_PATH_RATES  "./Atomic_data/rate_transitions_"
#define MCMD_PATH_COLL   "./Atomic_data/collisional_parameters_"
#define MCMD_PATH_WEIGHT "./Atomic_data/statistical_weight_"
#define MCMD_PATH_COLEX  "./Atomic_data/collisional_excitation_"
#define MCMD_PATH_CTAB_ION "./Atomic_data/collisional_ionization_table_"
#define MCMD_PATH_CTAB_EXC "./Atomic_data/collisional_excitation_table_"

/* Longest line accepted in a data file.  A single initial state can have many
 * final states on one line, so this needs headroom. */
#define MCMD_LINEBUF 1024

/* ================================================================ *
 * Element table
 *
 * Adding an element means adding a row below (at the end: the row is the
 * element index, which checkpoints store), and, above argon, supplying its
 * Atomic_data files.  An atom is identified by the nearest standard atomic
 * mass, which tells argon (39.948) from calcium (40.078); the rounded mass
 * stays the key of mcmd_mass2idx().
 * ================================================================ */

typedef struct
{
    int    mass;     /* rounded atomic mass */
    double std_mass; /* standard atomic weight, for mcmd_atom2idx() */
    char   symbol[3];
    int    ground_state[MCMD_NSHELL]; /* K, L, M occupations */
} t_mcmd_element;

static const t_mcmd_element mcmd_elements[MCMD_NUM_ELEMENTS] = {
    {  1,  1.008,  "H",  {1, 0, 0}  },
    { 12, 12.011,  "C",  {2, 4, 0}  },
    { 14, 14.007,  "N",  {2, 5, 0}  },
    { 16, 15.999,  "O",  {2, 6, 0}  },
    { 19, 18.998,  "F",  {2, 7, 0}  },
    { 24, 24.305,  "MG", {2, 8, 2}  },
    { 31, 30.974,  "P",  {2, 8, 5}  },
    { 32, 32.06,   "S",  {2, 8, 6}  },
    { 35, 35.45,   "CL", {2, 8, 7}  },
    { 40, 40.078,  "CA", {2, 8, 10} },
    { 56, 55.845,  "FE", {2, 8, 16} },
    { 59, 58.693,  "NI", {2, 8, 18} },
    { 28, 28.085,  "SI", {2, 8, 4}  },
    { 23, 22.990,  "NA", {2, 8, 1}  },
    /* Iodine's configuration is a placeholder copied from nickel: its bonded
     * behaviour is right but its ionization physics is not. */
    { 127, 126.904, "I", {2, 8, 18} },
    /* 2026-10-29: the rest of H to Ar, which the built-in data cover, except helium: a mass ~4 is far more
     * likely a heavy hydrogen (pdb2gmx -heavyh, mass repartitioning), which is refused rather than taken as
     * helium */
    {  7,  6.94,   "LI", {2, 1, 0}  },
    {  9,  9.0122, "BE", {2, 2, 0}  },
    { 11, 10.81,   "B",  {2, 3, 0}  },
    { 20, 20.180,  "NE", {2, 8, 0}  },
    { 27, 26.982,  "AL", {2, 8, 3}  },
    { 40, 39.948,  "AR", {2, 8, 8}  }
};

/* Deuterium.  Deliberately not a row of its own above: it is chemically
 * hydrogen and has hydrogen's electronic structure, so it shares hydrogen's
 * ground state and hydrogen's Atomic_data rate files.  A separate row would
 * mean a separate element index, hence a demand for energy_levels_D.txt that
 * does not exist, and it would no longer compare equal to mc->elem_H in the
 * charge-transfer pass. */
#define MCMD_MASS_D 2
#define MCMD_MASS_H 1

int mcmd_mass2idx(int mass)
{
    int i;

    /* Fold the isotope onto its element.  Only the species identity is
     * remapped; mdatoms->massT still carries the real mass 2, so a deuteron
     * recoils at its own mass. */
    if (mass == MCMD_MASS_D)
    {
        mass = MCMD_MASS_H;
    }

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

int mcmd_atom2idx(double mass)
{
    int    i, best = -1;
    double dbest = 0.5;                 /* within half a mass unit, as the rounded key was */

    if (mass >= 1.5 && mass < 2.5)
    {
        return mcmd_mass2idx(MCMD_MASS_D);
    }
    for (i = 0; i < MCMD_NUM_ELEMENTS; i++)
    {
        double d = fabs(mass - mcmd_elements[i].std_mass);

        if (d < dbest)
        {
            dbest = d;
            best  = i;
        }
    }
    return best;
}

const char *mcmd_idx2symbol(int idx)
{
    return (idx >= 0 && idx < MCMD_NUM_ELEMENTS) ? mcmd_elements[idx].symbol : "_";
}

const int *mcmd_idx_ground_state(int idx)
{
    return (idx >= 0 && idx < MCMD_NUM_ELEMENTS) ? mcmd_elements[idx].ground_state : NULL;
}

const char *mcmd_mass2symbol(int mass)
{
    return mcmd_idx2symbol(mcmd_mass2idx(mass));
}

const int *mcmd_ground_state(int mass)
{
    return mcmd_idx_ground_state(mcmd_mass2idx(mass));
}

/* Build "<prefix><SYMBOL>.txt".  Caller frees. */
static char *mcmd_data_path(const char *prefix, int elem)
{
    const char *symbol = mcmd_idx2symbol(elem);
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

static t_mcmd_dict *mcmd_dict_read(const char *path, gmx_bool bRequired)
{
    t_mcmd_dict *dict;
    FILE        *fp;
    int          state[MCMD_NSHELL];
    double       energy;
    int          i;

    fp = fopen(path, "r");
    if (fp == NULL)
    {
        if (!bRequired)
        {
            /* An optional table: absent is a supported configuration, not an
             * error.  The caller falls back to what it did before the file
             * existed. */
            return NULL;
        }
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

/* ---------------------------------------------------------------- *
 * Tabulated collisional rate coefficients (see atomic_data.h).  Rows hold one
 * value per grid temperature, so a line can be long: read it whole.
 * ---------------------------------------------------------------- */

/* One line of any length into *buf (grown as needed); FALSE at end of file. */
static gmx_bool mcmd_read_long_line(FILE *fp, char **buf, size_t *cap)
{
    size_t n = 0;
    int    ch;

    if (*buf == NULL)
    {
        *cap = 4096;
        snew(*buf, *cap);
    }
    while ((ch = fgetc(fp)) != EOF && ch != '\n')
    {
        if (n + 2 > *cap)
        {
            *cap *= 2;
            srenew(*buf, *cap);
        }
        (*buf)[n++] = (char)ch;
    }
    (*buf)[n] = '\0';
    return (ch != EOF || n > 0);
}

static gmx_bool mcmd_blank(const char *s)
{
    for (; *s != '\0'; s++)
    {
        if (*s != ' ' && *s != '\t' && *s != '\r')
        {
            return FALSE;
        }
    }
    return TRUE;
}

/* ln of a non-negative number written as text, any exponent ('3.1e-2650' is fine: the
 * coefficients include exp(-dE/T), far below the smallest double at low T_e).  -HUGE_VAL for 0;
 * *end == p if there is no number. */
static double mcmd_parse_ln(const char *p, char **end)
{
    const char *s = p, *e;
    char        mant[64];
    double      m;
    long        x = 0;
    size_t      n;

    while (*s == ' ' || *s == '\t') { s++; }
    e = s;
    while (*e != '\0' && *e != ' ' && *e != '\t' && *e != '\r' && *e != 'e' && *e != 'E') { e++; }
    n = (size_t)(e - s);
    if (n == 0 || n >= sizeof(mant))
    {
        *end = (char *)p;
        return 0;
    }
    memcpy(mant, s, n);
    mant[n] = '\0';
    m = strtod(mant, end);
    if (*end == mant)
    {
        *end = (char *)p;
        return 0;
    }
    *end = (char *)e;
    if (*e == 'e' || *e == 'E')
    {
        char *q;

        x = strtol(e + 1, &q, 10);
        *end = q;
    }
    if (m < 0 || !isfinite(m))
    {
        return NAN;
    }
    return (m == 0) ? -HUGE_VAL : log(m) + (double)x * M_LN10;
}

/* NULL if the file does not exist (the fit file is used instead). */
static t_mcmd_ctab *mcmd_ctab_read(const char *path)
{
    FILE        *fp = fopen(path, "r");
    t_mcmd_ctab *tab;
    char        *line = NULL, *p, *q;
    size_t       cap  = 0;
    int          cap_tr = 0, j;

    if (fp == NULL)
    {
        return NULL;
    }
    snew(tab, 1);
    /* header: 'T_e t1 t2 ... tn' */
    do
    {
        if (!mcmd_read_long_line(fp, &line, &cap))
        {
            gmx_fatal(FARGS, "%s: empty file (expected a 'T_e t1 ... tn' header line).", path);
        }
    }
    while (mcmd_blank(line));
    p = line;
    while (*p == ' ' || *p == '\t') { p++; }
    if (strncmp(p, "T_e", 3) != 0)
    {
        gmx_fatal(FARGS, "%s: the first line must be 'T_e t1 t2 ... tn' (eV).", path);
    }
    p += 3;
    for (;;)
    {
        double t = strtod(p, &q);

        if (q == p)
        {
            break;
        }
        if (!(t > 0) || (tab->nT > 0 && log(t) <= tab->lnT[tab->nT - 1]))
        {
            gmx_fatal(FARGS, "%s: the T_e grid must be positive and increasing (at %g).", path, t);
        }
        srenew(tab->lnT, tab->nT + 1);
        tab->lnT[tab->nT++] = log(t);
        p = q;
    }
    if (tab->nT < 2)
    {
        gmx_fatal(FARGS, "%s: the T_e grid needs at least two temperatures.", path);
    }
    /* rows: 'K L M K' L' M' k1 ... kn' */
    while (mcmd_read_long_line(fp, &line, &cap))
    {
        int st[2 * MCMD_NSHELL];

        if (mcmd_blank(line))
        {
            continue;
        }
        p = line;
        for (j = 0; j < 2 * MCMD_NSHELL; j++)
        {
            long v = strtol(p, &q, 10);

            if (q == p)
            {
                gmx_fatal(FARGS, "%s, transition %d: expected 'K L M K' L' M' k1 ... k%d'.", path, tab->ntr + 1,
                          tab->nT);
            }
            st[j] = (int)v;
            p     = q;
        }
        if (tab->ntr == cap_tr)
        {
            cap_tr = (cap_tr == 0) ? 64 : 2 * cap_tr;
            srenew(tab->ini, cap_tr);
            srenew(tab->fin, cap_tr);
            srenew(tab->lnk, (size_t)cap_tr * tab->nT);
        }
        for (j = 0; j < MCMD_NSHELL; j++)
        {
            tab->ini[tab->ntr][j] = st[j];
            tab->fin[tab->ntr][j] = st[MCMD_NSHELL + j];
        }
        for (j = 0; j < tab->nT; j++)
        {
            double lk = mcmd_parse_ln(p, &q);

            if (q == p || isnan(lk) || lk == HUGE_VAL)
            {
                gmx_fatal(FARGS, "%s, transition %d: needs %d rate coefficients >= 0 (cm^3/s), one per T_e.",
                          path, tab->ntr + 1, tab->nT);
            }
            tab->lnk[(size_t)tab->ntr * tab->nT + j] = lk;
            p = q;
        }
        while (*p == ' ' || *p == '\t' || *p == '\r') { p++; }
        if (*p != '\0')
        {
            gmx_fatal(FARGS, "%s, transition %d: more than %d rate coefficients.", path, tab->ntr + 1, tab->nT);
        }
        tab->ntr++;
    }
    fclose(fp);
    sfree(line);
    return tab;
}

static void mcmd_ctab_done(t_mcmd_ctab *tab)
{
    if (tab == NULL)
    {
        return;
    }
    sfree(tab->lnT); sfree(tab->ini); sfree(tab->fin); sfree(tab->lnk);
    sfree(tab);
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
 * Built-in screened-hydrogenic data (mcmd-atomic-data = hydrogenic)
 *
 * The same tables as Atomic_data/, computed at start-up instead of read: a
 * port of pqn_model.py (screened-hydrogenic shell model after the Appendix
 * of H.A. Scott & S.B. Hansen, High Energy Density Physics 6, 39 (2010), with
 * the screening constants of R. Marchand, S. Caille & Y.T. Lee, JQSRT 43, 149
 * (1990), level energies scaled to the NIST ionization energies) followed by
 * the conversion of make_atomic_data_cretin.py.  Every configuration of
 * shells n <= 3 (K <= 2, L <= 8, M <= 18) of every ionization stage is a
 * state, so Z <= 18 (no 4s; helium is not in the element table).  It
 * reproduces the CRETIN screened-hydrogenic models of the earlier production
 * data (rates within ~10 %; lysozyme at 1000 eV: charge and ion directions
 * unchanged, ion kinetic energy 3.5 % lower; manual, appendix C).
 *
 * Per state: the photoionization cross section out of each occupied shell at
 * the photon energy (Kramers, constant Gaunt factor, sigma ~ E^-3 up to the
 * fit limit of the models), the fluorescence rate of each radiative
 * transition into a lower shell, and the Auger rates summed per final state;
 * the outermost binding energy, the total binding energy and the statistical
 * weight.
 * ================================================================ */

#define MCMD_HYD_NMAX        3
#define MCMD_HYD_ZMAX        18
#define MCMD_HYD_RY          13.605693          /* eV */
#define MCMD_HYD_ALPHA       (1.0 / 137.035999)
#define MCMD_HYD_H_EV_S      4.135667696e-15    /* Planck constant, eV s */
#define MCMD_HYD_HC_EV_A     12398.419843       /* eV Angstrom */
#define MCMD_HYD_AUGER_GAUNT 0.12
#define MCMD_HYD_KRAMERS     7.907e-18          /* cm^2, hydrogen 1s at threshold */
#define MCMD_HYD_PHIS_EMAX   1.0e5              /* eV, upper end of the photoionization fits */
#define MCMD_HYD_C_CM_S      29979245800.0
#define MCMD_HYD_E_ESU       4.80320471e-10
#define MCMD_HYD_M_G         9.1093837e-28
#define MCMD_HYD_AUTOION     1e14               /* energy_levels placeholder for an autoionizing state */

static const int mcmd_hyd_cap[MCMD_HYD_NMAX] = { 2, 8, 18 };

/* Marchand, Caille & Lee (1990), Table 2: screening of shell n by one
 * electron in shell m <= n, [n-1][m-1]. */
static const double mcmd_hyd_sigma[MCMD_HYD_NMAX][MCMD_HYD_NMAX] = {
    { 0.5966, 0,      0      },
    { 0.8597, 0.6888, 0      },
    { 0.9923, 0.8877, 0.7322 }
};

/* Exact hydrogenic absorption oscillator strengths n -> m, [n-1][m-1]. */
static const double mcmd_hyd_fosc[MCMD_HYD_NMAX][MCMD_HYD_NMAX] = {
    { 0, 0.4162, 0.07910 },
    { 0, 0,      0.6407  },
    { 0, 0,      0       }
};

/* NIST ASD ground-state ionization energies (eV), [Z-1][nelec-1], nelec = bound electrons. */
static const double mcmd_hyd_ip[MCMD_HYD_ZMAX][MCMD_HYD_ZMAX] = {
    { 13.5984346 },
    { 54.41776553, 24.58738901 },
    { 122.4543591, 75.640097, 5.391714996 },
    { 217.7185846, 153.896205, 18.21115, 9.322699 },
    { 340.2260225, 259.374379, 37.93059, 25.15483, 8.298019 },
    { 489.9932078, 392.09056, 64.49352, 47.88778, 24.383143, 11.260288 },
    { 667.0461377, 552.06741, 97.8901, 77.4735, 47.4453, 29.60125, 14.53413 },
    { 871.4099138, 739.32697, 138.1189, 113.899, 77.4135, 54.93554, 35.12112, 13.618055 },
    { 1103.11753, 953.8983, 185.1868, 157.16311, 114.249, 87.175, 62.70798, 34.97081, 17.42282 },
    { 1362.199256, 1195.8082, 239.097, 207.271, 157.934, 126.247, 97.19, 63.4233, 40.96297, 21.564541 },
    { 1648.702285, 1465.0992, 299.856, 264.192, 208.504, 172.23, 138.404, 98.936, 71.62, 47.28636,
      5.13907696 },
    { 1962.663889, 1761.8049, 367.489, 327.99, 265.924, 225.02, 186.76, 141.33, 109.2654, 80.1436,
      15.035271, 7.646236 },
    { 2304.140359, 2085.97693, 442.005, 398.65, 330.21, 284.64, 241.76, 190.49, 153.8252, 119.9924,
      28.447642, 18.82855, 5.985769 },
    { 2673.177958, 2437.65805, 523.415, 476.273, 401.38, 351.28, 303.59, 246.57, 205.279, 166.767,
      45.14179, 33.493, 16.34585, 8.15168 },
    { 3069.842145, 2816.90868, 611.741, 560.62, 479.44, 424.4, 372.31, 309.6, 263.57, 220.43,
      65.02511, 51.44387, 30.20264, 19.76949, 10.486686 },
    { 3494.188518, 3223.78057, 706.994, 651.96, 564.41, 504.55, 447.7, 379.84, 328.794, 280.954,
      88.0529, 72.5945, 47.222, 34.86, 23.33788, 10.3600167 },
    { 3946.29179, 3658.34366, 809.198, 750.23, 656.3, 591.58, 530, 456.7, 400.851, 348.306,
      114.2013, 96.94, 67.68, 53.24, 39.8, 23.81364, 12.967633 },
    { 4426.22407, 4120.66559, 918.375, 855.5, 755.13, 685.5, 619, 540.4, 479.76, 422.6,
      143.4567, 124.41, 91.29, 74.84, 59.58, 40.735, 27.62967, 15.7596119 }
};

/* Nuclear charge of the elements in mcmd_elements[] the built-in data cover; 0 for the others. */
static int mcmd_hyd_Z(const char *symbol)
{
    static const char *sym[MCMD_HYD_ZMAX] = {
        "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR"
    };
    int z;

    for (z = 0; z < MCMD_HYD_ZMAX; z++)
    {
        if (strcmp(sym[z], symbol) == 0)
        {
            return z + 1;
        }
    }
    return 0;
}

static double mcmd_hyd_screened_charge(int Z, const int c[MCMD_NSHELL], int n)
{
    double q = Z + 0.5 * mcmd_hyd_sigma[n - 1][n - 1];
    int    m;

    for (m = 1; m <= n; m++)
    {
        if (c[m - 1] > 0)
        {
            q -= ((m == n) ? 0.5 : 1.0) * mcmd_hyd_sigma[n - 1][m - 1] * c[m - 1];
        }
    }
    return q;
}

/* Unscaled total binding energy (eV), with the relativistic corrections. */
static double mcmd_hyd_binding(int Z, const int c[MCMD_NSHELL])
{
    double total = 0;
    int    n;

    for (n = 1; n <= MCMD_HYD_NMAX; n++)
    {
        if (c[n - 1] > 0)
        {
            double q  = mcmd_hyd_screened_charge(Z, c, n) / n;
            double q2 = q * q;
            double rf = 2.0 * n / (n + 1) - 0.75;
            double a2 = MCMD_HYD_ALPHA * MCMD_HYD_ALPHA;

            total += c[n - 1] * MCMD_HYD_RY * q2 * (1.0 + a2 * q2 * rf + a2 * a2 * q2 * q2 / 6.0);
        }
    }
    return total;
}

static int mcmd_hyd_outermost(const int c[MCMD_NSHELL])
{
    int n;

    for (n = MCMD_HYD_NMAX; n > 1; n--)
    {
        if (c[n - 1] > 0)
        {
            return n;
        }
    }
    return 1;
}

static double mcmd_hyd_weight(const int c[MCMD_NSHELL])
{
    double g = 1;
    int    n, k;

    for (n = 0; n < MCMD_HYD_NMAX; n++)
    {
        for (k = 1; k <= c[n]; k++)                 /* binomial(cap, c) */
        {
            g = g * (mcmd_hyd_cap[n] - c[n] + k) / k;
        }
    }
    return floor(g + 0.5);
}

/* Absorption oscillator strength n -> m (n < m) from lower to upper, eq. (A6-7), without the occupation factors. */
static double mcmd_hyd_oscillator(int Z, const int lower[MCMD_NSHELL], const int upper[MCMD_NSHELL], int n, int m)
{
    double de = mcmd_hyd_binding(Z, lower) - mcmd_hyd_binding(Z, upper);
    double ql, qu, az2, rel;

    if (de <= 0)
    {
        return 0;
    }
    ql  = mcmd_hyd_screened_charge(Z, lower, n) / Z;
    qu  = mcmd_hyd_screened_charge(Z, upper, m) / Z;
    az2 = (MCMD_HYD_ALPHA * Z) * (MCMD_HYD_ALPHA * Z);
    rel = (1 + az2 / (n * n)) / (n * n) - (1 + az2 / (m * m)) / (m * m);
    return mcmd_hyd_fosc[n - 1][m - 1] * pow(ql, 4) * qu * qu * pow(Z * Z * MCMD_HYD_RY / de, 3) * pow(rel, 3);
}

/* The model of one element: ground configurations, scaling factors and energy offsets per stage. */
typedef struct
{
    int    Z;
    int    ground[MCMD_HYD_ZMAX + 1][MCMD_NSHELL];
    double eb_ground[MCMD_HYD_ZMAX + 1];
    double ip[MCMD_HYD_ZMAX + 1];        /* ionization energy of the stage with nelec electrons */
    double scale[MCMD_HYD_ZMAX + 1];
    double offset[MCMD_HYD_ZMAX + 1];    /* energy to strip the neutral atom down to nelec electrons */
} t_mcmd_hyd;

static void mcmd_hyd_init(t_mcmd_hyd *h, int Z)
{
    int ne, n, m;

    h->Z = Z;
    for (ne = 0; ne <= Z; ne++)
    {
        int left = ne;

        for (n = 0; n < MCMD_NSHELL; n++)        /* shells fill in order up to Z = 18 */
        {
            h->ground[ne][n] = (left < mcmd_hyd_cap[n]) ? left : mcmd_hyd_cap[n];
            left            -= h->ground[ne][n];
        }
        h->eb_ground[ne] = mcmd_hyd_binding(Z, h->ground[ne]);
    }
    h->ip[0]    = 0;
    h->scale[0] = 1;
    for (ne = 1; ne <= Z; ne++)
    {
        h->ip[ne]    = mcmd_hyd_ip[Z - 1][ne - 1];
        h->scale[ne] = h->ip[ne] / (h->eb_ground[ne] - h->eb_ground[ne - 1]);
    }
    for (ne = 0; ne <= Z; ne++)
    {
        h->offset[ne] = 0;
        for (m = ne + 1; m <= Z; m++)
        {
            h->offset[ne] += h->ip[m];
        }
    }
}

/* Excitation energy above the stage's ground configuration, eqs. (A4-5): scaled with the ionization energy of
 * the stage the lowest inner-shell hole belongs to. */
static double mcmd_hyd_excitation(const t_mcmd_hyd *h, const int c[MCMD_NSHELL])
{
    int ne = c[0] + c[1] + c[2];
    int valence = mcmd_hyd_outermost(h->ground[ne]);
    int hole = 0, n, k = ne;

    for (n = 1; n < valence; n++)
    {
        if (c[n - 1] < h->ground[ne][n - 1])
        {
            hole = n;
            break;
        }
    }
    if (hole == 1)
    {
        k = (ne < 2) ? ne : 2;
    }
    else if (hole == 2)
    {
        k = (ne < 10) ? ne : 10;
    }
    return h->scale[k] * (h->eb_ground[ne] - mcmd_hyd_binding(h->Z, c));
}

static double mcmd_hyd_absolute(const t_mcmd_hyd *h, const int c[MCMD_NSHELL])
{
    return h->offset[c[0] + c[1] + c[2]] + mcmd_hyd_excitation(h, c);
}

static gmx_bool mcmd_hyd_valid(const t_mcmd_hyd *h, const int c[MCMD_NSHELL])
{
    int n;

    for (n = 0; n < MCMD_NSHELL; n++)
    {
        if (c[n] < 0 || c[n] > mcmd_hyd_cap[n])
        {
            return FALSE;
        }
    }
    return (c[0] + c[1] + c[2] <= h->Z);
}

static void mcmd_hyd_add(t_mcmd_rate *r, const int f[MCMD_NSHELL], double rate, int type, gmx_bool bSum)
{
    int i = r->num_transitions;

    if (bSum)                                   /* Auger: one entry per final state */
    {
        for (i = 0; i < r->num_transitions; i++)
        {
            if (r->types[i] == type && mcmd_state_equal(r->final_states[i], f))
            {
                r->rates[i] += rate;
                return;
            }
        }
    }
    srenew(r->final_states, i + 1);
    srenew(r->rates, i + 1);
    srenew(r->types, i + 1);
    snew(r->final_states[i], MCMD_NSHELL);
    r->final_states[i][0] = f[0];
    r->final_states[i][1] = f[1];
    r->final_states[i][2] = f[2];
    r->rates[i]           = rate;
    r->types[i]           = type;
    r->num_transitions    = i + 1;
}

/* Photoionization, fluorescence and Auger out of state c. */
static void mcmd_hyd_transitions(const t_mcmd_hyd *h, const int c[MCMD_NSHELL], double energy, t_mcmd_rate *r)
{
    const int Z     = h->Z;
    const int ne    = c[0] + c[1] + c[2];
    double    start = mcmd_hyd_absolute(h, c);
    int       n, m, u, w, s;

    r->initial_state[0] = c[0];
    r->initial_state[1] = c[1];
    r->initial_state[2] = c[2];
    r->num_transitions  = 0;
    r->final_states     = NULL;
    r->rates            = NULL;
    r->types            = NULL;

    /* photoionization out of each occupied shell, eq. (A8); zero outside [threshold, fit limit] */
    for (n = 1; ne >= 1 && n <= MCMD_HYD_NMAX; n++)
    {
        int    f[MCMD_NSHELL] = { c[0], c[1], c[2] };
        double thr, emax, q;

        if (c[n - 1] == 0)
        {
            continue;
        }
        f[n - 1]--;
        thr = mcmd_hyd_absolute(h, f) - start;
        if (thr <= 0)
        {
            continue;
        }
        emax = (MCMD_HYD_PHIS_EMAX > 2.0 * thr) ? MCMD_HYD_PHIS_EMAX : 10.0 * thr;
        if (energy < thr || energy > emax)
        {
            continue;
        }
        q = mcmd_hyd_screened_charge(Z, c, n);
        if (q < 1e-3)
        {
            q = 1e-3;
        }
        mcmd_hyd_add(r, f, MCMD_HYD_KRAMERS * pow(q, 4) / pow(n, 5) * pow(MCMD_HYD_RY / energy, 3) * c[n - 1], 2,
                     FALSE);
    }

    /* fluorescence: c is the upper level of the absorption n -> m from lower = c with one electron m -> n;
     * A = 8 pi^2 e^2 nu^2 / (m_e c^3) (g_lower / g_upper) f */
    for (m = 2; m <= MCMD_HYD_NMAX; m++)
    {
        for (n = 1; n < m; n++)
        {
            int    lo[MCMD_NSHELL] = { c[0], c[1], c[2] };
            double photon, f, nu;

            if (c[m - 1] == 0)
            {
                continue;
            }
            lo[m - 1]--;
            lo[n - 1]++;
            if (!mcmd_hyd_valid(h, lo))
            {
                continue;
            }
            photon = start - mcmd_hyd_absolute(h, lo);
            if (photon <= 0)
            {
                continue;
            }
            f = mcmd_hyd_oscillator(Z, lo, c, n, m) * lo[n - 1] * (1.0 - (double)lo[m - 1] / mcmd_hyd_cap[m - 1]);
            if (f <= 0)
            {
                continue;
            }
            nu = MCMD_HYD_C_CM_S / (MCMD_HYD_HC_EV_A / photon * 1e-8);
            mcmd_hyd_add(r, lo, 8 * M_PI * M_PI * MCMD_HYD_E_ESU * MCMD_HYD_E_ESU * nu * nu /
                         (MCMD_HYD_M_G * MCMD_HYD_C_CM_S * MCMD_HYD_C_CM_S * MCMD_HYD_C_CM_S) *
                         mcmd_hyd_weight(lo) / mcmd_hyd_weight(c) * f, 1, FALSE);
        }
    }

    /* Auger, eq. (A12): the electron of shell s leaves while one of shell u drops into w < u;
     * A = 64 pi/sqrt3 q^2 Ry^2/h G g_f/g_a f_wu P_w (1 - P_u/2u^2)(1 - P_s/2s^2) / (E n_s^3), occupations in the
     * final state, q the charge of the final ion, E the energy of the u -> w jump */
    if (ne >= 2)
    {
        const double qf   = Z - (ne - 1);
        const double pref = 64 * M_PI / sqrt(3.0) * qf * qf * MCMD_HYD_RY * MCMD_HYD_RY / MCMD_HYD_H_EV_S *
            MCMD_HYD_AUGER_GAUNT / mcmd_hyd_weight(c);

        for (s = 1; s <= MCMD_HYD_NMAX; s++)
        {
            int    core[MCMD_NSHELL] = { c[0], c[1], c[2] };
            double ecore;

            if (c[s - 1] == 0)
            {
                continue;
            }
            core[s - 1]--;
            ecore = mcmd_hyd_absolute(h, core);
            for (u = 2; u <= MCMD_HYD_NMAX; u++)
            {
                for (w = 1; w < u; w++)
                {
                    int    fin[MCMD_NSHELL] = { core[0], core[1], core[2] };
                    double efin, jump, f, rate;

                    /* room in w is checked in the initial state, as in the models */
                    if (core[u - 1] == 0 || c[w - 1] >= mcmd_hyd_cap[w - 1])
                    {
                        continue;
                    }
                    fin[u - 1]--;
                    fin[w - 1]++;
                    efin = mcmd_hyd_absolute(h, fin);
                    jump = ecore - efin;
                    if (start - efin <= 0 || jump <= 0)
                    {
                        continue;
                    }
                    f    = mcmd_hyd_oscillator(Z, fin, core, w, u);
                    rate = pref * mcmd_hyd_weight(fin) / jump * f / (s * s * s) * fin[w - 1] *
                        (1.0 - (double)fin[u - 1] / mcmd_hyd_cap[u - 1]) *
                        (1.0 - (double)fin[s - 1] / mcmd_hyd_cap[s - 1]);
                    if (rate > 0)
                    {
                        mcmd_hyd_add(r, fin, rate, 0, TRUE);
                    }
                }
            }
        }
    }
}

static t_mcmd_dict *mcmd_dict_new(void)
{
    t_mcmd_dict *dict;
    int          i;

    snew(dict, 1);
    for (i = 0; i < MCMD_DICT_SIZE; i++)
    {
        dict->bucket[i] = NULL;
    }
    return dict;
}

/* Statistical weights of every initial state of rates[], from the shell combinatorics. */
static void mcmd_weights_compute(t_mcmd_weights *w, const t_mcmd_rate *rates, int num_rates)
{
    int i;

    w->num_states = num_rates;
    snew(w->states, num_rates);
    snew(w->weight, num_rates);
    for (i = 0; i < num_rates; i++)
    {
        snew(w->states[i], MCMD_NSHELL);
        w->states[i][0] = rates[i].initial_state[0];
        w->states[i][1] = rates[i].initial_state[1];
        w->states[i][2] = rates[i].initial_state[2];
        w->weight[i]    = mcmd_hyd_weight(rates[i].initial_state);
    }
}

/* A note on stderr where statistical_weight_X.txt departs from the shell combinatorics (it is used as given). */
static void mcmd_weights_compare(const t_mcmd_weights *w, const char *path)
{
    int i, n = 0, first = -1;

    for (i = 0; i < w->num_states; i++)
    {
        double g = mcmd_hyd_weight(w->states[i]);

        if (fabs(w->weight[i] - g) > 1e-5 * g)          /* 6-digit printing is fine */
        {
            if (first < 0)
            {
                first = i;
            }
            n++;
        }
    }
    if (n > 0)
    {
        fprintf(stderr, "\nNOTE: %s: %d of %d weights differ from C(2,K) C(8,L) C(18,M) (e.g. %d %d %d: %g, not %g); "
                "used as given.\n", path, n, w->num_states, w->states[first][0], w->states[first][1],
                w->states[first][2], w->weight[first], mcmd_hyd_weight(w->states[first]));
    }
}

/* Fill ad as mcmd_atomdata_read() does from files: energies, rates, weights. */
static void mcmd_hyd_fill(t_mcmd_atomdata *ad, double energy)
{
    const char *symbol = mcmd_idx2symbol(ad->elem);
    const int   Z      = mcmd_hyd_Z(symbol);
    t_mcmd_hyd  h;
    int         c[MCMD_NSHELL], ns = 0, n;

    if (Z == 0)
    {
        gmx_fatal(FARGS, "mcmd-atomic-data = hydrogenic covers H to Ar but He (shells K, L, M); %s needs "
                  "mcmd-atomic-data = files with its own Atomic_data.", symbol);
    }
    if (!(energy > 0))
    {
        gmx_fatal(FARGS, "mcmd-atomic-data = hydrogenic needs the photon energy (mcmd-pulse-photon-energy > 0).");
    }
    mcmd_hyd_init(&h, Z);

    ad->energy_levels  = mcmd_dict_new();
    ad->total_energies = mcmd_dict_new();
    for (n = 0; n < 2; n++)            /* pass 0 counts the states, pass 1 fills */
    {
        int i = 0;

        if (n == 1)
        {
            ad->num_rates          = ns;
            ad->weights.num_states = ns;
            snew(ad->rates, ns);
            snew(ad->weights.states, ns);
            snew(ad->weights.weight, ns);
        }
        for (c[0] = 0; c[0] <= mcmd_hyd_cap[0]; c[0]++)
        {
            for (c[1] = 0; c[1] <= mcmd_hyd_cap[1]; c[1]++)
            {
                for (c[2] = 0; c[2] <= mcmd_hyd_cap[2]; c[2]++)
                {
                    int    ne = c[0] + c[1] + c[2];
                    double exc, tg = 0, bd;
                    int    k;

                    if (ne > Z)
                    {
                        continue;
                    }
                    if (n == 0)
                    {
                        ns++;
                        continue;
                    }
                    exc = mcmd_hyd_excitation(&h, c);
                    for (k = 1; k <= ne; k++)
                    {
                        tg += h.ip[k];
                    }
                    bd = h.ip[ne] - exc;
                    mcmd_dict_insert(ad->energy_levels, c, (bd < 0) ? MCMD_HYD_AUTOION : bd);
                    mcmd_dict_insert(ad->total_energies, c, tg - exc);
                    snew(ad->weights.states[i], MCMD_NSHELL);
                    ad->weights.states[i][0] = c[0];
                    ad->weights.states[i][1] = c[1];
                    ad->weights.states[i][2] = c[2];
                    ad->weights.weight[i]    = mcmd_hyd_weight(c);
                    mcmd_hyd_transitions(&h, c, energy, &ad->rates[i]);
                    i++;
                }
            }
        }
    }
}

void mcmd_atomdata_write(const t_mcmd_atomdata *ad, const char *dir)
{
    const char *sym = mcmd_idx2symbol(ad->elem);
    const char *name[4] = { "energy_levels", "total_energies", "rate_transitions", "statistical_weight" };
    FILE       *fp[4];
    char        path[512];
    int         f, i, t;

    for (f = 0; f < 4; f++)
    {
        sprintf(path, "%s/%s_%s.txt", dir, name[f], sym);
        if ((fp[f] = fopen(path, "w")) == NULL)
        {
            gmx_fatal(FARGS, "Could not open '%s' for writing.", path);
        }
    }
    for (i = 0; i < ad->num_rates; i++)
    {
        const t_mcmd_rate *r = &ad->rates[i];
        const int         *s = r->initial_state;
        double             g = mcmd_weight_for_state(&ad->weights, s);

        fprintf(fp[0], "%d %d %d %.17g\n", s[0], s[1], s[2], mcmd_dict_get(ad->energy_levels, s));
        fprintf(fp[1], "%d %d %d %.17g\n", s[0], s[1], s[2], mcmd_dict_get(ad->total_energies, s));
        fprintf(fp[2], "%d %d %d ;", s[0], s[1], s[2]);
        for (t = 0; t < r->num_transitions; t++)
        {
            fprintf(fp[2], "%d %d %d %.17g %d ;", r->final_states[t][0], r->final_states[t][1],
                    r->final_states[t][2], r->rates[t], r->types[t]);
        }
        fprintf(fp[2], "\n");
        if (g != MCMD_NOT_FOUND)
        {
            fprintf(fp[3], "%d %d %d %.17g\n", s[0], s[1], s[2], g);
        }
    }
    for (f = 0; f < 4; f++)
    {
        fclose(fp[f]);
    }
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
                        int collmode, int source, double photon_energy)
{
    char *path;
    int   i;

    ad->mass          = mcmd_idx2mass(elem_idx);
    ad->elem          = elem_idx;
    ad->bWeightsFile  = FALSE;
    ad->bPresent      = TRUE;
    ad->energy_levels  = NULL;
    ad->total_energies = NULL;
    ad->rates          = NULL;
    ad->num_rates      = 0;
    ad->coll          = NULL;
    ad->num_coll      = 0;
    ad->weights.states     = NULL;
    ad->weights.weight     = NULL;
    ad->weights.num_states = 0;
    ad->max_transitions    = 0;
    ad->colex              = NULL;
    ad->num_colex          = 0;
    ad->coll_tab           = NULL;
    ad->colex_tab          = NULL;
    mcmd_state_index_clear(&ad->rate_index);
    mcmd_state_index_clear(&ad->coll_index);

    if (ad->mass < 0)
    {
        gmx_fatal(FARGS, "Invalid element index %d.", elem_idx);
    }

    if (source == 1)
    {
        /* built in: no file is read (mcmd-coll-data = files and excitation are refused by the caller) */
        mcmd_hyd_fill(ad, photon_energy);
    }
    else
    {
        path = mcmd_data_path(MCMD_PATH_ENERGY, ad->elem);
        ad->energy_levels = mcmd_dict_read(path, TRUE);
        sfree(path);

        path = mcmd_data_path(MCMD_PATH_TOTAL, ad->elem);
        ad->total_energies = mcmd_dict_read(path, FALSE);
        sfree(path);

        path = mcmd_data_path(MCMD_PATH_RATES, ad->elem);
        ad->rates = mcmd_rates_read(path, &ad->num_rates);
        sfree(path);

        if (collmode == 1)
        {
            /* a rate-coefficient table, where present, replaces the fit file */
            path = mcmd_data_path(MCMD_PATH_CTAB_ION, ad->elem);
            ad->coll_tab = mcmd_ctab_read(path);
            sfree(path);
            if (ad->coll_tab == NULL)
            {
                path = mcmd_data_path(MCMD_PATH_COLL, ad->elem);
                ad->coll = mcmd_coll_read(path, &ad->num_coll);
                sfree(path);
            }
        }
        if (collmode >= 1)
        {
            /* optional since 2026-10-29: the weight of a shell configuration is C(2,K) C(8,L) C(18,M) in any model;
             * a file, where present, is used as given */
            path = mcmd_data_path(MCMD_PATH_WEIGHT, ad->elem);
            if (mcmd_count_lines(path) >= 0)
            {
                mcmd_weights_read(&ad->weights, path);
                ad->bWeightsFile = TRUE;
                mcmd_weights_compare(&ad->weights, path);
            }
            else
            {
                mcmd_weights_compute(&ad->weights, ad->rates, ad->num_rates);
            }
            sfree(path);

            /* optional: older Atomic_data directories have no excitation file */
            ad->colex     = NULL;
            ad->num_colex = 0;
            path = mcmd_data_path(MCMD_PATH_CTAB_EXC, ad->elem);
            ad->colex_tab = mcmd_ctab_read(path);
            sfree(path);
            path = mcmd_data_path(MCMD_PATH_COLEX, ad->elem);
            if (ad->colex_tab == NULL && mcmd_count_lines(path) > 0)
            {
                ad->colex = mcmd_coll_read(path, &ad->num_colex);
            }
            sfree(path);
        }
    }

    /* Index the states once, so the per-event lookups are array reads rather
     * than scans over every state of the element (324 of them for sulphur). */
    if (ad->num_rates > 0)
    {
        mcmd_state_index_build(&ad->rate_index, ad->rates[0].initial_state,
                               sizeof(ad->rates[0]), ad->num_rates);
    }
    if (collmode >= 1 && ad->num_coll > 0)
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
    if (ad->total_energies != NULL)
    {
        mcmd_dict_done(ad->total_energies);
    }
    mcmd_rates_done(ad->rates, ad->num_rates);

    if (bCollisional)
    {
        mcmd_coll_done(ad->coll, ad->num_coll);
        mcmd_coll_done(ad->colex, ad->num_colex);
        mcmd_ctab_done(ad->coll_tab);
        mcmd_ctab_done(ad->colex_tab);
        ad->coll_tab = ad->colex_tab = NULL;
    }
    mcmd_weights_done(&ad->weights);       /* the built-in data carry weights without collisions too */

    sfree(ad->rate_index.idx);
    sfree(ad->coll_index.idx);
    mcmd_state_index_clear(&ad->rate_index);
    mcmd_state_index_clear(&ad->coll_index);

    ad->bPresent       = FALSE;
    ad->energy_levels  = NULL;
    ad->total_energies = NULL;
    ad->rates          = NULL;
    ad->coll           = NULL;
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
