/*
 * mcionize.h
 *
 * The MolDStruct hybrid Monte Carlo / MD ionization module: X-ray driven
 * photoionization, Auger-Meitner decay and fluorescence, a free-electron
 * plasma, and classical over-the-barrier charge transfer between atoms.
 *
 * Enabled by running mdrun with -ionize.  The mdp userint/userreal knobs it
 * reads are documented in mcionize.c above mcionize_init().
 *
 * This file is part of the MolDStruct modifications to GROMACS 4.5.4 and is
 * distributed under the same licence as GROMACS.
 */

#ifndef _mcionize_h
#define _mcionize_h

#include <stdio.h>

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Opaque module state.  One instance per mdrun. */
typedef struct t_mcionize t_mcionize;

/* Set up the module: read the atomic data for the elements present, allocate
 * the per-atom electronic state, create ./simulation_output/ and open the
 * output files.  Call once before the MD loop.  Returns NULL if -ionize was
 * not given.
 *
 * With bStartFromCpt the electronic state, the free-electron plasma and the
 * reference radius of gyration are restored from
 * ./simulation_output/mcionize_state.bin, which mcionize_done() writes at the
 * end of every run.  Restarting without that file is an error: re-deriving
 * the reference radius from an already-expanded sample would silently
 * rescale the plasma volume. */
t_mcionize *mcionize_init(FILE *fplog, const t_inputrec *ir,
                          t_mdatoms *mdatoms, const t_commrec *cr,
                          gmx_bool bStartFromCpt);

/* Advance the electronic state by one MD step: charge transfer, then the
 * kinetic Monte Carlo over electronic transitions, then the per-step output.
 * Updates mdatoms->chargeA, which must happen before do_force() so that this
 * step's forces see this step's charges.
 *
 * Under particle decomposition this gathers coordinates onto every rank,
 * runs the module on the master, and broadcasts the resulting charges and
 * configurations.  Domain decomposition is not supported and is refused. */
void mcionize_step(t_mcionize *mc, FILE *fplog, const t_inputrec *ir,
                   t_mdatoms *mdatoms, t_state *state,
                   double t, const t_commrec *cr, t_nrnb *nrnb);

/* Autostop test (mcmd-autostop): TRUE once E_kin/E_tot exceeds
 * mcmd-autostop-threshold.  Kept out of the step function because do_md only
 * knows the energies later in the step.
 *
 * bEnergiesGlobal must be do_md's bGStat: enerd->term[] is only reduced
 * across ranks on global-communication steps, and on every other step it
 * holds this rank's partial sums, so the ratio is not the system's.  The
 * decision is taken on one rank and broadcast, because every rank has to
 * reach the same answer - if one leaves the MD loop and the others do not,
 * they hang in the next collective. */
gmx_bool mcionize_autostop(const t_mcionize *mc, const t_inputrec *ir,
                           double E_kin, double E_tot, gmx_large_int_t step,
                           gmx_bool bEnergiesGlobal, const t_commrec *cr);

/* Write the end-of-run output (process statistics, final charges and
 * configurations) and release everything.  Safe to call with mc == NULL. */
void mcionize_done(t_mcionize *mc, t_mdatoms *mdatoms, const t_commrec *cr);

#ifdef __cplusplus
}
#endif

#endif /* _mcionize_h */
