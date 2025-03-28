/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "physics.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "index.h"
#include "symtab.h"
#include "string2.h"
#include "readinp.h"
#include "warninp.h"
#include "readir.h" 
#include "toputil.h"
#include "index.h"
#include "network.h"
#include "vec.h"
#include "pbc.h"
#include "mtop_util.h"
#include "chargegroup.h"
#include "inputrec.h"

#define MAXPTR 254
#define NOGID  255

/* Resource parameters 
 * Do not change any of these until you read the instruction
 * in readinp.h. Some cpp's do not take spaces after the backslash
 * (like the c-shell), which will give you a very weird compiler
 * message.
 */

static char tcgrps[STRLEN],tau_t[STRLEN],ref_t[STRLEN],
  acc[STRLEN],accgrps[STRLEN],freeze[STRLEN],frdim[STRLEN],
  energy[STRLEN],user1[STRLEN],user2[STRLEN],vcm[STRLEN],xtc_grps[STRLEN],
  couple_moltype[STRLEN],orirefitgrp[STRLEN],egptable[STRLEN],egpexcl[STRLEN],
  wall_atomtype[STRLEN],wall_density[STRLEN],deform[STRLEN],QMMM[STRLEN];
static char foreign_lambda[STRLEN];
static char **pull_grp;
static char anneal[STRLEN],anneal_npoints[STRLEN],
  anneal_time[STRLEN],anneal_temp[STRLEN];
static char QMmethod[STRLEN],QMbasis[STRLEN],QMcharge[STRLEN],QMmult[STRLEN],
  bSH[STRLEN],CASorbitals[STRLEN], CASelectrons[STRLEN],SAon[STRLEN],
  SAoff[STRLEN],SAsteps[STRLEN],bTS[STRLEN],bOPT[STRLEN]; 
static char efield_x[STRLEN],efield_xt[STRLEN],efield_y[STRLEN],
  efield_yt[STRLEN],efield_z[STRLEN],efield_zt[STRLEN];

enum {
    egrptpALL,         /* All particles have to be a member of a group.     */
    egrptpALL_GENREST, /* A rest group with name is generated for particles *
                        * that are not part of any group.                   */
    egrptpPART,        /* As egrptpALL_GENREST, but no name is generated    *
                        * for the rest group.                               */
    egrptpONE          /* Merge all selected groups into one group,         *
                        * make a rest group for the remaining particles.    */
};


void init_ir(t_inputrec *ir, t_gromppopts *opts)
{
  snew(opts->include,STRLEN); 
  snew(opts->define,STRLEN);
}

static void _low_check(gmx_bool b,char *s,warninp_t wi)
{
    if (b)
    {
        warning_error(wi,s);
    }
}

static void check_nst(const char *desc_nst,int nst,
                      const char *desc_p,int *p,
                      warninp_t wi)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        sprintf(buf,"%s should be a multiple of %s, changing %s to %d\n",
		desc_p,desc_nst,desc_p,*p);
        warning(wi,buf);
    }
}

static gmx_bool ir_NVE(const t_inputrec *ir)
{
    return ((ir->eI == eiMD || EI_VV(ir->eI)) && ir->etc == etcNO);
}

static int lcd(int n1,int n2)
{
    int d,i;
    
    d = 1;
    for(i=2; (i<=n1 && i<=n2); i++)
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }
    
  return d;
}

void check_ir(const char *mdparin,t_inputrec *ir, t_gromppopts *opts,
              warninp_t wi)
/* Check internal consistency */
{
    /* Strange macro: first one fills the err_buf, and then one can check 
     * the condition, which will print the message and increase the error
     * counter.
     */
#define CHECK(b) _low_check(b,err_buf,wi)
    char err_buf[256],warn_buf[STRLEN];
    int  ns_type=0;
    real dt_pcoupl;

  set_warning_line(wi,mdparin,-1);

  /* BASIC CUT-OFF STUFF */
  if (ir->rlist == 0 ||
      !((EEL_MIGHT_BE_ZERO_AT_CUTOFF(ir->coulombtype) && ir->rcoulomb > ir->rlist) ||
        (EVDW_MIGHT_BE_ZERO_AT_CUTOFF(ir->vdwtype)    && ir->rvdw     > ir->rlist))) {
    /* No switched potential and/or no twin-range:
     * we can set the long-range cut-off to the maximum of the other cut-offs.
     */
    ir->rlistlong = max_cutoff(ir->rlist,max_cutoff(ir->rvdw,ir->rcoulomb));
  } else if (ir->rlistlong < 0) {
    ir->rlistlong = max_cutoff(ir->rlist,max_cutoff(ir->rvdw,ir->rcoulomb));
    sprintf(warn_buf,"rlistlong was not set, setting it to %g (no buffer)",
	    ir->rlistlong);
    warning(wi,warn_buf);
  }
  if (ir->rlistlong == 0 && ir->ePBC != epbcNONE) {
      warning_error(wi,"Can not have an infinite cut-off with PBC");
  }
  if (ir->rlistlong > 0 && (ir->rlist == 0 || ir->rlistlong < ir->rlist)) {
      warning_error(wi,"rlistlong can not be shorter than rlist");
  }
  if (IR_TWINRANGE(*ir) && ir->nstlist <= 0) {
      warning_error(wi,"Can not have nstlist<=0 with twin-range interactions");
  }

    /* GENERAL INTEGRATOR STUFF */
    if (!(ir->eI == eiMD || EI_VV(ir->eI)))
    {
        ir->etc = etcNO;
    }
    if (!EI_DYNAMICS(ir->eI))
    {
        ir->epc = epcNO;
    }
    if (EI_DYNAMICS(ir->eI))
    {
        if (ir->nstcalcenergy < 0)
        {
            ir->nstcalcenergy = ir_optimal_nstcalcenergy(ir);
            if (ir->nstenergy != 0 && ir->nstenergy < ir->nstcalcenergy)
            {
                /* nstcalcenergy larger than nstener does not make sense.
                 * We ideally want nstcalcenergy=nstener.
                 */
                if (ir->nstlist > 0)
                {
                    ir->nstcalcenergy = lcd(ir->nstenergy,ir->nstlist);
                }
                else
                {
                    ir->nstcalcenergy = ir->nstenergy;
                }
            }
        }
        if (ir->epc != epcNO)
        {
            if (ir->nstpcouple < 0)
            {
                ir->nstpcouple = ir_optimal_nstpcouple(ir);
            }
        }
        if (IR_TWINRANGE(*ir))
        {
            check_nst("nstlist",ir->nstlist,
                      "nstcalcenergy",&ir->nstcalcenergy,wi);
            if (ir->epc != epcNO)
            {
                check_nst("nstlist",ir->nstlist,
                          "nstpcouple",&ir->nstpcouple,wi); 
            }
        }

        if (ir->nstcalcenergy > 1)
        {
            /* for storing exact averages nstenergy should be
             * a multiple of nstcalcenergy
             */
            check_nst("nstcalcenergy",ir->nstcalcenergy,
                      "nstenergy",&ir->nstenergy,wi);
            if (ir->efep != efepNO)
            {
                /* nstdhdl should be a multiple of nstcalcenergy */
                check_nst("nstcalcenergy",ir->nstcalcenergy,
                          "nstdhdl",&ir->nstdhdl,wi);
            }
        }
    }

  /* LD STUFF */
  if ((EI_SD(ir->eI) || ir->eI == eiBD) &&
      ir->bContinuation && ir->ld_seed != -1) {
      warning_note(wi,"You are doing a continuation with SD or BD, make sure that ld_seed is different from the previous run (using ld_seed=-1 will ensure this)");
  }

  /* TPI STUFF */
  if (EI_TPI(ir->eI)) {
    sprintf(err_buf,"TPI only works with pbc = %s",epbc_names[epbcXYZ]);
    CHECK(ir->ePBC != epbcXYZ);
    sprintf(err_buf,"TPI only works with ns = %s",ens_names[ensGRID]);
    CHECK(ir->ns_type != ensGRID);
    sprintf(err_buf,"with TPI nstlist should be larger than zero");
    CHECK(ir->nstlist <= 0);
    sprintf(err_buf,"TPI does not work with full electrostatics other than PME");
    CHECK(EEL_FULL(ir->coulombtype) && !EEL_PME(ir->coulombtype));
  }

  /* SHAKE / LINCS */
  if ( (opts->nshake > 0) && (opts->bMorse) ) {
    sprintf(warn_buf,
	    "Using morse bond-potentials while constraining bonds is useless");
    warning(wi,warn_buf);
  }
  
  sprintf(err_buf,"shake_tol must be > 0 instead of %g while using shake",
	  ir->shake_tol);
  CHECK(((ir->shake_tol <= 0.0) && (opts->nshake>0) && 
	 (ir->eConstrAlg == econtSHAKE)));
     
  /* PBC/WALLS */
  sprintf(err_buf,"walls only work with pbc=%s",epbc_names[epbcXY]);
  CHECK(ir->nwall && ir->ePBC!=epbcXY);

  /* VACUUM STUFF */
  if (ir->ePBC != epbcXYZ && ir->nwall != 2) {
    if (ir->ePBC == epbcNONE) {
      if (ir->epc != epcNO) {
          warning(wi,"Turning off pressure coupling for vacuum system");
          ir->epc = epcNO;
      }
    } else {
      sprintf(err_buf,"Can not have pressure coupling with pbc=%s",
	      epbc_names[ir->ePBC]);
      CHECK(ir->epc != epcNO);
    }
    sprintf(err_buf,"Can not have Ewald with pbc=%s",epbc_names[ir->ePBC]);
    CHECK(EEL_FULL(ir->coulombtype));
    
    sprintf(err_buf,"Can not have dispersion correction with pbc=%s",
	    epbc_names[ir->ePBC]);
    CHECK(ir->eDispCorr != edispcNO);
  }

  if (ir->rlist == 0.0) {
    sprintf(err_buf,"can only have neighborlist cut-off zero (=infinite)\n"
	    "with coulombtype = %s or coulombtype = %s\n"
	    "without periodic boundary conditions (pbc = %s) and\n"
	    "rcoulomb and rvdw set to zero",
	    eel_names[eelCUT],eel_names[eelUSER],epbc_names[epbcNONE]);
    CHECK(((ir->coulombtype != eelCUT) && (ir->coulombtype != eelUSER)) ||
	  (ir->ePBC     != epbcNONE) || 
	  (ir->rcoulomb != 0.0)      || (ir->rvdw != 0.0));

    if (ir->nstlist < 0) {
        warning_error(wi,"Can not have heuristic neighborlist updates without cut-off");
    }
    if (ir->nstlist > 0) {
        warning_note(wi,"Simulating without cut-offs is usually (slightly) faster with nstlist=0, nstype=simple and particle decomposition");
    }
  }

  /* COMM STUFF */
  if (ir->nstcomm == 0) {
    ir->comm_mode = ecmNO;
  }
  if (ir->comm_mode != ecmNO) {
    if (ir->nstcomm < 0) {
        warning(wi,"If you want to remove the rotation around the center of mass, you should set comm_mode = Angular instead of setting nstcomm < 0. nstcomm is modified to its absolute value");
      ir->nstcomm = abs(ir->nstcomm);
    }
    
    if (ir->nstcalcenergy > 0 && ir->nstcomm < ir->nstcalcenergy) {
        warning_note(wi,"nstcomm < nstcalcenergy defeats the purpose of nstcalcenergy, setting nstcomm to nstcalcenergy");
      ir->nstcomm = ir->nstcalcenergy;
    }

    if (ir->comm_mode == ecmANGULAR) {
      sprintf(err_buf,"Can not remove the rotation around the center of mass with periodic molecules");
      CHECK(ir->bPeriodicMols);
      if (ir->ePBC != epbcNONE)
          warning(wi,"Removing the rotation around the center of mass in a periodic system (this is not a problem when you have only one molecule).");
    }
  }
    
  if (EI_STATE_VELOCITY(ir->eI) && ir->ePBC == epbcNONE && ir->comm_mode != ecmANGULAR) {
      warning_note(wi,"Tumbling and or flying ice-cubes: We are not removing rotation around center of mass in a non-periodic system. You should probably set comm_mode = ANGULAR.");
  }
  
  sprintf(err_buf,"Free-energy not implemented for Ewald and PPPM");
  CHECK((ir->coulombtype==eelEWALD || ir->coulombtype==eelPPPM)
	&& (ir->efep!=efepNO));
  
  sprintf(err_buf,"Twin-range neighbour searching (NS) with simple NS"
	  " algorithm not implemented");
  CHECK(((ir->rcoulomb > ir->rlist) || (ir->rvdw > ir->rlist)) 
	&& (ir->ns_type == ensSIMPLE));
  
    /* TEMPERATURE COUPLING */
    if (ir->etc == etcYES)
    {
        ir->etc = etcBERENDSEN;
        warning_note(wi,"Old option for temperature coupling given: "
                     "changing \"yes\" to \"Berendsen\"\n");
    }
  
    if (ir->etc == etcNOSEHOOVER)
    {
        if (ir->opts.nhchainlength < 1) 
        {
            sprintf(warn_buf,"number of Nose-Hoover chains (currently %d) cannot be less than 1,reset to 1\n",ir->opts.nhchainlength);
            ir->opts.nhchainlength =1;
            warning(wi,warn_buf);
        }
        
        if (ir->etc==etcNOSEHOOVER && !EI_VV(ir->eI) && ir->opts.nhchainlength > 1)
        {
            warning_note(wi,"leapfrog does not yet support Nose-Hoover chains, nhchainlength reset to 1");
            ir->opts.nhchainlength = 1;
        }
    }
    else
    {
        ir->opts.nhchainlength = 0;
    }

    if (ir->etc == etcBERENDSEN)
    {
        sprintf(warn_buf,"The %s thermostat does not generate the correct kinetic energy distribution. You might want to consider using the %s thermostat.",
                ETCOUPLTYPE(ir->etc),ETCOUPLTYPE(etcVRESCALE));
        warning_note(wi,warn_buf);
    }

    if ((ir->etc==etcNOSEHOOVER || ir->etc==etcANDERSEN || ir->etc==etcANDERSENINTERVAL) 
        && ir->epc==epcBERENDSEN)
    {
        sprintf(warn_buf,"Using Berendsen pressure coupling invalidates the "
                "true ensemble for the thermostat");
        warning(wi,warn_buf);
    }

    /* PRESSURE COUPLING */
    if (ir->epc == epcISOTROPIC)
    {
        ir->epc = epcBERENDSEN;
        warning_note(wi,"Old option for pressure coupling given: "
                     "changing \"Isotropic\" to \"Berendsen\"\n"); 
    }

    if (ir->epc != epcNO)
    {
        dt_pcoupl = ir->nstpcouple*ir->delta_t;

        sprintf(err_buf,"tau_p must be > 0 instead of %g\n",ir->tau_p);
        CHECK(ir->tau_p <= 0);
        
        if (ir->tau_p/dt_pcoupl < pcouple_min_integration_steps(ir->epc))
        {
            sprintf(warn_buf,"For proper integration of the %s barostat, tau_p (%g) should be at least %d times larger than nstpcouple*dt (%g)",
                    EPCOUPLTYPE(ir->epc),ir->tau_p,pcouple_min_integration_steps(ir->epc),dt_pcoupl);
            warning(wi,warn_buf);
        }	
        
        sprintf(err_buf,"compressibility must be > 0 when using pressure" 
                " coupling %s\n",EPCOUPLTYPE(ir->epc));
        CHECK(ir->compress[XX][XX] < 0 || ir->compress[YY][YY] < 0 || 
              ir->compress[ZZ][ZZ] < 0 || 
              (trace(ir->compress) == 0 && ir->compress[YY][XX] <= 0 &&
               ir->compress[ZZ][XX] <= 0 && ir->compress[ZZ][YY] <= 0));
        
        sprintf(err_buf,"pressure coupling with PPPM not implemented, use PME");
        CHECK(ir->coulombtype == eelPPPM);
        
    }
    else if (ir->coulombtype == eelPPPM)
    {
        sprintf(warn_buf,"The pressure with PPPM is incorrect, if you need the pressure use PME");
        warning(wi,warn_buf);
    }
    
    if (EI_VV(ir->eI))
    {
        if (ir->epc > epcNO)
        {
            if (ir->epc!=epcMTTK)
            {
                warning_error(wi,"NPT only defined for vv using Martyna-Tuckerman-Tobias-Klein equations");	      
            }
        }
    }

  /* ELECTROSTATICS */
  /* More checks are in triple check (grompp.c) */
    if (ir->coulombtype == eelPPPM)
    {
        warning_error(wi,"PPPM is not functional in the current version, we plan to implement PPPM through a small modification of the PME code");
    }

  if (ir->coulombtype == eelSWITCH) {
    sprintf(warn_buf,"coulombtype = %s is only for testing purposes and can lead to serious artifacts, advice: use coulombtype = %s",
	    eel_names[ir->coulombtype],
	    eel_names[eelRF_ZERO]);
    warning(wi,warn_buf);
  }

  if (ir->epsilon_r!=1 && ir->implicit_solvent==eisGBSA) {
    sprintf(warn_buf,"epsilon_r = %g with GB implicit solvent, will use this value for inner dielectric",ir->epsilon_r);
    warning_note(wi,warn_buf);
  }

  if (EEL_RF(ir->coulombtype) && ir->epsilon_rf==1 && ir->epsilon_r!=1) {
    sprintf(warn_buf,"epsilon_r = %g and epsilon_rf = 1 with reaction field, assuming old format and exchanging epsilon_r and epsilon_rf",ir->epsilon_r);
    warning(wi,warn_buf);
    ir->epsilon_rf = ir->epsilon_r;
    ir->epsilon_r  = 1.0;
  }

  if (getenv("GALACTIC_DYNAMICS") == NULL) {  
    sprintf(err_buf,"epsilon_r must be >= 0 instead of %g\n",ir->epsilon_r);
    CHECK(ir->epsilon_r < 0);
  }
  
  if (EEL_RF(ir->coulombtype)) {
    /* reaction field (at the cut-off) */
    
    if (ir->coulombtype == eelRF_ZERO) {
       sprintf(err_buf,"With coulombtype = %s, epsilon_rf must be 0",
	       eel_names[ir->coulombtype]);
      CHECK(ir->epsilon_rf != 0);
    }

    sprintf(err_buf,"epsilon_rf must be >= epsilon_r");
    CHECK((ir->epsilon_rf < ir->epsilon_r && ir->epsilon_rf != 0) ||
	  (ir->epsilon_r == 0));
    if (ir->epsilon_rf == ir->epsilon_r) {
      sprintf(warn_buf,"Using epsilon_rf = epsilon_r with %s does not make sense",
	      eel_names[ir->coulombtype]);
      warning(wi,warn_buf);
    }
  }
  /* Allow rlist>rcoulomb for tabulated long range stuff. This just
   * means the interaction is zero outside rcoulomb, but it helps to
   * provide accurate energy conservation.
   */
  if (EEL_MIGHT_BE_ZERO_AT_CUTOFF(ir->coulombtype)) {
    if (EEL_SWITCHED(ir->coulombtype)) {
      sprintf(err_buf,
	      "With coulombtype = %s rcoulomb_switch must be < rcoulomb",
	      eel_names[ir->coulombtype]);
      CHECK(ir->rcoulomb_switch >= ir->rcoulomb);
    }
  } else if (ir->coulombtype == eelCUT || EEL_RF(ir->coulombtype)) {
    sprintf(err_buf,"With coulombtype = %s, rcoulomb must be >= rlist",
	    eel_names[ir->coulombtype]);
    CHECK(ir->rlist > ir->rcoulomb);
  }

  if (EEL_FULL(ir->coulombtype)) {
    if (ir->coulombtype==eelPMESWITCH || ir->coulombtype==eelPMEUSER ||
        ir->coulombtype==eelPMEUSERSWITCH) {
      sprintf(err_buf,"With coulombtype = %s, rcoulomb must be <= rlist",
	      eel_names[ir->coulombtype]);
      CHECK(ir->rcoulomb > ir->rlist);
    } else {
      if (ir->coulombtype == eelPME) {
	sprintf(err_buf,
		"With coulombtype = %s, rcoulomb must be equal to rlist\n"
		"If you want optimal energy conservation or exact integration use %s",
		eel_names[ir->coulombtype],eel_names[eelPMESWITCH]);
      } else { 
	sprintf(err_buf,
		"With coulombtype = %s, rcoulomb must be equal to rlist",
		eel_names[ir->coulombtype]);
      }
      CHECK(ir->rcoulomb != ir->rlist);
    }
  }

  if (EEL_PME(ir->coulombtype)) {
    if (ir->pme_order < 3) {
        warning_error(wi,"pme_order can not be smaller than 3");
    }
  }

  if (ir->nwall==2 && EEL_FULL(ir->coulombtype)) {
    if (ir->ewald_geometry == eewg3D) {
      sprintf(warn_buf,"With pbc=%s you should use ewald_geometry=%s",
	      epbc_names[ir->ePBC],eewg_names[eewg3DC]);
      warning(wi,warn_buf);
    }
    /* This check avoids extra pbc coding for exclusion corrections */
    sprintf(err_buf,"wall_ewald_zfac should be >= 2");
    CHECK(ir->wall_ewald_zfac < 2);
  }

  if (EVDW_SWITCHED(ir->vdwtype)) {
    sprintf(err_buf,"With vdwtype = %s rvdw_switch must be < rvdw",
	    evdw_names[ir->vdwtype]);
    CHECK(ir->rvdw_switch >= ir->rvdw);
  } else if (ir->vdwtype == evdwCUT) {
    sprintf(err_buf,"With vdwtype = %s, rvdw must be >= rlist",evdw_names[ir->vdwtype]);
    CHECK(ir->rlist > ir->rvdw);
  }
  if (EEL_IS_ZERO_AT_CUTOFF(ir->coulombtype)
      && (ir->rlistlong <= ir->rcoulomb)) {
    sprintf(warn_buf,"For energy conservation with switch/shift potentials, %s should be 0.1 to 0.3 nm larger than rcoulomb.",
	    IR_TWINRANGE(*ir) ? "rlistlong" : "rlist");
    warning_note(wi,warn_buf);
  }
  if (EVDW_SWITCHED(ir->vdwtype) && (ir->rlistlong <= ir->rvdw)) {
    sprintf(warn_buf,"For energy conservation with switch/shift potentials, %s should be 0.1 to 0.3 nm larger than rvdw.",
	    IR_TWINRANGE(*ir) ? "rlistlong" : "rlist");
    warning_note(wi,warn_buf);
  }

  if (ir->vdwtype == evdwUSER && ir->eDispCorr != edispcNO) {
      warning_note(wi,"You have selected user tables with dispersion correction, the dispersion will be corrected to -C6/r^6 beyond rvdw_switch (the tabulated interaction between rvdw_switch and rvdw will not be double counted). Make sure that you really want dispersion correction to -C6/r^6.");
  }

  if (ir->nstlist == -1) {
    sprintf(err_buf,
	    "nstlist=-1 only works with switched or shifted potentials,\n"
	    "suggestion: use vdw-type=%s and coulomb-type=%s",
	    evdw_names[evdwSHIFT],eel_names[eelPMESWITCH]);
    CHECK(!(EEL_MIGHT_BE_ZERO_AT_CUTOFF(ir->coulombtype) &&
            EVDW_MIGHT_BE_ZERO_AT_CUTOFF(ir->vdwtype)));

    sprintf(err_buf,"With nstlist=-1 rvdw and rcoulomb should be smaller than rlist to account for diffusion and possibly charge-group radii");
    CHECK(ir->rvdw >= ir->rlist || ir->rcoulomb >= ir->rlist);
  }
  sprintf(err_buf,"nstlist can not be smaller than -1");
  CHECK(ir->nstlist < -1);

  if (ir->eI == eiLBFGS && (ir->coulombtype==eelCUT || ir->vdwtype==evdwCUT)
     && ir->rvdw != 0) {
    warning(wi,"For efficient BFGS minimization, use switch/shift/pme instead of cut-off.");
  }

  if (ir->eI == eiLBFGS && ir->nbfgscorr <= 0) {
    warning(wi,"Using L-BFGS with nbfgscorr<=0 just gets you steepest descent.");
  }

  /* FREE ENERGY */
  if (ir->efep != efepNO) {
    sprintf(err_buf,"The soft-core power is %d and can only be 1 or 2",
	    ir->sc_power);
    CHECK(ir->sc_alpha!=0 && ir->sc_power!=1 && ir->sc_power!=2);
  }

    /* ENERGY CONSERVATION */
    if (ir_NVE(ir))
    {
        if (!EVDW_MIGHT_BE_ZERO_AT_CUTOFF(ir->vdwtype) && ir->rvdw > 0)
        {
            sprintf(warn_buf,"You are using a cut-off for VdW interactions with NVE, for good energy conservation use vdwtype = %s (possibly with DispCorr)",
                    evdw_names[evdwSHIFT]);
            warning_note(wi,warn_buf);
        }
        if (!EEL_MIGHT_BE_ZERO_AT_CUTOFF(ir->coulombtype) && ir->rcoulomb > 0)
        {
            sprintf(warn_buf,"You are using a cut-off for electrostatics with NVE, for good energy conservation use coulombtype = %s or %s",
                    eel_names[eelPMESWITCH],eel_names[eelRF_ZERO]);
            warning_note(wi,warn_buf);
        }
    }
  
  /* IMPLICIT SOLVENT */
  if(ir->coulombtype==eelGB_NOTUSED)
  {
    ir->coulombtype=eelCUT;
    ir->implicit_solvent=eisGBSA;
    fprintf(stderr,"Note: Old option for generalized born electrostatics given:\n"
	    "Changing coulombtype from \"generalized-born\" to \"cut-off\" and instead\n"
	    "setting implicit_solvent value to \"GBSA\" in input section.\n");
  }

  if(ir->sa_algorithm==esaSTILL)
  {
    sprintf(err_buf,"Still SA algorithm not available yet, use %s or %s instead\n",esa_names[esaAPPROX],esa_names[esaNO]);
    CHECK(ir->sa_algorithm == esaSTILL);
  }
  
  if(ir->implicit_solvent==eisGBSA)
  {
    sprintf(err_buf,"With GBSA implicit solvent, rgbradii must be equal to rlist.");
    CHECK(ir->rgbradii != ir->rlist);
	  
    if(ir->coulombtype!=eelCUT)
	  {
		  sprintf(err_buf,"With GBSA, coulombtype must be equal to %s\n",eel_names[eelCUT]);
		  CHECK(ir->coulombtype!=eelCUT);
	  }
	  if(ir->vdwtype!=evdwCUT)
	  {
		  sprintf(err_buf,"With GBSA, vdw-type must be equal to %s\n",evdw_names[evdwCUT]);
		  CHECK(ir->vdwtype!=evdwCUT);
	  }
    if(ir->nstgbradii<1)
    {
      sprintf(warn_buf,"Using GBSA with nstgbradii<1, setting nstgbradii=1");
      warning_note(wi,warn_buf);
      ir->nstgbradii=1;
    }
    if(ir->sa_algorithm==esaNO)
    {
      sprintf(warn_buf,"No SA (non-polar) calculation requested together with GB. Are you sure this is what you want?\n");
      warning_note(wi,warn_buf);
    }
    if(ir->sa_surface_tension<0 && ir->sa_algorithm!=esaNO)
    {
      sprintf(warn_buf,"Value of sa_surface_tension is < 0. Changing it to 2.05016 or 2.25936 kJ/nm^2/mol for Still and HCT/OBC respectively\n");
      warning_note(wi,warn_buf);
      
      if(ir->gb_algorithm==egbSTILL)
      {
        ir->sa_surface_tension = 0.0049 * CAL2JOULE * 100;
      }
      else
      {
        ir->sa_surface_tension = 0.0054 * CAL2JOULE * 100;
      }
    }
    if(ir->sa_surface_tension==0 && ir->sa_algorithm!=esaNO)
    {
      sprintf(err_buf, "Surface tension set to 0 while SA-calculation requested\n");
      CHECK(ir->sa_surface_tension==0 && ir->sa_algorithm!=esaNO);
    }
    
  }
}

static int str_nelem(const char *str,int maxptr,char *ptr[])
{
  int  np=0;
  char *copy0,*copy;
  
  copy0=strdup(str); 
  copy=copy0;
  ltrim(copy);
  while (*copy != '\0') {
    if (np >= maxptr)
      gmx_fatal(FARGS,"Too many groups on line: '%s' (max is %d)",
		  str,maxptr);
    if (ptr) 
      ptr[np]=copy;
    np++;
    while ((*copy != '\0') && !isspace(*copy))
      copy++;
    if (*copy != '\0') {
      *copy='\0';
      copy++;
    }
    ltrim(copy);
  }
  if (ptr == NULL)
    sfree(copy0);

  return np;
}

static void parse_n_double(char *str,int *n,double **r)
{
  char *ptr[MAXPTR];
  int  i;

  *n = str_nelem(str,MAXPTR,ptr);

  snew(*r,*n);
  for(i=0; i<*n; i++) {
    (*r)[i] = strtod(ptr[i],NULL);
  }
}

static void do_wall_params(t_inputrec *ir,
                           char *wall_atomtype, char *wall_density,
                           t_gromppopts *opts)
{
    int  nstr,i;
    char *names[MAXPTR];
    double dbl;

    opts->wall_atomtype[0] = NULL;
    opts->wall_atomtype[1] = NULL;

    ir->wall_atomtype[0] = -1;
    ir->wall_atomtype[1] = -1;
    ir->wall_density[0] = 0;
    ir->wall_density[1] = 0;
  
    if (ir->nwall > 0)
    {
        nstr = str_nelem(wall_atomtype,MAXPTR,names);
        if (nstr != ir->nwall)
        {
            gmx_fatal(FARGS,"Expected %d elements for wall_atomtype, found %d",
                      ir->nwall,nstr);
        }
        for(i=0; i<ir->nwall; i++)
        {
            opts->wall_atomtype[i] = strdup(names[i]);
        }
    
        if (ir->wall_type == ewt93 || ir->wall_type == ewt104) {
            nstr = str_nelem(wall_density,MAXPTR,names);
            if (nstr != ir->nwall)
            {
                gmx_fatal(FARGS,"Expected %d elements for wall_density, found %d",ir->nwall,nstr);
            }
            for(i=0; i<ir->nwall; i++)
            {
                sscanf(names[i],"%lf",&dbl);
                if (dbl <= 0)
                {
                    gmx_fatal(FARGS,"wall_density[%d] = %f\n",i,dbl);
                }
                ir->wall_density[i] = dbl;
            }
        }
    }
}

static void add_wall_energrps(gmx_groups_t *groups,int nwall,t_symtab *symtab)
{
  int  i;
  t_grps *grps;
  char str[STRLEN];
  
  if (nwall > 0) {
    srenew(groups->grpname,groups->ngrpname+nwall);
    grps = &(groups->grps[egcENER]);
    srenew(grps->nm_ind,grps->nr+nwall);
    for(i=0; i<nwall; i++) {
      sprintf(str,"wall%d",i);
      groups->grpname[groups->ngrpname] = put_symtab(symtab,str);
      grps->nm_ind[grps->nr++] = groups->ngrpname++;
    }
  }
}

void get_ir(const char *mdparin,const char *mdparout,
            t_inputrec *ir,t_gromppopts *opts,
            warninp_t wi)
{
  char      *dumstr[2];
  double    dumdub[2][6];
  t_inpfile *inp;
  const char *tmp;
  int       i,j,m,ninp;
  char      warn_buf[STRLEN];
  
  inp = read_inpfile(mdparin, &ninp, NULL, wi);

  snew(dumstr[0],STRLEN);
  snew(dumstr[1],STRLEN);

  REM_TYPE("title");
  REM_TYPE("cpp");
  REM_TYPE("domain-decomposition");
  REPL_TYPE("unconstrained-start","continuation");
  REM_TYPE("dihre-tau");
  REM_TYPE("nstdihreout");
  REM_TYPE("nstcheckpoint");

  CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
  CTYPE ("Preprocessor information: use cpp syntax.");
  CTYPE ("e.g.: -I/home/joe/doe -I/home/mary/roe");
  STYPE ("include",	opts->include,	NULL);
  CTYPE ("e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)");
  STYPE ("define",	opts->define,	NULL);
    
  CCTYPE ("RUN CONTROL PARAMETERS");
  EETYPE("integrator",  ir->eI,         ei_names);
  CTYPE ("Start time and timestep in ps");
  RTYPE ("tinit",	ir->init_t,	0.0);
  RTYPE ("dt",		ir->delta_t,	0.001);
  STEPTYPE ("nsteps",   ir->nsteps,     0);
  CTYPE ("For exact run continuation or redoing part of a run");
  STEPTYPE ("init_step",ir->init_step,  0);
  CTYPE ("Part index is updated automatically on checkpointing (keeps files separate)");
  ITYPE ("simulation_part", ir->simulation_part, 1);
  CTYPE ("mode for center of mass motion removal");
  EETYPE("comm-mode",   ir->comm_mode,  ecm_names);
  CTYPE ("number of steps for center of mass motion removal");
  ITYPE ("nstcomm",	ir->nstcomm,	10);
  CTYPE ("group(s) for center of mass motion removal");
  STYPE ("comm-grps",   vcm,            NULL);
  
  CCTYPE ("LANGEVIN DYNAMICS OPTIONS");
  CTYPE ("Friction coefficient (amu/ps) and random seed");
  RTYPE ("bd-fric",     ir->bd_fric,    0.0);
  ITYPE ("ld-seed",     ir->ld_seed,    1993);
  
  /* Em stuff */
  CCTYPE ("ENERGY MINIMIZATION OPTIONS");
  CTYPE ("Force tolerance and initial step-size");
  RTYPE ("emtol",       ir->em_tol,     10.0);
  RTYPE ("emstep",      ir->em_stepsize,0.01);
  CTYPE ("Max number of iterations in relax_shells");
  ITYPE ("niter",       ir->niter,      20);
  CTYPE ("Step size (ps^2) for minimization of flexible constraints");
  RTYPE ("fcstep",      ir->fc_stepsize, 0);
  CTYPE ("Frequency of steepest descents steps when doing CG");
  ITYPE ("nstcgsteep",	ir->nstcgsteep,	1000);
  ITYPE ("nbfgscorr",   ir->nbfgscorr,  10); 

  CCTYPE ("TEST PARTICLE INSERTION OPTIONS");
  RTYPE ("rtpi",	ir->rtpi,	0.05);

  /* Output options */
  CCTYPE ("OUTPUT CONTROL OPTIONS");
  CTYPE ("Output frequency for coords (x), velocities (v) and forces (f)");
  ITYPE ("nstxout",	ir->nstxout,	100);
  ITYPE ("nstvout",	ir->nstvout,	100);
  ITYPE ("nstfout",	ir->nstfout,	0);
  ir->nstcheckpoint = 1000;
  CTYPE ("Output frequency for energies to log file and energy file");
  ITYPE ("nstlog",	ir->nstlog,	100);
  ITYPE ("nstcalcenergy",ir->nstcalcenergy,	-1);
  ITYPE ("nstenergy",   ir->nstenergy,  100);
  CTYPE ("Output frequency and precision for .xtc file");
  ITYPE ("nstxtcout",   ir->nstxtcout,  0);
  RTYPE ("xtc-precision",ir->xtcprec,   1000.0);
  CTYPE ("This selects the subset of atoms for the .xtc file. You can");
  CTYPE ("select multiple groups. By default all atoms will be written.");
  STYPE ("xtc-grps",    xtc_grps,       NULL);
  CTYPE ("Selection of energy groups");
  STYPE ("energygrps",  energy,         NULL);

  /* Neighbor searching */  
  CCTYPE ("NEIGHBORSEARCHING PARAMETERS");
  CTYPE ("nblist update frequency");
  ITYPE ("nstlist",	ir->nstlist,	10);
  CTYPE ("ns algorithm (simple or grid)");
  EETYPE("ns-type",     ir->ns_type,    ens_names);
  /* set ndelta to the optimal value of 2 */
  ir->ndelta = 2;
  CTYPE ("Periodic boundary conditions: xyz, no, xy");
  EETYPE("pbc",         ir->ePBC,       epbc_names);
  EETYPE("periodic_molecules", ir->bPeriodicMols, yesno_names);
  CTYPE ("nblist cut-off");
  RTYPE ("rlist",	ir->rlist,	1.0);
  CTYPE ("long-range cut-off for switched potentials");
  RTYPE ("rlistlong",	ir->rlistlong,	-1);

  /* Electrostatics */
  CCTYPE ("OPTIONS FOR ELECTROSTATICS AND VDW");
  CTYPE ("Method for doing electrostatics");
  EETYPE("coulombtype",	ir->coulombtype,    eel_names);
  CTYPE ("cut-off lengths");
  RTYPE ("rcoulomb-switch",	ir->rcoulomb_switch,	0.0);
  RTYPE ("rcoulomb",	ir->rcoulomb,	1.0);
  CTYPE ("Relative dielectric constant for the medium and the reaction field");
  RTYPE ("epsilon_r",   ir->epsilon_r,  1.0);
  RTYPE ("epsilon_rf",  ir->epsilon_rf, 1.0);
  CTYPE ("Method for doing Van der Waals");
  EETYPE("vdw-type",	ir->vdwtype,    evdw_names);
  CTYPE ("cut-off lengths");
  RTYPE ("rvdw-switch",	ir->rvdw_switch,	0.0);
  RTYPE ("rvdw",	ir->rvdw,	1.0);
  CTYPE ("Apply long range dispersion corrections for Energy and Pressure");
  EETYPE("DispCorr",    ir->eDispCorr,  edispc_names);
  CTYPE ("Extension of the potential lookup tables beyond the cut-off");
  RTYPE ("table-extension", ir->tabext, 1.0);
  CTYPE ("Seperate tables between energy group pairs");
  STYPE ("energygrp_table", egptable,   NULL);
  CTYPE ("Spacing for the PME/PPPM FFT grid");
  RTYPE ("fourierspacing", opts->fourierspacing,0.12);
  CTYPE ("FFT grid size, when a value is 0 fourierspacing will be used");
  ITYPE ("fourier_nx",  ir->nkx,         0);
  ITYPE ("fourier_ny",  ir->nky,         0);
  ITYPE ("fourier_nz",  ir->nkz,         0);
  CTYPE ("EWALD/PME/PPPM parameters");
  ITYPE ("pme_order",   ir->pme_order,   4);
  RTYPE ("ewald_rtol",  ir->ewald_rtol, 0.00001);
  EETYPE("ewald_geometry", ir->ewald_geometry, eewg_names);
  RTYPE ("epsilon_surface", ir->epsilon_surface, 0.0);
  EETYPE("optimize_fft",ir->bOptFFT,  yesno_names);

  CCTYPE("IMPLICIT SOLVENT ALGORITHM");
  EETYPE("implicit_solvent", ir->implicit_solvent, eis_names);
	
  CCTYPE ("GENERALIZED BORN ELECTROSTATICS"); 
  CTYPE ("Algorithm for calculating Born radii");
  EETYPE("gb_algorithm", ir->gb_algorithm, egb_names);
  CTYPE ("Frequency of calculating the Born radii inside rlist");
  ITYPE ("nstgbradii", ir->nstgbradii, 1);
  CTYPE ("Cutoff for Born radii calculation; the contribution from atoms");
  CTYPE ("between rlist and rgbradii is updated every nstlist steps");
  RTYPE ("rgbradii",  ir->rgbradii, 1.0);
  CTYPE ("Dielectric coefficient of the implicit solvent");
  RTYPE ("gb_epsilon_solvent",ir->gb_epsilon_solvent, 80.0);	
  CTYPE ("Salt concentration in M for Generalized Born models");
  RTYPE ("gb_saltconc",  ir->gb_saltconc, 0.0); 
  CTYPE ("Scaling factors used in the OBC GB model. Default values are OBC(II)");
  RTYPE ("gb_obc_alpha", ir->gb_obc_alpha, 1.0);
  RTYPE ("gb_obc_beta", ir->gb_obc_beta, 0.8);
  RTYPE ("gb_obc_gamma", ir->gb_obc_gamma, 4.85);	
  RTYPE ("gb_dielectric_offset", ir->gb_dielectric_offset, 0.009);
  EETYPE("sa_algorithm", ir->sa_algorithm, esa_names);
  CTYPE ("Surface tension (kJ/mol/nm^2) for the SA (nonpolar surface) part of GBSA");
  CTYPE ("The value -1 will set default value for Still/HCT/OBC GB-models.");
  RTYPE ("sa_surface_tension", ir->sa_surface_tension, -1);
		 
  /* Coupling stuff */
  CCTYPE ("OPTIONS FOR WEAK COUPLING ALGORITHMS");
  CTYPE ("Temperature coupling");
  EETYPE("tcoupl",	ir->etc,        etcoupl_names);
  ITYPE ("nsttcouple", ir->nsttcouple,  -1);
  ITYPE("nh-chain-length",     ir->opts.nhchainlength, NHCHAINLENGTH);
  CTYPE ("Groups to couple separately");
  STYPE ("tc-grps",     tcgrps,         NULL);
  CTYPE ("Time constant (ps) and reference temperature (K)");
  STYPE ("tau-t",	tau_t,		NULL);
  STYPE ("ref-t",	ref_t,		NULL);
  CTYPE ("Pressure coupling");
  EETYPE("Pcoupl",	ir->epc,        epcoupl_names);
  EETYPE("Pcoupltype",	ir->epct,       epcoupltype_names);
  ITYPE ("nstpcouple", ir->nstpcouple,  -1);
  CTYPE ("Time constant (ps), compressibility (1/bar) and reference P (bar)");
  RTYPE ("tau-p",	ir->tau_p,	1.0);
  STYPE ("compressibility",	dumstr[0],	NULL);
  STYPE ("ref-p",       dumstr[1],      NULL);
  CTYPE ("Scaling of reference coordinates, No, All or COM");
  EETYPE ("refcoord_scaling",ir->refcoord_scaling,erefscaling_names);

  CTYPE ("Random seed for Andersen thermostat");
  ITYPE ("andersen_seed", ir->andersen_seed, 815131);

  /* QMMM */
  CCTYPE ("OPTIONS FOR QMMM calculations");
  EETYPE("QMMM", ir->bQMMM, yesno_names);
  CTYPE ("Groups treated Quantum Mechanically");
  STYPE ("QMMM-grps",  QMMM,          NULL);
  CTYPE ("QM method");
  STYPE("QMmethod",     QMmethod, NULL);
  CTYPE ("QMMM scheme");
  EETYPE("QMMMscheme",  ir->QMMMscheme,    eQMMMscheme_names);
  CTYPE ("QM basisset");
  STYPE("QMbasis",      QMbasis, NULL);
  CTYPE ("QM charge");
  STYPE ("QMcharge",    QMcharge,NULL);
  CTYPE ("QM multiplicity");
  STYPE ("QMmult",      QMmult,NULL);
  CTYPE ("Surface Hopping");
  STYPE ("SH",          bSH, NULL);
  CTYPE ("CAS space options");
  STYPE ("CASorbitals",      CASorbitals,   NULL);
  STYPE ("CASelectrons",     CASelectrons,  NULL);
  STYPE ("SAon", SAon, NULL);
  STYPE ("SAoff",SAoff,NULL);
  STYPE ("SAsteps",  SAsteps, NULL);
  CTYPE ("Scale factor for MM charges");
  RTYPE ("MMChargeScaleFactor", ir->scalefactor, 1.0);
  CTYPE ("Optimization of QM subsystem");
  STYPE ("bOPT",          bOPT, NULL);
  STYPE ("bTS",          bTS, NULL);

  /* Simulated annealing */
  CCTYPE("SIMULATED ANNEALING");
  CTYPE ("Type of annealing for each temperature group (no/single/periodic)");
  STYPE ("annealing",   anneal,      NULL);
  CTYPE ("Number of time points to use for specifying annealing in each group");
  STYPE ("annealing_npoints", anneal_npoints, NULL);
  CTYPE ("List of times at the annealing points for each group");
  STYPE ("annealing_time",       anneal_time,       NULL);
  CTYPE ("Temp. at each annealing point, for each group.");
  STYPE ("annealing_temp",  anneal_temp,  NULL);
  
  /* Startup run */
  CCTYPE ("GENERATE VELOCITIES FOR STARTUP RUN");
  EETYPE("gen-vel",     opts->bGenVel,  yesno_names);
  RTYPE ("gen-temp",    opts->tempi,    300.0);
  ITYPE ("gen-seed",    opts->seed,     173529);
  
  /* Shake stuff */
  CCTYPE ("OPTIONS FOR BONDS");
  EETYPE("constraints",	opts->nshake,	constraints);
  CTYPE ("Type of constraint algorithm");
  EETYPE("constraint-algorithm",  ir->eConstrAlg, econstr_names);
  CTYPE ("Do not constrain the start configuration");
  EETYPE("continuation", ir->bContinuation, yesno_names);
  CTYPE ("Use successive overrelaxation to reduce the number of shake iterations");
  EETYPE("Shake-SOR", ir->bShakeSOR, yesno_names);
  CTYPE ("Relative tolerance of shake");
  RTYPE ("shake-tol", ir->shake_tol, 0.0001);
  CTYPE ("Highest order in the expansion of the constraint coupling matrix");
  ITYPE ("lincs-order", ir->nProjOrder, 4);
  CTYPE ("Number of iterations in the final step of LINCS. 1 is fine for");
  CTYPE ("normal simulations, but use 2 to conserve energy in NVE runs.");
  CTYPE ("For energy minimization with constraints it should be 4 to 8.");
  ITYPE ("lincs-iter", ir->nLincsIter, 1);
  CTYPE ("Lincs will write a warning to the stderr if in one step a bond"); 
  CTYPE ("rotates over more degrees than");
  RTYPE ("lincs-warnangle", ir->LincsWarnAngle, 30.0);
  CTYPE ("Convert harmonic bonds to morse potentials");
  EETYPE("morse",       opts->bMorse,yesno_names);

  /* Energy group exclusions */
  CCTYPE ("ENERGY GROUP EXCLUSIONS");
  CTYPE ("Pairs of energy groups for which all non-bonded interactions are excluded");
  STYPE ("energygrp_excl", egpexcl,     NULL);
  
  /* Walls */
  CCTYPE ("WALLS");
  CTYPE ("Number of walls, type, atom types, densities and box-z scale factor for Ewald");
  ITYPE ("nwall", ir->nwall, 0);
  EETYPE("wall_type",     ir->wall_type,   ewt_names);
  RTYPE ("wall_r_linpot", ir->wall_r_linpot, -1);
  STYPE ("wall_atomtype", wall_atomtype, NULL);
  STYPE ("wall_density",  wall_density,  NULL);
  RTYPE ("wall_ewald_zfac", ir->wall_ewald_zfac, 3);
  
  /* COM pulling */
  CCTYPE("COM PULLING");
  CTYPE("Pull type: no, umbrella, constraint or constant_force");
  EETYPE("pull",          ir->ePull, epull_names);
  if (ir->ePull != epullNO) {
    snew(ir->pull,1);
    pull_grp = read_pullparams(&ninp,&inp,ir->pull,&opts->pull_start,wi);
  }

  /* Refinement */
  CCTYPE("NMR refinement stuff");
  CTYPE ("Distance restraints type: No, Simple or Ensemble");
  EETYPE("disre",       ir->eDisre,     edisre_names);
  CTYPE ("Force weighting of pairs in one distance restraint: Conservative or Equal");
  EETYPE("disre-weighting", ir->eDisreWeighting, edisreweighting_names);
  CTYPE ("Use sqrt of the time averaged times the instantaneous violation");
  EETYPE("disre-mixed", ir->bDisreMixed, yesno_names);
  RTYPE ("disre-fc",	ir->dr_fc,	1000.0);
  RTYPE ("disre-tau",	ir->dr_tau,	0.0);
  CTYPE ("Output frequency for pair distances to energy file");
  ITYPE ("nstdisreout", ir->nstdisreout, 100);
  CTYPE ("Orientation restraints: No or Yes");
  EETYPE("orire",       opts->bOrire,   yesno_names);
  CTYPE ("Orientation restraints force constant and tau for time averaging");
  RTYPE ("orire-fc",	ir->orires_fc,	0.0);
  RTYPE ("orire-tau",	ir->orires_tau,	0.0);
  STYPE ("orire-fitgrp",orirefitgrp,    NULL);
  CTYPE ("Output frequency for trace(SD) and S to energy file");
  ITYPE ("nstorireout", ir->nstorireout, 100);
  CTYPE ("Dihedral angle restraints: No or Yes");
  EETYPE("dihre",       opts->bDihre,   yesno_names);
  RTYPE ("dihre-fc",	ir->dihre_fc,	1000.0);

  /* Free energy stuff */
  CCTYPE ("Free energy control stuff");
  EETYPE("free-energy",	ir->efep, efep_names);
  RTYPE ("init-lambda",	ir->init_lambda,0.0);
  RTYPE ("delta-lambda",ir->delta_lambda,0.0);
  STYPE ("foreign_lambda", foreign_lambda, NULL);
  RTYPE ("sc-alpha",ir->sc_alpha,0.0);
  ITYPE ("sc-power",ir->sc_power,0);
  RTYPE ("sc-sigma",ir->sc_sigma,0.3);
  ITYPE ("nstdhdl",     ir->nstdhdl, 10);
  EETYPE("separate-dhdl-file", ir->separate_dhdl_file, 
                               separate_dhdl_file_names);
  EETYPE("dhdl-derivatives", ir->dhdl_derivatives, dhdl_derivatives_names);
  ITYPE ("dh_hist_size", ir->dh_hist_size, 0);
  RTYPE ("dh_hist_spacing", ir->dh_hist_spacing, 0.1);
  STYPE ("couple-moltype",  couple_moltype,  NULL);
  EETYPE("couple-lambda0", opts->couple_lam0, couple_lam);
  EETYPE("couple-lambda1", opts->couple_lam1, couple_lam);
  EETYPE("couple-intramol", opts->bCoupleIntra, yesno_names);

  /* Non-equilibrium MD stuff */  
  CCTYPE("Non-equilibrium MD stuff");
  STYPE ("acc-grps",    accgrps,        NULL);
  STYPE ("accelerate",  acc,            NULL);
  STYPE ("freezegrps",  freeze,         NULL);
  STYPE ("freezedim",   frdim,          NULL);
  RTYPE ("cos-acceleration", ir->cos_accel, 0);
  STYPE ("deform",      deform,         NULL);

  /* Electric fields */
  CCTYPE("Electric fields");
  CTYPE ("Format is number of terms (int) and for all terms an amplitude (real)");
  CTYPE ("and a phase angle (real)");
  STYPE ("E-x",   	efield_x,	NULL);
  STYPE ("E-xt",	efield_xt,	NULL);
  STYPE ("E-y",   	efield_y,	NULL);
  STYPE ("E-yt",	efield_yt,	NULL);
  STYPE ("E-z",   	efield_z,	NULL);
  STYPE ("E-zt",	efield_zt,	NULL);
  
  /* User defined thingies */
  CCTYPE ("User defined thingies");
  STYPE ("user1-grps",  user1,          NULL);
  STYPE ("user2-grps",  user2,          NULL);

  ITYPE ("userint1",    ir->userint1,   1); // Enable altererd forcefields 
  ITYPE ("userint2",    ir->userint2,   1); // Do charge transfer
  ITYPE ("userint3",    ir->userint3,   0); // Enable stopping when threshgold is reached (userreal6)
  ITYPE ("userint4",    ir->userint4,   0); // Read states from previous sim
  ITYPE ("userint5",    ir->userint5,   0); // Enable logging
  ITYPE ("userint6",    ir->userint6,   0); // Enable collsional ionization (Not implemented)
  ITYPE ("userint7",    ir->userint7,   0); // FREE
  ITYPE ("userint8",    ir->userint8,   0); // FREE
  ITYPE ("userint9",    ir->userint9,   0); // Read charges from file


  RTYPE ("userreal1",   ir->userreal1,  0); // Gaussian peak in ps
  RTYPE ("userreal2",   ir->userreal2,  0); // Number of photons
  RTYPE ("userreal3",   ir->userreal3,  0); // Sigma value of gaussian
  RTYPE ("userreal4",   ir->userreal4,  0); // Focal diamater
  RTYPE ("userreal5",   ir->userreal5,  0); // Photon energy
  RTYPE ("userreal6",   ir->userreal6,  0.99); // Threshold for stopping sim. userint3 must be turned on.
  RTYPE ("userreal7",   ir->userreal7,  0);
  RTYPE ("userreal8",   ir->userreal8,  0);
  RTYPE ("userreal9",   ir->userreal9,  0);


#undef CTYPE

  write_inpfile(mdparout,ninp,inp,FALSE,wi);
  for (i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  /* Process options if necessary */
  for(m=0; m<2; m++) {
    for(i=0; i<2*DIM; i++)
      dumdub[m][i]=0.0;
    if(ir->epc) {
      switch (ir->epct) {
      case epctISOTROPIC:
	if (sscanf(dumstr[m],"%lf",&(dumdub[m][XX]))!=1) {
        warning_error(wi,"Pressure coupling not enough values (I need 1)");
	}
	dumdub[m][YY]=dumdub[m][ZZ]=dumdub[m][XX];
	break;
      case epctSEMIISOTROPIC:
      case epctSURFACETENSION:
	if (sscanf(dumstr[m],"%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][ZZ]))!=2) {
        warning_error(wi,"Pressure coupling not enough values (I need 2)");
	}
	dumdub[m][YY]=dumdub[m][XX];
	break;
      case epctANISOTROPIC:
	if (sscanf(dumstr[m],"%lf%lf%lf%lf%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][YY]),&(dumdub[m][ZZ]),
		   &(dumdub[m][3]),&(dumdub[m][4]),&(dumdub[m][5]))!=6) {
        warning_error(wi,"Pressure coupling not enough values (I need 6)");
	}
	break;
      default:
	gmx_fatal(FARGS,"Pressure coupling type %s not implemented yet",
		    epcoupltype_names[ir->epct]);
      }
    }
  }
  clear_mat(ir->ref_p);
  clear_mat(ir->compress);
  for(i=0; i<DIM; i++) {
    ir->ref_p[i][i]    = dumdub[1][i];
    ir->compress[i][i] = dumdub[0][i];
  }
  if (ir->epct == epctANISOTROPIC) {
    ir->ref_p[XX][YY] = dumdub[1][3];
    ir->ref_p[XX][ZZ] = dumdub[1][4];
    ir->ref_p[YY][ZZ] = dumdub[1][5];
    if (ir->ref_p[XX][YY]!=0 && ir->ref_p[XX][ZZ]!=0 && ir->ref_p[YY][ZZ]!=0) {
      warning(wi,"All off-diagonal reference pressures are non-zero. Are you sure you want to apply a threefold shear stress?\n");
    }
    ir->compress[XX][YY] = dumdub[0][3];
    ir->compress[XX][ZZ] = dumdub[0][4];
    ir->compress[YY][ZZ] = dumdub[0][5];
    for(i=0; i<DIM; i++) {
      for(m=0; m<i; m++) {
	ir->ref_p[i][m] = ir->ref_p[m][i];
	ir->compress[i][m] = ir->compress[m][i];
      }
    }
  } 
  
  if (ir->comm_mode == ecmNO)
    ir->nstcomm = 0;

  opts->couple_moltype = NULL;
  if (strlen(couple_moltype) > 0) {
    if (ir->efep != efepNO) {
      opts->couple_moltype = strdup(couple_moltype);
      if (opts->couple_lam0 == opts->couple_lam1)
	warning(wi,"The lambda=0 and lambda=1 states for coupling are identical");
      if (ir->eI == eiMD && (opts->couple_lam0 == ecouplamNONE ||
			     opts->couple_lam1 == ecouplamNONE)) {
	warning(wi,"For proper sampling of the (nearly) decoupled state, stochastic dynamics should be used");
      }
    } else {
      warning(wi,"Can not couple a molecule with free_energy = no");
    }
  }

  do_wall_params(ir,wall_atomtype,wall_density,opts);
  
  if (opts->bOrire && str_nelem(orirefitgrp,MAXPTR,NULL)!=1) {
      warning_error(wi,"ERROR: Need one orientation restraint fit group\n");
  }

  clear_mat(ir->deform);
  for(i=0; i<6; i++)
    dumdub[0][i] = 0;
  m = sscanf(deform,"%lf %lf %lf %lf %lf %lf",
	     &(dumdub[0][0]),&(dumdub[0][1]),&(dumdub[0][2]),
	     &(dumdub[0][3]),&(dumdub[0][4]),&(dumdub[0][5]));
  for(i=0; i<3; i++)
    ir->deform[i][i] = dumdub[0][i];
  ir->deform[YY][XX] = dumdub[0][3];
  ir->deform[ZZ][XX] = dumdub[0][4];
  ir->deform[ZZ][YY] = dumdub[0][5];
  if (ir->epc != epcNO) {
    for(i=0; i<3; i++)
      for(j=0; j<=i; j++)
	if (ir->deform[i][j]!=0 && ir->compress[i][j]!=0) {
        warning_error(wi,"A box element has deform set and compressibility > 0");
	}
    for(i=0; i<3; i++)
      for(j=0; j<i; j++)
	if (ir->deform[i][j]!=0) {
	  for(m=j; m<DIM; m++)
	    if (ir->compress[m][j]!=0) {
	      sprintf(warn_buf,"An off-diagonal box element has deform set while compressibility > 0 for the same component of another box vector, this might lead to spurious periodicity effects.");
	      warning(wi,warn_buf);
	    }
	}
  }

  if (ir->efep != efepNO) {
    parse_n_double(foreign_lambda,&ir->n_flambda,&ir->flambda);
    if (ir->n_flambda > 0 && ir->rlist < max(ir->rvdw,ir->rcoulomb)) {
      warning_note(wi,"For foreign lambda free energy differences it is assumed that the soft-core interactions have no effect beyond the neighborlist cut-off");
    }
  } else {
    ir->n_flambda = 0;
  }

  sfree(dumstr[0]);
  sfree(dumstr[1]);
}

static int search_QMstring(char *s,int ng,const char *gn[])
{
  /* same as normal search_string, but this one searches QM strings */
  int i;

  for(i=0; (i<ng); i++)
    if (gmx_strcasecmp(s,gn[i]) == 0)
      return i;

  gmx_fatal(FARGS,"this QM method or basisset (%s) is not implemented\n!",s);

  return -1;

} /* search_QMstring */


int search_string(char *s,int ng,char *gn[])
{
  int i;
  
  for(i=0; (i<ng); i++)
    if (gmx_strcasecmp(s,gn[i]) == 0)
      return i;
      
  gmx_fatal(FARGS,"Group %s not found in indexfile.\nMaybe you have non-default goups in your .mdp file, while not using the '-n' option of grompp.\nIn that case use the '-n' option.\n",s);
  
  return -1;
}

static gmx_bool do_numbering(int natoms,gmx_groups_t *groups,int ng,char *ptrs[],
                         t_blocka *block,char *gnames[],
                         int gtype,int restnm,
                         int grptp,gmx_bool bVerbose,
                         warninp_t wi)
{
    unsigned short *cbuf;
    t_grps *grps=&(groups->grps[gtype]);
    int    i,j,gid,aj,ognr,ntot=0;
    const char *title;
    gmx_bool   bRest;
    char   warn_buf[STRLEN];

    if (debug)
    {
        fprintf(debug,"Starting numbering %d groups of type %d\n",ng,gtype);
    }
  
    title = gtypes[gtype];
    
    snew(cbuf,natoms);
    /* Mark all id's as not set */
    for(i=0; (i<natoms); i++)
    {
        cbuf[i] = NOGID;
    }
  
    snew(grps->nm_ind,ng+1); /* +1 for possible rest group */
    for(i=0; (i<ng); i++)
    {
        /* Lookup the group name in the block structure */
        gid = search_string(ptrs[i],block->nr,gnames);
        if ((grptp != egrptpONE) || (i == 0))
        {
            grps->nm_ind[grps->nr++]=gid;
        }
        if (debug) 
        {
            fprintf(debug,"Found gid %d for group %s\n",gid,ptrs[i]);
        }
    
        /* Now go over the atoms in the group */
        for(j=block->index[gid]; (j<block->index[gid+1]); j++)
        {

            aj=block->a[j];
      
            /* Range checking */
            if ((aj < 0) || (aj >= natoms)) 
            {
                gmx_fatal(FARGS,"Invalid atom number %d in indexfile",aj);
            }
            /* Lookup up the old group number */
            ognr = cbuf[aj];
            if (ognr != NOGID)
            {
                gmx_fatal(FARGS,"Atom %d in multiple %s groups (%d and %d)",
                          aj+1,title,ognr+1,i+1);
            }
            else
            {
                /* Store the group number in buffer */
                if (grptp == egrptpONE)
                {
                    cbuf[aj] = 0;
                }
                else
                {
                    cbuf[aj] = i;
                }
                ntot++;
            }
        }
    }
    
    /* Now check whether we have done all atoms */
    bRest = FALSE;
    if (ntot != natoms)
    {
        if (grptp == egrptpALL)
        {
            gmx_fatal(FARGS,"%d atoms are not part of any of the %s groups",
                      natoms-ntot,title);
        }
        else if (grptp == egrptpPART)
        {
            sprintf(warn_buf,"%d atoms are not part of any of the %s groups",
                    natoms-ntot,title);
            warning_note(wi,warn_buf);
        }
        /* Assign all atoms currently unassigned to a rest group */
        for(j=0; (j<natoms); j++)
        {
            if (cbuf[j] == NOGID)
            {
                cbuf[j] = grps->nr;
                bRest = TRUE;
            }
        }
        if (grptp != egrptpPART)
        {
            if (bVerbose)
            {
                fprintf(stderr,
                        "Making dummy/rest group for %s containing %d elements\n",
                        title,natoms-ntot);
            }
            /* Add group name "rest" */ 
            grps->nm_ind[grps->nr] = restnm;
            
            /* Assign the rest name to all atoms not currently assigned to a group */
            for(j=0; (j<natoms); j++)
            {
                if (cbuf[j] == NOGID)
                {
                    cbuf[j] = grps->nr;
                }
            }
            grps->nr++;
        }
    }
    
    if (grps->nr == 1)
    {
        groups->ngrpnr[gtype] = 0;
        groups->grpnr[gtype]  = NULL;
    }
    else
    {
        groups->ngrpnr[gtype] = natoms;
        snew(groups->grpnr[gtype],natoms);
        for(j=0; (j<natoms); j++)
        {
            groups->grpnr[gtype][j] = cbuf[j];
        }
    }
    
    sfree(cbuf);

    return (bRest && grptp == egrptpPART);
}

static void calc_nrdf(gmx_mtop_t *mtop,t_inputrec *ir,char **gnames)
{
  t_grpopts *opts;
  gmx_groups_t *groups;
  t_pull  *pull;
  int     natoms,ai,aj,i,j,d,g,imin,jmin,nc;
  t_iatom *ia;
  int     *nrdf2,*na_vcm,na_tot;
  double  *nrdf_tc,*nrdf_vcm,nrdf_uc,n_sub=0;
  gmx_mtop_atomloop_all_t aloop;
  t_atom  *atom;
  int     mb,mol,ftype,as;
  gmx_molblock_t *molb;
  gmx_moltype_t *molt;

  /* Calculate nrdf. 
   * First calc 3xnr-atoms for each group
   * then subtract half a degree of freedom for each constraint
   *
   * Only atoms and nuclei contribute to the degrees of freedom...
   */

  opts = &ir->opts;
  
  groups = &mtop->groups;
  natoms = mtop->natoms;

  /* Allocate one more for a possible rest group */
  /* We need to sum degrees of freedom into doubles,
   * since floats give too low nrdf's above 3 million atoms.
   */
  snew(nrdf_tc,groups->grps[egcTC].nr+1);
  snew(nrdf_vcm,groups->grps[egcVCM].nr+1);
  snew(na_vcm,groups->grps[egcVCM].nr+1);
  
  for(i=0; i<groups->grps[egcTC].nr; i++)
    nrdf_tc[i] = 0;
  for(i=0; i<groups->grps[egcVCM].nr+1; i++)
    nrdf_vcm[i] = 0;

  snew(nrdf2,natoms);
  aloop = gmx_mtop_atomloop_all_init(mtop);
  while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) {
    nrdf2[i] = 0;
    if (atom->ptype == eptAtom || atom->ptype == eptNucleus) {
      g = ggrpnr(groups,egcFREEZE,i);
      /* Double count nrdf for particle i */
      for(d=0; d<DIM; d++) {
	if (opts->nFreeze[g][d] == 0) {
	  nrdf2[i] += 2;
	}
      }
      nrdf_tc [ggrpnr(groups,egcTC ,i)] += 0.5*nrdf2[i];
      nrdf_vcm[ggrpnr(groups,egcVCM,i)] += 0.5*nrdf2[i];
    }
  }

  as = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molb = &mtop->molblock[mb];
    molt = &mtop->moltype[molb->type];
    atom = molt->atoms.atom;
    for(mol=0; mol<molb->nmol; mol++) {
      for (ftype=F_CONSTR; ftype<=F_CONSTRNC; ftype++) {
	ia = molt->ilist[ftype].iatoms;
	for(i=0; i<molt->ilist[ftype].nr; ) {
	  /* Subtract degrees of freedom for the constraints,
	   * if the particles still have degrees of freedom left.
	   * If one of the particles is a vsite or a shell, then all
	   * constraint motion will go there, but since they do not
	   * contribute to the constraints the degrees of freedom do not
	   * change.
	   */
	  ai = as + ia[1];
	  aj = as + ia[2];
	  if (((atom[ia[1]].ptype == eptNucleus) ||
	       (atom[ia[1]].ptype == eptAtom)) &&
	      ((atom[ia[2]].ptype == eptNucleus) ||
	       (atom[ia[2]].ptype == eptAtom))) {
	    if (nrdf2[ai] > 0) 
	      jmin = 1;
	    else
	      jmin = 2;
	    if (nrdf2[aj] > 0)
	      imin = 1;
	    else
	      imin = 2;
	    imin = min(imin,nrdf2[ai]);
	    jmin = min(jmin,nrdf2[aj]);
	    nrdf2[ai] -= imin;
	    nrdf2[aj] -= jmin;
	    nrdf_tc [ggrpnr(groups,egcTC ,ai)] -= 0.5*imin;
	    nrdf_tc [ggrpnr(groups,egcTC ,aj)] -= 0.5*jmin;
	    nrdf_vcm[ggrpnr(groups,egcVCM,ai)] -= 0.5*imin;
	    nrdf_vcm[ggrpnr(groups,egcVCM,aj)] -= 0.5*jmin;
	  }
	  ia += interaction_function[ftype].nratoms+1;
	  i  += interaction_function[ftype].nratoms+1;
	}
      }
      ia = molt->ilist[F_SETTLE].iatoms;
      for(i=0; i<molt->ilist[F_SETTLE].nr; ) {
	/* Subtract 1 dof from every atom in the SETTLE */
	for(ai=as+ia[1]; ai<as+ia[1]+3; ai++) {
	  imin = min(2,nrdf2[ai]);
	  nrdf2[ai] -= imin;
	  nrdf_tc [ggrpnr(groups,egcTC ,ai)] -= 0.5*imin;
	  nrdf_vcm[ggrpnr(groups,egcVCM,ai)] -= 0.5*imin;
	}
	ia += 2;
	i  += 2;
      }
      as += molt->atoms.nr;
    }
  }

  if (ir->ePull == epullCONSTRAINT) {
    /* Correct nrdf for the COM constraints.
     * We correct using the TC and VCM group of the first atom
     * in the reference and pull group. If atoms in one pull group
     * belong to different TC or VCM groups it is anyhow difficult
     * to determine the optimal nrdf assignment.
     */
    pull = ir->pull;
    if (pull->eGeom == epullgPOS) {
      nc = 0;
      for(i=0; i<DIM; i++) {
	if (pull->dim[i])
	  nc++;
      }
    } else {
      nc = 1;
    }
    for(i=0; i<pull->ngrp; i++) {
      imin = 2*nc;
      if (pull->grp[0].nat > 0) {
	/* Subtract 1/2 dof from the reference group */
	ai = pull->grp[0].ind[0];
	if (nrdf_tc[ggrpnr(groups,egcTC,ai)] > 1) {
	  nrdf_tc [ggrpnr(groups,egcTC ,ai)] -= 0.5;
	  nrdf_vcm[ggrpnr(groups,egcVCM,ai)] -= 0.5;
	  imin--;
	}
      }
      /* Subtract 1/2 dof from the pulled group */
      ai = pull->grp[1+i].ind[0];
      nrdf_tc [ggrpnr(groups,egcTC ,ai)] -= 0.5*imin;
      nrdf_vcm[ggrpnr(groups,egcVCM,ai)] -= 0.5*imin;
      if (nrdf_tc[ggrpnr(groups,egcTC,ai)] < 0)
	gmx_fatal(FARGS,"Center of mass pulling constraints caused the number of degrees of freedom for temperature coupling group %s to be negative",gnames[groups->grps[egcTC].nm_ind[ggrpnr(groups,egcTC,ai)]]);
    }
  }
  
  if (ir->nstcomm != 0) {
    /* Subtract 3 from the number of degrees of freedom in each vcm group
     * when com translation is removed and 6 when rotation is removed
     * as well.
     */
    switch (ir->comm_mode) {
    case ecmLINEAR:
      n_sub = ndof_com(ir);
      break;
    case ecmANGULAR:
      n_sub = 6;
      break;
    default:
      n_sub = 0;
      gmx_incons("Checking comm_mode");
    }
    
    for(i=0; i<groups->grps[egcTC].nr; i++) {
      /* Count the number of atoms of TC group i for every VCM group */
      for(j=0; j<groups->grps[egcVCM].nr+1; j++)
	na_vcm[j] = 0;
      na_tot = 0;
      for(ai=0; ai<natoms; ai++)
	if (ggrpnr(groups,egcTC,ai) == i) {
	  na_vcm[ggrpnr(groups,egcVCM,ai)]++;
	  na_tot++;
	}
      /* Correct for VCM removal according to the fraction of each VCM
       * group present in this TC group.
       */
      nrdf_uc = nrdf_tc[i];
      if (debug) {
	fprintf(debug,"T-group[%d] nrdf_uc = %g, n_sub = %g\n",
		i,nrdf_uc,n_sub);
      }
      nrdf_tc[i] = 0;
      for(j=0; j<groups->grps[egcVCM].nr+1; j++) {
	if (nrdf_vcm[j] > n_sub) {
	  nrdf_tc[i] += nrdf_uc*((double)na_vcm[j]/(double)na_tot)*
	    (nrdf_vcm[j] - n_sub)/nrdf_vcm[j];
	}
	if (debug) {
	  fprintf(debug,"  nrdf_vcm[%d] = %g, nrdf = %g\n",
		  j,nrdf_vcm[j],nrdf_tc[i]);
	}
      }
    }
  }
  for(i=0; (i<groups->grps[egcTC].nr); i++) {
    opts->nrdf[i] = nrdf_tc[i];
    if (opts->nrdf[i] < 0)
      opts->nrdf[i] = 0;
    fprintf(stderr,
	    "Number of degrees of freedom in T-Coupling group %s is %.2f\n",
	    gnames[groups->grps[egcTC].nm_ind[i]],opts->nrdf[i]);
  }
  
  sfree(nrdf2);
  sfree(nrdf_tc);
  sfree(nrdf_vcm);
  sfree(na_vcm);
}

static void decode_cos(char *s,t_cosines *cosine,gmx_bool bTime)
{
  char   *t;
  char   format[STRLEN],f1[STRLEN];
  double a,phi;
  int    i;
  
  t=strdup(s);
  trim(t);
  
  cosine->n=0;
  cosine->a=NULL;
  cosine->phi=NULL;
  if (strlen(t)) {
    sscanf(t,"%d",&(cosine->n));
    if (cosine->n <= 0) {
      cosine->n=0;
    } else {
      snew(cosine->a,cosine->n);
      snew(cosine->phi,cosine->n);
      
      sprintf(format,"%%*d");
      for(i=0; (i<cosine->n); i++) {
	strcpy(f1,format);
	strcat(f1,"%lf%lf");
	if (sscanf(t,f1,&a,&phi) < 2)
	  gmx_fatal(FARGS,"Invalid input for electric field shift: '%s'",t);
	cosine->a[i]=a;
	cosine->phi[i]=phi;
	strcat(format,"%*lf%*lf");
      }
    }
  }
  sfree(t);
}

static gmx_bool do_egp_flag(t_inputrec *ir,gmx_groups_t *groups,
			const char *option,const char *val,int flag)
{
  /* The maximum number of energy group pairs would be MAXPTR*(MAXPTR+1)/2.
   * But since this is much larger than STRLEN, such a line can not be parsed.
   * The real maximum is the number of names that fit in a string: STRLEN/2.
   */
#define EGP_MAX (STRLEN/2)
  int  nelem,i,j,k,nr;
  char *names[EGP_MAX];
  char ***gnames;
  gmx_bool bSet;

  gnames = groups->grpname;

  nelem = str_nelem(val,EGP_MAX,names);
  if (nelem % 2 != 0)
    gmx_fatal(FARGS,"The number of groups for %s is odd",option);
  nr = groups->grps[egcENER].nr;
  bSet = FALSE;
  for(i=0; i<nelem/2; i++) {
    j = 0;
    while ((j < nr) &&
	   gmx_strcasecmp(names[2*i],*(gnames[groups->grps[egcENER].nm_ind[j]])))
      j++;
    if (j == nr)
      gmx_fatal(FARGS,"%s in %s is not an energy group\n",
		  names[2*i],option);
    k = 0;
    while ((k < nr) &&
	   gmx_strcasecmp(names[2*i+1],*(gnames[groups->grps[egcENER].nm_ind[k]])))
      k++;
    if (k==nr)
      gmx_fatal(FARGS,"%s in %s is not an energy group\n",
	      names[2*i+1],option);
    if ((j < nr) && (k < nr)) {
      ir->opts.egp_flags[nr*j+k] |= flag;
      ir->opts.egp_flags[nr*k+j] |= flag;
      bSet = TRUE;
    }
  }

  return bSet;
}

void do_index(const char* mdparin, const char *ndx,
              gmx_mtop_t *mtop,
              gmx_bool bVerbose,
              t_inputrec *ir,rvec *v,
              warninp_t wi)
{
  t_blocka *grps;
  gmx_groups_t *groups;
  int     natoms;
  t_symtab *symtab;
  t_atoms atoms_all;
  char    warnbuf[STRLEN],**gnames;
  int     nr,ntcg,ntau_t,nref_t,nacc,nofg,nSA,nSA_points,nSA_time,nSA_temp;
  real    tau_min;
  int     nstcmin;
  int     nacg,nfreeze,nfrdim,nenergy,nvcm,nuser;
  char    *ptr1[MAXPTR],*ptr2[MAXPTR],*ptr3[MAXPTR];
  int     i,j,k,restnm;
  real    SAtime;
  gmx_bool    bExcl,bTable,bSetTCpar,bAnneal,bRest;
  int     nQMmethod,nQMbasis,nQMcharge,nQMmult,nbSH,nCASorb,nCASelec,
    nSAon,nSAoff,nSAsteps,nQMg,nbOPT,nbTS;
  char    warn_buf[STRLEN];

  if (bVerbose)
    fprintf(stderr,"processing index file...\n");
  debug_gmx();
  if (ndx == NULL) {
    snew(grps,1);
    snew(grps->index,1);
    snew(gnames,1);
    atoms_all = gmx_mtop_global_atoms(mtop);
    analyse(&atoms_all,grps,&gnames,FALSE,TRUE);
    free_t_atoms(&atoms_all,FALSE);
  } else {
    grps = init_index(ndx,&gnames);
  }

  groups = &mtop->groups;
  natoms = mtop->natoms;
  symtab = &mtop->symtab;

  snew(groups->grpname,grps->nr+1);
  
  for(i=0; (i<grps->nr); i++) {
    groups->grpname[i] = put_symtab(symtab,gnames[i]);
  }
  groups->grpname[i] = put_symtab(symtab,"rest");
  restnm=i;
  srenew(gnames,grps->nr+1);
  gnames[restnm] = *(groups->grpname[i]);
  groups->ngrpname = grps->nr+1;

  set_warning_line(wi,mdparin,-1);

  ntau_t = str_nelem(tau_t,MAXPTR,ptr1);
  nref_t = str_nelem(ref_t,MAXPTR,ptr2);
  ntcg   = str_nelem(tcgrps,MAXPTR,ptr3);
  if ((ntau_t != ntcg) || (nref_t != ntcg)) {
    gmx_fatal(FARGS,"Invalid T coupling input: %d groups, %d ref_t values and "
		"%d tau_t values",ntcg,nref_t,ntau_t);
  }

  bSetTCpar = (ir->etc || EI_SD(ir->eI) || ir->eI==eiBD || EI_TPI(ir->eI));
  do_numbering(natoms,groups,ntcg,ptr3,grps,gnames,egcTC,
               restnm,bSetTCpar ? egrptpALL : egrptpALL_GENREST,bVerbose,wi);
  nr = groups->grps[egcTC].nr;
  ir->opts.ngtc = nr;
  snew(ir->opts.nrdf,nr);
  snew(ir->opts.tau_t,nr);
  snew(ir->opts.ref_t,nr);
  if (ir->eI==eiBD && ir->bd_fric==0) {
    fprintf(stderr,"bd_fric=0, so tau_t will be used as the inverse friction constant(s)\n"); 
  }

  if (bSetTCpar)
  {
      if (nr != nref_t)
      {
          gmx_fatal(FARGS,"Not enough ref_t and tau_t values!");
      }
      
      tau_min = 1e20;
      for(i=0; (i<nr); i++)
      {
          ir->opts.tau_t[i] = strtod(ptr1[i],NULL);
          if ((ir->eI == eiBD || ir->eI == eiSD2) && ir->opts.tau_t[i] <= 0)
          {
              sprintf(warn_buf,"With integrator %s tau_t should be larger than 0",ei_names[ir->eI]);
              warning_error(wi,warn_buf);
          }
          if ((ir->etc == etcVRESCALE && ir->opts.tau_t[i] >= 0) || 
              (ir->etc != etcVRESCALE && ir->opts.tau_t[i] >  0))
          {
              tau_min = min(tau_min,ir->opts.tau_t[i]);
          }
      }
      if (ir->etc != etcNO && ir->nsttcouple == -1)
      {
            ir->nsttcouple = ir_optimal_nsttcouple(ir);
      }
      if (EI_VV(ir->eI)) 
      {
          if ((ir->epc==epcMTTK) && (ir->etc>etcNO))
          {
              int mincouple;
              mincouple = ir->nsttcouple;
              if (ir->nstpcouple < mincouple)
              {
                  mincouple = ir->nstpcouple;
              }
              ir->nstpcouple = mincouple;
              ir->nsttcouple = mincouple;
              warning_note(wi,"for current Trotter decomposition methods with vv, nsttcouple and nstpcouple must be equal.  Both have been reset to min(nsttcouple,nstpcouple)");
          }
      }
      nstcmin = tcouple_min_integration_steps(ir->etc);
      if (nstcmin > 1)
      {
          if (tau_min/(ir->delta_t*ir->nsttcouple) < nstcmin)
          {
              sprintf(warn_buf,"For proper integration of the %s thermostat, tau_t (%g) should be at least %d times larger than nsttcouple*dt (%g)",
                      ETCOUPLTYPE(ir->etc),
                      tau_min,nstcmin,
                      ir->nsttcouple*ir->delta_t);
              warning(wi,warn_buf);
          }
      }
      for(i=0; (i<nr); i++)
      {
          ir->opts.ref_t[i] = strtod(ptr2[i],NULL);
          if (ir->opts.ref_t[i] < 0)
          {
              gmx_fatal(FARGS,"ref_t for group %d negative",i);
          }
      }
  }
    
  /* Simulated annealing for each group. There are nr groups */
  nSA = str_nelem(anneal,MAXPTR,ptr1);
  if (nSA == 1 && (ptr1[0][0]=='n' || ptr1[0][0]=='N'))
     nSA = 0;
  if(nSA>0 && nSA != nr) 
    gmx_fatal(FARGS,"Not enough annealing values: %d (for %d groups)\n",nSA,nr);
  else {
    snew(ir->opts.annealing,nr);
    snew(ir->opts.anneal_npoints,nr);
    snew(ir->opts.anneal_time,nr);
    snew(ir->opts.anneal_temp,nr);
    for(i=0;i<nr;i++) {
      ir->opts.annealing[i]=eannNO;
      ir->opts.anneal_npoints[i]=0;
      ir->opts.anneal_time[i]=NULL;
      ir->opts.anneal_temp[i]=NULL;
    }
    if (nSA > 0) {
      bAnneal=FALSE;
      for(i=0;i<nr;i++) { 
	if(ptr1[i][0]=='n' || ptr1[i][0]=='N') {
	  ir->opts.annealing[i]=eannNO;
	} else if(ptr1[i][0]=='s'|| ptr1[i][0]=='S') {
	  ir->opts.annealing[i]=eannSINGLE;
	  bAnneal=TRUE;
	} else if(ptr1[i][0]=='p'|| ptr1[i][0]=='P') {
	  ir->opts.annealing[i]=eannPERIODIC;
	  bAnneal=TRUE;
	} 
      } 
      if(bAnneal) {
	/* Read the other fields too */
	nSA_points = str_nelem(anneal_npoints,MAXPTR,ptr1);
	if(nSA_points!=nSA) 
	  gmx_fatal(FARGS,"Found %d annealing_npoints values for %d groups\n",nSA_points,nSA);
	for(k=0,i=0;i<nr;i++) {
	  ir->opts.anneal_npoints[i]=strtol(ptr1[i],NULL,10);
	  if(ir->opts.anneal_npoints[i]==1)
	    gmx_fatal(FARGS,"Please specify at least a start and an end point for annealing\n");
	  snew(ir->opts.anneal_time[i],ir->opts.anneal_npoints[i]);
	  snew(ir->opts.anneal_temp[i],ir->opts.anneal_npoints[i]);
	  k += ir->opts.anneal_npoints[i];
	}

	nSA_time = str_nelem(anneal_time,MAXPTR,ptr1);
	if(nSA_time!=k) 
	  gmx_fatal(FARGS,"Found %d annealing_time values, wanter %d\n",nSA_time,k);
	nSA_temp = str_nelem(anneal_temp,MAXPTR,ptr2);
	if(nSA_temp!=k) 
	  gmx_fatal(FARGS,"Found %d annealing_temp values, wanted %d\n",nSA_temp,k);

	for(i=0,k=0;i<nr;i++) {
	  
	  for(j=0;j<ir->opts.anneal_npoints[i];j++) {
	    ir->opts.anneal_time[i][j]=strtod(ptr1[k],NULL);
	    ir->opts.anneal_temp[i][j]=strtod(ptr2[k],NULL);
	    if(j==0) {
	      if(ir->opts.anneal_time[i][0] > (ir->init_t+GMX_REAL_EPS))
		gmx_fatal(FARGS,"First time point for annealing > init_t.\n");      
	    } else { 
	      /* j>0 */
	      if(ir->opts.anneal_time[i][j]<ir->opts.anneal_time[i][j-1])
		gmx_fatal(FARGS,"Annealing timepoints out of order: t=%f comes after t=%f\n",
			    ir->opts.anneal_time[i][j],ir->opts.anneal_time[i][j-1]);
	    }
	    if(ir->opts.anneal_temp[i][j]<0) 
	      gmx_fatal(FARGS,"Found negative temperature in annealing: %f\n",ir->opts.anneal_temp[i][j]);    
	    k++;
	  }
	}
	/* Print out some summary information, to make sure we got it right */
	for(i=0,k=0;i<nr;i++) {
	  if(ir->opts.annealing[i]!=eannNO) {
	    j = groups->grps[egcTC].nm_ind[i];
	    fprintf(stderr,"Simulated annealing for group %s: %s, %d timepoints\n",
		    *(groups->grpname[j]),eann_names[ir->opts.annealing[i]],
		    ir->opts.anneal_npoints[i]);
	    fprintf(stderr,"Time (ps)   Temperature (K)\n");
	    /* All terms except the last one */
	    for(j=0;j<(ir->opts.anneal_npoints[i]-1);j++) 
		fprintf(stderr,"%9.1f      %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	    
	    /* Finally the last one */
	    j = ir->opts.anneal_npoints[i]-1;
	    if(ir->opts.annealing[i]==eannSINGLE)
	      fprintf(stderr,"%9.1f-     %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	    else {
	      fprintf(stderr,"%9.1f      %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	      if(fabs(ir->opts.anneal_temp[i][j]-ir->opts.anneal_temp[i][0])>GMX_REAL_EPS)
		warning_note(wi,"There is a temperature jump when your annealing loops back.\n");
	    }
	  }
	} 
      }
    }
  }	

  if (ir->ePull != epullNO) {
    make_pull_groups(ir->pull,pull_grp,grps,gnames);
  }

  nacc = str_nelem(acc,MAXPTR,ptr1);
  nacg = str_nelem(accgrps,MAXPTR,ptr2);
  if (nacg*DIM != nacc)
    gmx_fatal(FARGS,"Invalid Acceleration input: %d groups and %d acc. values",
		nacg,nacc);
  do_numbering(natoms,groups,nacg,ptr2,grps,gnames,egcACC,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  nr = groups->grps[egcACC].nr;
  snew(ir->opts.acc,nr);
  ir->opts.ngacc=nr;
  
  for(i=k=0; (i<nacg); i++)
    for(j=0; (j<DIM); j++,k++)
      ir->opts.acc[i][j]=strtod(ptr1[k],NULL);
  for( ;(i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.acc[i][j]=0;
  
  nfrdim  = str_nelem(frdim,MAXPTR,ptr1);
  nfreeze = str_nelem(freeze,MAXPTR,ptr2);
  if (nfrdim != DIM*nfreeze)
    gmx_fatal(FARGS,"Invalid Freezing input: %d groups and %d freeze values",
		nfreeze,nfrdim);
  do_numbering(natoms,groups,nfreeze,ptr2,grps,gnames,egcFREEZE,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  nr = groups->grps[egcFREEZE].nr;
  ir->opts.ngfrz=nr;
  snew(ir->opts.nFreeze,nr);
  for(i=k=0; (i<nfreeze); i++)
    for(j=0; (j<DIM); j++,k++) {
      ir->opts.nFreeze[i][j]=(gmx_strncasecmp(ptr1[k],"Y",1)==0);
      if (!ir->opts.nFreeze[i][j]) {
	if (gmx_strncasecmp(ptr1[k],"N",1) != 0) {
	  sprintf(warnbuf,"Please use Y(ES) or N(O) for freezedim only "
		  "(not %s)", ptr1[k]);
	  warning(wi,warn_buf);
	}
      }
    }
  for( ; (i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.nFreeze[i][j]=0;
  
  nenergy=str_nelem(energy,MAXPTR,ptr1);
  do_numbering(natoms,groups,nenergy,ptr1,grps,gnames,egcENER,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  add_wall_energrps(groups,ir->nwall,symtab);
  ir->opts.ngener = groups->grps[egcENER].nr;
  nvcm=str_nelem(vcm,MAXPTR,ptr1);
  bRest =
    do_numbering(natoms,groups,nvcm,ptr1,grps,gnames,egcVCM,
                 restnm,nvcm==0 ? egrptpALL_GENREST : egrptpPART,bVerbose,wi);
  if (bRest) {
    warning(wi,"Some atoms are not part of any center of mass motion removal group.\n"
	    "This may lead to artifacts.\n"
	    "In most cases one should use one group for the whole system.");
  }

  /* Now we have filled the freeze struct, so we can calculate NRDF */ 
  calc_nrdf(mtop,ir,gnames);

  if (v && NULL) {
    real fac,ntot=0;
    
    /* Must check per group! */
    for(i=0; (i<ir->opts.ngtc); i++) 
      ntot += ir->opts.nrdf[i];
    if (ntot != (DIM*natoms)) {
      fac = sqrt(ntot/(DIM*natoms));
      if (bVerbose)
	fprintf(stderr,"Scaling velocities by a factor of %.3f to account for constraints\n"
		"and removal of center of mass motion\n",fac);
      for(i=0; (i<natoms); i++)
	svmul(fac,v[i],v[i]);
    }
  }
  
  nuser=str_nelem(user1,MAXPTR,ptr1);
  do_numbering(natoms,groups,nuser,ptr1,grps,gnames,egcUser1,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  nuser=str_nelem(user2,MAXPTR,ptr1);
  do_numbering(natoms,groups,nuser,ptr1,grps,gnames,egcUser2,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  nuser=str_nelem(xtc_grps,MAXPTR,ptr1);
  do_numbering(natoms,groups,nuser,ptr1,grps,gnames,egcXTC,
               restnm,egrptpONE,bVerbose,wi);
  nofg = str_nelem(orirefitgrp,MAXPTR,ptr1);
  do_numbering(natoms,groups,nofg,ptr1,grps,gnames,egcORFIT,
               restnm,egrptpALL_GENREST,bVerbose,wi);

  /* QMMM input processing */
  nQMg          = str_nelem(QMMM,MAXPTR,ptr1);
  nQMmethod     = str_nelem(QMmethod,MAXPTR,ptr2);
  nQMbasis      = str_nelem(QMbasis,MAXPTR,ptr3);
  if((nQMmethod != nQMg)||(nQMbasis != nQMg)){
    gmx_fatal(FARGS,"Invalid QMMM input: %d groups %d basissets"
	      " and %d methods\n",nQMg,nQMbasis,nQMmethod);
  }
  /* group rest, if any, is always MM! */
  do_numbering(natoms,groups,nQMg,ptr1,grps,gnames,egcQMMM,
               restnm,egrptpALL_GENREST,bVerbose,wi);
  nr = nQMg; /*atoms->grps[egcQMMM].nr;*/
  ir->opts.ngQM = nQMg;
  snew(ir->opts.QMmethod,nr);
  snew(ir->opts.QMbasis,nr);
  for(i=0;i<nr;i++){
    /* input consists of strings: RHF CASSCF PM3 .. These need to be
     * converted to the corresponding enum in names.c
     */
    ir->opts.QMmethod[i] = search_QMstring(ptr2[i],eQMmethodNR,
                                           eQMmethod_names);
    ir->opts.QMbasis[i]  = search_QMstring(ptr3[i],eQMbasisNR,
                                           eQMbasis_names);

  }
  nQMmult   = str_nelem(QMmult,MAXPTR,ptr1);
  nQMcharge = str_nelem(QMcharge,MAXPTR,ptr2);
  nbSH      = str_nelem(bSH,MAXPTR,ptr3);
  snew(ir->opts.QMmult,nr);
  snew(ir->opts.QMcharge,nr);
  snew(ir->opts.bSH,nr);

  for(i=0;i<nr;i++){
    ir->opts.QMmult[i]   = strtol(ptr1[i],NULL,10);
    ir->opts.QMcharge[i] = strtol(ptr2[i],NULL,10);
    ir->opts.bSH[i]      = (gmx_strncasecmp(ptr3[i],"Y",1)==0);
  }

  nCASelec  = str_nelem(CASelectrons,MAXPTR,ptr1);
  nCASorb   = str_nelem(CASorbitals,MAXPTR,ptr2);
  snew(ir->opts.CASelectrons,nr);
  snew(ir->opts.CASorbitals,nr);
  for(i=0;i<nr;i++){
    ir->opts.CASelectrons[i]= strtol(ptr1[i],NULL,10);
    ir->opts.CASorbitals[i] = strtol(ptr2[i],NULL,10);
  }
  /* special optimization options */

  nbOPT = str_nelem(bOPT,MAXPTR,ptr1);
  nbTS = str_nelem(bTS,MAXPTR,ptr2);
  snew(ir->opts.bOPT,nr);
  snew(ir->opts.bTS,nr);
  for(i=0;i<nr;i++){
    ir->opts.bOPT[i] = (gmx_strncasecmp(ptr1[i],"Y",1)==0);
    ir->opts.bTS[i]  = (gmx_strncasecmp(ptr2[i],"Y",1)==0);
  }
  nSAon     = str_nelem(SAon,MAXPTR,ptr1);
  nSAoff    = str_nelem(SAoff,MAXPTR,ptr2);
  nSAsteps  = str_nelem(SAsteps,MAXPTR,ptr3);
  snew(ir->opts.SAon,nr);
  snew(ir->opts.SAoff,nr);
  snew(ir->opts.SAsteps,nr);

  for(i=0;i<nr;i++){
    ir->opts.SAon[i]    = strtod(ptr1[i],NULL);
    ir->opts.SAoff[i]   = strtod(ptr2[i],NULL);
    ir->opts.SAsteps[i] = strtol(ptr3[i],NULL,10);
  }
  /* end of QMMM input */

  if (bVerbose)
    for(i=0; (i<egcNR); i++) {
      fprintf(stderr,"%-16s has %d element(s):",gtypes[i],groups->grps[i].nr); 
      for(j=0; (j<groups->grps[i].nr); j++)
	fprintf(stderr," %s",*(groups->grpname[groups->grps[i].nm_ind[j]]));
      fprintf(stderr,"\n");
    }

  nr = groups->grps[egcENER].nr;
  snew(ir->opts.egp_flags,nr*nr);

  bExcl = do_egp_flag(ir,groups,"energygrp_excl",egpexcl,EGP_EXCL);
  if (bExcl && EEL_FULL(ir->coulombtype))
    warning(wi,"Can not exclude the lattice Coulomb energy between energy groups");

  bTable = do_egp_flag(ir,groups,"energygrp_table",egptable,EGP_TABLE);
  if (bTable && !(ir->vdwtype == evdwUSER) && 
      !(ir->coulombtype == eelUSER) && !(ir->coulombtype == eelPMEUSER) &&
      !(ir->coulombtype == eelPMEUSERSWITCH))
    gmx_fatal(FARGS,"Can only have energy group pair tables in combination with user tables for VdW and/or Coulomb");

  decode_cos(efield_x,&(ir->ex[XX]),FALSE);
  decode_cos(efield_xt,&(ir->et[XX]),TRUE);
  decode_cos(efield_y,&(ir->ex[YY]),FALSE);
  decode_cos(efield_yt,&(ir->et[YY]),TRUE);
  decode_cos(efield_z,&(ir->ex[ZZ]),FALSE);
  decode_cos(efield_zt,&(ir->et[ZZ]),TRUE);
  
  for(i=0; (i<grps->nr); i++)
    sfree(gnames[i]);
  sfree(gnames);
  done_blocka(grps);
  sfree(grps);

}



static void check_disre(gmx_mtop_t *mtop)
{
  gmx_ffparams_t *ffparams;
  t_functype *functype;
  t_iparams  *ip;
  int i,ndouble,ftype;
  int label,old_label;
  
  if (gmx_mtop_ftype_count(mtop,F_DISRES) > 0) {
    ffparams  = &mtop->ffparams;
    functype  = ffparams->functype;
    ip        = ffparams->iparams;
    ndouble   = 0;
    old_label = -1;
    for(i=0; i<ffparams->ntypes; i++) {
      ftype = functype[i];
      if (ftype == F_DISRES) {
	label = ip[i].disres.label;
	if (label == old_label) {
	  fprintf(stderr,"Distance restraint index %d occurs twice\n",label);
	  ndouble++;
	}
	old_label = label;
      }
    }
    if (ndouble>0)
      gmx_fatal(FARGS,"Found %d double distance restraint indices,\n"
		"probably the parameters for multiple pairs in one restraint "
		"are not identical\n",ndouble);
  }
}

static gmx_bool absolute_reference(t_inputrec *ir,gmx_mtop_t *sys,ivec AbsRef)
{
  int d,g,i;
  gmx_mtop_ilistloop_t iloop;
  t_ilist *ilist;
  int nmol;
  t_iparams *pr;

  /* Check the COM */
  for(d=0; d<DIM; d++) {
    AbsRef[d] = (d < ndof_com(ir) ? 0 : 1);
  }
  /* Check for freeze groups */
  for(g=0; g<ir->opts.ngfrz; g++) {
    for(d=0; d<DIM; d++) {
      if (ir->opts.nFreeze[g][d] != 0) {
	AbsRef[d] = 1;
      }
    }
  }
  /* Check for position restraints */
  iloop = gmx_mtop_ilistloop_init(sys);
  while (gmx_mtop_ilistloop_next(iloop,&ilist,&nmol)) {
    if (nmol > 0) {
      for(i=0; i<ilist[F_POSRES].nr; i+=2) {
	pr = &sys->ffparams.iparams[ilist[F_POSRES].iatoms[i]];
	for(d=0; d<DIM; d++) {
	  if (pr->posres.fcA[d] != 0) {
	    AbsRef[d] = 1;
	  }
	}
      }
    }
  }

  return (AbsRef[XX] != 0 && AbsRef[YY] != 0 && AbsRef[ZZ] != 0);
}

void triple_check(const char *mdparin,t_inputrec *ir,gmx_mtop_t *sys,
                  warninp_t wi)
{
  char err_buf[256];
  int  i,m,g,nmol,npct;
  gmx_bool bCharge,bAcc;
  real gdt_max,*mgrp,mt;
  rvec acc;
  gmx_mtop_atomloop_block_t aloopb;
  gmx_mtop_atomloop_all_t aloop;
  t_atom *atom;
  ivec AbsRef;
  char warn_buf[STRLEN];

  set_warning_line(wi,mdparin,-1);

  if (EI_DYNAMICS(ir->eI) && !EI_SD(ir->eI) && ir->eI != eiBD &&
      ir->comm_mode == ecmNO &&
      !(absolute_reference(ir,sys,AbsRef) || ir->nsteps <= 10)) {
    warning(wi,"You are not using center of mass motion removal (mdp option comm-mode), numerical rounding errors can lead to build up of kinetic energy of the center of mass");
  }
  
  bCharge = FALSE;
  aloopb = gmx_mtop_atomloop_block_init(sys);
  while (gmx_mtop_atomloop_block_next(aloopb,&atom,&nmol)) {
    if (atom->q != 0 || atom->qB != 0) {
      bCharge = TRUE;
    }
  }
  
  if (!bCharge) {
    if (EEL_FULL(ir->coulombtype)) {
      sprintf(err_buf,
	      "You are using full electrostatics treatment %s for a system without charges.\n"
	      "This costs a lot of performance for just processing zeros, consider using %s instead.\n",
	      EELTYPE(ir->coulombtype),EELTYPE(eelCUT));
      warning(wi,err_buf);
    }
  } else {
    if (ir->coulombtype == eelCUT && ir->rcoulomb > 0 && !ir->implicit_solvent) {
      sprintf(err_buf,
	      "You are using a plain Coulomb cut-off, which might produce artifacts.\n"
	      "You might want to consider using %s electrostatics.\n",
	      EELTYPE(eelPME));
      warning_note(wi,err_buf);
    }
  }

  /* Generalized reaction field */  
  if (ir->opts.ngtc == 0) {
    sprintf(err_buf,"No temperature coupling while using coulombtype %s",
	    eel_names[eelGRF]);
    CHECK(ir->coulombtype == eelGRF);
  }
  else {
    sprintf(err_buf,"When using coulombtype = %s"
	    " ref_t for temperature coupling should be > 0",
	    eel_names[eelGRF]);
    CHECK((ir->coulombtype == eelGRF) && (ir->opts.ref_t[0] <= 0));
  }
    
  if (ir->eI == eiSD1) {
    gdt_max = 0;
    for(i=0; (i<ir->opts.ngtc); i++)
      gdt_max = max(gdt_max,ir->delta_t/ir->opts.tau_t[i]);
    if (0.5*gdt_max > 0.0015) {
      sprintf(warn_buf,"The relative error with integrator %s is 0.5*delta_t/tau_t = %g, you might want to switch to integrator %s\n",
	      ei_names[ir->eI],0.5*gdt_max,ei_names[eiSD2]);
      warning_note(wi,warn_buf);
    }
  }

  bAcc = FALSE;
  for(i=0; (i<sys->groups.grps[egcACC].nr); i++) {
    for(m=0; (m<DIM); m++) {
      if (fabs(ir->opts.acc[i][m]) > 1e-6) {
	bAcc = TRUE;
      }
    }
  }
  if (bAcc) {
    clear_rvec(acc);
    snew(mgrp,sys->groups.grps[egcACC].nr);
    aloop = gmx_mtop_atomloop_all_init(sys);
    while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) {
      mgrp[ggrpnr(&sys->groups,egcACC,i)] += atom->m;
    }
    mt = 0.0;
    for(i=0; (i<sys->groups.grps[egcACC].nr); i++) {
      for(m=0; (m<DIM); m++)
	acc[m] += ir->opts.acc[i][m]*mgrp[i];
      mt += mgrp[i];
    }
    for(m=0; (m<DIM); m++) {
      if (fabs(acc[m]) > 1e-6) {
	const char *dim[DIM] = { "X", "Y", "Z" };
	fprintf(stderr,
		"Net Acceleration in %s direction, will %s be corrected\n",
		dim[m],ir->nstcomm != 0 ? "" : "not");
	if (ir->nstcomm != 0 && m < ndof_com(ir)) {
	  acc[m] /= mt;
	  for (i=0; (i<sys->groups.grps[egcACC].nr); i++)
	    ir->opts.acc[i][m] -= acc[m];
	}
      }
    }
    sfree(mgrp);
  }

  if (ir->efep != efepNO && ir->sc_alpha != 0 &&
      !gmx_within_tol(sys->ffparams.reppow,12.0,10*GMX_DOUBLE_EPS)) {
    gmx_fatal(FARGS,"Soft-core interactions are only supported with VdW repulsion power 12");
  }

  if (ir->ePull != epullNO) {
    if (ir->pull->grp[0].nat == 0) {
      absolute_reference(ir,sys,AbsRef);
      for(m=0; m<DIM; m++) {
	if (ir->pull->dim[m] && !AbsRef[m]) {
	  warning(wi,"You are using an absolute reference for pulling, but the rest of the system does not have an absolute reference. This will lead to artifacts.");
	  break;
	}
      }
    }

    if (ir->pull->eGeom == epullgDIRPBC) {
      for(i=0; i<3; i++) {
	for(m=0; m<=i; m++) {
	  if ((ir->epc != epcNO && ir->compress[i][m] != 0) ||
	      ir->deform[i][m] != 0) {
	    for(g=1; g<ir->pull->ngrp; g++) {
	      if (ir->pull->grp[g].vec[m] != 0) {
		gmx_fatal(FARGS,"Can not have dynamic box while using pull geometry '%s' (dim %c)",EPULLGEOM(ir->pull->eGeom),'x'+m);
	      }
	    }
	  }
	}
      }
    }
  }

  check_disre(sys);
}

void double_check(t_inputrec *ir,matrix box,gmx_bool bConstr,warninp_t wi)
{
  real min_size;
  gmx_bool bTWIN;
  char warn_buf[STRLEN];
  const char *ptr;
  
  ptr = check_box(ir->ePBC,box);
  if (ptr) {
      warning_error(wi,ptr);
  }  

  if (bConstr && ir->eConstrAlg == econtSHAKE) {
    if (ir->shake_tol <= 0.0) {
      sprintf(warn_buf,"ERROR: shake_tol must be > 0 instead of %g\n",
              ir->shake_tol);
      warning_error(wi,warn_buf);
    }

    if (IR_TWINRANGE(*ir) && ir->nstlist > 1) {
      sprintf(warn_buf,"With twin-range cut-off's and SHAKE the virial and the pressure are incorrect.");
      if (ir->epc == epcNO) {
	warning(wi,warn_buf);
      } else {
          warning_error(wi,warn_buf);
      }
    }
  }

  if( (ir->eConstrAlg == econtLINCS) && bConstr) {
    /* If we have Lincs constraints: */
    if(ir->eI==eiMD && ir->etc==etcNO &&
       ir->eConstrAlg==econtLINCS && ir->nLincsIter==1) {
      sprintf(warn_buf,"For energy conservation with LINCS, lincs_iter should be 2 or larger.\n");
      warning_note(wi,warn_buf);
    }
    
    if ((ir->eI == eiCG || ir->eI == eiLBFGS) && (ir->nProjOrder<8)) {
      sprintf(warn_buf,"For accurate %s with LINCS constraints, lincs_order should be 8 or more.",ei_names[ir->eI]);
      warning_note(wi,warn_buf);
    }
    if (ir->epc==epcMTTK) {
        warning_error(wi,"MTTK not compatible with lincs -- use shake instead.");
    }
  }

  if (ir->LincsWarnAngle > 90.0) {
    sprintf(warn_buf,"lincs-warnangle can not be larger than 90 degrees, setting it to 90.\n");
    warning(wi,warn_buf);
    ir->LincsWarnAngle = 90.0;
  }

  if (ir->ePBC != epbcNONE) {
    if (ir->nstlist == 0) {
      warning(wi,"With nstlist=0 atoms are only put into the box at step 0, therefore drifting atoms might cause the simulation to crash.");
    }
    bTWIN = (ir->rlistlong > ir->rlist);
    if (ir->ns_type == ensGRID) {
      if (sqr(ir->rlistlong) >= max_cutoff2(ir->ePBC,box)) {
          sprintf(warn_buf,"ERROR: The cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element. Increase the box size or decrease %s.\n",
		bTWIN ? (ir->rcoulomb==ir->rlistlong ? "rcoulomb" : "rvdw"):"rlist");
          warning_error(wi,warn_buf);
      }
    } else {
      min_size = min(box[XX][XX],min(box[YY][YY],box[ZZ][ZZ]));
      if (2*ir->rlistlong >= min_size) {
          sprintf(warn_buf,"ERROR: One of the box lengths is smaller than twice the cut-off length. Increase the box size or decrease rlist.");
          warning_error(wi,warn_buf);
	if (TRICLINIC(box))
	  fprintf(stderr,"Grid search might allow larger cut-off's than simple search with triclinic boxes.");
      }
    }
  }
}

void check_chargegroup_radii(const gmx_mtop_t *mtop,const t_inputrec *ir,
                             rvec *x,
                             warninp_t wi)
{
    real rvdw1,rvdw2,rcoul1,rcoul2;
    char warn_buf[STRLEN];

    calc_chargegroup_radii(mtop,x,&rvdw1,&rvdw2,&rcoul1,&rcoul2);

    if (rvdw1 > 0)
    {
        printf("Largest charge group radii for Van der Waals: %5.3f, %5.3f nm\n",
               rvdw1,rvdw2);
    }
    if (rcoul1 > 0)
    {
        printf("Largest charge group radii for Coulomb:       %5.3f, %5.3f nm\n",
               rcoul1,rcoul2);
    }

    if (ir->rlist > 0)
    {
        if (rvdw1  + rvdw2  > ir->rlist ||
            rcoul1 + rcoul2 > ir->rlist)
        {
            sprintf(warn_buf,"The sum of the two largest charge group radii (%f) is larger than rlist (%f)\n",max(rvdw1+rvdw2,rcoul1+rcoul2),ir->rlist);
            warning(wi,warn_buf);
        }
        else
        {
            /* Here we do not use the zero at cut-off macro,
             * since user defined interactions might purposely
             * not be zero at the cut-off.
             */
            if (EVDW_IS_ZERO_AT_CUTOFF(ir->vdwtype) &&
                rvdw1 + rvdw2 > ir->rlist - ir->rvdw)
            {
                sprintf(warn_buf,"The sum of the two largest charge group radii (%f) is larger than rlist (%f) - rvdw (%f)\n",
                        rvdw1+rvdw2,
                        ir->rlist,ir->rvdw);
                if (ir_NVE(ir))
                {
                    warning(wi,warn_buf);
                }
                else
                {
                    warning_note(wi,warn_buf);
                }
            }
            if (EEL_IS_ZERO_AT_CUTOFF(ir->coulombtype) &&
                rcoul1 + rcoul2 > ir->rlistlong - ir->rcoulomb)
            {
                sprintf(warn_buf,"The sum of the two largest charge group radii (%f) is larger than %s (%f) - rcoulomb (%f)\n",
                        rcoul1+rcoul2,
                        ir->rlistlong > ir->rlist ? "rlistlong" : "rlist",
                        ir->rlistlong,ir->rcoulomb);
                if (ir_NVE(ir))
                {
                    warning(wi,warn_buf);
                }
                else
                {
                    warning_note(wi,warn_buf);
                }
            }
        }
    }
}
