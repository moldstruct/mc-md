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

#include <signal.h>
#include <stdlib.h>

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
#include "ionize.h"
#include "disre.h"
#include "orires.h"
#include "dihre.h"
#include "pppm.h"
#include "pme.h"
#include "mdatoms.h"
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "string2.h"
#include <time.h>
#include <omp.h>



#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

int USERINT1;
double debye_length;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//////////////////////////////////////////////////////////////////////////////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                 This is dictionary town                    //////////////
////////////         This town also houses other structures             //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////
//////////////////////////////////////////////////////////////////////////////////////
// struct for transition

struct transition {
    int* final_state;
    double rate;
                   };



// STRUCTURE FOR RATES

struct Rate {
int   mass;
int*  initial_state;
int   num_transitions;
int** final_states;
double* rates;          
};


// Function to count the number of lines in a file
static int countLinesInFile(const char* filePath) {
    FILE* file = fopen(filePath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening the file in '%s' in 'countLinesInFile'.\n",filePath);
        return -1;
    }

    int count = 0;
    char ch;
    int lastCharWasNewline = 0;  // To track if the last character was a newline

    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            count++;
            lastCharWasNewline = 1;
        } else {
            lastCharWasNewline = 0;
        }
    }

    // Check if the last character was not a newline
    if (lastCharWasNewline == 0) {
        count++;
    }

    // Reset the file pointer to the beginning
    rewind(file);
    fclose(file);

    return count;
}

// Function to initialize a struct Rate from file data
static struct Rate initializeRateFromFile(FILE* file, int mass) {
    struct Rate rate;
    rate.mass = mass;

    // Read the line containing transitions
    char line[1024]; // Adjust the buffer size as needed
    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error reading transitions from the file.\n");
        exit(1);
    }

    // Read the initial state
    rate.initial_state = (int*)malloc(3 * sizeof(int));
    if (sscanf(line, "%d %d %d;", &rate.initial_state[0], &rate.initial_state[1], &rate.initial_state[2]) != 3) {
        fprintf(stderr, "Error reading initial state from the file.\n");
        exit(1);
    }

    // Split the line and determine the number of transitions
    // Count the number of semicolons to calculate the number of transitions
    int num_transitions = 0;
    int i;
    for (i = 0; line[i] != '\0'; i++) {
        if (line[i] == ';') {
            num_transitions++;
        }
    }

    rate.num_transitions = num_transitions-1;

    // Allocate memory for final_states and rates
    rate.final_states = (int**)malloc(rate.num_transitions * sizeof(int*));
    rate.rates = (double*)malloc(rate.num_transitions * sizeof(double));
    for (i = 0; i < rate.num_transitions; i++) {
        rate.final_states[i] = (int*)malloc(3 * sizeof(int));
    }
    // Split the line into individual transitions and process them
    char* token = strtok(line, ";");
    token = strtok(NULL, ";");
    int transitionIndex = 0;
    while (token != NULL) {
        sscanf(token, "%d %d %d %lf", &rate.final_states[transitionIndex][0], &rate.final_states[transitionIndex][1], &rate.final_states[transitionIndex][2], &rate.rates[transitionIndex]);
        transitionIndex++;

        if (transitionIndex == rate.num_transitions) {
            break; 
        }
        token = strtok(NULL, ";");
    }

    return rate;
}



//initialze an array of struct Rate from a file
static struct Rate* initializeRatesArrayFromFile(const char* filePath,int mass) {
    FILE* file = fopen(filePath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening the file in 'initializeRatesArrayFromFile'.\n");
        return NULL;
    }

    // Determine the number of Rate instances based on the number of lines in the file
    int numRates = countLinesInFile(filePath);

    // Create an array of struct Rate
    struct Rate* rates = (struct Rate*)malloc(numRates * sizeof(struct Rate));

    // Initialize each Rate instance from the file
    int i;
    for (i = 0; i < numRates; i++) {
        // Position the file pointer to the start of the current Rate data
        if (i > 0) {
            if (fscanf(file, "%*[\n]")) {
                fprintf(stderr, "Error reading newline character.\n");
                exit(1);
            }
        }

        rates[i] = initializeRateFromFile(file,mass);
    }

    // Close the file
    fclose(file);

    return rates;
}



static void freeRatesArray(struct Rate* rates, int numRates) {
    if (rates == NULL) {
        return; // Nothing to free
    }
    int i,j;
    for (i = 0; i < numRates; i++) {
        // Free dynamically allocated memory for each Rate instance
        free(rates[i].initial_state);
        for (j = 0; j < rates[i].num_transitions; j++) {
            free(rates[i].final_states[j]);
        }
        free(rates[i].final_states);
        free(rates[i].rates);
    }

    // Free memory for the array of Rate instances
    free(rates);
}


// Function to check if a state is in an array of Rate
static int RatesStateIndex(struct Rate* rates, int numRates, int* state) {
    int i;
    for (i = 0; i < numRates; i++) {
        if (rates[i].initial_state[0] == state[0] &&
            rates[i].initial_state[1] == state[1] &&
            rates[i].initial_state[2] == state[2]) {
            return i;  // return the index of matching state. Will always evaluate to true
        }
    }
    return -1;  // The state is not found in the array, false
}

////  STRUCTURE FOR COLLISIONAL STUFF


// Very similar to Rates but this has 5 rates per transition instead of 1
struct Coll {
int   mass;
int*  initial_state;
int   num_transitions;
int** final_states;
double** coll_rates;         // 5 rates for each transition

};


// Function to initialize a struct Coll from file data
static struct Coll initializeCollFromFile(FILE* file, int mass) {
    struct Coll coll;
    coll.mass = mass;

    // Read the line containing transitions
    char line[1024]; // Adjust the buffer size as needed
    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error reading transitions from the file.\n");
        exit(1);
    }

    // Read the initial state
    coll.initial_state = (int*)malloc(3 * sizeof(int));
    if (sscanf(line, "%d %d %d;", &coll.initial_state[0], &coll.initial_state[1], &coll.initial_state[2]) != 3) {
        fprintf(stderr, "Error reading initial state from the file.\n");
        exit(1);
    }

    // Split the line and determine the number of transitions
    // Count the number of semicolons to calculate the number of transitions
    int num_transitions = 0;
    int i;
    for (i = 0; line[i] != '\0'; i++) {
        if (line[i] == ';') {
            num_transitions++;
        }
    }

    coll.num_transitions = num_transitions-1;

    // Allocate memory for final_states and rates
    coll.final_states = (int**)malloc(coll.num_transitions * sizeof(int*));
    coll.coll_rates = (double**)malloc(coll.num_transitions * sizeof(double*));
    for (i = 0; i < coll.num_transitions; i++) {
        coll.final_states[i] = (int*)malloc(3 * sizeof(int));
        coll.coll_rates[i] = (double*)malloc(5 * sizeof(double));
    }
    // Split the line into individual transitions and process them
    char* token = strtok(line, ";");
    token = strtok(NULL, ";"); // Skip the initial state
    int transitionIndex = 0;
    while (token != NULL) {
        sscanf(token, "%d %d %d %lf %lf %lf %lf %lf",
                &coll.final_states[transitionIndex][0],
                &coll.final_states[transitionIndex][1],
                 &coll.final_states[transitionIndex][2],
                  &coll.coll_rates[transitionIndex][0],
                  &coll.coll_rates[transitionIndex][1],
                  &coll.coll_rates[transitionIndex][2],
                  &coll.coll_rates[transitionIndex][3],
                  &coll.coll_rates[transitionIndex][4]);
        transitionIndex++;

        if (transitionIndex == coll.num_transitions) {
            break; 
        }
        token = strtok(NULL, ";");
    }

    return coll;
}

//initialze an array of struct Rate from a file
static struct Coll* initializeCollArrayFromFile(const char* filePath,int mass) {
    FILE* file = fopen(filePath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening the file '%s' in 'initializeCollArrayFromFile'.\n",filePath);
        return NULL;
    }

    // Determine the number of Rate instances based on the number of lines in the file
    int numColl = countLinesInFile(filePath);

    // Create an array of struct Rate
    struct Coll* coll = (struct Coll*)malloc(numColl * sizeof(struct Rate));

    // Initialize each Rate instance from the file
    int i;
    for (i = 0; i < numColl; i++) {
        // Position the file pointer to the start of the current Rate data
        if (i > 0) {
            if (fscanf(file, "%*[\n]")) {
                fprintf(stderr, "Error reading newline character.\n");
                exit(1);
            }
        }
        coll[i] = initializeCollFromFile(file,mass);
    }

    // Close the file
    fclose(file);

    return coll;
}


static void freeCollArray(struct Coll* coll, int numColl) {
    if (coll == NULL) {
        return; // Nothing to free
    }
    int i,j;
    for (i = 0; i < numColl; i++) {
        // Free dynamically allocated memory for each Rate instance
        free(coll[i].initial_state);
        for (j = 0; j < coll[i].num_transitions; j++) {
            free(coll[i].final_states[j]);
            free(coll[i].coll_rates[j]);
        }
        free(coll[i].final_states);
        free(coll[i].coll_rates);
    }

    // Free memory for the array of Rate instances
    free(coll);
}

// Function to check if a state is in an array of Coll
static int CollStateIndex(struct Coll* coll, int numColl, int* state) {
    int i;                                           
    for (i = 0; i < numColl; i++) {
        if (coll[i].initial_state[0] == state[0] &&
            coll[i].initial_state[1] == state[1] &&
            coll[i].initial_state[2] == state[2]) {
            return i;  // return the index of matching state. Will always evaluate to true
        }
    }
    return -1.0;  // The state is not found in the array, false
}

/// Statistical weight

struct Weights {
int mass;
int num_states;
int** states;
double* weight;
                };

// Initalize statistical weights from a file
static struct Weights InitializeWeightsFromFile(const char* filePath, int mass) {
    struct Weights weights;
    weights.mass = mass;
    weights.num_states = countLinesInFile(filePath);

    weights.states = (int**)malloc(weights.num_states * sizeof(int*));
    weights.weight = (double*)malloc(weights.num_states * sizeof(double));

    FILE* file = fopen(filePath, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening the file '%s' in 'InitializeWeightsFromFile'.\n",filePath);
        exit(1);
    }
    int i;
    for (i = 0; i < weights.num_states; i++) {
        weights.states[i] = (int*)malloc(3 * sizeof(int));
        if (fscanf(file, "%d %d %d %lf", &weights.states[i][0], &weights.states[i][1], &weights.states[i][2], &weights.weight[i]) != 4) {
            fprintf(stderr, "Error reading data from the file.\n");
            exit(1);
        }
    }
    fclose(file);
    return weights;
}

// Function to free the memory of the struct Weights
static void freeWeights(struct Weights weights) {
    int i;
    for (i = 0; i < weights.num_states; i++) {
        free(weights.states[i]);
    }
    free(weights.states);
    free(weights.weight);
}

static double getWeightForState(struct Weights weights, int state[3]) {
    int i;
    for (i = 0; i < weights.num_states; i++) {
        if (weights.states[i][0] == state[0] &&
            weights.states[i][1] == state[1] &&
            weights.states[i][2] == state[2]) {
            return weights.weight[i];
        }
    }
    // If the state is not found, return a default value (you can change this as needed)
    return -1.0;
}





///// DICTIONARIES FOR ENERGIES


#define TABLE_SIZE 100

// Define a structure for key-value pairs
struct KeyValuePair {
    char* key;
    double value;
    struct KeyValuePair* next;
};

// Define a dictionary structure
struct Dictionary {
    struct KeyValuePair* table[TABLE_SIZE];
};

// Hash function to map keys to an index in the table
static int hash(char* key) {
    int hash = 0;
    int i;
    for (i = 0; key[i] != '\0'; i++) {
        hash = (hash + key[i]) % TABLE_SIZE;
    }
    return hash;
}

// Initialize a dictionary
static struct Dictionary* createDictionary() {
    struct Dictionary* dict = (struct Dictionary*)malloc(sizeof(struct Dictionary));
    int i;
    for (i = 0; i < TABLE_SIZE; i++) {
        dict->table[i] = NULL;
    }
    return dict;
}

// Insert a key-value pair into the dictionary
static void insert(struct Dictionary* dict, char* key, double value) {
    int index = hash(key);
    struct KeyValuePair* newPair = (struct KeyValuePair*)malloc(sizeof(struct KeyValuePair));
    newPair->key = strdup(key);
    newPair->value = value;
    newPair->next = dict->table[index];
    dict->table[index] = newPair;
}

// Get the value associated with a key
static double get(struct Dictionary* dict, char* key) {
    int index = hash(key);
    struct KeyValuePair* current = dict->table[index];
    while (current != NULL) {
        if (strcmp(current->key, key) == 0) {
            return current->value;
        }
        current = current->next;
    }
    // Key not found
    return -1.0;
}

static void printDictionary(struct Dictionary* dict) {
   int i;
    for (i = 0; i < TABLE_SIZE; i++) {
        struct KeyValuePair* current = dict->table[i];
        while (current != NULL) {
            // Process the key-value pair
            printf("%s -> %f\n",current->key, current->value );
            current = current->next;
        }
    }
}

// Function to read the file and populate the dictionary
static struct Dictionary* readFileAndCreateEnergyDictionary(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Failed to open the file '%s' in 'readFileAndCreateEnergyDictionary'.\n",filename);
        return NULL;
    }

    struct Dictionary* dict = createDictionary();
    double energy;
    int a,b,c;
    char key[20]; 

    while (fscanf(file, "%d %d %d %lf", &a, &b, &c, &energy) == 4) {
        snprintf(key, sizeof(key), "%d %d %d", a, b, c);
        insert(dict, key, energy);
    }

    fclose(file);
    return dict;
}


// Free the memory used by a dictionary
static void freeDictionary(struct Dictionary* dict) {
    int i;
    for (i = 0; i < TABLE_SIZE; i++) {
        struct KeyValuePair* current = dict->table[i];
        while (current != NULL) {
            struct KeyValuePair* temp = current;
            current = current->next;
            free(temp->key);
            free(temp);
        }
    }
    free(dict);
}




//////////////////////////////////////////////////////////////////////////////////////
////////////                                                            //////////////                                                          //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                 Here ends dictionary town                  //////////////
////////////                                                            //////////////
////////////                                                            //////////////
////////////                                                            //////////////                                                //////////////
//////////////////////////////////////////////////////////////////////////////////////




/* simulation conditions to transmit. Keep in mind that they are 
   transmitted to other nodes through an MPI_Reduce after
   casting them to a real (so the signals can be sent together with other 
   data). This means that the only meaningful values are positive, 
   negative or zero. */
enum { eglsNABNSB, eglsCHKPT, eglsSTOPCOND, eglsRESETCOUNTERS, eglsNR };
/* Is the signal in one simulation independent of other simulations? */
gmx_bool gs_simlocal[eglsNR] = { TRUE, FALSE, FALSE, TRUE };

typedef struct {
    int nstms;       /* The frequency for intersimulation communication */
    int sig[eglsNR]; /* The signal set by one process in do_md */
    int set[eglsNR]; /* The communicated signal, equal for all processes */
} globsig_t;


//// EI ////

////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_Ei.c                                            //
// Routine(s):                                                                //
//    Exponential_Integral_Ei                                                 //
//    xExponential_Integral_Ei                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for fabsl(), expl() and logl()        
#include <float.h>          // required for LDBL_EPSILON, DBL_MAX

//                         Internally Defined Routines                        //
double      Exponential_Integral_Ei( double x );
long double xExponential_Integral_Ei( long double x );

static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);


//                         Internally Defined Constants                       //
static const long double epsilon = 10.0 * LDBL_EPSILON;

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Ei( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Exponential_Integral_Ei( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Exponential_Integral_Ei( double x )
{
   return (double) xExponential_Integral_Ei( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
   if ( x < -5.0L ) return Continued_Fraction_Ei(x);
   if ( x == 0.0L ) return -DBL_MAX;
   if ( x < 6.8L )  return Power_Series_Ei(x);
   if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
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

static long double Continued_Fraction_Ei( long double x )
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
   while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {
      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0L;
      } else {
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

static long double Power_Series_Ei( long double x )
{ 
   long double xn = -x;
   long double Sn = -x;
   long double Sm1 = 0.0L;
   long double hsum = 1.0L;
   long double g = 0.5772156649015328606065121L;
   long double y = 1.0L;
   long double factorial = 1.0L;
  
   if ( x == 0.0L ) return (long double) -DBL_MAX;
 
   while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
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
      1.915047433355013959531e2L,  4.403798995348382689974e2L,
      1.037878290717089587658e3L,  2.492228976241877759138e3L,
      6.071406374098611507965e3L,  1.495953266639752885229e4L,
      3.719768849068903560439e4L,  9.319251363396537129882e4L,
      2.349558524907683035782e5L,  5.955609986708370018502e5L,
      1.516637894042516884433e6L,  3.877904330597443502996e6L,
      9.950907251046844760026e6L,  2.561565266405658882048e7L,
      6.612718635548492136250e7L,  1.711446713003636684975e8L,
      4.439663698302712208698e8L,  1.154115391849182948287e9L,
      3.005950906525548689841e9L,  7.842940991898186370453e9L,
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
      3.975442747903744836007e19L, 1.058563689713169096306e20L
   };
   int  k = (int) (x + 0.5);
   int  j = 0;
   long double xx = (long double) k;
   long double dx = x - xx;
   long double xxj = xx;
   long double edx = expl(dx);
   long double Sm = 1.0L;
   long double Sn = (edx - 1.0L) / xxj;
   long double term = DBL_MAX;
   long double factorial = 1.0L;
   long double dxj = 1.0L;

   while (fabsl(term) > epsilon * fabsl(Sn) ) {
      j++;
      factorial *= (long double) j;
      xxj *= xx;
      dxj *= (-dx);
      Sm += (dxj / factorial);
      term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
      Sn += term;
   }
   
   return ei[k-7] + Sn * expl(xx); 
}



/* check which of the multisim simulations has the shortest number of
   steps and return that number of nsteps */
static gmx_large_int_t get_multisim_nsteps(const t_commrec *cr,
                                           gmx_large_int_t nsteps)
{
    gmx_large_int_t steps_out;

    if MASTER(cr)
    {
        gmx_large_int_t *buf;
        int s;

        snew(buf,cr->ms->nsim);

        buf[cr->ms->sim] = nsteps;
        gmx_sumli_sim(cr->ms->nsim, buf, cr->ms);

        steps_out=-1;
        for(s=0; s<cr->ms->nsim; s++)
        {
            /* find the smallest positive number */
            if (buf[s]>= 0 && ((steps_out < 0) || (buf[s]<steps_out)) )
            {
                steps_out=buf[s];
            }
        }
        sfree(buf);

        /* if we're the limiting simulation, don't do anything */
        if (steps_out>=0 && steps_out<nsteps) 
        {
            char strbuf[255];
            snprintf(strbuf, 255, "Will stop simulation %%d after %s steps (another simulation will end then).\n", gmx_large_int_pfmt);
            fprintf(stderr, strbuf, cr->ms->sim, steps_out);
        }
    }
    /* broadcast to non-masters */
    gmx_bcast(sizeof(gmx_large_int_t), &steps_out, cr);
    return steps_out;
}

static int multisim_min(const gmx_multisim_t *ms,int nmin,int n)
{
    int  *buf;
    gmx_bool bPos,bEqual;
    int  s,d;

    snew(buf,ms->nsim);
    buf[ms->sim] = n;
    gmx_sumi_sim(ms->nsim,buf,ms);
    bPos   = TRUE;
    bEqual = TRUE;
    for(s=0; s<ms->nsim; s++)
    {
        bPos   = bPos   && (buf[s] > 0);
        bEqual = bEqual && (buf[s] == buf[0]);
    }
    if (bPos)
    {
        if (bEqual)
        {
            nmin = min(nmin,buf[0]);
        }
        else
        {
            /* Find the least common multiple */
            for(d=2; d<nmin; d++)
            {
                s = 0;
                while (s < ms->nsim && d % buf[s] == 0)
                {
                    s++;
                }
                if (s == ms->nsim)
                {
                    /* We found the LCM and it is less than nmin */
                    nmin = d;
                    break;
                }
            }
        }
    }
    sfree(buf);

    return nmin;
}

static int multisim_nstsimsync(const t_commrec *cr,
                               const t_inputrec *ir,int repl_ex_nst)
{
    int nmin;

    if (MASTER(cr))
    {
        nmin = INT_MAX;
        nmin = multisim_min(cr->ms,nmin,ir->nstlist);
        nmin = multisim_min(cr->ms,nmin,ir->nstcalcenergy);
        nmin = multisim_min(cr->ms,nmin,repl_ex_nst);
        if (nmin == INT_MAX)
        {
            gmx_fatal(FARGS,"Can not find an appropriate interval for inter-simulation communication, since nstlist, nstcalcenergy and -replex are all <= 0");
        }
        /* Avoid inter-simulation communication at every (second) step */
        if (nmin <= 2)
        {
            nmin = 10;
        }
    }

    gmx_bcast(sizeof(int),&nmin,cr);

    return nmin;
}

static void init_global_signals(globsig_t *gs,const t_commrec *cr,
                                const t_inputrec *ir,int repl_ex_nst)
{
    int i;

    if (MULTISIM(cr))
    {
        gs->nstms = multisim_nstsimsync(cr,ir,repl_ex_nst);
        if (debug)
        {
            fprintf(debug,"Syncing simulations for checkpointing and termination every %d steps\n",gs->nstms);
        }
    }
    else
    {
        gs->nstms = 1;
    }

    for(i=0; i<eglsNR; i++)
    {
        gs->sig[i] = 0;
        gs->set[i] = 0;
    }
}

static void copy_coupling_state(t_state *statea,t_state *stateb, 
                                gmx_ekindata_t *ekinda,gmx_ekindata_t *ekindb, t_grpopts* opts) 
{
    
    /* MRS note -- might be able to get rid of some of the arguments.  Look over it when it's all debugged */
    
    int i,j,nc;

    /* Make sure we have enough space for x and v */
    if (statea->nalloc > stateb->nalloc)
    {
        stateb->nalloc = statea->nalloc;
        srenew(stateb->x,stateb->nalloc);
        srenew(stateb->v,stateb->nalloc);
    }

    stateb->natoms     = statea->natoms;
    stateb->ngtc       = statea->ngtc;
    stateb->nnhpres    = statea->nnhpres;
    stateb->veta       = statea->veta;
    if (ekinda) 
    {
        copy_mat(ekinda->ekin,ekindb->ekin);
        for (i=0; i<stateb->ngtc; i++) 
        {
            ekindb->tcstat[i].T = ekinda->tcstat[i].T;
            ekindb->tcstat[i].Th = ekinda->tcstat[i].Th;
            copy_mat(ekinda->tcstat[i].ekinh,ekindb->tcstat[i].ekinh);
            copy_mat(ekinda->tcstat[i].ekinf,ekindb->tcstat[i].ekinf);
            ekindb->tcstat[i].ekinscalef_nhc =  ekinda->tcstat[i].ekinscalef_nhc;
            ekindb->tcstat[i].ekinscaleh_nhc =  ekinda->tcstat[i].ekinscaleh_nhc;
            ekindb->tcstat[i].vscale_nhc =  ekinda->tcstat[i].vscale_nhc;
        }
    }
    copy_rvecn(statea->x,stateb->x,0,stateb->natoms);
    copy_rvecn(statea->v,stateb->v,0,stateb->natoms);
    copy_mat(statea->box,stateb->box);
    copy_mat(statea->box_rel,stateb->box_rel);
    copy_mat(statea->boxv,stateb->boxv);

    for (i = 0; i<stateb->ngtc; i++) 
    { 
        nc = i*opts->nhchainlength;
        for (j=0; j<opts->nhchainlength; j++) 
        {
            stateb->nosehoover_xi[nc+j]  = statea->nosehoover_xi[nc+j];
            stateb->nosehoover_vxi[nc+j] = statea->nosehoover_vxi[nc+j];
        }
    }
    if (stateb->nhpres_xi != NULL)
    {
        for (i = 0; i<stateb->nnhpres; i++) 
        {
            nc = i*opts->nhchainlength;
            for (j=0; j<opts->nhchainlength; j++) 
            {
                stateb->nhpres_xi[nc+j]  = statea->nhpres_xi[nc+j];
                stateb->nhpres_vxi[nc+j] = statea->nhpres_vxi[nc+j];
            }
        }
    }
}

static real compute_conserved_from_auxiliary(t_inputrec *ir, t_state *state, t_extmass *MassQ)
{
    real quantity = 0;
    switch (ir->etc) 
    {
    case etcNO:
        break;
    case etcBERENDSEN:
        break;
    case etcNOSEHOOVER:
        quantity = NPT_energy(ir,state,MassQ);                
        break;
    case etcVRESCALE:
        quantity = vrescale_energy(&(ir->opts),state->therm_integral);
        break;
    default:
        break;
    }
    return quantity;
}

static void compute_globals(FILE *fplog, gmx_global_stat_t gstat, t_commrec *cr, t_inputrec *ir, 
                            t_forcerec *fr, gmx_ekindata_t *ekind, 
                            t_state *state, t_state *state_global, t_mdatoms *mdatoms, 
                            t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                            gmx_enerdata_t *enerd,tensor force_vir, tensor shake_vir, tensor total_vir, 
                            tensor pres, rvec mu_tot, gmx_constr_t constr, 
                            globsig_t *gs,gmx_bool bInterSimGS,
                            matrix box, gmx_mtop_t *top_global, real *pcurr, 
                            int natoms, gmx_bool *bSumEkinhOld, int flags)
{
    int  i,gsi;
    real gs_buf[eglsNR];
    tensor corr_vir,corr_pres,shakeall_vir;
    gmx_bool bEner,bPres,bTemp, bVV;
    gmx_bool bRerunMD, bStopCM, bGStat, bIterate, 
        bFirstIterate,bReadEkin,bEkinAveVel,bScaleEkin, bConstrain;
    real ekin,temp,prescorr,enercorr,dvdlcorr;
    
    /* translate CGLO flags to gmx_booleans */
    bRerunMD = flags & CGLO_RERUNMD;
    bStopCM = flags & CGLO_STOPCM;
    bGStat = flags & CGLO_GSTAT;

    bReadEkin = (flags & CGLO_READEKIN);
    bScaleEkin = (flags & CGLO_SCALEEKIN);
    bEner = flags & CGLO_ENERGY;
    bTemp = flags & CGLO_TEMPERATURE;
    bPres  = (flags & CGLO_PRESSURE);
    bConstrain = (flags & CGLO_CONSTRAINT);
    bIterate = (flags & CGLO_ITERATE);
    bFirstIterate = (flags & CGLO_FIRSTITERATE);

    /* we calculate a full state kinetic energy either with full-step velocity verlet
       or half step where we need the pressure */
    
    bEkinAveVel = (ir->eI==eiVV || (ir->eI==eiVVAK && bPres) || bReadEkin);
    
    /* in initalization, it sums the shake virial in vv, and to 
       sums ekinh_old in leapfrog (or if we are calculating ekinh_old) for other reasons */

    /* ########## Kinetic energy  ############## */
    
    if (bTemp) 
    {
        /* Non-equilibrium MD: this is parallellized, but only does communication
         * when there really is NEMD.
         */
        
        if (PAR(cr) && (ekind->bNEMD)) 
        {
            accumulate_u(cr,&(ir->opts),ekind);
        }
        debug_gmx();
        if (bReadEkin)
        {
            restore_ekinstate_from_state(cr,ekind,&state_global->ekinstate);
        }
        else 
        {

            calc_ke_part(state,&(ir->opts),mdatoms,ekind,nrnb,bEkinAveVel,bIterate);
        }
        
        debug_gmx();
        
        /* Calculate center of mass velocity if necessary, also parallellized */
        if (bStopCM && !bRerunMD && bEner) 
        {
            calc_vcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms,
                         state->x,state->v,vcm);
        }
    }

    if (bTemp || bPres || bEner || bConstrain) 
    {
        if (!bGStat)
        {
            /* We will not sum ekinh_old,                                                            
             * so signal that we still have to do it.                                                
             */
            *bSumEkinhOld = TRUE;

        }
        else
        {
            if (gs != NULL)
            {
                for(i=0; i<eglsNR; i++)
                {
                    gs_buf[i] = gs->sig[i];
                }
            }
            if (PAR(cr)) 
            {
                wallcycle_start(wcycle,ewcMoveE);
                GMX_MPE_LOG(ev_global_stat_start);
                global_stat(fplog,gstat,cr,enerd,force_vir,shake_vir,mu_tot,
                            ir,ekind,constr,vcm,
                            gs != NULL ? eglsNR : 0,gs_buf,
                            top_global,state,
                            *bSumEkinhOld,flags);
                GMX_MPE_LOG(ev_global_stat_finish);
                wallcycle_stop(wcycle,ewcMoveE);
            }
            if (gs != NULL)
            {
                if (MULTISIM(cr) && bInterSimGS)
                {
                    if (MASTER(cr))
                    {
                        /* Communicate the signals between the simulations */
                        gmx_sum_sim(eglsNR,gs_buf,cr->ms);
                    }
                    /* Communicate the signals form the master to the others */
                    gmx_bcast(eglsNR*sizeof(gs_buf[0]),gs_buf,cr);
                }
                for(i=0; i<eglsNR; i++)
                {
                    if (bInterSimGS || gs_simlocal[i])
                    {
                        /* Set the communicated signal only when it is non-zero,
                         * since signals might not be processed at each MD step.
                         */
                        gsi = (gs_buf[i] >= 0 ?
                               (int)(gs_buf[i] + 0.5) :
                               (int)(gs_buf[i] - 0.5));
                        if (gsi != 0)
                        {
                            gs->set[i] = gsi;
                        }
                        /* Turn off the local signal */
                        gs->sig[i] = 0;
                    }
                }
            }
            *bSumEkinhOld = FALSE;
        }
    }
    
    if (!ekind->bNEMD && debug && bTemp && (vcm->nr > 0))
    {
        correct_ekin(debug,
                     mdatoms->start,mdatoms->start+mdatoms->homenr,
                     state->v,vcm->group_p[0],
                     mdatoms->massT,mdatoms->tmass,ekind->ekin);
    }
    
    if (bEner) {
        /* Do center of mass motion removal */
        if (bStopCM && !bRerunMD) /* is this correct?  Does it get called too often with this logic? */
        {
            check_cm_grp(fplog,vcm,ir,1);
            do_stopcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms->cVCM,
                          state->x,state->v,vcm);
            inc_nrnb(nrnb,eNR_STOPCM,mdatoms->homenr);
        }

        /* Calculate the amplitude of the cosine velocity profile */
        ekind->cosacc.vcos = ekind->cosacc.mvcos/mdatoms->tmass;
    }

    if (bTemp) 
    {
        /* Sum the kinetic energies of the groups & calc temp */
        /* compute full step kinetic energies if vv, or if vv-avek and we are computing the pressure with IR_NPT_TROTTER */
        /* three maincase:  VV with AveVel (md-vv), vv with AveEkin (md-vv-avek), leap with AveEkin (md).  
           Leap with AveVel is not supported; it's not clear that it will actually work.  
           bEkinAveVel: If TRUE, we simply multiply ekin by ekinscale to get a full step kinetic energy. 
           If FALSE, we average ekinh_old and ekinh*ekinscale_nhc to get an averaged half step kinetic energy.
           bSaveEkinOld: If TRUE (in the case of iteration = bIterate is TRUE), we don't reset the ekinscale_nhc.  
           If FALSE, we go ahead and erase over it.
        */ 
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,&(enerd->term[F_DKDL]),
                                       bEkinAveVel,bIterate,bScaleEkin);
 
        enerd->term[F_EKIN] = trace(ekind->ekin);
    }
    
    /* ##########  Long range energy information ###### */
    
    if (bEner || bPres || bConstrain) 
    {
        calc_dispcorr(fplog,ir,fr,0,top_global->natoms,box,state->lambda,
                      corr_pres,corr_vir,&prescorr,&enercorr,&dvdlcorr);
    }
    
    if (bEner && bFirstIterate) 
    {
        enerd->term[F_DISPCORR] = enercorr;
        enerd->term[F_EPOT] += enercorr;
        enerd->term[F_DVDL] += dvdlcorr;
        if (fr->efep != efepNO) {
            enerd->dvdl_lin += dvdlcorr;
        }
    }
    
    /* ########## Now pressure ############## */
    if (bPres || bConstrain) 
    {
        
        m_add(force_vir,shake_vir,total_vir);
        
        /* Calculate pressure and apply LR correction if PPPM is used.
         * Use the box from last timestep since we already called update().
         */
        
        enerd->term[F_PRES] = calc_pres(fr->ePBC,ir->nwall,box,ekind->ekin,total_vir,pres,
                                        (fr->eeltype==eelPPPM)?enerd->term[F_COUL_RECIP]:0.0);
        
        /* Calculate long range corrections to pressure and energy */
        /* this adds to enerd->term[F_PRES] and enerd->term[F_ETOT], 
           and computes enerd->term[F_DISPCORR].  Also modifies the 
           total_vir and pres tesors */
        
        m_add(total_vir,corr_vir,total_vir);
        m_add(pres,corr_pres,pres);
        enerd->term[F_PDISPCORR] = prescorr;
        enerd->term[F_PRES] += prescorr;
        *pcurr = enerd->term[F_PRES];
        /* calculate temperature using virial */
        enerd->term[F_VTEMP] = calc_temp(trace(total_vir),ir->opts.nrdf[0]);
        
    }    
}


/* Definitions for convergence of iterated constraints */

/* iterate constraints up to 50 times  */
#define MAXITERCONST       50

/* data type */
typedef struct
{
    real f,fprev,x,xprev;  
    int iter_i;
    gmx_bool bIterate;
    real allrelerr[MAXITERCONST+2];
    int num_close; /* number of "close" violations, caused by limited precision. */
} gmx_iterate_t;
  
#ifdef GMX_DOUBLE
#define CONVERGEITER  0.000000001
#define CLOSE_ENOUGH  0.000001000
#else
#define CONVERGEITER  0.0001
#define CLOSE_ENOUGH  0.0050
#endif

/* we want to keep track of the close calls.  If there are too many, there might be some other issues.
   so we make sure that it's either less than some predetermined number, or if more than that number,
   only some small fraction of the total. */
#define MAX_NUMBER_CLOSE        50
#define FRACTION_CLOSE       0.001
  
/* maximum length of cyclic traps to check, emerging from limited numerical precision  */
#define CYCLEMAX            20

static void gmx_iterate_init(gmx_iterate_t *iterate,gmx_bool bIterate)
{
    int i;

    iterate->iter_i = 0;
    iterate->bIterate = bIterate;
    iterate->num_close = 0;
    for (i=0;i<MAXITERCONST+2;i++) 
    {
        iterate->allrelerr[i] = 0;
    }
}

static gmx_bool done_iterating(const t_commrec *cr,FILE *fplog, int nsteps, gmx_iterate_t *iterate, gmx_bool bFirstIterate, real fom, real *newf) 
{    
    /* monitor convergence, and use a secant search to propose new
       values.  
                                                                  x_{i} - x_{i-1}
       The secant method computes x_{i+1} = x_{i} - f(x_{i}) * ---------------------
                                                                f(x_{i}) - f(x_{i-1})
       
       The function we are trying to zero is fom-x, where fom is the
       "figure of merit" which is the pressure (or the veta value) we
       would get by putting in an old value of the pressure or veta into
       the incrementor function for the step or half step.  I have
       verified that this gives the same answer as self consistent
       iteration, usually in many fewer steps, especially for small tau_p.
       
       We could possibly eliminate an iteration with proper use
       of the value from the previous step, but that would take a bit
       more bookkeeping, especially for veta, since tests indicate the
       function of veta on the last step is not sufficiently close to
       guarantee convergence this step. This is
       good enough for now.  On my tests, I could use tau_p down to
       0.02, which is smaller that would ever be necessary in
       practice. Generally, 3-5 iterations will be sufficient */

    real relerr,err,xmin;
    char buf[256];
    int i;
    gmx_bool incycle;
    
    if (bFirstIterate) 
    {
        iterate->x = fom;
        iterate->f = fom-iterate->x;
        iterate->xprev = 0;
        iterate->fprev = 0;
        *newf = fom;
    } 
    else 
    {
        iterate->f = fom-iterate->x; /* we want to zero this difference */
        if ((iterate->iter_i > 1) && (iterate->iter_i < MAXITERCONST)) 
        {
            if (iterate->f==iterate->fprev) 
            {
                *newf = iterate->f;
            } 
            else 
            {
                *newf = iterate->x - (iterate->x-iterate->xprev)*(iterate->f)/(iterate->f-iterate->fprev); 
            }
        } 
        else 
        {
            /* just use self-consistent iteration the first step to initialize, or 
               if it's not converging (which happens occasionally -- need to investigate why) */
            *newf = fom; 
        }
    }
    /* Consider a slight shortcut allowing us to exit one sooner -- we check the
       difference between the closest of x and xprev to the new
       value. To be 100% certain, we should check the difference between
       the last result, and the previous result, or
       
       relerr = (fabs((x-xprev)/fom));
       
       but this is pretty much never necessary under typical conditions.
       Checking numerically, it seems to lead to almost exactly the same
       trajectories, but there are small differences out a few decimal
       places in the pressure, and eventually in the v_eta, but it could
       save an interation.
       
       if (fabs(*newf-x) < fabs(*newf - xprev)) { xmin = x;} else { xmin = xprev;}
       relerr = (fabs((*newf-xmin) / *newf));
    */
    
    err = fabs((iterate->f-iterate->fprev));
    relerr = fabs(err/fom);

    iterate->allrelerr[iterate->iter_i] = relerr;
    
    if (iterate->iter_i > 0) 
    {
        if (debug) 
        {
            fprintf(debug,"Iterating NPT constraints: %6i %20.12f%14.6g%20.12f\n",
                    iterate->iter_i,fom,relerr,*newf);
        }
        
        if ((relerr < CONVERGEITER) || (err < CONVERGEITER) || (fom==0) || ((iterate->x == iterate->xprev) && iterate->iter_i > 1))
        {
            iterate->bIterate = FALSE;
            if (debug) 
            {
                fprintf(debug,"Iterating NPT constraints: CONVERGED\n");
            }
            return TRUE;
        }
        if (iterate->iter_i > MAXITERCONST)
        {
            if (relerr < CLOSE_ENOUGH)
            {
                incycle = FALSE;
                for (i=1;i<CYCLEMAX;i++) {
                    if ((iterate->allrelerr[iterate->iter_i-(1+i)] == iterate->allrelerr[iterate->iter_i-1]) &&
                        (iterate->allrelerr[iterate->iter_i-(1+i)] == iterate->allrelerr[iterate->iter_i-(1+2*i)])) {
                        incycle = TRUE;
                        if (debug) 
                        {
                            fprintf(debug,"Exiting from an NPT iterating cycle of length %d\n",i);
                        }
                        break;
                    }
                }
                
                if (incycle) {
                    /* step 1: trapped in a numerical attractor */
                    /* we are trapped in a numerical attractor, and can't converge any more, and are close to the final result.
                       Better to give up convergence here than have the simulation die.
                    */
                    iterate->num_close++;
                    return TRUE;
                } 
                else 
                {
                    /* Step #2: test if we are reasonably close for other reasons, then monitor the number.  If not, die */
                    
                    /* how many close calls have we had?  If less than a few, we're OK */
                    if (iterate->num_close < MAX_NUMBER_CLOSE) 
                    {
                        sprintf(buf,"Slight numerical convergence deviation with NPT at step %d, relative error only %10.5g, likely not a problem, continuing\n",nsteps,relerr);
                        md_print_warning(cr,fplog,buf);
                        iterate->num_close++;
                        return TRUE;
                        /* if more than a few, check the total fraction.  If too high, die. */
                    } else if (iterate->num_close/(double)nsteps > FRACTION_CLOSE) {
                        gmx_fatal(FARGS,"Could not converge NPT constraints, too many exceptions (%d%%\n",iterate->num_close/(double)nsteps);
                    } 
                }
            }
            else 
            {
                gmx_fatal(FARGS,"Could not converge NPT constraints\n");
            }
        }
    }
    
    iterate->xprev = iterate->x;
    iterate->x = *newf;
    iterate->fprev = iterate->f;
    iterate->iter_i++;
    
    return FALSE;
}

static void check_nst_param(FILE *fplog,t_commrec *cr,
                            const char *desc_nst,int nst,
                            const char *desc_p,int *p)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        sprintf(buf,"NOTE: %s changes %s to %d\n",desc_nst,desc_p,*p);
        md_print_warning(cr,fplog,buf);
    }
}

static void reset_all_counters(FILE *fplog,t_commrec *cr,
                               gmx_large_int_t step,
                               gmx_large_int_t *step_rel,t_inputrec *ir,
                               gmx_wallcycle_t wcycle,t_nrnb *nrnb,
                               gmx_runtime_t *runtime)
{
    char buf[STRLEN],sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    sprintf(buf,"Step %s: resetting all time and cycle counters\n",
            gmx_step_str(step,sbuf));
    md_print_warning(cr,fplog,buf);

    wallcycle_stop(wcycle,ewcRUN);
    wallcycle_reset_all(wcycle);
    if (DOMAINDECOMP(cr))
    {
        reset_dd_statistics_counters(cr->dd);
    }
    init_nrnb(nrnb);
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel = 0;
    wallcycle_start(wcycle,ewcRUN);
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Restarted time",runtime);
}

static void min_zero(int *n,int i)
{
    if (i > 0 && (*n == 0 || i < *n))
    {
        *n = i;
    }
}

static int lcd4(int i1,int i2,int i3,int i4)
{
    int nst;

    nst = 0;
    min_zero(&nst,i1);
    min_zero(&nst,i2);
    min_zero(&nst,i3);
    min_zero(&nst,i4);
    if (nst == 0)
    {
        gmx_incons("All 4 inputs for determininig nstglobalcomm are <= 0");
    }
    
    while (nst > 1 && ((i1 > 0 && i1 % nst != 0)  ||
                       (i2 > 0 && i2 % nst != 0)  ||
                       (i3 > 0 && i3 % nst != 0)  ||
                       (i4 > 0 && i4 % nst != 0)))
    {
        nst--;
    }

    return nst;
}

static int check_nstglobalcomm(FILE *fplog,t_commrec *cr,
                               int nstglobalcomm,t_inputrec *ir)
{
    char buf[STRLEN];

    if (!EI_DYNAMICS(ir->eI))
    {
        nstglobalcomm = 1;
    }

    if (nstglobalcomm == -1)
    {
        if (!(ir->nstcalcenergy > 0 ||
              ir->nstlist > 0 ||
              ir->etc != etcNO ||
              ir->epc != epcNO))
        {
            nstglobalcomm = 10;
            if (ir->nstenergy > 0 && ir->nstenergy < nstglobalcomm)
            {
                nstglobalcomm = ir->nstenergy;
            }
        }
        else
        {
            /* Ensure that we do timely global communication for
             * (possibly) each of the four following options.
             */
            nstglobalcomm = lcd4(ir->nstcalcenergy,
                                 ir->nstlist,
                                 ir->etc != etcNO ? ir->nsttcouple : 0,
                                 ir->epc != epcNO ? ir->nstpcouple : 0);
        }
    }
    else
    {
        if (ir->nstlist > 0 &&
            nstglobalcomm > ir->nstlist && nstglobalcomm % ir->nstlist != 0)
        {
            nstglobalcomm = (nstglobalcomm / ir->nstlist)*ir->nstlist;
            sprintf(buf,"WARNING: nstglobalcomm is larger than nstlist, but not a multiple, setting it to %d\n",nstglobalcomm);
            md_print_warning(cr,fplog,buf);
        }
        if (ir->nstcalcenergy > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstcalcenergy",&ir->nstcalcenergy);
        }
        if (ir->etc != etcNO && ir->nsttcouple > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nsttcouple",&ir->nsttcouple);
        }
        if (ir->epc != epcNO && ir->nstpcouple > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstpcouple",&ir->nstpcouple);
        }

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstenergy",&ir->nstenergy);

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstlog",&ir->nstlog);
    }

    if (ir->comm_mode != ecmNO && ir->nstcomm < nstglobalcomm)
    {
        sprintf(buf,"WARNING: Changing nstcomm from %d to %d\n",
                ir->nstcomm,nstglobalcomm);
        md_print_warning(cr,fplog,buf);
        ir->nstcomm = nstglobalcomm;
    }

    return nstglobalcomm;
}

void check_ir_old_tpx_versions(t_commrec *cr,FILE *fplog,
                               t_inputrec *ir,gmx_mtop_t *mtop)
{
    /* Check required for old tpx files */
    if (IR_TWINRANGE(*ir) && ir->nstlist > 1 &&
        ir->nstcalcenergy % ir->nstlist != 0)
    {
        md_print_warning(cr,fplog,"Old tpr file with twin-range settings: modifying energy calculation and/or T/P-coupling frequencies");

        if (gmx_mtop_ftype_count(mtop,F_CONSTR) +
            gmx_mtop_ftype_count(mtop,F_CONSTRNC) > 0 &&
            ir->eConstrAlg == econtSHAKE)
        {
            md_print_warning(cr,fplog,"With twin-range cut-off's and SHAKE the virial and pressure are incorrect");
            if (ir->epc != epcNO)
            {
                gmx_fatal(FARGS,"Can not do pressure coupling with twin-range cut-off's and SHAKE");
            }
        }
        check_nst_param(fplog,cr,"nstlist",ir->nstlist,
                        "nstcalcenergy",&ir->nstcalcenergy);
        if (ir->epc != epcNO)
        {
            check_nst_param(fplog,cr,"nstlist",ir->nstlist,
                            "nstpcouple",&ir->nstpcouple);
        }
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "nstenergy",&ir->nstenergy);
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "nstlog",&ir->nstlog);
        if (ir->efep != efepNO)
        {
            check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                            "nstdhdl",&ir->nstdhdl);
        }
    }
}

typedef struct {
    gmx_bool       bGStatEveryStep;
    gmx_large_int_t step_ns;
    gmx_large_int_t step_nscheck;
    gmx_large_int_t nns;
    matrix     scale_tot;
    int        nabnsb;
    double     s1;
    double     s2;
    double     ab;
    double     lt_runav;
    double     lt_runav2;
} gmx_nlheur_t;

static void reset_nlistheuristics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    nlh->lt_runav  = 0;
    nlh->lt_runav2 = 0;
    nlh->step_nscheck = step;
}

static void init_nlistheuristics(gmx_nlheur_t *nlh,
                                 gmx_bool bGStatEveryStep,gmx_large_int_t step)
{
    nlh->bGStatEveryStep = bGStatEveryStep;
    nlh->nns       = 0;
    nlh->nabnsb    = 0;
    nlh->s1        = 0;
    nlh->s2        = 0;
    nlh->ab        = 0;

    reset_nlistheuristics(nlh,step);
}

static void update_nliststatistics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    gmx_large_int_t nl_lt;
    char sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];

    /* Determine the neighbor list life time */
    nl_lt = step - nlh->step_ns;
    if (debug)
    {
        fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %s steps\n",nlh->nabnsb,gmx_step_str(nl_lt,sbuf));
    }
    nlh->nns++;
    nlh->s1 += nl_lt;
    nlh->s2 += nl_lt*nl_lt;
    nlh->ab += nlh->nabnsb;
    if (nlh->lt_runav == 0)
    {
        nlh->lt_runav  = nl_lt;
        /* Initialize the fluctuation average
         * such that at startup we check after 0 steps.
         */
        nlh->lt_runav2 = sqr(nl_lt/2.0);
    }
    /* Running average with 0.9 gives an exp. history of 9.5 */
    nlh->lt_runav2 = 0.9*nlh->lt_runav2 + 0.1*sqr(nlh->lt_runav - nl_lt);
    nlh->lt_runav  = 0.9*nlh->lt_runav  + 0.1*nl_lt;
    if (nlh->bGStatEveryStep)
    {
        /* Always check the nlist validity */
        nlh->step_nscheck = step;
    }
    else
    {
        /* We check after:  <life time> - 2*sigma
         * The factor 2 is quite conservative,
         * but we assume that with nstlist=-1 the user
         * prefers exact integration over performance.
         */
        nlh->step_nscheck = step
                  + (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)) - 1;
    }
    if (debug)
    {
        fprintf(debug,"nlist life time %s run av. %4.1f sig %3.1f check %s check with -gcom %d\n",
                gmx_step_str(nl_lt,sbuf),nlh->lt_runav,sqrt(nlh->lt_runav2),
                gmx_step_str(nlh->step_nscheck-step+1,sbuf2),
                (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)));
    }
}

static void set_nlistheuristics(gmx_nlheur_t *nlh,gmx_bool bReset,gmx_large_int_t step)
{
    int d;

    if (bReset)
    {
        reset_nlistheuristics(nlh,step);
    }
    else
    {
        update_nliststatistics(nlh,step);
    }

    nlh->step_ns = step;
    /* Initialize the cumulative coordinate scaling matrix */
    clear_mat(nlh->scale_tot);
    for(d=0; d<DIM; d++)
    {
        nlh->scale_tot[d][d] = 1.0;
    }
}

static void rerun_parallel_comm(t_commrec *cr,t_trxframe *fr,
                                gmx_bool *bNotLastFrame)
{
    gmx_bool bAlloc;
    rvec *xp,*vp;

    bAlloc = (fr->natoms == 0);

    if (MASTER(cr) && !*bNotLastFrame)
    {
        fr->natoms = -1;
    }
    xp = fr->x;
    vp = fr->v;
    gmx_bcast(sizeof(*fr),fr,cr);
    fr->x = xp;
    fr->v = vp;

    *bNotLastFrame = (fr->natoms >= 0);

    if (*bNotLastFrame && PARTDECOMP(cr))
    {
        /* x and v are the only variable size quantities stored in trr
         * that are required for rerun (f is not needed).
         */
        if (bAlloc)
        {
            snew(fr->x,fr->natoms);
            snew(fr->v,fr->natoms);
        }
        if (fr->bX)
        {
            gmx_bcast(fr->natoms*sizeof(fr->x[0]),fr->x[0],cr);
        }
        if (fr->bV)
        {
            gmx_bcast(fr->natoms*sizeof(fr->v[0]),fr->v[0],cr);
        }
    }
}

double do_md(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose,gmx_bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite,gmx_constr_t constr,
             int stepout,t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb,gmx_wallcycle_t wcycle,
             gmx_edsam_t ed,t_forcerec *fr,
             int repl_ex_nst,int repl_ex_seed,
             real cpt_period,real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime)
{
    gmx_mdoutf_t *outf;
    gmx_large_int_t step,step_rel;
    double     run_time;
    double     t,t0,lam0;
    gmx_bool       bGStatEveryStep,bGStat,bNstEner,bCalcEnerPres;
    gmx_bool       bNS,bNStList,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
               bFirstStep,bStateFromTPX,bInitStep,bLastStep,
               bBornRadii,bStartingFromCpt;
    gmx_bool       bDoDHDL=FALSE;
    gmx_bool       do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
               bForceUpdate=FALSE,bCPT;
    int        mdof_flags;
    gmx_bool       bMasterState;
    int        force_flags,cglo_flags;
    tensor     force_vir,shake_vir,total_vir,tmp_vir,pres;
    int        i,m;
    t_trxstatus *status;
    rvec       mu_tot;
    t_vcm      *vcm;
    t_state    *bufstate=NULL;   
    matrix     *scale_tot,pcoupl_mu,M,ebox;
    gmx_nlheur_t nlh;
    t_trxframe rerun_fr;
    gmx_repl_ex_t repl_ex=NULL;
    int        nchkpt=1;

    gmx_localtop_t *top;    
    t_mdebin *mdebin=NULL;
    t_state    *state=NULL;
    rvec       *f_global=NULL;
    int        n_xtc=-1;
    rvec       *x_xtc=NULL;
    gmx_enerdata_t *enerd;
    rvec       *f=NULL;
    gmx_global_stat_t gstat;
    gmx_update_t upd=NULL;
    t_graph    *graph=NULL;
    globsig_t   gs;

    gmx_bool        bFFscan;
    gmx_groups_t *groups;
    gmx_ekindata_t *ekind, *ekind_save;
    gmx_shellfc_t shellfc;
    int         count,nconverged=0;
    real        timestep=0;
    double      tcount=0;
    gmx_bool        bIonize=FALSE;
    gmx_bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
    gmx_bool        bAppend;
    gmx_bool        bResetCountersHalfMaxH=FALSE;
    gmx_bool        bVV,bIterations,bFirstIterate,bTemp,bPres,bTrotter;
    real        temp0,mu_aver=0,dvdl;
    int         a0,a1,gnx=0,ii;
    atom_id     *grpindex=NULL;
    char        *grpname;
    t_coupl_rec *tcr=NULL;
    rvec        *xcopy=NULL,*vcopy=NULL,*cbuf=NULL;
    matrix      boxcopy={{0}},lastbox;
    tensor      tmpvir;
    real        fom,oldfom,veta_save,pcurr,scalevir,tracevir;
    real        vetanew = 0;
    double      cycles;
    real        saved_conserved_quantity = 0;
    real        last_ekin = 0;
    int         iter_i;
    t_extmass   MassQ;
    int         **trotter_seq; 
    char        sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
    int         handled_stop_condition=gmx_stop_cond_none; /* compare to get_stop_condition*/
    gmx_iterate_t iterate;
    gmx_large_int_t multisim_nsteps=-1; /* number of steps to do  before first multisim 
                                          simulation stops. If equal to zero, don't
                                          communicate any more between multisims.*/
#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////                                                                                     /////////
/////////                                                                                     /////////
/////////                               Initialization                                        /////////
/////////                                                                                     /////////
/////////                                                                                     /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////    
// (approximative) Masses, used as atomic species identifiers 
#define MASS_H 1
#define MASS_C 12
#define MASS_N 14
#define MASS_O 16
#define MASS_P 31
#define MASS_S 32
#define MASS_FE 56

// Energy levels
#define ENERGY_H "./Atomic_data/energy_levels_H.txt"
#define ENERGY_C "./Atomic_data/energy_levels_C.txt"
#define ENERGY_N "./Atomic_data/energy_levels_N.txt"
#define ENERGY_O "./Atomic_data/energy_levels_O.txt"
#define ENERGY_P "./Atomic_data/energy_levels_P.txt"
#define ENERGY_S "./Atomic_data/energy_levels_S.txt"
#define ENERGY_FE "./Atomic_data/energy_levels_F.txt"  

// Transition rates
#define RATES_H "./Atomic_data/rate_transitions_to_gromacs_H.txt"
#define RATES_C "./Atomic_data/rate_transitions_to_gromacs_C.txt"
#define RATES_N "./Atomic_data/rate_transitions_to_gromacs_N.txt"
#define RATES_O "./Atomic_data/rate_transitions_to_gromacs_O.txt"
#define RATES_P "./Atomic_data/rate_transitions_to_gromacs_P.txt"
#define RATES_S "./Atomic_data/rate_transitions_to_gromacs_S.txt"
#define RATES_FE "./Atomic_data/rate_transitions_to_gromacs_F.txt"   

// Collisional parameters
#define COLL_H "./Atomic_data/collisional_parameters_H.txt"
#define COLL_C "./Atomic_data/collisional_parameters_C.txt"
#define COLL_N "./Atomic_data/collisional_parameters_N.txt"
#define COLL_O "./Atomic_data/collisional_parameters_O.txt"
#define COLL_P "./Atomic_data/collisional_parameters_P.txt"
#define COLL_S "./Atomic_data/collisional_parameters_S.txt"
#define COLL_FE "./Atomic_data/collisional_parameters_F.txt"

// Statistical weight 
#define WEIGHT_H "./Atomic_data/statistical_weight_H.txt"
#define WEIGHT_C "./Atomic_data/statistical_weight_C.txt"
#define WEIGHT_N "./Atomic_data/statistical_weight_N.txt"
#define WEIGHT_O "./Atomic_data/statistical_weight_O.txt"
#define WEIGHT_P "./Atomic_data/statistical_weight_P.txt"
#define WEIGHT_S "./Atomic_data/statistical_weight_S.txt"
#define WEIGHT_FE "./Atomic_data/statistical_weight_F.txt"



// Create dictionaries for energy levels
struct Dictionary* Energylevels_H;
struct Dictionary* Energylevels_C;
struct Dictionary* Energylevels_N;
struct Dictionary* Energylevels_O;
struct Dictionary* Energylevels_P;
struct Dictionary* Energylevels_S;
struct Dictionary* Energylevels_FE;

// Construct "rate arrays" 
// These consists of each initial states and all possible final states together with the rate
struct Rate* transitionRates_H;
struct Rate* transitionRates_C;
struct Rate* transitionRates_N;
struct Rate* transitionRates_O;
struct Rate* transitionRates_P;
struct Rate* transitionRates_S;
struct Rate* transitionRates_FE;;

// Construct "coll arrays" 
// These consists of each initial states and all possible final states together with their parameters
struct Coll* collisions_H;
struct Coll* collisions_C;
struct Coll* collisions_N;
struct Coll* collisions_O;
struct Coll* collisions_P;
struct Coll* collisions_S;
struct Coll* collisions_FE;

// Construct weights
struct Weights Weigths_H;
struct Weights Weigths_C;
struct Weights Weigths_N;
struct Weights Weigths_O;
struct Weights Weigths_P;
struct Weights Weigths_S;
struct Weights Weigths_FE;

int numColl_H;
int numColl_C;
int numColl_N;
int numColl_O;
int numColl_P;
int numColl_S;
int numColl_FE;

int numRates_H;
int numRates_C;
int numRates_N;
int numRates_O;
int numRates_P;
int numRates_S;
int numRates_FE;


int** atom_configurations;
int** atomic_transitions;
int** GS_configurations;


const double boltzmann = 1.380649e-23; 
const double elementary_charge = 1.60217663e-19;
const double ev_to_joule = 1.60218e-19;
const double kelvin_to_ev = 8.61732814974493e-5;
const double pi_const = 3.1415;



    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
    bIonize  = (Flags & MD_IONIZE);
    bFFscan  = (Flags & MD_FFSCAN);
    bAppend  = (Flags & MD_APPENDFILES);
    if (Flags & MD_RESETCOUNTERSHALFWAY)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle,ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control 
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf,top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */ 
    bIterations = ((IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || (IR_NVT_TROTTER(ir))));        
    
    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

    nstglobalcomm = check_nstglobalcomm(fplog,cr,nstglobalcomm,ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    if (bRerunMD || bFFscan)
    {
        ir->nstxtcout = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&outf,&mdebin,
            force_vir,shake_vir,mu_tot,&bSimAnn,&vcm,state_global,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f,top_global->natoms);
    }

    /* Kinetic energy data */
    snew(ekind,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind);
    /* needed for iteration of constraints */
    snew(ekind_save,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind_save);
    /* Copy the cos acceleration to the groups struct */    
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,n_flexible_constraints(constr),
                                 (ir->bContinuation || 
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);

    if (DEFORM(*ir))
    {
#ifdef GMX_THREADS
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }

    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }

    if (DOMAINDECOMP(cr)) {
        top = dd_init_local_top(top_global);

        snew(state,1);
        dd_init_local_state(cr->dd,state_global,state);

        if (DDMASTER(cr->dd) && ir->nstfout) {
            snew(f_global,state_global->natoms);
        }
    } else {
        if (PAR(cr)) {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog,top_global,ir,cr);

            pd_cg_range(cr,&fr->cg0,&fr->hcg);
            pd_at_range(cr,&a0,&a1);
        } else {
            top = gmx_mtop_generate_local_top(top_global,ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        state = partdec_init_local_state(cr,state_global);
        f_global = f;

        atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

        if (vsite) {
            set_vsite_top(vsite,top,mdatoms,cr);
        }

        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
            graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
        }

        if (shellfc) {
            make_local_shells(cr,mdatoms,shellfc);
        }

        if (ir->pull && PAR(cr)) {
            dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog,ir->init_step,cr,TRUE,1,
                            state_global,top_global,ir,
                          
  state,&f,mdatoms,top,fr,
                            vsite,shellfc,constr,
                            nrnb,wcycle,FALSE);
    }

    update_mdatoms(mdatoms,state->lambda);

    if (MASTER(cr))
    {
        if (opt2bSet("-cpi",nfile,fnm))
        {
            /* Update mdebin with energy history if appending to output files */
            if ( Flags & MD_APPENDFILES )
            {
                restore_energyhistory_from_state(mdebin,&state_global->enerhist);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                done_energyhistory(&state_global->enerhist);
                init_energyhistory(&state_global->enerhist);
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(&state_global->enerhist,mdebin);
    }   

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG)) {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd,state);
    }

    /* Initialize constraints */
    if (constr) {
        if (!DOMAINDECOMP(cr))
            set_constraints(constr,top,ir,mdatoms,cr);
    }

    /* Check whether we have to GCT stuff */
    bTCR = ftp2bSet(efGCT,nfile,fnm);
    if (bTCR) {
        if (MASTER(cr)) {
            fprintf(stderr,"Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex,gnx);
        for(i=0; (i<gnx); i++) {
            grpindex[i] = i;
        }
    }

    if (repl_ex_nst > 0)
    {
        /* We need to be sure replica exchange can only occur
         * when the energies are current */
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "repl_ex_nst",&repl_ex_nst);
        /* This check needs to happen before inter-simulation
         * signals are initialized, too */
    }
    if (repl_ex_nst > 0 && MASTER(cr))
        repl_ex = init_replica_exchange(fplog,cr->ms,state_global,ir,
                                        repl_ex_nst,repl_ex_seed);

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
            {
                for(m=0; m<DIM; m++)
                {
                    if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                    {
                        state->v[i][m] = 0;
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog,constr,ir,mdatoms,state,f,
                               graph,cr,nrnb,fr,top,shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,NULL,
                             top->idef.iparams,top->idef.il,
                             fr->ePBC,fr->bMolPBC,graph,cr,state->box);
        }
    }

    debug_gmx();
  
    /* I'm assuming we need global communication the first time! MRS */
    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | (bVV ? CGLO_PRESSURE:0)
                  | (bVV ? CGLO_CONSTRAINT:0)
                  | (bRerunMD ? CGLO_RERUNMD:0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN:0));
    
    bSumEkinhOld = FALSE;
    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                    constr,NULL,FALSE,state->box,
                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,cglo_flags);
    if (ir->eI == eiVVAK) {
        /* a second call to get the half step temperature initialized as well */ 
        /* we do the same call as above, but turn the pressure off -- internally to 
           compute_globals, this is recognized as a velocity verlet half-step 
           kinetic energy calculation.  This minimized excess variables, but 
           perhaps loses some logic?*/
        
        compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                        wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                        constr,NULL,FALSE,state->box,
                        top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                        cglo_flags &~ CGLO_PRESSURE);
    }
    
    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT)) 
    {
        for(i=0; (i<ir->opts.ngtc); i++) 
        {
            copy_mat(ekind->tcstat[i].ekinh,ekind->tcstat[i].ekinh_old);
        } 
    }
    if (ir->eI != eiVV) 
    {
        enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                     and there is no previous step */
    }
    temp0 = enerd->term[F_TEMP];
    
    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterations) 
    {
            bufstate = init_bufstate(state);
    }
    if (bFFscan) 
    {
        snew(xcopy,state->natoms);
        snew(vcopy,state->natoms);
        copy_rvecn(state->x,xcopy,0,state->natoms);
        copy_rvecn(state->v,vcopy,0,state->natoms);
        copy_mat(state->box,boxcopy);
    } 
    
    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir,state,&MassQ,bTrotter);
    
    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr,FALSE));
        }
        fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
        if (bRerunMD)
        {
            fprintf(stderr,"starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name),opt2fn("-rerun",nfile,fnm));
            if (bVerbose)
            {
                fprintf(stderr,"Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
        {
            char tbuf[20];
            fprintf(stderr,"starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->nsteps >= 0)
            {
                sprintf(tbuf,"%8.1f",(ir->init_step+ir->nsteps)*ir->delta_t);
            }
            else
            {
                sprintf(tbuf,"%s","infinite");
            }
            if (ir->init_step > 0)
            {
                fprintf(stderr,"%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps,sbuf),tbuf,
                        gmx_step_str(ir->init_step,sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr,"%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps,sbuf),tbuf);
            }
        }
        fprintf(fplog,"\n");
    }

    /* Set and write start time */
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Started mdrun",runtime);
    wallcycle_start(wcycle,ewcRUN);
    if (fplog)
        fprintf(fplog,"\n");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
    chkpt_ret=fcCheckPointParallel( cr->nodeid,
                                    NULL,0);
    if ( chkpt_ret == 0 ) 
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", 0 );
#endif

    debug_gmx();

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////                                                                                     /////////
/////////                                                                                     /////////
/////////                           Hydrogen Index initialization                             /////////
/////////                                                                                     /////////
/////////                                                                                     /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////  



// Work out the hydrogen index
   int  num_hydrogen = 0;
   int* hydrogen_idx;
   int counter = 0;

   for(i=mdatoms->start; (i<mdatoms->nr); i++) {
      if (round(mdatoms->massT[i])==1) {
         num_hydrogen +=1;
      }
   }


   hydrogen_idx = (int*)malloc(num_hydrogen*sizeof(int));
   for(i=mdatoms->start; (i<mdatoms->nr); i++) {
      if ((int)round(mdatoms->massT[i])==1) {
         hydrogen_idx[counter] = i;
         counter +=1;
      }
   }



int count_H = 0;
int count_C = 0;
int count_N = 0;
int count_O = 0;
int count_P = 0;
int count_S = 0;
int count_FE = 0;
int approx_mass;


// Count atoms
for(i=mdatoms->start; (i<mdatoms->nr); i++) {
    approx_mass = (int)(round(mdatoms->massT[i]));
    switch(approx_mass) {
        case MASS_H: // Hydrogen
            count_H +=1;
            break;

        case MASS_C:
            count_C +=1;
            break;
            
                    
        case MASS_N:
            count_N +=1;
            break;

        case MASS_O:
            count_O +=1;
            break;

        case MASS_P:
            count_P +=1;
            break;

        case MASS_S:
            count_S +=1;
            break;

        case MASS_FE:
            count_FE +=1;
            break;
    }
}
if (count_H) {
    Energylevels_H    = readFileAndCreateEnergyDictionary(ENERGY_H);
    transitionRates_H = initializeRatesArrayFromFile(RATES_H,MASS_H);
    numRates_H        = countLinesInFile(RATES_H);
    collisions_H      = initializeCollArrayFromFile(COLL_H,MASS_H);
    numColl_H         = countLinesInFile(COLL_H);
    Weigths_H         = InitializeWeightsFromFile(WEIGHT_H,MASS_H);          
}
if (count_C) {
    Energylevels_C    = readFileAndCreateEnergyDictionary(ENERGY_C);
    transitionRates_C = initializeRatesArrayFromFile(RATES_C,MASS_C);
    numRates_C        = countLinesInFile(RATES_C);
    collisions_C      = initializeCollArrayFromFile(COLL_C,MASS_C);
    numColl_C         = countLinesInFile(COLL_C);
    Weigths_C         = InitializeWeightsFromFile(WEIGHT_C,MASS_C);         
}

if (count_N) {
    Energylevels_N    = readFileAndCreateEnergyDictionary(ENERGY_N);
    transitionRates_N = initializeRatesArrayFromFile(RATES_N,MASS_N);
    numRates_N        = countLinesInFile(RATES_N);
    collisions_N      = initializeCollArrayFromFile(COLL_N,MASS_N);
    numColl_N         = countLinesInFile(COLL_N);
    Weigths_N         = InitializeWeightsFromFile(WEIGHT_N,MASS_N);
}

if (count_O) {
    Energylevels_O    = readFileAndCreateEnergyDictionary(ENERGY_O);
    transitionRates_O = initializeRatesArrayFromFile(RATES_O,MASS_O);
    numRates_O        = countLinesInFile(RATES_O);
    collisions_O      = initializeCollArrayFromFile(COLL_O,MASS_O);
    numColl_O         = countLinesInFile(COLL_O);
    Weigths_O         = InitializeWeightsFromFile(WEIGHT_O,MASS_O);
}

if (count_P) {
    Energylevels_P    = readFileAndCreateEnergyDictionary(ENERGY_P);
    transitionRates_P = initializeRatesArrayFromFile(RATES_P,MASS_P);
    numRates_P        = countLinesInFile(RATES_P);
    collisions_P      = initializeCollArrayFromFile(COLL_P,MASS_P);
    numColl_P         = countLinesInFile(COLL_P);
    Weigths_P         = InitializeWeightsFromFile(WEIGHT_P,MASS_P);
}

if (count_S) {
    Energylevels_S    = readFileAndCreateEnergyDictionary(ENERGY_S);
    transitionRates_S = initializeRatesArrayFromFile(RATES_S,MASS_S);
    numRates_S        = countLinesInFile(RATES_S);
    collisions_S      = initializeCollArrayFromFile(COLL_S,MASS_S);
    numColl_S         = countLinesInFile(COLL_S);
    Weigths_S         = InitializeWeightsFromFile(WEIGHT_S,MASS_S);            
}

if (count_FE) {
    Energylevels_FE    = readFileAndCreateEnergyDictionary(ENERGY_FE);
    transitionRates_FE = initializeRatesArrayFromFile(RATES_FE,MASS_FE);
    numRates_FE        = countLinesInFile(RATES_FE);
    collisions_FE      = initializeCollArrayFromFile(COLL_FE,MASS_FE);
    numColl_FE         = countLinesInFile(COLL_FE);
    Weigths_FE         = InitializeWeightsFromFile(WEIGHT_FE,MASS_FE);
}


    /***********************************************************
     *
     *             Loop over MD steps 
     *
     ************************************************************/

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        rerun_fr.natoms = 0;
        if (MASTER(cr))
        {
            bNotLastFrame = read_first_frame(oenv,&status,
                                             opt2fn("-rerun",nfile,fnm),
                                             &rerun_fr,TRX_NEED_X | TRX_READ_V);
            if (rerun_fr.natoms != top_global->natoms)
            {
                gmx_fatal(FARGS,
                          "Number of atoms in trajectory (%d) does not match the "
                          "run input file (%d)\n",
                          rerun_fr.natoms,top_global->natoms);
            }
            if (ir->ePBC != epbcNONE)
            {
                if (!rerun_fr.bBox)
                {
                    gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
                }
                if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
                {
                    gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
                }
            }
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr,&rerun_fr,&bNotLastFrame);
        }

        if (ir->ePBC != epbcNONE)
        {
            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box,fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
    bInitStep = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep    = FALSE;
    bSumEkinhOld = FALSE;
    bExchanged   = FALSE;

    init_global_signals(&gs,cr,ir,repl_ex_nst);

    step = ir->init_step;
    step_rel = 0;




    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh,bGStatEveryStep,step);
    }

    if (MULTISIM(cr) && (repl_ex_nst <=0 ))
    {
        /* check how many steps are left in other sims */
        multisim_nsteps=get_multisim_nsteps(cr, ir->nsteps);
    }


    /* and stop now if we should */
    bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps) ||
                 ((multisim_nsteps >= 0) && (step_rel >= multisim_nsteps )));
    while (!bLastStep || (bRerunMD && bNotLastFrame)) {  // THIS IS THE MAIN MD LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        if (bRerunMD) {
            if (rerun_fr.bStep) {
                step = rerun_fr.step;
                step_rel = step - ir->init_step;
            }
            if (rerun_fr.bTime) {
                t = rerun_fr.time;
            }
            else
            {
                t = step;
            }
        } 
        else 
        {
            bLastStep = (step_rel == ir->nsteps);
            t = t0 + step*ir->delta_t;
        }

        if (ir->efep != efepNO)
        {
            if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
            {
                state_global->lambda = rerun_fr.lambda;
            }
            else
            {
                state_global->lambda = lam0 + step*ir->delta_lambda;
            }
            state->lambda = state_global->lambda;
            bDoDHDL = do_per_step(step,ir->nstdhdl);
        }

        if (bSimAnn) 
        {
            update_annealing_target_temp(&(ir->opts),t);
        }

        if (bRerunMD)
        {
            if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
            {
                for(i=0; i<state_global->natoms; i++)
                {
                    copy_rvec(rerun_fr.x[i],state_global->x[i]);
                }
                if (rerun_fr.bV)
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        copy_rvec(rerun_fr.v[i],state_global->v[i]);
                    }
                }
                else
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        clear_rvec(state_global->v[i]);
                    }
                    if (bRerunWarnNoV)
                    {
                        fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                                "         Ekin, temperature and pressure are incorrect,\n"
                                "         the virial will be incorrect when constraints are present.\n"
                                "\n");
                        bRerunWarnNoV = FALSE;
                    }
                }
            }
            copy_mat(rerun_fr.box,state_global->box);
            copy_mat(state_global->box,state->box);

            if (vsite && (Flags & MD_RERUN_VSITE))
            {
                if (DOMAINDECOMP(cr))
                {
                    gmx_fatal(FARGS,"Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
                }
                if (graph)
                {
                    /* Following is necessary because the graph may get out of sync
                     * with the coordinates if we only have every N'th coordinate set
                     */
                    mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
                if (graph)
                {
                    unshift_self(graph,state->box,state->x);
                }
            }
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step,ir->nstcomm));

        /* Copy back starting coordinates in case we're doing a forcefield scan */
        if (bFFscan)
        {
            for(ii=0; (ii<state->natoms); ii++)
            {
                copy_rvec(xcopy[ii],state->x[ii]);
                copy_rvec(vcopy[ii],state->v[ii]);
            }
            copy_mat(boxcopy,state->box);
        }

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS = (bFirstStep || ir->nstlist != 0);
            bNStList = bNS;
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);
            
            bNS = (bFirstStep || bExchanged || bNStList ||
                   (ir->nstlist == -1 && nlh.nabnsb > 0));

            if (bNS && ir->nstlist == -1)
            {
                set_nlistheuristics(&nlh,bFirstStep || bExchanged,step);
            }
        } 

        /* check whether we should stop because another simulation has 
           stopped. */
        if (MULTISIM(cr))
        {
            if ( (multisim_nsteps >= 0) &&  (step_rel >= multisim_nsteps)  &&  
                 (multisim_nsteps != ir->nsteps) )  
            {
                if (bNS)
                {
                    if (MASTER(cr))
                    {
                        fprintf(stderr, 
                                "Stopping simulation %d because another one has finished\n",
                                cr->ms->sim);
                    }
                    bLastStep=TRUE;
                    gs.sig[eglsCHKPT] = 1;
                }
            }
        }

        /* < 0 means stop at next step, > 0 means stop at next NS step */
        if ( (gs.set[eglsSTOPCOND] < 0 ) ||
             ( (gs.set[eglsSTOPCOND] > 0 ) && ( bNS || ir->nstlist==0)) )
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii=bFirstStep;
        if (ir->implicit_solvent && (step % ir->nstgbradii==0))
        {
            bBornRadii=TRUE;
        }
        
        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
                  (step % stepout == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
        {
            if (bRerunMD)
            {
                bMasterState = TRUE;
            }
            else
            {
                bMasterState = FALSE;
                /* Correct the new box if it is too skewed */
                if (DYNAMIC_BOX(*ir))
                {
                    if (correct_box(fplog,step,state->box,graph))
                    {
                        bMasterState = TRUE;
                    }
                }
                if (DOMAINDECOMP(cr) && bMasterState)
                {
                    dd_collect_state(cr->dd,state,state_global);
                }
            }

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                wallcycle_start(wcycle,ewcDOMDEC);
                dd_partition_system(fplog,step,cr,
                                    bMasterState,nstglobalcomm,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,do_verbose);
                wallcycle_stop(wcycle,ewcDOMDEC);
                /* If using an iterative integrator, reallocate space to match the decomposition */
            }
        }

        if (MASTER(cr) && do_log && !bFFscan)
        {
            print_ebin_header(fplog,step,t,state->lambda);
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms,state->lambda); 
        }

        if (bRerunMD && rerun_fr.bV)
        {
            
            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,NULL,NULL,NULL,NULL,mu_tot,
                            constr,NULL,FALSE,state->box,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
        }
        clear_mat(force_vir);
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////                                                ////////
////////                  IONIZE MODULE                 ////////
////////                                                ////////
////////                                                ////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

        /* Ionize the atoms if necessary */
        if (bIonize){       
            

                        
        int j;
        int n = top_global->natoms;  // the size of the array

        FILE *fp;
   
      static double thermalized_free_electron_energy;
      static double num_free_electrons;
      static double radius_gyration_equilibrium;
      double electron_temperature;
      double electron_density = 0;
      static double escaped_photoelectrons; 
      double current_num_free_electrons;

      double radius_of_gyration = 0;
      static double total_mass = 0;
      double rg_factor;
      double CURRENT_VOLUME;
      static double sample_volume; //in nm^3
      static double radius_of_sample;
      double T_MEAN,imax,width,rho,nphot;
      double INTENSITY;
      double max_x = 0;
      double max_y = 0;
      double max_z = 0;
      int logging = ir->userint5;
      int read_states = ir->userint4; 
      /* if 0 the simulation will initiate all charges and electronic states 
      in the ground state configuration. If 1 it will try to read 
      states and charges from ___. and ___.
        */
    /*
    LIST OF USERINT AND USERREAL THAT ARE IMPORTANT
    userint1 - Alter forcefield. If set to zero the everything should run as unmodifed gromacs 4.5.4 (hopefully) (default = 1)
    userint2 - Do charge transfer. Enables the charge transfer module (default = 1)
    userint3 - Use screened potential. Enables Debye shielding (default = 1)
    userint4 - Read electronic states from file. Reads electronic states from file. Useful for continued sims (default = 0)
    userint5 - Enable logging of electronic dynamics,writes a bunch of useful information, big performance drop due to I/O. (default = 0)

    The userreal are mostly the FEL parameters
    userreal1 - Peak of the gaussian pulse in ps
    userreal2 - Total number of photons in the pulse
    userreal3 - Width of the peak in ps. This is the sigma value of the guassian.
    userreal4 - Diameter of the focal spot (nm) 
    */


      USERINT1 = ir->userint1; // alter forcefield


      T_MEAN   = ir->userreal1*1e-12;  /* Peak of the gaussian pulse in ps           */
      nphot    = (double)ir->userreal2;  /* Number of photons                      */
      width    = ir->userreal3*1e-12;  /* Width of the peak (in time, in ps) SIGMA           */
      rho      = ((double)ir->userreal4)*1e-9*100;  /* Diameter of the focal spot (nm)       */

     
    if (nphot > 0) {
        imax  = (nphot/(pi_const*sqr(rho/2.0)))/(width*sqrt(2.0*pi_const)); // this quantity should have units cm^-2s^-1 (hence multiplication by 1e-12 ps -> s)
    } else {
        imax = 0;
    }

    INTENSITY = imax*exp(-0.5*sqr((t*1e-12-T_MEAN)/width));


            if (t==0) {
                    atom_configurations = (int**)malloc(n*sizeof(int*));  // allocate memory for the array of arrays
                    atomic_transitions= (int**)malloc(n*sizeof(int*));  // allocate memory for the array of arrays
                    GS_configurations= (int**)malloc(n*sizeof(int*));  // allocate memory for the array of arrays

            if (atom_configurations == NULL) {
                // Memory allocation failed
                printf("Failed to allocate memory\n");
                exit(0);
            }

            // init files
             fp = fopen("./simulation_output/mean_charge_vs_time.txt", "w");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fclose(fp);  // close the file

            fp = fopen("./simulation_output/pulse_profile.txt", "w");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fclose(fp);  // close the file


            fp = fopen("./simulation_output/electron_data.txt", "w");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fclose(fp);  // close the file

            fp = fopen("./simulation_output/electronic_transition_log.txt", "w");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fclose(fp);  // close the file


            electron_temperature = 0.0;
            num_free_electrons = 0.0;
            thermalized_free_electron_energy=0.0;
            escaped_photoelectrons = 0.0;

            double y;   
            double x=5.0;                                                                                                                 
            y = Exponential_Integral_Ei( x ); 
           // printf("y: %lf and x=%lf", y,x); 
           // exit(0); 

           // determine the radius of the sampl, estimated by finding the atoms the maximum coordintes of x, y and z.

            for (i = 0; i<mdatoms->nr; i++) {
                    if (state->x[i][0] > max_x) {
                        max_x = state->x[i][0] ;
                    }

                    if (state->x[i][1] > max_y) {
                        max_y = state->x[i][1] ;
                    }

                    if (state->x[i][2] > max_z) {
                        max_z = state->x[i][2] ;
                    }

                }

            radius_of_sample = sqrt(max_x*max_x + max_y*max_y + max_z*max_z);
            //printf("Radius of sample: %lf", radius_of_sample);

             for (i=mdatoms->start; i <mdatoms->nr; i++) {
                radius_of_gyration += mdatoms->massT[i]*(state->x[i][0]*state->x[i][0] + state->x[i][1]*state->x[i][1] + state->x[i][2]*state->x[i][2]);
                total_mass += mdatoms->massT[i];
                mdatoms->chargeA[i] = 0.0;
                
             }
             radius_of_gyration = sqrt(radius_of_gyration/total_mass);
             //printf("Radis of gyration at timestep 0: %lf\n", radius_of_gyration);
             radius_gyration_equilibrium = radius_of_gyration;
    
            srand(time(NULL));
            num_free_electrons=0;
            thermalized_free_electron_energy = 0;
            sample_volume= 100*100*100*1e-9*1e-9*1e-9*radius_of_sample*radius_of_sample*radius_of_sample;


            for (i = 0; i < n; i++) { // works in serial
                atom_configurations[i] = (int*)malloc(3*sizeof(int));  // allocate memory for each subarray
                GS_configurations[i] = (int*)malloc(3*sizeof(int)); 
                atomic_transitions[i] = (int*)malloc(3*sizeof(int));  // allocate memory for each subarray and set initial values to 0

             if (atomic_transitions[i] == NULL) {
                        // Memory allocation failed for subarray
                        printf("Failed to allocate memory for subarray %d\n", i);
                        // You may need to free the memory allocated so far and exit the program
                    }

             if (atom_configurations[i] == NULL) {
                        // Memory allocation failed for subarray
                        printf("Failed to allocate memory for subarray %d\n", i);
                        exit(0);
                        // You may need to free the memory allocated so far and exit the program
                    }
             if (GS_configurations[i] == NULL) {
                        // Memory allocation failed for subarray
                        printf("Failed to allocate memory for subarray %d\n", i);
                        exit(0);
                        // You may need to free the memory allocated so far and exit the program
                    }
            }


            // TODO Initialize groundstates outside MD loop
            for(i=mdatoms->start; (i<mdatoms->nr); i++) {

                switch((int)round(mdatoms->massT[i])) {
                    case MASS_H: // Hydrogen
                        atom_configurations[i][0] = 1;// Hydrogen
                        atom_configurations[i][1] = 0;
                        atom_configurations[i][2] = 0;
                        GS_configurations[i][0] = 1;
                        GS_configurations[i][1] = 0;
                        GS_configurations[i][2] = 0; 
                        break;

                    case MASS_C:
                        atom_configurations[i][0] = 2;// Carbon
                        atom_configurations[i][1] = 4;
                        atom_configurations[i][2] = 0;
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 4;
                        GS_configurations[i][2] = 0; 
                        break;
  
                    
                    case MASS_N:
                        atom_configurations[i][0] = 2;// Nitrogen
                        atom_configurations[i][1] = 5;
                        atom_configurations[i][2] = 0;
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 5;
                        GS_configurations[i][2] = 0; 
                        break;

                    case MASS_O:
                        atom_configurations[i][0] = 2;// Oxygen
                        atom_configurations[i][1] = 6;
                        atom_configurations[i][2] = 0;
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 6;
                        GS_configurations[i][2] = 0; 
                        break;

                    case MASS_P:
                        atom_configurations[i][0] = 2;// Phosphorus
                        atom_configurations[i][1] = 8;
                        atom_configurations[i][2] = 5;
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 8;
                        GS_configurations[i][2] = 5; 
                        break;

                    case MASS_S:
                        atom_configurations[i][0] = 2;// sulfur
                        atom_configurations[i][1] = 8;
                        atom_configurations[i][2] = 6;
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 8;
                        GS_configurations[i][2] = 6; 
                        break;

                    case MASS_FE:
                        atom_configurations[i][0] = 2;// iron
                        atom_configurations[i][1] = 8;
                        atom_configurations[i][2] = 16; 
                        GS_configurations[i][0] = 2;
                        GS_configurations[i][1] = 8;
                        GS_configurations[i][2] = 16; 
                        break;
                }
            }
            // We read charges and electronic states from a previous simulation.
            if (read_states) {
                int n = top_global->natoms; 

                fp = fopen("configurations.bin", "rb");
                if (fp == NULL) {
                    printf("Error: could not open file for reading.\n");
                    exit(1);
                }

                for (i = 0; i < n; i++) {
                    fread(atom_configurations[i], sizeof(int), 3, fp);
                }

                fclose(fp);

                fp = fopen("charges.bin", "rb");
                if (fp == NULL) {
                    printf("Error: could not open file for reading.\n");
                    exit(1);
                }

                for (i = 0; i < n; i++) {
                    fread(&mdatoms->chargeA[i], sizeof(double), 1, fp);
                }

                fclose(fp);


            }


        } // IF t !=0!! 

        
            int match[3];
            current_num_free_electrons  = num_free_electrons;
   
            double photon_energy = 600;  // in eV
            double critical_potential;
            double ionization_energy;
            int found_match, found_match_collisional;
            double min_value;
            int min_value_index;
            double net_charge;
            double Q_D, Q_A, V_B, V_D, R_cob,R_crit;
            double DT_current;
            int charge_transfer = ir->userint2; // 1;
            double tev;
            double deltaE;
            double B;
            double sigma_collisional;
            double sigma_recombination;
            double qb;
            double gi, gf; 
          
///////////////////////////////////////////////////
////                                           ////
////                                           ////
////     Charge Transfer Module START          ////
////                                           ////
////                                           ////
///////////////////////////////////////////////////

        int i;
        if (charge_transfer && (t!=0.0)) {
           
            for (i = 0; i<n; i++){ //set initial values to 0
                atomic_transitions[i][0] = 0;
                atomic_transitions[i][1] = 0;
                atomic_transitions[i][2] = 0;
            }


      // Charge transfer module
      double E_donor;
      int indx;
      int number_of_charge_transfers = 0;
      int donor_idx;
      int acceptor_idx;
      int j;
      double R_min = 1e7;
      double R_crit = 1e7;
      int R_min_idx;
      int R_crit_idx;
      int* idx_map = (int*)malloc((mdatoms->nr)*sizeof(int)); 

      // Populate map
      for(j=mdatoms->start; (j<mdatoms->nr); j++) {
         idx_map[j] = j;
      }

// shuffles/randomizes an array
void shuffle(int arr[], int size) {
    // Initialize the random number generator with the current time
    srand(time(NULL));

    // Perform Fisher-Yates shuffle
    for (i = size - 1; i > 0; i--) {
        // Generate a random index between 0 and i (inclusive)
      j = rand() % (i + 1);

        // Swap the elements at indices i and j
        int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}
            // Shuffe map
            shuffle(idx_map,mdatoms->nr);
         int j2,k;
         for(j2=mdatoms->start; (j2<mdatoms->nr); j2++) {
            j = idx_map[j2]; // Use randomized index
            R_min = 1e7;
            R_min_idx = -1; // Reset index, if no suitible transfer is found we check if == -1 
            for(k=0; k<num_hydrogen; k++) {
                  i = hydrogen_idx[k];
                  if (i == j) {
                     continue;
                  }

                  if (atom_configurations[i][0] == 0 && 
                      atom_configurations[i][1] == 0 && 
                      atom_configurations[i][2] == 0
                      ) {
                     continue;
                  }


                  if (abs((mdatoms->chargeA[j])-mdatoms->chargeA[i])<1e-5 ) {
                      continue;
                  }
                                R_cob = sqrt((state->x[i][0]-state->x[j][0])*(state->x[i][0]-state->x[j][0]) + (state->x[i][1]-state->x[j][1])*(state->x[i][1]-state->x[j][1])+ (state->x[i][2]-state->x[j][2])*(state->x[i][2]-state->x[j][2]));

                  if (R_cob > 1) {   
                     continue;
                  }
 
                  if (mdatoms->chargeA[i] > mdatoms->chargeA[j]) {
                     // the one with more net charge is the acceptor
                     Q_D = mdatoms->chargeA[j];        
                     Q_A = mdatoms->chargeA[i]; 
                     // j donor
                     donor_idx = j;
                     acceptor_idx = i;

                  }
                  else {
                     Q_D = mdatoms->chargeA[i];
                     Q_A = mdatoms->chargeA[j]; 
                     // i donor
                     donor_idx = i;
                     acceptor_idx = j;
                  }

                              // Find state of donor 
                                 char state_str[20];
                                 switch((int)round(mdatoms->massT[donor_idx])) { // Match the mass and lookup the energy of the state.
                                    case MASS_H: // Hydrogen
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_H,state_str);
                                       break;
                                    case MASS_C: // Carbon
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_C,state_str);
                                       break;
                                    case MASS_N: // Nitrogen
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_N,state_str);
                                       break;
                                    case MASS_O: // Oxygen
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_O,state_str);
                                       break;
                                    case MASS_P: // Phosphor
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_P,state_str);
                                       break;
                                    case MASS_S: // Sulfur
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_S,state_str);
                                       break;
                                    case MASS_FE: // IRON
                                       sprintf(state_str,"%d %d %d",atom_configurations[donor_idx][0],atom_configurations[donor_idx][1],atom_configurations[donor_idx][2]);
                                       E_donor = get(Energylevels_FE,state_str);
                                       break;
                                                                                }

                               if ((int)E_donor ==-1) {
                                        printf("No matching state for atom %d with mass: %lf \n",donor_idx,round(mdatoms->massT[donor_idx]));
                                        exit(0);
                                       }

                              R_crit = (Q_D + 1 + 2*sqrt((Q_D+1)*Q_A))/E_donor;

                              if (R_crit < R_cob) {
                                 continue;
                                 }

                              if (R_cob < R_min) {
                                 R_min = R_cob;
                                 R_min_idx = acceptor_idx;
                                 } 
                              } // End of hydrogen loop

                              // If no suitable transfer is found.
                              // Then continue and don't do one.
                              if (R_min_idx == -1) { 
                                 continue;
                              }

                              int INDEX1;
                              int INDEX2;

                              for(indx=2; (indx>-1); indx--) {
                                 if ((double)atom_configurations[donor_idx][indx]<1e-5) {
                                    if (indx ==0){
                                       INDEX1 = 100;
                                    }
                                 } else {
                                    INDEX1 = indx;
                                    break;
                                 }
                              }
                              for(indx=2; (indx>-1); indx--) {
                                 if ((double)atom_configurations[R_min_idx][indx]<1e-5) {
                                    if (indx ==0) {
                                       INDEX2 = 100;
                                    }   
                                 } else {
                                    INDEX2 = indx;
                                    break;
                                 }
                              }

                              if (INDEX1 == 100 || INDEX2 == 100 ){
                                 continue;
                              }
                                    if (mdatoms->chargeA[R_min_idx] > mdatoms->chargeA[donor_idx]) {

                                        if ((atom_configurations[donor_idx][INDEX1] + atomic_transitions[donor_idx][INDEX1] -1 ) < 0  || 
                                            (atom_configurations[R_min_idx][INDEX2] + atomic_transitions[R_min_idx][INDEX2] + 1 ) > GS_configurations[R_min_idx][INDEX2]
                                            ) {
                                                    // Make sure that we are obeying how many electrons this orbital can occupy
                                                    break;

                                             }
                                        else {
                                            atomic_transitions[donor_idx][INDEX1]-=1; 
                                            atomic_transitions[R_min_idx][INDEX2]+=1; 

                                        }
                                    }


                                    if (mdatoms->chargeA[j] > mdatoms->chargeA[R_min_idx]) {

   if ((atom_configurations[R_min_idx][INDEX2] + atomic_transitions[R_min_idx][INDEX2] -1 ) < 0  || 
       (atom_configurations[donor_idx][INDEX1] + atomic_transitions[donor_idx][INDEX1] + 1 ) > GS_configurations[donor_idx][INDEX1]
       ) {
                                                    // Make sure that we are obeying how many electrons this orbital can occupy
                                                    break; 
                                                 }
                                        else {
                                            atomic_transitions[R_min_idx][INDEX2]-=1; 
                                            atomic_transitions[donor_idx][INDEX1]+=1; 

                                        }
                                    }
                           number_of_charge_transfers+=1;


                // update all atoms' configurations due to charge transfer
                if (atomic_transitions[R_min_idx][0]>0 ||
                    atomic_transitions[R_min_idx][1]>0 ||
                    atomic_transitions[R_min_idx][2]>0 ||
                    atomic_transitions[R_min_idx][0]<0 ||
                    atomic_transitions[R_min_idx][1]<0 ||
                    atomic_transitions[R_min_idx][2]<0) {
                  

                    if (logging) {
                        fp = fopen("./simulation_output/electronic_transition_log.txt", "a");  // open the file in write mode
                        if (fp == NULL) {
                            printf("Error: could not open file.\n");
                            return 1;
                        }
                        fprintf(fp, "Charge transfer at t: %lf, Previous state: [%lf, %lf, %lf]. Charge init: %lf ", t,atom_configurations[R_min_idx][0], atom_configurations[R_min_idx][1],atom_configurations[R_min_idx][2], mdatoms->chargeA[R_min_idx]);  // write transition to file
                    }

                    // acceptor
                    atom_configurations[R_min_idx][0] = atom_configurations[R_min_idx][0] + atomic_transitions[R_min_idx][0];
                    atom_configurations[R_min_idx][1] = atom_configurations[R_min_idx][1] + atomic_transitions[R_min_idx][1];
                    atom_configurations[R_min_idx][2] = atom_configurations[R_min_idx][2] + atomic_transitions[R_min_idx][2];

                    mdatoms->chargeA[R_min_idx] =  (double)((GS_configurations[R_min_idx][0] - atom_configurations[R_min_idx][0]) + 
                                                            (GS_configurations[R_min_idx][1] - atom_configurations[R_min_idx][1]) + 
                                                            (GS_configurations[R_min_idx][2] - atom_configurations[R_min_idx][2])); // minus means electron lost, + means got one 
                    // Donator 
                    atom_configurations[donor_idx][0] = atom_configurations[donor_idx][0] + atomic_transitions[donor_idx][0];
                    atom_configurations[donor_idx][1] = atom_configurations[donor_idx][1] + atomic_transitions[donor_idx][1];
                    atom_configurations[donor_idx][2] = atom_configurations[donor_idx][2] + atomic_transitions[donor_idx][2];

                    mdatoms->chargeA[donor_idx] =  (double)((GS_configurations[donor_idx][0] - atom_configurations[donor_idx][0]) + 
                                                            (GS_configurations[donor_idx][1] - atom_configurations[donor_idx][1]) + 
                                                            (GS_configurations[donor_idx][2] - atom_configurations[donor_idx][2])); // minus means electron lost, + means got one
                    if (logging) {
                        fprintf(fp, " New state: [%lf, %lf, %lf]. With %i charge transfers occuring. Charge after: %lf \n", atom_configurations[R_min_idx][0], atom_configurations[R_min_idx][1],atom_configurations[R_min_idx][2],number_of_charge_transfers, mdatoms->chargeA[R_min_idx]);
                        fclose(fp);  // close the file
                    }

                } else {
                atom_configurations[R_min_idx][0] = atom_configurations[R_min_idx][0] + atomic_transitions[R_min_idx][0];
                atom_configurations[R_min_idx][1] = atom_configurations[R_min_idx][1] + atomic_transitions[R_min_idx][1];
                atom_configurations[R_min_idx][2] = atom_configurations[R_min_idx][2] + atomic_transitions[R_min_idx][2];
                }
               } // ATOM LOOP
               free(idx_map);
            } // CHARGE TRANSFER IF


///////////////////////////////////////////////////
////                                           ////
////       Charge Transfer Module END          ////
////  Photon matter interaction Module START   ////
////                                           ////
///////////////////////////////////////////////////
            
            int currentNumRates;
            int currentNumColl;
            struct Rate* currentRate;
            struct Dictionary* currentDict;
            struct Coll* currentColl;
            struct Weights currentWeights;

            radius_of_gyration = 0;
            for (i=mdatoms->start; i <mdatoms->nr; i++) {
                // calculate radius of gyration for current timestep to estimate the expansion of the system
                radius_of_gyration += mdatoms->massT[i]*(state->x[i][0]*state->x[i][0] + state->x[i][1]*state->x[i][1] + state->x[i][2]*state->x[i][2]);
             }
             radius_of_gyration = sqrt(radius_of_gyration/total_mass);
             rg_factor = radius_of_gyration/radius_gyration_equilibrium;
             CURRENT_VOLUME = sample_volume*rg_factor*rg_factor*rg_factor; // multiply by rg_factor^3 since sample volume is the volume of a sphere

         int k;
         for(k=mdatoms->start; (k<mdatoms->nr); k++) { // This is the big loop which handles the electronic states for each atom in the system  
            DT_current = 0.0;
            //printf("Current atom %d \n\n\n\n\n",(int)round(mdatoms->massT[k]));
            // Crude way to handle the different atom species.
            switch ((int)round(mdatoms->massT[k])) {
               case MASS_H:
                  currentDict = Energylevels_H;
                  currentRate = transitionRates_H;
                  currentNumRates = numRates_H;
                  currentColl = collisions_H;
                  currentNumColl = numColl_H;
                  currentWeights = Weigths_H;
                  break;

               case MASS_C:
                  currentDict = Energylevels_C;
                  currentRate = transitionRates_C;
                  currentNumRates = numRates_C;
                  currentColl = collisions_C;
                  currentNumColl = numColl_C;
                  currentWeights = Weigths_C;
                  break;

               case MASS_N:
                  currentDict = Energylevels_N;
                  currentRate = transitionRates_N;
                  currentNumRates = numRates_N;
                  currentColl = collisions_N;
                  currentNumColl = numColl_N;
                  currentWeights = Weigths_N;
                  break;

               case MASS_O:
                  currentDict = Energylevels_O;
                  currentRate = transitionRates_O;
                  currentNumRates = numRates_O;
                  currentColl = collisions_O;
                  currentNumColl = numColl_O;
                  currentWeights = Weigths_O;
                  break;

               case MASS_P:
                  currentDict = Energylevels_P;
                  currentRate = transitionRates_P;
                  currentNumRates = numRates_P;
                  currentColl = collisions_P;
                  currentNumColl = numColl_P;
                  currentWeights = Weigths_P;
                  break;

               case MASS_S:
                  currentDict = Energylevels_S;
                  currentRate = transitionRates_S;
                  currentNumRates = numRates_S;
                  currentColl = collisions_S;
                  currentNumColl = numColl_S;
                  currentWeights = Weigths_S;
                  break;

               case MASS_FE:
                  currentDict = Energylevels_FE;
                  currentRate = transitionRates_FE;
                  currentNumRates = numRates_FE;
                  currentColl = collisions_FE;
                  currentNumColl = numColl_FE;
                  currentWeights = Weigths_FE;
                  break;
               }

            double DT_scale=0.0;
            // do electronic stuff untill we reach timescales bigger than MD timestep
            while (DT_current<ir->delta_t*1e-12) {


            match[0] = atom_configurations[k][0];
            match[1] = atom_configurations[k][1];
            match[2] = atom_configurations[k][2];

            found_match = 0;      
            found_match_collisional =0 ;     
            char state_str[20];
            int* state_diff;
            int* final_state;
            double transition_rate;
            double* collision_params; 
            struct transition* possible_transitions;
            int total_transitions;
            int b;

            int match_index = RatesStateIndex(currentRate, currentNumRates, match);
            if (match_index != -1) found_match = 1;
            int match_index_coll = CollStateIndex(currentColl, currentNumColl, match);
            if (match_index_coll != -1) found_match_collisional = 1;



            
            // If there is any match
            if (found_match && found_match_collisional) {  // we could change this to a && and not do 3 checks

                total_transitions = (currentRate[match_index].num_transitions + currentColl[match_index_coll].num_transitions);
                possible_transitions = (struct transition*)malloc(total_transitions*sizeof(struct transition)); 
                for (i = 0; i<total_transitions; i++) {
                    possible_transitions[i].final_state = (int*)malloc(3*sizeof(int));
                }

                state_diff = (int*)malloc(3*sizeof(int));
                final_state = (int*)malloc(3*sizeof(int));

                for (i = 0; i < currentRate[match_index].num_transitions; i++) { // This loops over the "normal" possible transitions

                    // Store current final state
                    final_state[0] = currentRate[match_index].final_states[i][0];
                    final_state[1] = currentRate[match_index].final_states[i][1];
                    final_state[2] = currentRate[match_index].final_states[i][2];
                    transition_rate = currentRate[match_index].rates[i];

                    // Calculate the difference between the states
                    state_diff[0] = (int)abs(match[0] - final_state[0]);
                    state_diff[1] = (int)abs(match[1] - final_state[1]);
                    state_diff[2] = (int)abs(match[2] - final_state[2]);

                    // Check for auger decay
                    if ((state_diff[0] == 2) ||
                        (state_diff[1] == 2) ||
                        (state_diff[2] == 2)) {

                        // The transition is Auger decay since the occupations in one shell differs by 2, so no need to scale by pulse profile //
                        possible_transitions[i].final_state[0] = final_state[0];
                        possible_transitions[i].final_state[1] = final_state[1];
                        possible_transitions[i].final_state[2] = final_state[2];
                        possible_transitions[i].rate = transition_rate;

                    } // Check for fluoresence 
                    else if ((state_diff[0] + state_diff[1] + state_diff[2]) == 0) {

                        possible_transitions[i].final_state[0] = final_state[0];
                        possible_transitions[i].final_state[1] = final_state[1];
                        possible_transitions[i].final_state[2] = final_state[2];
                        possible_transitions[i].rate = transition_rate;

                    } // If none of fluoresence or auger decay it is photoionization
                    else {

                        // This is a photoionization transition, since the difference is 1 in netcharge, we scale the rate by the pulse 
                        possible_transitions[i].final_state[0] = final_state[0];
                        possible_transitions[i].final_state[1] = final_state[1];
                        possible_transitions[i].final_state[2] = final_state[2];
                        possible_transitions[i].rate = INTENSITY*transition_rate;
                    }
                }

                free(state_diff);


                deltaE = 200.0; // energy difference between the two levels
                tev = thermalized_free_electron_energy; //50.0;
                B = deltaE/tev;
                electron_density = num_free_electrons/(CURRENT_VOLUME); //current_num_free_electrons/sample_volume;

            int offset = currentRate[match_index].num_transitions;
            for (i = offset; i < total_transitions; i++) { // This loops over 

                    // Assign for ease of access
                    final_state[0] = currentColl[match_index_coll].final_states[i-offset][0];
                    final_state[1] = currentColl[match_index_coll].final_states[i-offset][1];
                    final_state[2] = currentColl[match_index_coll].final_states[i-offset][2];
                    collision_params = currentColl[match_index_coll].coll_rates;

                    // Save the state of the transtion
                    possible_transitions[i].final_state[0] = final_state[0];
                    possible_transitions[i].final_state[1] = final_state[1];
                    possible_transitions[i].final_state[2] = final_state[2];

                    // This would be equivalent
                    qb = (collision_params[0] + 
                          collision_params[1]*sqrt(pi_const*B)*erfc(B)*exp(B*B) +
                          collision_params[2]*B + 
                          collision_params[3]*Exponential_Integral_Ei(B)*exp(B) * 
                          collision_params[4]
                          );

                    // Weight initial state
                    gi = getWeightForState(currentWeights,match);

                    // Weight final state 
                    gf = getWeightForState(currentWeights,currentColl[match_index_coll].final_states[i-offset]);
            

                    // Calculate the collisional transition rate based on the collisional parameters
                    if (abs(current_num_free_electrons-0)<1e-5 || isnan(qb) || qb < 0.0) {
                        sigma_collisional = 0.00001;
                        sigma_recombination = 0.00001;
                    } else {
                        // This is probably the main case
                        sigma_collisional = electron_density*(1.09*1e-6*qb*exp(-b))/(deltaE*sqrt(tev));
                        sigma_recombination = sigma_collisional*1.66*1e-22*(gi/gf)*(electron_density/pow(tev, 1.5))*exp(-(deltaE)/(boltzmann*tev));
                  
                        if (sigma_collisional < 0 || abs(sigma_collisional)<1e-5) {
                            sigma_collisional = 0.00001; 
                        }
                        if (sigma_recombination < 0 || abs(sigma_recombination)<1e-5) {
                            sigma_recombination = 0.0001; 
                        }
                        if (sigma_collisional > 1e-3 || sigma_recombination> 1e-3) {
                        // Nothing happens here?
                        }
                    }

                    // save the final collsional rate
                    if ((match[0] - final_state[0] +
                        match[1] - final_state[1] + 
                        match[2] - final_state[2]) > 0) {
                        // Collision
                        possible_transitions[i].rate = sigma_collisional;
                    } else {
                        // Recombination
                        possible_transitions[i].rate = sigma_recombination;
                    }
                    // Handle exception
                    if (isinf(sigma_collisional) || isnan(sigma_collisional)) {
                        possible_transitions[i].rate;
                    }

            } // For-loop ends here
            free(final_state);
            // Now read collisional ionization data ///

            } else { // I no match is found we go here

                if (atom_configurations[k][0]==0 && 
                    atom_configurations[k][1]==0 && 
                    atom_configurations[k][2]==0) {
                    printf("All is zero.. but still possible to get recombination.. keep this state (do not continue)\n");
                } else {
                    printf("No match found for match[0]=%d, match[1]=%d, match[2]=%d, mass= %lf, k=%i\n",
                            atom_configurations[k][0], 
                            atom_configurations[k][1],
                            atom_configurations[k][2], 
                            mdatoms->massT[k],
                            k
                            );
                    exit(0);
                }
            }

            ionization_energy = 871.0;

            // Function to generate random array
            double* generate_random_array(int size) {
                double* arr_rand = malloc(size * sizeof(double));
               // srand(time(NULL));  // Seed the random number generator
                for (i = 0; i < size; i++) {
                   // srand(time(NULL)+1);  // Seed the random number generator
                    arr_rand[i] = (double) rand() / (RAND_MAX+1.0);
                }
                return arr_rand;
            }


            double* arr_rand = generate_random_array(total_transitions); // Generate num_transitions random numbers.
            double* dt_processes = (double*)malloc(total_transitions*sizeof(double));  

            min_value= 1.1* ir->delta_t*1e-12; // Set a bit bigger than MD timestep

            for (i = 0; i < total_transitions; i++) { // Loop over possible transitions
                if (possible_transitions[i].rate<1e-7) {
                // if rate is zero, we set the DT to be large.. no transition should occur
                dt_processes[i]= 1000;
                } else { 
                dt_processes[i]= (double)(-1.0*(log(arr_rand[i])/possible_transitions[i].rate)); // This is the magical MC-step
                }

                if (dt_processes[i] < min_value) {
                    min_value = dt_processes[i];
                    min_value_index = i; // Index to transtion with minimum value
                }
                if (isinf(dt_processes[i])) {     
                    printf("DT is infinite. Exiting...\n");
                    exit(0); // If DT is infinite we exit.
                }
            }
            //printf("DT_current is %e and md timestep is %e\n",DT_current,ir->delta_t*1e-12);
            DT_current +=min_value;
            if (DT_current>ir->delta_t*1e-12) { // Break if we go over md timestep.
                free(dt_processes);
                free(arr_rand);
                for (i = 0; i<total_transitions; i++) {
                    free(possible_transitions[i].final_state);
                }
                free(possible_transitions);
                break;
        
             } else {
                if (logging) {
                    fp = fopen("./simulation_output/electronic_transition_log.txt", "a");  // open the file in write mode
                    if (fp == NULL) {
                        printf("Error: could not open file.\n");
                        return 1;
                    }
                    fprintf(fp, "Making a transition, the time is %lf and current Monte-Carlo DT is: %e.\n", t, DT_current);  // write transition to file
                    fclose(fp);  // close the file
                }

                if (isinf(DT_current)) {
                    printf("DT_current is infinite, this should not happen. Exiting...\n");
                    exit(0);
                }
            }
            net_charge = (double)(match[0] - possible_transitions[min_value_index].final_state[0] +
                          match[1] - possible_transitions[min_value_index].final_state[1] +
                          match[2] - possible_transitions[min_value_index].final_state[2]);  

            critical_potential = 100*(elementary_charge*elementary_charge*4*pi_const*9e9*1e9*(escaped_photoelectrons/CURRENT_VOLUME)*state->box[XX][XX]*state->box[XX][XX])/3; // box XX is temporary radius
            if ((photon_energy-ionization_energy)*1.60218e-19 > critical_potential) {  // check if kinetic energy is enough for it to escape
                escaped_photoelectrons+=1;
            } else {

                num_free_electrons+= net_charge;
                thermalized_free_electron_energy += net_charge*25; // add the energy of the trapped electron to the total
            }

            electron_density = num_free_electrons/(CURRENT_VOLUME); // multiply by rg_factor^3, since the volume is for a sphere
            electron_temperature = kelvin_to_ev*((2.0/(3.0*boltzmann))*((ev_to_joule*thermalized_free_electron_energy)/num_free_electrons));

            if (logging) {
                fp = fopen("./simulation_output/electronic_transition_log.txt", "a");  // open the file in write mode
                    if (fp == NULL) {
                        printf("Error: could not open file.\n");
                        return 1;
                    }
                fprintf(fp, "New transition: time is t: %lf, netcharge is: %lf. Previous state: [%d, %d, %d]. New state: [%d, %d, %d].\n",
                        t,
                        net_charge,
                        match[0],
                        match[1],
                        match[2],
                        possible_transitions[min_value_index].final_state[0], 
                        possible_transitions[min_value_index].final_state[1],
                        possible_transitions[min_value_index].final_state[2]
                        );  // write transition to file

                fprintf(fp, "Number of free electrons: %lf, electron_density: %lf and escaped electrons: %lf and electron temperature: %lf \n", 
                        num_free_electrons, 
                        electron_density,
                        escaped_photoelectrons,
                        electron_temperature
                        );  // write current state of system to file

                fclose(fp);  // close the file
            }
                match[0] = possible_transitions[min_value_index].final_state[0];
                match[1] = possible_transitions[min_value_index].final_state[1];
                match[2] = possible_transitions[min_value_index].final_state[2];

                atom_configurations[k][0] = match[0];
                atom_configurations[k][1] = match[1];
                atom_configurations[k][2] = match[2];

                mdatoms->chargeA[k] += net_charge;

                free(dt_processes);
                free(arr_rand);
                for (i = 0; i<total_transitions; i++) {
                    free(possible_transitions[i].final_state);
                }
                free(possible_transitions);

                
            } // THIS IS THE DT WHILE LOOP 
        } // THIS IS THE FOR LOOP
     
    
       double mean_charge = 0;
       double N_atoms;
       for (i = 0; i <mdatoms->nr; i++) {
            mean_charge+=mdatoms->chargeA[i];            
        }
    
        N_atoms = (double)mdatoms->nr; 
        mean_charge = mean_charge/N_atoms; 

        if (electron_density > 0 && electron_temperature >0) {
            debye_length = 1e6*69*sqrt((11606*electron_temperature)/electron_density); // calculate debye length 
        }
        else {
            debye_length = 100;
        }

        if (ir->userint3 == 0) { 
            // dont use screened coulomb potentials, we get back regular coulomb by settting the debye length a very large number. exp(-x/debye) goes to 1
            debye_length = 1e9;
        }


        if (logging) {
            fp = fopen("./simulation_output/mean_charge_vs_time.txt", "a");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fprintf(fp, "%lf %lf\n", t, mean_charge);  // write each element to the file
            fclose(fp);  // close the file

            fp = fopen("./simulation_output/pulse_profile.txt", "a");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }
            fprintf(fp, "%lf %lf\n", t, INTENSITY);  // write each element to the file
            fclose(fp);  // close the file

            fp = fopen("./simulation_output/electron_data.txt", "a");  // open the file in write mode
            if (fp == NULL) {
                printf("Error: could not open file.\n");
                return 1;
            }

            fprintf(fp,"Time: %lf. e- density: %lf. e- temperature: %lf. debye_length: %lf. rg_factor: %lf.\n",
                t,
                electron_density,
                electron_temperature,
                debye_length,
                rg_factor
                );  // write each element to the file

            fclose(fp);  // close the file
        }
    } // IONIZE IF STATEMENT

///////////////////////////////////
////                           ////
////  End of electron dynamics ////
////                           ////
///////////////////////////////////        
        


        /* Update force field in ffscan program */
        if (bFFscan)
        {
            if (update_forcefield(fplog,
                                  nfile,fnm,fr,
                                  mdatoms->nr,state->x,state->box)) {
                if (gmx_parallel_env_initialized())
                {
                    gmx_finalize();
                }
                exit(0);
            }
        }

        GMX_MPE_LOG(ev_timestep2);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
        bCPT = (((gs.set[eglsCHKPT] && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step && !bRerunMD);
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
        }

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        bNstEner = do_per_step(step,ir->nstcalcenergy);
        bCalcEnerPres =
            (bNstEner ||
             (ir->epc != epcNO && do_per_step(step,ir->nstpcouple)));

        /* Do we need global communication ? */
        bGStat = (bCalcEnerPres || bStopCM ||
                  do_per_step(step,nstglobalcomm) ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcEnerPres = TRUE;
            bGStat        = TRUE;
        }
        
        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                      (bStopCM ? CGLO_STOPCM : 0) |
                      (bGStat ? CGLO_GSTAT : 0)
            );
        
        force_flags = (GMX_FORCE_STATECHANGED |
                       ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bNStList ? GMX_FORCE_DOLR : 0) |
                       GMX_FORCE_SEPLRF |
                       (bCalcEnerPres ? GMX_FORCE_VIRIAL : 0) |
                       (bDoDHDL ? GMX_FORCE_DHDL : 0)
            );
        
        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count=relax_shell_flexcon(fplog,cr,bVerbose,bFFscan ? step+1 : step,
                                      ir,bNS,force_flags,
                                      bStopCM,top,top_global,
                                      constr,enerd,fcd,
                                      state,f,force_vir,mdatoms,
                                      nrnb,wcycle,graph,groups,
                                      shellfc,fr,bBornRadii,t,mu_tot,
                                      state->natoms,&bConverged,vsite,
                                      outf->fp_field);
            tcount+=count;

            if (bConverged)
            {
                nconverged++;
            }
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too. 
             * Check comments in sim_util.c
             */
        
            do_force(fplog,cr,ir,step,nrnb,wcycle,top,top_global,groups,
                     state->box,state->x,&state->hist,
                     f,force_vir,mdatoms,enerd,fcd,
                     state->lambda,graph,
                     fr,vsite,mu_tot,t,outf->fp_field,ed,bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }
    
        GMX_BARRIER(cr->mpi_comm_mygroup);
        
        if (bTCR)
        {
            mu_aver = calc_mu_aver(cr,state->x,mdatoms->chargeA,
                                   mu_tot,&top_global->mols,mdatoms,gnx,grpindex);
        }
        
        if (bTCR && bFirstStep)
        {
            tcr=init_coupling(fplog,nfile,fnm,cr,fr,mdatoms,&(top->idef));
            fprintf(fplog,"Done init_coupling\n"); 
            fflush(fplog);
        }
        
        if (bVV && !bStartingFromCpt && !bRerunMD)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            if (ir->eI==eiVV && bInitStep) 
            {
                /* if using velocity verlet with full time step Ekin,
                 * take the first half step only to compute the 
                 * virial for the first step. From there,
                 * revert back to the initial coordinates
                 * so that the input is actually the initial step.
                 */
                copy_rvecn(state->v,cbuf,0,state->natoms); /* should make this better for parallelizing? */
            } else {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ1);            
            }

            update_coords(fplog,step,ir,mdatoms,state,
                          f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                          ekind,M,wcycle,upd,bInitStep,etrtVELOCITY1,
                          cr,nrnb,constr,&top->idef);
            
            if (bIterations)
            {
                gmx_iterate_init(&iterate,bIterations && !bInitStep);
            }
            /* for iterations, we save these vectors, as we will be self-consistently iterating
               the calculations */

            /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */
            
            /* save the state */
            if (bIterations && iterate.bIterate) { 
                copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
            }
            
            bFirstIterate = TRUE;
            while (bFirstIterate || (bIterations && iterate.bIterate))
            {
                if (bIterations && iterate.bIterate) 
                {
                    copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
                    if (bFirstIterate && bTrotter) 
                    {
                        /* The first time through, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */
                        
                        veta_save = state->veta;
                        trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ0);
                        vetanew = state->veta;
                        state->veta = veta_save;
                    } 
                } 
                
                bOK = TRUE;
                if ( !bRerunMD || rerun_fr.bV || bForceUpdate) {  /* Why is rerun_fr.bV here?  Unclear. */
                    dvdl = 0;
                    
                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,shake_vir,NULL,
                                       cr,nrnb,wcycle,upd,constr,
                                       bInitStep,TRUE,bCalcEnerPres,vetanew);
                    
                    if (!bOK && !bFFscan)
                    {
                        gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }
                    
                } 
                else if (graph)
                { /* Need to unshift here if a do_force has been
                     called in the previous step */
                    unshift_self(graph,state->box,state->x);
                }
                
                
                /* if VV, compute the pressure and constraints */
                /* For VV2, we strictly only need this if using pressure
                 * control, but we really would like to have accurate pressures
                 * printed out.
                 * Think about ways around this in the future?
                 * For now, keep this choice in comments.
                 */
                /*bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir)); */
                    /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));*/
                bPres = TRUE;
                bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK));
                compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                constr,NULL,FALSE,state->box,
                                top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                cglo_flags 
                                | CGLO_ENERGY 
                                | (bTemp ? CGLO_TEMPERATURE:0) 
                                | (bPres ? CGLO_PRESSURE : 0) 
                                | (bPres ? CGLO_CONSTRAINT : 0)
                                | ((bIterations && iterate.bIterate) ? CGLO_ITERATE : 0)  
                                | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                | CGLO_SCALEEKIN 
                    );
                /* explanation of above: 
                   a) We compute Ekin at the full time step
                   if 1) we are using the AveVel Ekin, and it's not the
                   initial step, or 2) if we are using AveEkin, but need the full
                   time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                   b) If we are using EkinAveEkin for the kinetic energy for the temperture control, we still feed in 
                   EkinAveVel because it's needed for the pressure */
                
                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (!bInitStep) 
                {
                    if (bTrotter)
                    {
                        trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ2);
                    } 
                    else 
                    {
                        update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
                    }
                }
                
                if (bIterations &&
                    done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                                   state->veta,&vetanew)) 
                {
                    break;
                }
                bFirstIterate = FALSE;
            }

            if (bTrotter && !bInitStep) {
                copy_mat(shake_vir,state->svir_prev);
                copy_mat(force_vir,state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI==eiVV) {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,NULL,(ir->eI==eiVV),FALSE,FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (bInitStep && ir->eI==eiVV) {
                copy_rvecn(cbuf,state->v,0,state->natoms);
            }
            
            if (fr->bSepDVDL && fplog && do_log) 
            {
                fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
            }
            enerd->term[F_DHDL_CON] += dvdl;
            
            GMX_MPE_LOG(ev_timestep1);
        }
    
        /* MRS -- now done iterating -- compute the conserved quantity */
        if (bVV) {
            saved_conserved_quantity = compute_conserved_from_auxiliary(ir,state,&MassQ);
            if (ir->eI==eiVV) 
            {
                last_ekin = enerd->term[F_EKIN]; /* does this get preserved through checkpointing? */
            }
            if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres)) 
            {
                saved_conserved_quantity -= enerd->term[F_DISPCORR];
            }
        }
        
        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        
        /* ################## START TRAJECTORY OUTPUT ################# */
        
        /* Now we have the energies and forces corresponding to the 
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        mdof_flags = 0;
        if (do_per_step(step,ir->nstxout)) { mdof_flags |= MDOF_X; }
        if (do_per_step(step,ir->nstvout)) { mdof_flags |= MDOF_V; }
        if (do_per_step(step,ir->nstfout)) { mdof_flags |= MDOF_F; }
        if (do_per_step(step,ir->nstxtcout)) { mdof_flags |= MDOF_XTC; }
        if (bCPT) { mdof_flags |= MDOF_CPT; };

#if defined(GMX_FAHCORE) || defined(GMX_WRITELASTSTEP)
        if (bLastStep)
        {
            /* Enforce writing positions and velocities at end of run */
            mdof_flags |= (MDOF_X | MDOF_V);
        }
#endif
#ifdef GMX_FAHCORE
        if (MASTER(cr))
            fcReportProgress( ir->nsteps, step );

        /* sync bCPT and fc record-keeping */
        if (bCPT && MASTER(cr))
            fcRequestCheckPoint();
#endif
        
        if (mdof_flags != 0)
        {
            wallcycle_start(wcycle,ewcTRAJ);
            if (bCPT)
            {
                if (state->flags & (1<<estLD_RNG))
                {
                    get_stochd_state(upd,state);
                }
                if (MASTER(cr))
                {
                    if (bSumEkinhOld)
                    {
                        state_global->ekinstate.bUpToDate = FALSE;
                    }
                    else
                    {
                        update_ekinstate(&state_global->ekinstate,ekind);
                        state_global->ekinstate.bUpToDate = TRUE;
                    }
                    update_energyhistory(&state_global->enerhist,mdebin);
                }
            }
            write_traj(fplog,cr,outf,mdof_flags,top_global,
                       step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
            if (bCPT)
            {
                nchkpt++;
                bCPT = FALSE;
            }
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj,
                 * because a checkpoint file will always be written
                 * at the last step.
                 */
                fprintf(stderr,"\nWriting final coordinates.\n");
                if (ir->ePBC != epbcNONE && !ir->bPeriodicMols &&
                    DOMAINDECOMP(cr))
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog,ir->ePBC,state->box,top_global,state_global->x);
                }
                write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                                    *top_global->name,top_global,
                                    state_global->x,state_global->v,
                                    ir->ePBC,state->box);
                debug_gmx();
            }
            wallcycle_stop(wcycle,ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);
        
        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV) 
        {
            copy_mat(state->svir_prev,shake_vir);
            copy_mat(state->fvir_prev,force_vir);
        }
        /*  ################## END TRAJECTORY OUTPUT ################ */
        
        /* Determine the wallclock run time up till now */
        run_time = gmx_gettime() - (double)runtime->real;

        /* Check whether everything is still allright */    
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREADS
            && MASTER(cr)
#endif
            )
        {
            /* this is just make gs.sig compatible with the hack 
               of sending signals around by MPI_Reduce with together with
               other floats */
            if ( gmx_get_stop_condition() == gmx_stop_cond_next_ns )
                gs.sig[eglsSTOPCOND]=1;
            if ( gmx_get_stop_condition() == gmx_stop_cond_next )
                gs.sig[eglsSTOPCOND]=-1;
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition=(int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.sig[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            run_time > max_hours*60.0*60.0*0.495)
        {
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        if (ir->nstlist == -1 && !bRerunMD)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,
                                                     nlh.scale_tot,state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 || 
                            run_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }
  
        if (bIterations)
        {
            gmx_iterate_init(&iterate,bIterations);
        }
    
        /* for iterations, we save these vectors, as we will be redoing the calculations */
        if (bIterations && iterate.bIterate) 
        {
            copy_coupling_state(state,bufstate,ekind,ekind_save,&(ir->opts));
        }
        bFirstIterate = TRUE;
        while (bFirstIterate || (bIterations && iterate.bIterate))
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */    
            if (bIterations) 
            { 
                copy_coupling_state(bufstate,state,ekind_save,ekind,&(ir->opts));
            }

            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */
            
            /* #########   START SECOND UPDATE STEP ################# */
            GMX_MPE_LOG(ev_update_start);
            /* Box is changed in update() when we do pressure coupling,
             * but we should still use the old box for energy corrections and when
             * writing it to the energy file, so it matches the trajectory files for
             * the same timestep above. Make a copy in a separate array.
             */
            copy_mat(state->box,lastbox);

            bOK = TRUE;
            if (!(bRerunMD && !rerun_fr.bV && !bForceUpdate))
            {
                wallcycle_start(wcycle,ewcUPDATE);
                dvdl = 0;
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
                if (bTrotter) 
                {
                    if (bIterations && iterate.bIterate) 
                    {
                        if (bFirstIterate) 
                        {
                            scalevir = 1;
                        }
                        else 
                        {
                            /* we use a new value of scalevir to converge the iterations faster */
                            scalevir = tracevir/trace(shake_vir);
                        }
                        msmul(shake_vir,scalevir,shake_vir); 
                        m_add(force_vir,shake_vir,total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ3);
                /* We can only do Berendsen coupling after we have summed
                 * the kinetic energy or virial. Since the happens
                 * in global_state after update, we should only do it at
                 * step % nstlist = 1 with bGStatEveryStep=FALSE.
                 */
                }
                else 
                {
                    update_tcouple(fplog,step,ir,state,ekind,wcycle,upd,&MassQ,mdatoms);
                    update_pcouple(fplog,step,ir,state,pcoupl_mu,M,wcycle,
                                   upd,bInitStep);
                }

                if (bVV)
                {
                    /* velocity half-step update */
                    update_coords(fplog,step,ir,mdatoms,state,f,
                                  fr->bTwinRange && bNStList,fr->f_twin,fcd,
                                  ekind,M,wcycle,upd,FALSE,etrtVELOCITY2,
                                  cr,nrnb,constr,&top->idef);
                }

                /* Above, initialize just copies ekinh into ekin,
                 * it doesn't copy position (for VV),
                 * and entire integrator for MD.
                 */
                
                if (ir->eI==eiVVAK) 
                {
                    copy_rvecn(state->x,cbuf,0,state->natoms);
                }
                
                update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                              ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                wallcycle_stop(wcycle,ewcUPDATE);

                update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                   &top->idef,shake_vir,force_vir,
                                   cr,nrnb,wcycle,upd,constr,
                                   bInitStep,FALSE,bCalcEnerPres,state->veta);  
                
                if (ir->eI==eiVVAK) 
                {
                    /* erase F_EKIN and F_TEMP here? */
                    /* just compute the kinetic energy at the half step to perform a trotter step */
                    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                    constr,NULL,FALSE,lastbox,
                                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                    cglo_flags | CGLO_TEMPERATURE    
                        );
                    wallcycle_start(wcycle,ewcUPDATE);
                    trotter_update(ir,step,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq,ettTSEQ4);            
                    /* now we know the scaling, we can compute the positions again again */
                    copy_rvecn(cbuf,state->x,0,state->natoms);

                    update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                                  ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                    wallcycle_stop(wcycle,ewcUPDATE);

                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure. 
                     * For now, will call without actually constraining, constr=NULL*/
                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,tmp_vir,force_vir,
                                       cr,nrnb,wcycle,upd,NULL,
                                       bInitStep,FALSE,bCalcEnerPres,
                                       state->veta);  
                }
                if (!bOK && !bFFscan) 
                {
                    gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                }
                
                if (fr->bSepDVDL && fplog && do_log) 
                {
                    fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
                }
                enerd->term[F_DHDL_CON] += dvdl;
            } 
            else if (graph) 
            {
                /* Need to unshift here */
                unshift_self(graph,state->box,state->x);
            }
            
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_update_finish);

            if (vsite != NULL) 
            {
                wallcycle_start(wcycle,ewcVSITECONSTR);
                if (graph != NULL) 
                {
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
                
                if (graph != NULL) 
                {
                    unshift_self(graph,state->box,state->x);
                }
                wallcycle_stop(wcycle,ewcVSITECONSTR);
            }
            
            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints ############ */
            if (ir->nstlist == -1 && bFirstIterate)
            {
                gs.sig[eglsNABNSB] = nlh.nabnsb;
            }
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                            constr,
                            bFirstIterate ? &gs : NULL, 
                            (step_rel % gs.nstms == 0) && 
                                (multisim_nsteps<0 || (step_rel<multisim_nsteps)),
                            lastbox,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            cglo_flags 
                            | (!EI_VV(ir->eI) ? CGLO_ENERGY : 0) 
                            | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0) 
                            | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0) 
                            | (bIterations && iterate.bIterate ? CGLO_ITERATE : 0) 
                            | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                            | CGLO_CONSTRAINT 
                );
            if (ir->nstlist == -1 && bFirstIterate)
            {
                nlh.nabnsb = gs.set[eglsNABNSB];
                gs.set[eglsNABNSB] = 0;
            }
            /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
            /* #############  END CALC EKIN AND PRESSURE ################# */
        
            /* Note: this is OK, but there are some numerical precision issues with using the convergence of
               the virial that should probably be addressed eventually. state->veta has better properies,
               but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
               generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

            if (bIterations && 
                done_iterating(cr,fplog,step,&iterate,bFirstIterate,
                               trace(shake_vir),&tracevir)) 
            {
                break;
            }
            bFirstIterate = FALSE;
        }

        update_box(fplog,step,ir,mdatoms,state,graph,f,
                   ir->nstlist==-1 ? &nlh.scale_tot : NULL,pcoupl_mu,nrnb,wcycle,upd,bInitStep,FALSE);
        
        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */
    
        /* The coordinates (x) were unshifted in update */
        if (bFFscan && (shellfc==NULL || bConverged))
        {
            if (print_forcefield(fplog,enerd->term,mdatoms->homenr,
                                 f,NULL,xcopy,
                                 &(top_global->mols),mdatoms->massT,pres))
            {
                if (gmx_parallel_env_initialized())
                {
                    gmx_finalize();
                }
                fprintf(stderr,"\n");
                exit(0);
            }
        }
        if (!bGStat)
        {
            /* We will not sum ekinh_old,                                                            
             * so signal that we still have to do it.                                                
             */
            bSumEkinhOld = TRUE;
        }
        
        if (bTCR)
        {
            /* Only do GCT when the relaxation of shells (minimization) has converged,
             * otherwise we might be coupling to bogus energies. 
             * In parallel we must always do this, because the other sims might
             * update the FF.
             */

            /* Since this is called with the new coordinates state->x, I assume
             * we want the new box state->box too. / EL 20040121
             */
            do_coupling(fplog,oenv,nfile,fnm,tcr,t,step,enerd->term,fr,
                        ir,MASTER(cr),
                        mdatoms,&(top->idef),mu_aver,
                        top_global->mols.nr,cr,
                        state->box,total_vir,pres,
                        mu_tot,state->x,f,bConverged);
            debug_gmx();
        }

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */
        
        /* sum up the foreign energy and dhdl terms */
        sum_dhdl(enerd,state->lambda,ir);

        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && ir->eI==eiVV) 
        {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];
        
        if (bVV)
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
        }
        else 
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + compute_conserved_from_auxiliary(ir,state,&MassQ);
        }
        /* Check for excessively large energies */
        if (bIonize) 
        {
#ifdef GMX_DOUBLE
            real etot_max = 1e200;
#else
            real etot_max = 1e30;
#endif
            if (fabs(enerd->term[F_ETOT]) > etot_max) 
            {
                fprintf(stderr,"Energy too large (%g), giving up\n",
                        enerd->term[F_ETOT]);
            }
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */
        
        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep) 
        {
            runtime_upd_proc(runtime);
        }
        
        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool do_dr,do_or;
            
            if (!(bStartingFromCpt && (EI_VV(ir->eI)))) 
            {
                if (bNstEner)
                {
                    upd_mdebin(mdebin,bDoDHDL, TRUE,
                               t,mdatoms->tmass,enerd,state,lastbox,
                               shake_vir,force_vir,total_vir,pres,
                               ekind,mu_tot,constr);
                }
                else
                {
                    upd_mdebin_step(mdebin);
                }
                
                do_dr  = do_per_step(step,ir->nstdisreout);
                do_or  = do_per_step(step,ir->nstorireout);
                
                print_ebin(outf->fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,
                           step,t,
                           eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
            }
            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull,step,t);
            }
            
            if (do_per_step(step,ir->nstlog))
            {
                if(fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
                }
            }
        }


        /* Remaining runtime */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal() ))
        {
            if (shellfc) 
            {
                fprintf(stderr,"\n");
            }
            print_time(stderr,runtime,step,ir,cr);
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step,repl_ex_nst)) 
        {
            bExchanged = replica_exchange(fplog,cr,repl_ex,
                                          state_global,enerd->term,
                                          state,step,t);

            if (bExchanged && DOMAINDECOMP(cr)) 
            {
                dd_partition_system(fplog,step,cr,TRUE,1,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,FALSE);
            }
        }
        
        bFirstStep = FALSE;
        bInitStep = FALSE;
        bStartingFromCpt = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & (1<<estPRES_PREV)) &&
            (bGStatEveryStep ||
             (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres,state->pres_prev);
        }
        
        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */
        
        if (bRerunMD) 
        {
            if (MASTER(cr))
            {
                /* read next frame from input trajectory */
                bNotLastFrame = read_next_frame(oenv,status,&rerun_fr);
            }

            if (PAR(cr))
            {
                rerun_parallel_comm(cr,&rerun_fr,&bNotLastFrame);
            }
        }
        
        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }
        
        cycles = wallcycle_stop(wcycle,ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd,cycles,ddCyclStep);
        }
        
        if (step_rel == wcycle_get_reset_counters(wcycle) ||
            gs.set[eglsRESETCOUNTERS] != 0)
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog,cr,step,&step_rel,ir,wcycle,nrnb,runtime);
            wcycle_set_reset_counters(wcycle,-1);
            /* Correct max_hours for the elapsed time */
            max_hours -= run_time/(60.0*60.0);
            bResetCountersHalfMaxH = FALSE;
            gs.set[eglsRESETCOUNTERS] = 0;
        }


    }
    /* End of main MD loop */
    if (bIonize) {
        int n = top_global->natoms; 

        // Write final atom configurations
        FILE* fp = fopen("configurations.bin", "wb");
        if (fp == NULL) {
            printf("Error: could not open file for writing.\n");
            exit(1);
        }

        for (i = 0; i < n; i++) {
            fwrite(atom_configurations[i], sizeof(int), 3, fp);  // write each subarray
        }
        fclose(fp);

        // Write the charges of the atoms
        fp = fopen("charges.bin", "wb");
        if (fp == NULL) {
            printf("Error: could not open file for writing.\n");
            exit(1);
        }

        for (i = 0; i < n; i++) {
            fwrite(&mdatoms->chargeA[i], sizeof(double), 1, fp);  // write each charge
        }
        fclose(fp);


        // Free each sub-array
        for (i = 0; i<top_global->natoms; i++){
            free(atom_configurations[i]);
            free(GS_configurations[i]);
            free(atomic_transitions[i]);
        }
        free(atom_configurations);
        free(GS_configurations);
        free(atomic_transitions);

    }




    debug_gmx();
    
    /* Stop the time */
    runtime_end(runtime);
    
    if (bRerunMD && MASTER(cr))
    {
        close_trj(status);
    }
    
    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_finish(cr);
    }
    
    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0 && !bRerunMD) 
        {
            print_ebin(outf->fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
        }
    }

    done_mdoutf(outf);

    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",nlh.s1/nlh.nns,sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",nlh.ab/nlh.nns);
    }
    
    if (shellfc && fplog)
    {
        fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }
    
    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog,repl_ex);
    }
    
    runtime->nsteps_done = step_rel;

    // FREE STUFF
    free(hydrogen_idx);


    // FREE DICTIONARIES 
    freeDictionary(Energylevels_H);
    freeDictionary(Energylevels_C);
    freeDictionary(Energylevels_N);
    freeDictionary(Energylevels_O);
    freeDictionary(Energylevels_S);
    freeDictionary(Energylevels_FE);

    // FREE RATE ARRAYS
    freeRatesArray(transitionRates_H,numRates_H);
    freeRatesArray(transitionRates_C,numRates_C);
    freeRatesArray(transitionRates_N,numRates_N);
    freeRatesArray(transitionRates_O,numRates_O);
    freeRatesArray(transitionRates_S,numRates_S);
    freeRatesArray(transitionRates_FE,numRates_FE);

    // Free Coll_rates
    freeCollArray(collisions_H,numColl_H);
    freeCollArray(collisions_C,numColl_C);
    freeCollArray(collisions_N,numColl_N);
    freeCollArray(collisions_O,numColl_O);
    freeCollArray(collisions_S,numColl_S);
    freeCollArray(collisions_FE,numColl_FE);


    // Free Weights
    freeWeights(Weigths_H);
    freeWeights(Weigths_C);
    freeWeights(Weigths_N);
    freeWeights(Weigths_O);
    freeWeights(Weigths_S);
    freeWeights(Weigths_FE);


    return 0;
}
