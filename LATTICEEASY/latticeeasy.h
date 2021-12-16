/*
This file contains the global variable declarations, function declarations, 
and some definitions used in many of the routines. The global variables are 
defined in the file latticeeasy.cpp.
*/

#ifndef _LATTICEEASYHEADER_
#define _LATTICEEASYHEADER_

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

const float pi = (float)(2.*asin(1.));
inline float pw2(float x) {return x*x;} // Useful macro for squaring floats

/////////////////////////////////INCLUDE ADJUSTABLE PARAMETERS///////////////////
#include "parameters.h"

/////////////////////////////////GLOBAL DYNAMIC VARIABLES////////////////////////
extern float t,t0; // Current time and initial time (t0=0 unless the run is a continuation of a previous one)
extern float a,ad,ad2,aterm; // Scale factor and its derivatives (aterm is a combination of the others used in the equations of motion)
extern float hubble_init; // Initial value of the Hubble constant
extern int run_number; // 0 for a first run, 1 for a continuation of a "0" run, etc.. Stored in the grid image (see checkpoint() function).
extern int no_initialization; // If this variable is set to 1 by the model file then the fields will not be initialized in the normal way.
extern float rescaling; // Rescaling for output. This is left as 1 unless the model file modifies it.
extern char ext_[500]; // Extension for filenames - set once and used by all functions
extern int nfldsout; // Number of fields to output
extern char mode_[]; // Mode in which to open files, i.e. write ("w") or append ("a+"). Depends on the variable continue_run and on whether a previous grid image was found.

/////////////////////////////////NON-ADJUSTABLE VARIABLES////////////////////////
const float dx=L/(float)N; // Distance between adjacent gridpoints

/////////////////////////////////DIMENSIONAL SPECIFICATIONS//////////////////////
#if NDIMS==1
extern float f[nflds][N],fd[nflds][N]; // Field values and derivatives
const int gridsize=N; // Number of spatial points in the grid
#define FIELD(fld) f[fld][i]
#define FIELDD(fld) fd[fld][i]
#define FIELDPOINT(fld,i,j,k) f[fld][k]
#define LOOP for(i=0;i<N;i++)
#define INDEXLIST int i, ...
#define DECLARE_INDICES int i;
#elif NDIMS==2
extern float f[nflds][N][N],fd[nflds][N][N]; // Field values and derivatives
const int gridsize=N*N; // Number of spatial points in the grid
#define FIELD(fld) f[fld][i][j]
#define FIELDD(fld) fd[fld][i][j]
#define FIELDPOINT(fld,i,j,k) f[fld][j][k]
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++)
#define INDEXLIST int i, int j, ...
#define DECLARE_INDICES int i,j;
#elif NDIMS==3
extern float f[nflds][N][N][N],fd[nflds][N][N][N]; // Field values and derivatives
const int gridsize=N*N*N; // Number of spatial points in the grid
#define FIELD(fld) f[fld][i][j][k]
#define FIELDD(fld) fd[fld][i][j][k]
#define FIELDPOINT(fld,i,j,k) f[fld][i][j][k]
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
#define INDEXLIST int i, int j, int k
#define DECLARE_INDICES int i,j,k;
#endif

/////////////////////////////////INCLUDE SPECIFIC MODEL//////////////////////////
#include "model.h"

/////////////////////////////////FUNCTION DECLARATIONS///////////////////////////
// initialize.cpp
void initialize(); // Set initial parameters and field values
// evolution.cpp
float gradient_energy(int fld); // Calculate the gradient energy, <|Grad(f)|^2>=<-f Lapl(f)>, of a field
void evolve_scale(float d); // Calculate the scale factor and its derivatives
void evolve_fields(float d); // Advance the field values and scale factor using the first derivatives
void evolve_derivs(float d); // Calculate second derivatives of fields and use them to advance first derivatives. Also calls evolve_scale().
// output.cpp
void output_parameters(); // Output information about the run parameters
void save(int force); // Calculate and save quantities (means, variances, spectra, etc.)
// ffteasy.cpp
void fftr1(float f[], int N, int forward); // Do a Fourier transform of a 1D array of real numbers. Used when NDIMS=1.
void fftrn(float f[], float fnyquist[], int ndims, int size[], int forward); // Do a Fourier transform of an ndims dimensional array of real numbers

#endif // End of conditional for definition of _LATTICEEASYHEADER_ macro






