/*
MODEL.H file for the TWOFLDLAMBDA model.
Written by Gary Felder (gfelder@email.smith.edu)
Last Modified 3/7/01

This file contains the model specific functions and definitions for the model:
V = 1/4 lambda phi^4 + 1/2 g^2 phi^2 chi^2
Note that this model is simply a special case of the NFLDLAMBDA model. Choosing nflds=2 for that model is equivalent to running this one.
Comments: The inflaton potential is taken to be pure lambda phi^4. The inflaton decays (generally resonantly) to the single field chi. Neither of the fields have bare masses and the theory is therefore conformal. In practice this means that after appropriate rescaling of the time, space, and field variables the only effect of the expansion of the universe on the field equations is via the term a''/a, which rapidly becomes negligible after inflation.

We use the default rescalings (see below), which give A=1/f0, B=sqrt(lambda) f0, r=1, s=-1
The potential in program units is given by V_{pr} = 1/4 phi_{pr}^4 + 1/2 gl phi_{pr}^2 chi_{pr}^2
*/

/*
General comments about the model.h file:
This file contains the following functions - all called externally
modelinfo(FILE *info_) outputs information about the model and model-specific parameters to a file.
modelinitialize() performs any model-specific initialization
potential_energy(int term, float *field_values) calculates the average potential energy density term by term. The variable num_potential_terms (just above this function) specifies how many separate potential terms are used in this model.
dvdf(int fld, int i, int j, int k) calculates the potential term in the equation of motion, dV/dfield, for the field fld at the lattice point (i,j,k)
effective_mass(float mass_sq[], float *field_values) calculates the square masses of the fields and puts them in the array mass_sq. The parameter beginning tells the function to use initial field values - if this parameter is zero then the field quantities will be calculated dynamically.
model_output(int flush, char *ext_) allows each model to include its own specialized output function(s). The parameter flush is set to 1 when infrequent calculations are being performed and 0 otherwise. The string ext_ gives the extension for output filenames.
*/

/* - This section should be copied into parameters.h when using this model
// ---Adjustable parameters for TWOFLDLAMBDA model--- //
const float lambda=9.e-14; // Self-coupling of inflaton
const float gl=200; // Resonance parameter g^2/lambda
*/

// Rescaling parameters.
// The program variables ("_pr") are defined as
//   f_pr = rescale_A a^rescale_r f (f=field value)
//   x_pr = rescale_B x (x=distance)
//   dt_pr = rescale_B a^rescale_s dt (t=time)
// The constants beta, cpl, and f0 are used to set the variable rescalings rescale_A, rescale_B, rescale_r, and rescale_s.
// These rescaling constants may be reset independently; their settings in terms of beta, cpl, and f0 are suggestions. See the documentation for more details.
// These rescalings are intrinsic to the model and probably shouldn't be changed for individual runs. Adjustable parameters are stored in "parameters.h"
const float beta=4.; // Exponent of the dominant term in the potential
const float cpl=lambda; // Coefficient of the dominant term in the potential (up to numerical factors - see documentation)
const float f0=0.342; // Initial value of phi in Planck units, typically the point at which phi'=0
// By default these are automatically set to A=1/f0, B=sqrt(cpl) f0^(-1+beta/2), R=6/(2+beta), S=3(2-beta)/(2+beta). They may be adjusted to different values, but the relationship S=2R-3 must be maintained for the program equations to remain correct.
const float rescale_A=1./f0;
const float rescale_B=sqrt(cpl)*pow(f0,-1.+beta/2.);
const float rescale_r=6./(2.+beta);
// The value of S in terms of R SHOULD NOT be changed.
const float rescale_s=2.*rescale_r-3.;

// Other global variables
// The array model_vars is intended to hold any model-specific, non-constant global variables.
// These variables should be initialized in modelinitialize() below
// Even if you're not using any, num_model_vars should be at least 1 to keep some compilers happy.
const int num_model_vars=1;
// Model specific variables: None are defined for this model.
extern float model_vars[num_model_vars];

// Macros to make the equations more readable: The values of fld are 0=Phi,1=Chi
#define PHI FIELD(0)
#define CHI FIELD(1)

// Model specific details about the run to be output to an information file
inline void modelinfo(FILE *info_)
{
  // Name and description of model
  fprintf(info_,"Two Field Lambda Model\n");
  fprintf(info_,"V = 1/4 lambda phi^4 + 1/2 g^2 phi^2 chi^2\n\n");

  // Model specific parameter values
  fprintf(info_,"gl (g^2/lambda) = %f\n",gl);
  fprintf(info_,"lambda = %e\n",lambda);

  // Rescalings
  fprintf(info_,"rescale_A = %e\n",rescale_A);
  fprintf(info_,"rescale_B = %e\n",rescale_B);
  fprintf(info_,"rescale_r = %e\n",rescale_r);
  fprintf(info_,"rescale_s = %e\n",rescale_s);
}

// Perform any model specific initialization
// This function is called twice, once before initializing the fields (which_call=1) and once after (which_call=2)
inline void modelinitialize(int which_call)
{
  if(which_call==1)
  {
    if(nflds!=2)
    {
      printf("Number of fields for TWOFLDLAMBDA model must be 2. Exiting.\n");
      exit(1);
    }
  }
}

// The constant num_potential_terms must be defined for use by outside functions
// Terms: term=0: 1/4 lambda phi^4 --- term=1: 1/2 g^2 phi^2 chi^2
const int num_potential_terms=2; // Number of terms that are calculated separately in the potential
// Potential energy terms
// See documentation for normalization of these terms.
// When setting initial conditions field values will be supplied in the array field_values. Otherwise this function will calculate them on the lattice.
inline float potential_energy(int term, float *field_values)
{
  DECLARE_INDICES
  float potential=0.;

  if(field_values==NULL) // If no values are given calculate averages on the lattice
  {
    // Loop over grid to calculate potential term
    LOOP
    {
      if(term==0)
        potential += pw2(pw2(PHI));
      else if(term==1)
        potential += pw2(PHI*CHI);
    }

    // Convert sum to average
    potential /= gridsize;
  }
  else // If field values are given then use them instead
  {
    if(term==0)
      potential = pw2(pw2(field_values[0]));
    else if(term==1)
      potential = pw2(field_values[0]*field_values[1]);
  }

  // Include numerical coefficients
  if(term==0) // 1/4 lambda phi^4
    potential *= .25;
  else if(term==1) // 1/2 g^2 phi^2 chi^2 
    potential *= .5*gl;

  return (potential);
}

// Potential terms in the equations of motion, dV/dfield, evaluated at point (i,j,k)
// See documentation for details on the normalization of these terms
inline float dvdf(int fld, INDEXLIST)
{
  if(fld==0) // Phi
    return( (pw2(PHI) + gl*pw2(CHI))*PHI );
  else // Chi
    return( gl*pw2(PHI)*CHI );
}

// Calculate effective mass squared and put it into the array mass_sq[] (used for initial conditions and power spectra)
// See documentation for normalization of these terms
// When setting initial conditions field values will be supplied in the array field_values. Otherwise this function will calculate them on the lattice.
inline void effective_mass(float mass_sq[], float *field_values)
{
  DECLARE_INDICES
  int fld;
  float fldsqrd[nflds]; // Square value of field
  float correction; // Used to adjust masses by the appropriate power of the scale factor. (See documentation.)

  // Loop over fields to find mean-square value
  if(field_values==NULL) // If no values are given calculate averages on the lattice
  {
    for(fld=0;fld<nflds;fld++)
    {
      fldsqrd[fld]=0.;
      LOOP
        fldsqrd[fld]+=pw2(FIELD(fld));
      fldsqrd[fld] /= (float)gridsize;
    }
  }
  else // If field values are given then use them instead
    for(fld=0;fld<nflds;fld++)
      fldsqrd[fld]=pw2(field_values[fld]);

  mass_sq[0] = 3.*fldsqrd[0] + gl*fldsqrd[1];
  mass_sq[1] = gl*fldsqrd[0];

  // Put in scale factor correction. This calculation should be the same for all models.
  if(expansion>0) // If there's no expansion don't bother with this.
  {
    correction = pow(a,2.*rescale_s+2.);
    for(fld=0;fld<nflds;fld++)
      mass_sq[fld] *= correction;
  }
}

// Model-specific output functions
inline void model_output(int flush,char *ext_){}

#undef PHI
#undef CHI
