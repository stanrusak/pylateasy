/*
This file contains all the adjustable parameters for running the latticeeasy program. Once a particular model has been defined all runs can be done by changing the parameters in this file and recompiling. See the documentation for how to define new models.
*/

// Model specific adjustable parameters - This section should be replaced by the appropriate one for the model being used.
// (These definitions can be left in model.h but by putting them here you have all the adjustable parameters in one file.)
// ---Adjustable parameters for TWOFLDLAMBDA model--- //
const float lambda=9.e-14; // Self-coupling of inflaton
const float gl=300.0; // Resonance parameter g^2/lambda

// Adjustable run parameters
#define NDIMS 2
const int N=256; // Number of points along each edge of the cubical lattice
const int nflds=2;  //Number of fields
const float L=20.; // Size of box (i.e. length of each edge) in rescaled distance units
const float dt=.02; // Size of time step
const float tf=500.0; // Final time
const int seed=1; // Random number seed. Should be a positive integer
const float initfield[]={1.}; // Initial values of the fields in program units. All nonspecified values are taken to be zero.
const float initderivs[]={0.}; // Initial values of the field derivatives in program units. All nonspecified values are taken to be zero.
const int expansion=2; // Whether to use no expansion (0), power-law expansion (1), or self-consistent expansion (2)
const float expansion_power=.5; // Power of t in power law expansion. Only used when expansion=1.. Set to .5 for radiation or .67 for matter domination.
const float kcutoff=10.0; // Momentum for initial lowpass filter. Set to 0 to not filter

// If and how to continue previous runs.
// If no grid image is available in the run directory then a new run will be started irrespective of continue_run.
// 0=Start a new run at t=0. 1=Continue old run, appending new data to old data files. 2=Continue old run, creating new data files for all output. (Old ones will not in general be overwritten.)
const int continue_run=0;

// Variables controlling output
const int noutput_flds=0; // Number of fields to output information about. All fields will be output nflds if noutput=0 or noutput>nflds
const char alt_extension[]=""; // Optional alternative extension for output files (default is the run number)
const int noutput_times=20000; // Number of times at which to calculate and save output variables
const int print_interval=1; // Interval in seconds between successive outputs of time
const int screen_updates=1; // Set to 1 for time to be periodically output to screen (0 otherwise)
const float checkpoint_interval=0.5; // How often to output a grid image and perform infrequent calculations (see list below). Only done at end if checkpoint_interval=0.
const float store_lattice_times[]={0.}; // An optional list of times at which to close the grid image file and open a new one.
// The variables s<name> control what will be saved (1=save, 0=don't save)
const int smeansvars=1; // Output means and variances. This function is also used to check for exponential instability, so it is generally wise to leave it on.
const int sexpansion=1; // Output scale factor, Hubble constant, and a'' (This is ignored except in self-consistent expansion.)
const int smodel=1; // Call model-specific output functions. This must be on for the model file to set the rescaling for the other output functions.
// The following calculations are performed at intervals given by checkpoint_interval
const float t_start_output=0.; // Time to start doing these calculations. This can be reset for any individual calculation.
const int scheckpoint=1; // Save an image of the grid
  const float tcheckpoint=t_start_output;
const int sspectra=1; // Output power spectra
  const float tspectra=t_start_output;
const int senergy=1; // Output components of energy density
  const float tenergy=t_start_output;
const int shistograms=1; // Output histograms of fields using nbins as the number of bins
  const float thistograms=t_start_output;
  const int nbins=256; // Number of bins
  const float histogram_min=0., histogram_max=0.; // Upper and lower limits of the histograms. To use all current field values set these two to be equal.
const int shistograms2d=1; // Output two dimensional histograms of fields
  const float thistograms2d=t_start_output;
  const int nbins2d=10, nbins0=nbins2d, nbins1=nbins2d; // Number of bins in each field direction
  const float histogram2d_min=0., histogram2d_max=0.; // Upper and lower limits of the histograms. To use all current field values set these two to be equal.
  const int hist2dflds[]={0,1}; // Pairs of fields to be evaluated. This array should always contain an even number of integers.
const int sslices=0; // Output the field values on a slice through the lattice
  const float tslices=t_start_output;
  const int slicedim=2; // Dimensions of slice to be output. (If slicedim>=NDIMS the whole lattice will be output.) Warning: If slicedim=3 the resulting file may be very large.
  const int slicelength=N, sliceskip=2; // The slices will use every <sliceskip> point up to a total of <slicelength>.  Set length=N and skip=1 to output all points in the slice.
  const int sliceaverage=1; // If sliceskip>1 and sliceaverage=1 the slices will contain averages over all the field values in between sampled points.  Otherwise they will just use the sampled values.

