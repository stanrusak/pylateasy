/*
LATTICEEASY consists of the C++ files ``latticeeasy.cpp,''
``initialize.cpp,'' ``evolution.cpp,'' ``output.cpp,''
``latticeeasy.h,'' ``parameters.h,''. (The distribution also includes
the file ffteasy.cpp but this file is distributed separately and
therefore not considered part of the LATTICEEASY distribution in what
follows.) LATTICEEASY is free. We are not in any way, shape, or form
expecting to make money off of these routines. We wrote them for the
sake of doing good science and we're putting them out on the Internet
in case other people might find them useful. Feel free to download
them, incorporate them into your code, modify them, translate the
comment lines into Swahili, or whatever else you want. What we do want
is the following:
1) Leave this notice (i.e. this entire paragraph beginning with
``LATTICEEASY consists of...'' and ending with our email addresses) in
with the code wherever you put it. Even if you're just using it
in-house in your department, business, or wherever else we would like
these credits to remain with it. This is partly so that people can...
2) Give us feedback. Did LATTICEEASY work great for you and help
your work?  Did you hate it? Did you find a way to improve it, or
translate it into another programming language? Whatever the case
might be, we would love to hear about it. Please let us know at the
email address below.
3) Finally, insofar as we have the legal right to do so we forbid
you to make money off of this code without our consent. In other words
if you want to publish these functions in a book or bundle them into
commercial software or anything like that contact us about it
first. We'll probably say yes, but we would like to reserve that
right.

For any comments or questions you can reach us at
gfelder@email.smith.edu
Igor.Tkachev@cern.ch

Enjoy LATTICEEASY!

Gary Felder and Igor Tkachev
*/

#include "latticeeasy.h"

// Global variables used by various functions. They are all declared as extern in latticeeasy.h
// There are many more global constants declared as const in parameters.h;
//   these are the global variables that are dynamically changed by the program.
// Field values and derivatives
#if NDIMS==1
float f[nflds][N],fd[nflds][N];
#elif NDIMS==2
float f[nflds][N][N],fd[nflds][N][N];
#else
float f[nflds][N][N][N],fd[nflds][N][N][N];
#endif
float t,t0; // Current time and initial time (t0=0 unless the run is a continuation of a previous one)
float a=1.,ad=0.,ad2=0.,aterm=0.; // Scale factor and its derivatives (aterm is a combination of the others used in the equations of motion). Values are initialized to their defaults for the case of no expansion.
float hubble_init=0.; // Initial value of the Hubble constant
int run_number; // 0 for a first run, 1 for a continuation of a "0" run, etc.. Stored in the grid image (see checkpoint() function).
int no_initialization=0; // If this variable is set to 1 by the model file then the fields will not be initialized in the normal way.
char mode_[10]="w"; // Mode in which to open files, i.e. write ("w") or append ("a+"). Depends on the variable continue_run and on whether a previous grid image was found.
float rescaling=1.; // Rescaling for output. This is left as 1 unless the model file modifies it.
char ext_[500]="_0.dat"; // Extension for filenames - set once and used by all output functions
int nfldsout; // Number of fields to output
float model_vars[num_model_vars]; // Model-specific variables

int main()
{
  int numsteps=0,output_interval; // Quantities for counting how often to calculate and output derived quantities
  FILE *output_=fopen("output.txt","w"); // Outputs time. Used to remotely monitor progress
  int update_time; // Controls when to output time to output file and screen 

  if(seed<1) // The use of seed<1 turns off certain functions (random numbers, fourier transforms, gradients, and potential energy) and should only be used for debugging
    printf("Warning: The parameter seed has been set to %d, which will result in incorrect output. For correct output set seed to a positive integer.",seed);
  initialize(); // Set parameter values and initial conditions
  t=t0;
  output_interval = (int)((tf-t0)/dt)/noutput_times + 1; // Set the interval between saves

  // Take Initial Half Time Step if this is a new run
  if(run_number==0)
    evolve_fields(.5*dt);

  update_time=time(NULL)+print_interval; // Set initial time for update
  while(t<=tf) // Main time evolution loop
  {
    evolve_derivs(dt);
    evolve_fields(dt);
    numsteps++;
    if(numsteps%output_interval == 0 && t<tf)
      save(0); // Calculate and output grid-averaged quantities (means, variances, etc.)
    if(time(NULL) >= update_time) // Print an update whenever elapsed time exceeds print_interval
    {
      if(screen_updates) // This option determines whether or not to update progress on the screen
        printf("%f\n",t);
      fprintf(output_,"%f\n",t); // Output progress to a file for monitoring progress
      fflush(output_); // Make sure output file is always up to date
      update_time += print_interval; // Set time for next update
    }
  }

  printf("Saving final data\n");
  save(1); // Calculate and save quantities. Force infrequently calculated quantities to be calculated.
  output_parameters(); // Save run parameters and elapsed time
  fprintf(output_,"LATTICEEASY program finished\n");
  printf("LATTICEEASY program finished\n");

  return(0);
}
