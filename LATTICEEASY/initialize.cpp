/*
This file contains the initialization function for setting field and parameter values.

The only function here that should be called externally is initialize(). It sets field values in momentum space and then transforms them to position space. It also sets parameters such as the scale factor and its derivatives.
*/

#include "latticeeasy.h"

// Generate a uniform deviate between 0 and 1 using the Park-Miller minimum standard algorithm
#define randa 16807
#define randm 2147483647
#define randq 127773
#define randr 2836
float rand_uniform(void)
{
  if(seed<1) return(0.33); // *DEBUG* This is used to avoid randomness, for debugging purposes only.
  static int i=0;
  static int next=seed;
  if(!(next>0)) // Guard against 0, negative, or other invalid seeds
  {
    printf("Invalid seed used in random number function. Using seed=1\n");
    next=1;
  }
  if(i==0) // On the first call run through 100 calls. This allows small seeds without the first results necessarily being small.
    for(i=1;i<100;i++)
      rand_uniform();
  next = randa*(next%randq) - randr*(next/randq);
  if(next<0) next += randm;
  return ((float)next/(float)randm);
}
#undef randa
#undef randm
#undef randq
#undef randr

// Set the amplitude and phase of a mode of vacuum fluctuations
// Phase is set randomly (except for modes that must be real).
// Amplitude is set with a Rayleigh distribution about an rms value of 1/sqrt(omega) (times normalization terms).
void set_mode(float p2, float m2, float *field, float *deriv, int real)
{
  float phase,amplitude,rms_amplitude,omega;
  float re_f_left,im_f_left,re_f_right,im_f_right;
#if NDIMS==1
  static float norm = rescale_A*rescale_B*pow(L/pw2(dx),.5)/sqrt(4.*pi);
#elif NDIMS==2
  static float norm = rescale_A*rescale_B*L/pw2(dx)/sqrt(2.*pi);
#elif NDIMS==3
  static float norm = rescale_A*rescale_B*pow(L/pw2(dx),1.5)/sqrt(2.);
#endif
  static float hbterm = hubble_init*(rescale_r-1.);
  static int tachyonic=0; // Avoid printing the same error repeatedly

  // Momentum cutoff. If kcutoff!=0 then eliminate all initial modes with k>kcutoff.
  static float k2cutoff = (kcutoff<2.*pi*(float)N/L ? pw2(kcutoff) : 0.);
  if(k2cutoff>0. && p2>k2cutoff)
  {
    field[0]=0.;
    field[1]=0.;
    deriv[0]=0.;
    deriv[1]=0.;
    return;
  }

  if(p2+m2>0.) // Check to avoid floating point errors
    omega=sqrt(p2+m2); // Omega = Sqrt(p^2 + m^2)
  else
  {
    if(tachyonic==0)
      printf("Warning: Tachyonic mode(s) may be initialized inaccurately\n");
    omega=sqrt(p2); // If p^2 + m^2 < 0 use m^2=0
    tachyonic=1;
  }
  if(omega>0.) // Avoid dividing by zero
    rms_amplitude=norm/sqrt(omega)*pow(p2,.75-(float)NDIMS/4.);
  else
    rms_amplitude=0.;

  // Amplitude = RMS amplitude x Rayleigh distributed random number
  // The same amplitude is used for left and right moving waves to generate standing waves. The extra 1/sqrt(2) normalizes the initial occupation number correctly.
  amplitude=rms_amplitude/sqrt(2.)*sqrt(log(1./rand_uniform()));
  // Left moving component
  phase=2.*pi*rand_uniform(); // Set phase randomly
  re_f_left = amplitude*cos(phase);
  im_f_left = amplitude*sin(phase);
  // Right moving component
  phase=2.*pi*rand_uniform(); // Set phase randomly
  re_f_right = amplitude*cos(phase);
  im_f_right = amplitude*sin(phase);

  field[0] = re_f_left + re_f_right; // Re(field)
  field[1] = im_f_left + im_f_right; // Im(field)
  deriv[0] = omega*(im_f_left - im_f_right) + hbterm*field[0]; // Field derivative
  deriv[1] = -omega*(re_f_left - re_f_right) + hbterm*field[1];
  if(real==1) // For real modes set the imaginary parts to zero
  {
    field[1] = 0.;
    deriv[1] = 0.;
  }

  return;
}

/////////////////////////////////////////////////////
// Externally called function(s)
/////////////////////////////////////////////////////

// Set initial parameters and field values
void initialize()
{
  int fld;
  float p2; // Total squared momentum
  float dp2=pw2(2.*pi/L); // Square of grid spacing in momentum space
  float mass_sq[nflds]; // Effective mass squared of fields
  float initial_field_values[nflds]; // Initial value of fields (set to zero unless specified in parameters.h)
  float initial_field_derivs[nflds]; // Initial value of field derivatives (set to zero unless specified in parameters.h)
  float fsquared=0.,fdsquared=0.,ffd=0.,pot_energy=0.; // Sum(field^2), Sum(field_dot^2), Sum(f*f_dot), and potential energy - used for calculating the Hubble constant
  FILE *old_grid_; // Used for reading in previously generated data as initial conditions
#if NDIMS==1
  int i,k;
  float pz; // Component of momentum in units of grid spacing
#elif NDIMS==2
  int i,j,k,jconj; // jconj is the index where conjugate modes to j are stored
  float py,pz; // Components of momentum in units of grid spacing
  float fnyquist[2*N],fdnyquist[2*N]; // Used by FFT routine to store modes with k=nyquist
  int arraysize[]={N,N}; // Array of grid size in all dimensions - used by FFT routine
#elif NDIMS==3
  int i,j,k,iconj,jconj; // iconj and jconj are the indices where conjugate modes to i and j are stored
  float px,py,pz; // Components of momentum in units of grid spacing
  float fnyquist[N][2*N],fdnyquist[N][2*N]; // Used by FFT routine to store modes with k=nyquist
  int arraysize[]={N,N,N}; // Array of grid size in all dimensions - used by FFT routine
#endif
  // Check to make sure time step is small enough to satisfy Courant condition, dt/dx < 1/Sqrt(ndims)
  if(dt>dx/sqrt((double)NDIMS))
  {
    printf("Time step too large. The ratio dt/dx is currently %f but for stability should never exceed 1/sqrt(%d) (%f)\n",dt/dx,NDIMS,1./sqrt((double)NDIMS));
    printf("Adjust dt to AT MOST %e, and preferably somewhat smaller than that.\n",dx/sqrt((double)NDIMS));
    exit(1);
  }

  modelinitialize(1); // This allows specific models to perform any needed initialization

  // Output initializations - Set values of nfldsout and ext_
  nfldsout = (noutput_flds==0 || noutput_flds>nflds ? nflds : noutput_flds); // Set number of files to output
  if(alt_extension[0]!='\0') // If an alternate extension was given use that instead of the default "_<run_number>.dat"
    sprintf(ext_,"%s",alt_extension);

  // If an old grid image can be opened use it to read in initial conditions
  if(continue_run>0 && (old_grid_=fopen("grid.img","rb")))
  {
    printf("Previously generated grid image found. Reading in data...\n");
    size_t fread_result = fread(&run_number,sizeof(run_number),1,old_grid_);
    run_number++;
    fread_result = fread(&t0,sizeof(t0),1,old_grid_);
    if(t0>=tf) // Check to make sure that the time in the old grid is earlier than the final time for this run
    {
      printf("A grid image file was found in this directory with values stored at t=%f. To continue that run set tf to a later time. To start a new run move or rename the file grid.img.\n",t0);
      exit(1);
    }
    fread_result = fread(&a,sizeof(a),1,old_grid_);
    fread_result = fread(&ad,sizeof(ad),1,old_grid_);
    fread_result = fread(f,sizeof(float),nflds*gridsize,old_grid_);
    fread_result = fread(fd,sizeof(float),nflds*gridsize,old_grid_);
    fclose(old_grid_);
    if(continue_run==1) // Option to append new data to old data files
      sprintf(mode_,"a+");
    else if(alt_extension[0]=='\0') // If no alternate extension was given set filename extension to indicate run number
      sprintf(ext_,"_%d.dat",run_number);
    printf("Data read. Resuming run at t=%f\n",t0);
    output_parameters(); // Save information about the model and parameters and start the clock for timing the run
    return;
  }

  // If the variable no_initialization is set to 1 by the model file then don't initialize the modes
  if(no_initialization==1)
  {
    save(0); // Save field values and derived quantities
    output_parameters(); // Save information about the model and parameters and start the clock for timing the run
    return;
  }

  // If no old grid image is found generate new initial conditions and set run_number=0
  printf("Generating initial conditions for new run at t=0\n");
  t0=0;
  run_number=0;

  for(fld=0;fld<nflds;fld++)
  {
    // Set initial field values
    if(fld<(int)(sizeof initfield/sizeof(float))) // Use preset values for initial field values if given
      initial_field_values[fld]=initfield[fld];
    else // Otherwise initialize field to zero
      initial_field_values[fld]=0.;

    // Set initial field derivative values
    if(fld<(int)(sizeof initderivs/sizeof(float))) // Use preset values for initial field derivative values if given
      initial_field_derivs[fld]=initderivs[fld];
    else // Otherwise initialize derivatives to zero
      initial_field_derivs[fld]=0.;
  }

  // Set initial values of effective mass.
  effective_mass(mass_sq,initial_field_values);

  // Set initial value of Hubble constant - See the documentation for an explanation of the formulas used
  if(expansion>0)
  {
    if(expansion==1)
      printf("The initial value of the fields is being used to determine the initial Hubble constant.\nFrom here on power law expansion will be used\n");
    for(fld=0;fld<nflds;fld++) // Find sum of squares of fields and derivatives
    {
      fsquared += pw2(initial_field_values[fld]);
      fdsquared += pw2(initial_field_derivs[fld]);
      ffd += initial_field_values[fld]*initial_field_derivs[fld];
    }
    for(i=0;i<num_potential_terms;i++) // Find potential energy
      pot_energy += potential_energy(i,initial_field_values);
    hubble_init = sqrt( 3.*pw2(rescale_A)/(4.*pi)*fdsquared + 2.*pot_energy*(3.*pw2(rescale_A)/(4.*pi)-pw2(rescale_r)*fsquared) );
    hubble_init -= rescale_r*ffd;
    hubble_init /= 3.*pw2(rescale_A)/(4.*pi)-pw2(rescale_r)*fsquared;
    if(!(hubble_init>=0.)) // Make sure Hubble isn't negative or undefined
    {
      printf("Error in calculating initial Hubble constant. Exiting.\n");
      exit(1);
    }
    ad=hubble_init;
  }

  for(fld=0;fld<nflds;fld++) // Set initial conditions for each field
  {
#if NDIMS==1
    // Set mode with k=N/2. This is a real number stored in the position f[fld][1].
    // This is done before the loop because this would overwrite the real part of k=1 otherwise.
    p2=dp2*pw2(N/2); // Total momentum squared for k=N/2
    set_mode(p2,mass_sq[fld],&f[fld][1],&fd[fld][1],1); // The last argument specifies that a real value should be set

    // Loop over gridpoints.
    for(k=1;k<N/2;k++)
    {
      pz=k; // z-component of momentum of modes at z=k
      p2=dp2*pw2(pz); // Total momentum squared
      set_mode(p2,mass_sq[fld],&f[fld][2*k],&fd[fld][2*k],0); // Set mode
    }

    f[fld][0]=0.; // Set zeromode of field and derivative to zero (it gets added in position space)
    fd[fld][0]=0.;

    // *DEBUG* The option to not use the FFTs is for debugging purposes only. For actual simulations seed should always be positive
    if(seed>=0)
    {
      fftr1((float *)f[fld],N,-1); // Inverse Fourier transform of field
      fftr1((float *)fd[fld],N,-1); // Inverse Fourier transform of field derivatives
    }
#elif NDIMS==2
    // Loop over gridpoints.
    // py goes over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
    for(j=0;j<N;j++)
    {
      py=(j<=N/2 ? j : j-N); // y-component of momentum of modes at y=j

      // Set all modes with 0<k<N/2. The complex conjugates of these modes do not need to be set.
      for(k=1;k<N/2;k++)
      {
        pz=k; // z-component of momentum of modes at z=k
        p2=dp2*(pw2(py)+pw2(pz)); // Total momentum squared
        set_mode(p2,mass_sq[fld],&f[fld][j][2*k],&fd[fld][j][2*k],0); // Set mode
      }

      // Set modes with k=0 or N/2.
      if(j>N/2) // The complex conjugates of these modes appear explicitly on the lattice and must be set to satisfy f(-p)=f(p)*
      {
        jconj=N-j; // Index where complex conjugates of modes at y=j are stored
        // k=0
        p2=dp2*pw2(py); // Total momentum squared
        set_mode(p2,mass_sq[fld],&f[fld][j][0],&fd[fld][j][0],0); // Set mode
        f[fld][jconj][0]=f[fld][j][0]; // Set complex conjugate mode
        f[fld][jconj][1]=-f[fld][j][1];
        fd[fld][jconj][0]=fd[fld][j][0];
        fd[fld][jconj][1]=-fd[fld][j][1];
        // k=N/2
        p2=dp2*(pw2(py)+pw2(N/2)); // Total momentum squared
        set_mode(p2,mass_sq[fld],&fnyquist[2*j],&fdnyquist[2*j],0); // Set mode
        fnyquist[2*jconj]=fnyquist[2*j]; // Set complex conjugate mode
        fnyquist[2*jconj+1]=-fnyquist[2*j+1];
        fdnyquist[2*jconj]=fdnyquist[2*j];
        fdnyquist[2*jconj+1]=-fdnyquist[2*j+1];
      }
      else if(j==0 || j==N/2) // The 4 "corners" of the lattice are set to real values
      {
        p2=dp2*pw2(py); // Total momentum squared for k=0
        if(p2>0.) // Don't set the zeromode here (see below)
          set_mode(p2,mass_sq[fld],&f[fld][j][0],&fd[fld][j][0],1); // The last argument specifies that a real value should be set
        p2=dp2*(pw2(py)+pw2(N/2)); // Total momentum squared for k=N/2
        set_mode(p2,mass_sq[fld],&fnyquist[2*j],&fdnyquist[2*j],1); // The last argument specifies that a real value should be set
      }
    } // End of loop over j (y-index on lattice)
    f[fld][0][0]=0.; // Set zeromode of field and derivative to zero (it gets added in position space)
    f[fld][0][1]=0.;
    fd[fld][0][0]=0.;
    fd[fld][0][1]=0.;

    // *DEBUG* The option to not use the FFTs is for debugging purposes only. For actual simulations seed should always be positive
    if(seed>=0)
    {
      fftrn((float *)f[fld],(float *)fnyquist,NDIMS,arraysize,-1); // Inverse Fourier transform of field
      fftrn((float *)fd[fld],(float *)fdnyquist,NDIMS,arraysize,-1); // Inverse Fourier transform of field derivatives
    }
#elif NDIMS==3
    // Loop over gridpoints.
    // px and py go over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
    for(i=0;i<N;i++)
    {
      px=(i<=N/2 ? i : i-N); // x-component of momentum of modes at x=i
      iconj=(i==0 ? 0 : N-i); // Index where complex conjugates of modes at x=i are stored (only used for k=0 or k=N/2)
      for(j=0;j<N;j++)
      {
        py=(j<=N/2 ? j : j-N); // y-component of momentum of modes at y=j

        // Set all modes with 0<k<N/2. The complex conjugates of these modes do not need to be set.
        for(k=1;k<N/2;k++)
        {
          pz=k; // z-component of momentum of modes at z=k
          p2=dp2*(pw2(px)+pw2(py)+pw2(pz)); // Total momentum squared
          set_mode(p2,mass_sq[fld],&f[fld][i][j][2*k],&fd[fld][i][j][2*k],0); // Set mode
        }

        // Set modes with k=0 or N/2.
        if(j>N/2 || (i>N/2 && (j==0 || j==N/2))) // The complex conjugates of these modes appear explicitly on the lattice and must be set to satisfy f(-p)=f(p)*
        {
          jconj=(j==0 ? 0 : N-j); // Index where complex conjugates of modes at y=j are stored
          // k=0
          p2=dp2*(pw2(px)+pw2(py)); // Total momentum squared
          set_mode(p2,mass_sq[fld],&f[fld][i][j][0],&fd[fld][i][j][0],0); // Set mode
          f[fld][iconj][jconj][0]=f[fld][i][j][0]; // Set complex conjugate mode
          f[fld][iconj][jconj][1]=-f[fld][i][j][1];
          fd[fld][iconj][jconj][0]=fd[fld][i][j][0];
          fd[fld][iconj][jconj][1]=-fd[fld][i][j][1];
          // k=N/2
          p2=dp2*(pw2(px)+pw2(py)+pw2(N/2)); // Total momentum squared
          set_mode(p2,mass_sq[fld],&fnyquist[i][2*j],&fdnyquist[i][2*j],0); // Set mode
          fnyquist[iconj][2*jconj]=fnyquist[i][2*j]; // Set complex conjugate mode
          fnyquist[iconj][2*jconj+1]=-fnyquist[i][2*j+1];
          fdnyquist[iconj][2*jconj]=fdnyquist[i][2*j];
          fdnyquist[iconj][2*jconj+1]=-fdnyquist[i][2*j+1];
        }
        else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
        {
          p2=dp2*(pw2(px)+pw2(py)); // Total momentum squared for k=0
          if(p2>0.) // Don't set the zeromode here (see below)
            set_mode(p2,mass_sq[fld],&f[fld][i][j][0],&fd[fld][i][j][0],1); // The last argument specifies that a real value should be set
          p2=dp2*(pw2(px)+pw2(py)+pw2(N/2)); // Total momentum squared for k=N/2
          set_mode(p2,mass_sq[fld],&fnyquist[i][2*j],&fdnyquist[i][2*j],1); // The last argument specifies that a real value should be set
        }
      } // End of loop over j (y-index on lattice)
    } // End of loop over i (x-index on lattice)
    f[fld][0][0][0]=0.; // Set zeromode of field and derivative to zero (it gets added in position space)
    f[fld][0][0][1]=0.;
    fd[fld][0][0][0]=0.;
    fd[fld][0][0][1]=0.;

    // *DEBUG* The option to not use the FFTs is for debugging purposes only. For actual simulations seed should always be positive
    if(seed>=0)
    {
      fftrn((float *)f[fld],(float *)fnyquist,NDIMS,arraysize,-1); // Inverse Fourier transform of field
      fftrn((float *)fd[fld],(float *)fdnyquist,NDIMS,arraysize,-1); // Inverse Fourier transform of field derivatives
    }
#endif
    // Add zeromode
    LOOP
    {
      FIELD(fld) += initial_field_values[fld];
      FIELDD(fld) += initial_field_derivs[fld];
    }
  } // End loop over fields

  modelinitialize(2); // This allows specific models to perform any needed initialization
  save(0); // Save field values and derived quantities
  output_parameters(); // Save information about the model and parameters and start the clock for timing the run

  printf("Finished initial conditions\n");
  
  return;
}







