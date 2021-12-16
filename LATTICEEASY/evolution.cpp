/*
This file contains the functions for taking time steps to solve differential equations on the lattice.

The method used is a staggered leapfrog on a spatial grid with fixed time step. In other words the field values on the grid are stored in an array f at some time t and the first derivatives are stored in fd at a time t +/- 1/2 dt.

The functions here which should be called externally are gradient_energy, evolve_scale, evolve_fields and evolve_derivs.
(Only evolve_fields() and evolve_derivs() are called externally by the main function, but certain output functions call gradient_energy() and evolve_scale().)
The function gradient_energy(int fld) calculates the average gradient energy density of the field fld on the lattice.
The function evolve_scale(float d) advances the value of the scale factor. The behavior of this function depends on the type of expansion being used. For no expansion (0) the function does nothing. For power-law expansion (1) the function calculates a and its derivatives at time t according to a preset formula. For consistent expansion (2) the function calculates the second derivative of a and uses it to advance adot by an interval d. In this case the scale factor itself is advanced in the function evolve_fields().
The function evolve_fields(float d) advances the field values (f) and the scale factor (a) by a time interval d assuming that the field derivatives (fd) and the scale factor derivative (ad) are known at a time in the middle of that interval.
The function evolve_derivs(float d) advances the field derivatives by a time interval d. It assumes f is known at the midpoint of that interval, uses it to calculate the second derivative, and uses that to advance the first derivatives. It also calls evolve_fields() to advance the scale factor.
*/

#include "latticeeasy.h"

// Increments a grid location accounting for periodic wrapping
inline int INCREMENT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}

// Decrements a grid location accounting for periodic wrapping
inline int DECREMENT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}

// Calculate the Laplacian of a field point in the bulk. (The result must be divided by dx^2 to give a Laplacian.)
inline float lapl(int fld, INDEXLIST)
{
#if NDIMS==1
  return (f[fld][i+1] + f[fld][i-1] - 2.*f[fld][i]);
#elif NDIMS==2
  return (f[fld][i][j+1] + f[fld][i][j-1]
         +f[fld][i+1][j] + f[fld][i-1][j]
         -4.*f[fld][i][j]);
#elif NDIMS==3
  return (f[fld][i][j][k+1] + f[fld][i][j][k-1]
         +f[fld][i][j+1][k] + f[fld][i][j-1][k]
         +f[fld][i+1][j][k] + f[fld][i-1][j][k]
         -6.*f[fld][i][j][k]);
#endif
}

// Calculate the Laplacian of a field point on the boundary. (The result must be divided by dx^2 to give a Laplacian.)
inline float laplb(int fld, INDEXLIST)
{
#if NDIMS==1
  return (f[fld][INCREMENT(i)] + f[fld][DECREMENT(i)] - 2.*f[fld][i]);
#elif NDIMS==2
  return (f[fld][i][INCREMENT(j)] + f[fld][i][DECREMENT(j)]
         +f[fld][INCREMENT(i)][j] + f[fld][DECREMENT(i)][j]
         -4.*f[fld][i][j]);
#elif NDIMS==3
  return (f[fld][i][j][INCREMENT(k)] + f[fld][i][j][DECREMENT(k)]
         +f[fld][i][INCREMENT(j)][k] + f[fld][i][DECREMENT(j)][k]
         +f[fld][INCREMENT(i)][j][k] + f[fld][DECREMENT(i)][j][k]
         -6.*f[fld][i][j][k]);
#endif
}

/////////////////////////////////////////////////////
// Externally called function(s)
/////////////////////////////////////////////////////

// Calculate the gradient energy, 1/2 <|Grad(f)|^2> = 1/2 <-f Lapl(f)>, of a field
float gradient_energy(int fld)
{
  int i=0,j=0,k=0; // The initial values are just put in to make the compiler happy
  double gradient=0.;
  float norm=1./pw2(dx); // Converts the output of lapl() to an actual Laplacian

  // Iterate over the bulk - Using separate loops for the bulk and boundary avoids unnecessary calls to INCREMENT and DECREMENT.
  for(i=1;i<N-1;i++)
#if NDIMS>1
    for(j=1;j<N-1;j++)
#endif
#if NDIMS>2
      for(k=1;k<N-1;k++)
#endif
        gradient -= FIELD(fld)*lapl(fld,i,j,k); // Extra indices are ignored
  
  // Iterate over the boundary
#if NDIMS==1
  gradient -= f[fld][0]*laplb(fld,0) + f[fld][N-1]*laplb(fld,N-1); // Right and left endpoints (x=0,N-1)
#elif NDIMS==2
  for(i=0;i<N;i++) // Index in x direction
  {
    gradient -= f[fld][i][0]*laplb(fld,i,0) + f[fld][i][N-1]*laplb(fld,i,N-1); // Front and back surfaces (y=0,N-1)
    if(i==0 || i==N-1) continue; // Reinterpret i as a y index. Don't double-count front and back points.
    gradient -= f[fld][0][i]*laplb(fld,0,i) + f[fld][N-1][i]*laplb(fld,N-1,i); // Right and left surfaces (x=0,N-1)
  }
#elif NDIMS==3
  for(i=0;i<N;i++) // Index in x direction
    for(j=0;j<N;j++) // Index in y direction
    {
      gradient -= f[fld][i][j][0]*laplb(fld,i,j,0) + f[fld][i][j][N-1]*laplb(fld,i,j,N-1); // Top and bottom surfaces (z=0,N-1)
      if(j==0 || j==N-1) continue; // Reinterpret j as a z index. Don't double-count top and bottom points.
      gradient -= f[fld][i][0][j]*laplb(fld,i,0,j) + f[fld][i][N-1][j]*laplb(fld,i,N-1,j); // Front and back surfaces (y=0,N-1)
      if(i==0 || i==N-1) continue; // Reinterpret i as a y index. Don't double-count front and back points.
      gradient -= f[fld][0][i][j]*laplb(fld,0,i,j) + f[fld][N-1][i][j]*laplb(fld,N-1,i,j); // Right and left surfaces (x=0,N-1)
    }
#endif

  // norm converts the results of lapl() to an actual laplacian and gridsize converts the sum over the lattice to an average
  return(.5*gradient*norm/(float)gridsize);
}

// Calculate the scale factor and its derivatives
// Use d=0 to indicate that all quantities are known at the same time. Otherwise it's assumed that they are known at staggered times with time step d.
void evolve_scale(float d)
{
  int term ,fld;
  float grad_energy=0.,pot_energy=0.; // Gradient and potential energies of fields
  float sfexponent,sfbase; // Model-dependent terms in power-law expansion
  float sfev1,sfev2,sfev3; // Model-dependent terms in self-consistent expansion

  if(expansion==0) // In case of no expansion do nothing
    return;
  else if(expansion==1) // In case of power-law expansion set all scale factor variables at time t using the model dependent parameter sfexponent.
  {
    sfexponent = expansion_power/(expansion_power*rescale_s+1.); // Exponent in power-law expansion expression for scale factor evolution
    sfbase = t*hubble_init/sfexponent + 1.; // Base of the exponent in power-law expansion expression for scale factor evolution
    a = pow(sfbase,sfexponent); // Scale factor
    ad = hubble_init/sfbase*a; // First derivative of scale factor
    ad2 = (sfexponent-1.)/sfexponent*pw2(hubble_init)/pw2(sfbase)*a; // Second derivative of scale factor
    aterm = rescale_r*(rescale_s-rescale_r+2.)*pw2(ad/a) + rescale_r*ad2/a; // Term used in evolving fields
  }
  else if(expansion==2) // In case of self-consistent expansion calculate sfdd
  {
    sfev1=rescale_s+2.; // See documentation for an explanation of these terms in the evolution equation for a
    sfev2=2.*(rescale_r+rescale_s+1.);
    sfev3=2.*(rescale_s+1.);
    for(fld=0;fld<nflds;fld++) // Sum gradient energy over all fields
      grad_energy += gradient_energy(fld);
    for(term=0;term<num_potential_terms;term++)
      pot_energy += potential_energy(term,NULL);
    if(d==0) // If all quantities are known at the same time calculate ad2 directly. (This option is called by the output routine while synchronizing values.)
      ad2 = -sfev1*pw2(ad)/a + 8.*pi/pw2(rescale_A)/pow(a,sfev2-1.)*(2.*grad_energy/3.+pow(a,sfev3)*pot_energy);
    else // Otherwise use the leapfrog correction, and use ad2 to calculate aterm and advance ad.
    {
      ad2 = (-2.*ad - 2.*a/d/sfev1*(1.-sqrt(1.+2.*d*sfev1*ad/a+8.*pi*pw2(d/rescale_A)*sfev1/pow(a,sfev2)*(2.*grad_energy/3.+pow(a,sfev3)*pot_energy))))/d;
      ad += .5*d*ad2; // Advance ad to time t for field evolution equations
      aterm = rescale_r*(rescale_s-rescale_r+2.)*pw2(ad/a) + rescale_r*ad2/a; // Term used in evolving fields
      ad += .5*d*ad2; // Advance ad to time t+dt/2 for evolving the scale factor
    }
  }
}

// Advance the field values and scale factor using the first derivatives
void evolve_fields(float d)
{
  DECLARE_INDICES
  int fld;

  // Advance time
  t += d;

  // Advance field values
  for(fld=0;fld<nflds;fld++)
    LOOP // Loop over all gridpoints
      FIELD(fld) += d*FIELDD(fld);

  // In case of self-consistent expansion advance scale factor
  if(expansion==2) 
    a += d*ad;
  
  return;
}


// Calculate second derivatives of fields and use them to advance first derivatives, and call evolve_scale to advance the scale factor
void evolve_derivs(float d)
{
  int i=0,j=0,k=0,fld; // The initial values are just put in to make the compiler happy
  float laplnorm = 1./pw2(dx)/pow(a,2.*rescale_s+2.); // Set coefficient for laplacian term in equations of motion. The dx^2 converts the output of lapl() to a laplacian and the scale factor term accounts for model dependent rescalings of the equations of motion.

  evolve_scale(d); // Calculate the scale factor and its derivatives

  for(fld=0;fld<nflds;fld++)
  {
    // Iterate over the bulk - Using separate loops for the bulk and boundary avoids unnecessary calls to INCREMENT and DECREMENT.
    for(i=1;i<N-1;i++)
#if NDIMS>1
      for(j=1;j<N-1;j++)
#endif
#if NDIMS>2
        for(k=1;k<N-1;k++)
#endif
          FIELDD(fld) += d*(laplnorm*lapl(fld,i,j,k)+aterm*FIELD(fld)-dvdf(fld,i,j,k)); // Extra indices are ignored

    // Iterate over the boundary
#if NDIMS==1
    fd[fld][0] += d*(laplnorm*laplb(fld,0)+aterm*f[fld][0]-dvdf(fld,0)); // Right and left surfaces (x=0,N-1)
    fd[fld][N-1] += d*(laplnorm*laplb(fld,N-1)+aterm*f[fld][N-1]-dvdf(fld,N-1));
#elif NDIMS==2
    for(i=0;i<N;i++) // Index in x direction
    {
      fd[fld][i][0] += d*(laplnorm*laplb(fld,i,0)+aterm*f[fld][i][0]-dvdf(fld,i,0)); // Front and back surfaces (y=0,N-1)
      fd[fld][i][N-1] += d*(laplnorm*laplb(fld,i,N-1)+aterm*f[fld][i][N-1]-dvdf(fld,i,N-1));
      if(i==0 || i==N-1) continue; // Reinterpret i as a y index. Don't double-count front and back points.
      fd[fld][0][i] += d*(laplnorm*laplb(fld,0,i)+aterm*f[fld][0][i]-dvdf(fld,0,i)); // Right and left surfaces (x=0,N-1)
      fd[fld][N-1][i] += d*(laplnorm*laplb(fld,N-1,i)+aterm*f[fld][N-1][i]-dvdf(fld,N-1,i));
    }
#elif NDIMS==3
    for(i=0;i<N;i++) // Index in x direction
      for(j=0;j<N;j++) // Index in y direction
      {
        fd[fld][i][j][0] += d*(laplnorm*laplb(fld,i,j,0)+aterm*f[fld][i][j][0]-dvdf(fld,i,j,0)); // Top and bottom surfaces (z=0,N-1)
        fd[fld][i][j][N-1] += d*(laplnorm*laplb(fld,i,j,N-1)+aterm*f[fld][i][j][N-1]-dvdf(fld,i,j,N-1));
        if(j==0 || j==N-1) continue; // Reinterpret j as a z index. Don't double-count top and bottom points.
        fd[fld][i][0][j] += d*(laplnorm*laplb(fld,i,0,j)+aterm*f[fld][i][0][j]-dvdf(fld,i,0,j)); // Front and back surfaces (y=0,N-1)
        fd[fld][i][N-1][j] += d*(laplnorm*laplb(fld,i,N-1,j)+aterm*f[fld][i][N-1][j]-dvdf(fld,i,N-1,j));
        if(i==0 || i==N-1) continue; // Reinterpret i as a y index. Don't double-count front and back points.
        fd[fld][0][i][j] += d*(laplnorm*laplb(fld,0,i,j)+aterm*f[fld][0][i][j]-dvdf(fld,0,i,j)); // Right and left surfaces (x=0,N-1)
        fd[fld][N-1][i][j] += d*(laplnorm*laplb(fld,N-1,i,j)+aterm*f[fld][N-1][i][j]-dvdf(fld,N-1,i,j));
      }
#endif
  } // End loop over fields
}
