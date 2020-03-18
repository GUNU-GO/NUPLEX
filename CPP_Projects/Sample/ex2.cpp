// 1-d PIC code to solve plasma two-stream instability problem.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <blitz/array.h>
#include <fftw.h>

using namespace blitz;

void Output (char* fn1, char* fn2, double t, 
	     Array<double,1> r, Array<double,1> v);
void Density (Array<double,1> r, Array<double,1>& n);
void Electric (Array<double,1> phi, Array<double,1>& E);
void Poisson1D (Array<double,1>& u, Array<double,1> v, double kappa);
void rk4_fixed (double& x, Array<double,1>& y, 
              void (*rhs_eval)(double, Array<double,1>, Array<double,1>&), 
              double h);
void rhs_eval (double t, Array<double,1> y, Array<double,1>& dydt);
void Load (Array<double,1> r, Array<double,1> v, Array<double,1>& y);
void UnLoad (Array<double,1> y, Array<double,1>& r, Array<double,1>& v);
double distribution (double vb);

double L; int N, J;

int main()
{
  // Parameters
  L;            // Domain of solution 0 <= x <= L (in Debye lengths)
  N;            // Number of electrons
  J;            // Number of grid points
  double vb;    // Beam velocity
  double dt;    // Time-step (in inverse plasma frequencies)
  double tmax;  // Simulation run from t = 0. to t = tmax

  // Get parameters
  printf ("Please input N:  "); scanf ("%d", &N);
  printf ("Please input vb:  "); scanf ("%lf", &vb); 
  printf ("Please input L:  "); scanf ("%lf", &L);
  printf ("Please input J:  "); scanf ("%d", &J);
  printf ("Please input dt:  "); scanf ("%lf", &dt); 
  printf ("Please input tmax:  "); scanf ("%lf", &tmax);
  int skip = int (tmax / dt) / 10;
  if ((N < 1) || (J < 2) || (L <= 0.) || (vb <= 0.) 
      || (dt <= 0.) || (tmax <= 0.) || (skip < 1))
    {
      printf ("Error - invalid input parameters\n");
      exit (1);
    }

  // Set names of output files
  char* phase[11]; char* data[11];
  phase[0] = "phase0.out";phase[1] = "phase1.out";phase[2] = "phase2.out"; 
  phase[3] = "phase3.out";phase[4] = "phase4.out";phase[5] = "phase5.out"; 
  phase[6] = "phase6.out";phase[7] = "phase7.out";phase[8] = "phase8.out";
  phase[9] = "phase9.out";phase[10] = "phase10.out";data[0] = "data0.out";  
  data[1] = "data1.out"; data[2] = "data2.out"; data[3] = "data3.out"; 
  data[4] = "data4.out"; data[5] = "data5.out"; data[6] = "data6.out"; 
  data[7] = "data7.out"; data[8] = "data8.out"; data[9] = "data9.out"; 
  data[10] = "data10.out";

  // Initialize solution
  double t = 0.;
  int seed = time (NULL); srand (seed);
  Array<double,1> r(N), v(N);
  for (int i = 0; i < N; i++)
    {
      r(i) = L * double (rand ()) / double (RAND_MAX);
      v(i) = distribution (vb);
    }
  Output (phase[0], data[0], t, r, v);

  // Evolve solution
  Array<double,1> y(2*N);
  Load (r, v, y);
  for (int k = 1; k <= 10; k++)
    {
      for (int kk = 0; kk < skip; kk++)
        {
           // Take time-step
           rk4_fixed (t, y, rhs_eval, dt);

           // Make sure all coordinates in range 0 to L.
           for (int i = 0; i < N; i++)
             {
               if (y(i) < 0.) y(i) += L;
               if (y(i) > L) y(i) -= L;
             }
	  
           printf ("t = %11.4e\n", t);
        }
      printf ("Plot %3d\n", k);

      // Output data
      UnLoad (y, r, v);
      Output(phase[k], data[k], t, r, v);
    }

  return 0;
}
The following routine outputs the simulation data to various data-files.

// Write data to output files

void Output (char* fn1, char* fn2, double t,
	     Array<double,1> r, Array<double,1> v)
{  
  // Write phase-space data
  FILE* file = fopen (fn1, "w");
  for (int i = 0; i < N; i++)
    fprintf (file, "%e %e\n", r(i), v(i));
  fclose (file);

  // Write electric field data
  Array<double,1> ne(J), n(J), phi(J), E(J);  
  Density (r, ne);
  for (int j = 0; j < J; j++)
    n(j) = double (J) * ne(j) / double (N) - 1.;
  double kappa = 2. * M_PI / L; 
  Poisson1D (phi, n, kappa);
  Electric (phi, E);

  file = fopen (fn2, "w");
  for (int j = 0; j < J; j++)
    {
      double x = double (j) * L / double (J);
      fprintf (file, "%e %e %e %e\n", x, ne(j), n(j), E(j));
    }
  double x = L;
  fprintf (file, "%e %e %e %e\n", x, ne(0), n(0), E(0));
  fclose (file);
}
The following routine returns a random velocity distributed on a double Maxwellian distribution function corresponding to two counter-streaming beams. The algorithm used to achieve this is called the rejection method, and will be discussed later in this course.

// Function to distribute electron velocities randomly so as 
// to generate two counter propagating warm beams of thermal
// velocities unity and mean velocities +/- vb.
// Uses rejection method.

double distribution (double vb)
{ 
  // Initialize random number generator
  static int flag = 0;
  if (flag == 0)
    {
      int seed = time (NULL);
      srand (seed);
      flag = 1;
    }

  // Generate random v value
  double fmax = 0.5 * (1. + exp (-2. * vb * vb));
  double vmin = - 5. * vb;
  double vmax = + 5. * vb;
  double v = vmin + (vmax - vmin) * double (rand ()) / double (RAND_MAX);

  // Accept/reject value
  double f = 0.5 * (exp (-(v - vb) * (v - vb) / 2.) +
		    exp (-(v + vb) * (v + vb) / 2.));
  double x = fmax * double (rand ()) / double (RAND_MAX);
  if (x > f) return distribution (vb);
  else return v;
}
The routine below evaluates the electron number density on an evenly spaced mesh given the instantaneous electron coordinates.

// Evaluates electron number density n(0:J-1) from 
// array r(0:N-1) of electron coordinates.

void Density (Array<double,1> r, Array<double,1>& n)
{
  // Initialize 
  double dx = L / double (J);
  n = 0.;

  // Evaluate number density.
  for (int i = 0; i < N; i++)
    {
      int j = int (r(i) / dx);
      double y = r(i) / dx - double (j);
      n(j) += (1. - y) / dx;
      if (j+1 == J) n(0) += y / dx;
      else n(j+1) += y / dx;
    }
}
The following functions are wrapper routines for using the fftw library with periodic functions.

// Functions to calculate Fourier transforms of real data 
// using fftw Fast-Fourier-Transform routine.
// Input/ouput arrays are assumed to be of extent J.

// Calculates Fourier transform of array f in arrays Fr and Fi
void fft_forward (Array<double,1>f, Array<double,1>&Fr, 
      Array<double,1>& Fi)
{
  fftw_complex ff[J], FF[J];

  // Load data
  for (int j = 0; j < J; j++)
    {
      c_re (ff[j]) = f(j); c_im (ff[j]) = 0.;
    }

  // Call fftw routine
  fftw_plan p = fftw_create_plan (J, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_one (p, ff, FF);
  fftw_destroy_plan (p); 

  // Unload data
  for (int j = 0; j < J; j++)
    {
      Fr(j) = c_re (FF[j]); Fi(j) = c_im (FF[j]);
    }

  // Normalize data
  Fr /= double (J);
  Fi /= double (J);
}

// Calculates inverse Fourier transform of arrays Fr and Fi in array f
void fft_backward (Array<double,1> Fr, Array<double,1> Fi, 
      Array<double,1>& f)
{
  fftw_complex ff[J], FF[J];

  // Load data
  for (int j = 0; j < J; j++)
    {
      c_re (FF[j]) = Fr(j); c_im (FF[j]) = Fi(j);
    }

  // Call fftw routine
  fftw_plan p = fftw_create_plan (J, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_one (p, FF, ff);
  fftw_destroy_plan (p); 

  // Unload data
  for (int j = 0; j < J; j++)
      f(j) = c_re (ff[j]); 
}
The following routine solves Poisson's equation in 1-D to find the instantaneous electric potential on a uniform grid.

// Solves 1-d Poisson equation:
//    d^u / dx^2 = v   for  0 <= x <= L
// Periodic boundary conditions:
//    u(x + L) = u(x),  v(x + L) = v(x)
// Arrays u and v assumed to be of length J.
// Now, jth grid point corresponds to
//    x_j = j dx  for j = 0,J-1
// where dx = L / J.
// Also,
//    kappa = 2 pi / L

void Poisson1D (Array<double,1>& u, Array<double,1> v, double kappa)
{
  // Declare local arrays.
  Array<double,1> Vr(J), Vi(J), Ur(J), Ui(J);

  // Fourier transform source term
  fft_forward (v, Vr, Vi);

  // Calculate Fourier transform of u
  Ur(0) = Ui(0) = 0.;
  for (int j = 1; j <= J/2; j++)
    {
      Ur(j) = - Vr(j) / double (j * j) / kappa / kappa;
      Ui(j) = - Vi(j) / double (j * j) / kappa / kappa;
    } 
  for (int j = J/2; j < J; j++)
    {
      Ur(j) = Ur(J-j);
      Ui(j) = - Ui(J-j);
    }

  // Inverse Fourier transform to obtain u
  fft_backward (Ur, Ui, u);
}
The following function evaluates the electric field on a uniform grid from the electric potential.

// Calculate electric field from potential

void Electric (Array<double,1> phi, Array<double,1>& E)
{
  double dx = L / double (J);

  for (int j = 1; j < J-1; j++)
    E(j) = (phi(j-1) - phi(j+1)) / 2. / dx;
  E(0) = (phi(J-1) - phi(1)) / 2. / dx;
  E(J-1) = (phi(J-2) - phi(0)) / 2. / dx;
}
The following routine is the right-hand side routine for the electron equations of motion. Is is designed to be used with the fixed-step RK4 solver described earlier in this course.

// Electron equations of motion:
//    y(0:N-1)  = r_i
//    y(N:2N-1) = dr_i/dt

void rhs_eval (double t, Array<double,1> y, Array<double,1>& dydt)
{
  // Declare local arrays
  Array<double,1> r(N), v(N), rdot(N), vdot(N), r0(N);
  Array<double,1> ne(J), rho(J), phi(J), E(J);

  // Unload data from y
  UnLoad (y, r, v);

  // Make sure all coordinates in range 0 to L
  r0 = r;
  for (int i = 0; i < N; i++)
    {
      if (r0(i) < 0.) r0(i) += L;
      if (r0(i) > L) r0(i) -= L;
    }

  // Calculate electron number density
  Density (r0, ne);

  // Solve Poisson's equation
  double n0 = double (N) / L;
  for (int j = 0; j < J; j++)
    rho(j) = ne(j) / n0 - 1.;
  double kappa = 2. * M_PI / L; 
  Poisson1D (phi, rho, kappa);

  // Calculate electric field
  Electric (phi, E);

  // Equations of motion
  for (int i = 0; i < N; i++)
    {
      double dx = L / double (J);
      int j = int (r0(i) / dx);
      double y = r0(i) / dx - double (j);
      
      double Efield;
      if (j+1 == J)
         Efield = E(j) * (1. - y) + E(0) * y;
      else
         Efield = E(j) * (1. - y) + E(j+1) * y;

      rdot(i) = v(i);
      vdot(i) = - Efield;
    }

  // Load data into dydt
  Load (rdot, vdot, dydt);
}
The following functions load and unload the electron phase-space coordinates into the solution vector y used by the RK4 routine.

// Load particle coordinates into solution vector

void Load (Array<double,1> r, Array<double,1> v, Array<double,1>& y)
{
  for (int i = 0; i < N; i++)
    {
      y(i) = r(i);
      y(N+i) = v(i);
    }
}

// Unload particle coordinates from solution vector

void UnLoad (Array<double,1> y, Array<double,1>& r, Array<double,1>& v)
{
  for (int i = 0; i < N; i++)
    {
      r(i) = y(i);
      v(i) = y(N+i);
    }
}