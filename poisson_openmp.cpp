# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <omp.h>

using namespace std;

# define NX 161
# define NY 161

int main ( int argc, char *argv[] );
double r8mat_rms ( int m, int n, double a[NX][NY] );
void rhs ( int nx, int ny, double f[NX][NY] );
void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY], 
  int itold, int itnew, double u[NX][NY], double unew[NX][NY] );
void timestamp ( );
double u_exact ( double x, double y );
double uxxyy_exact ( double x, double y );

//****************************************************************************80

int main ( int argc, char *argv[] )

{
  bool converged;
  double diff;
  double dx;
  double dy;
  double error;
  double f[NX][NY];
  int i;
  int id;
  int itnew;
  int itold;
  int j;
  int jt;
  int jt_max = 20;
  int nx = NX;
  int ny = NY;
  double tolerance = 0.000001;
  double u[NX][NY];
  double u_norm;
  double udiff[NX][NY];
  double uexact[NX][NY];
  double unew[NX][NY];
  double unew_norm;
  double wtime;
  double x;
  double y;

  dx = 1.0 / ( double ) ( nx - 1 );
  dy = 1.0 / ( double ) ( ny - 1 );
//
//  Print a message.
//
  timestamp ( );
  cout << "\n";
  cout << "POISSON_OPENMP:\n";
  cout << "  C++ version\n";
  cout << "  A program for solving the Poisson equation.\n";
  cout << "\n";
  cout << "  Use OpenMP for parallel execution.\n";
  cout << "  The number of processors is " << omp_get_num_procs ( ) << "\n";
# pragma omp parallel
{
  id = omp_get_thread_num ( );
  if ( id == 0 )
  {
    cout << "  The maximum number of threads is " << omp_get_num_threads ( ) << "\n"; 
  }
}
  cout << "\n";
  cout << "  -DEL^2 U = F(X,Y)\n";
  cout << "\n";
  cout << "  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.\n";
  cout << "\n";
  cout << "  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )\n";
  cout << "\n";
  cout << "  The number of interior X grid points is " << nx << "\n";
  cout << "  The number of interior Y grid points is " << ny << "\n";
  cout << "  The X grid spacing is " << dx << "\n";
  cout << "  The Y grid spacing is " << dy << "\n";
//
//  Set the right hand side array F.
//
  rhs ( nx, ny, f );
//
//  Set the initial solution estimate UNEW.
//  We are "allowed" to pick up the boundary conditions exactly.
//
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        unew[i][j] = f[i][j];
      }
      else
      {
        unew[i][j] = 0.0;
      }
    }
  }
  unew_norm = r8mat_rms ( nx, ny, unew );
//
//  Set up the exact solution UEXACT.
//
  for ( j = 0; j < ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
      uexact[i][j] = u_exact ( x, y );
    }
  }
  u_norm = r8mat_rms ( nx, ny, uexact );
  cout << "  RMS of exact solution = " << u_norm << "\n";
//
//  Do the iteration.
//
  converged = false;

  cout << "\n";
  cout << "  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||\n";
  cout << "\n";

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      udiff[i][j] = unew[i][j] - uexact[i][j];
    }
  }
  error = r8mat_rms ( nx, ny, udiff );
  cout << "  " << setw(4) << 0
       << "  " << setw(14) << unew_norm
       << "  " << "              "
       << "  " << setw(14) << error << "\n";

  wtime = omp_get_wtime ( );

  itnew = 0;

  for ( ; ; )
  {
    itold = itnew;
    itnew = itold + 500;
//
//  SWEEP carries out 500 Jacobi steps in parallel before we come
//  back to check for convergence.
//
    sweep ( nx, ny, dx, dy, f, itold, itnew, u, unew );
//
//  Check for convergence.
//
    u_norm = unew_norm;
    unew_norm = r8mat_rms ( nx, ny, unew );

    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        udiff[i][j] = unew[i][j] - u[i][j];
      }
    }
    diff = r8mat_rms ( nx, ny, udiff );

    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        udiff[i][j] = unew[i][j] - uexact[i][j];
      }
    }
    error = r8mat_rms ( nx, ny, udiff );

    cout << "  " << setw(4)  << itnew
         << "  " << setw(14) << unew_norm
         << "  " << setw(14) << diff
         << "  " << setw(14) << error << "\n";

    if ( diff <= tolerance )
    {
      converged = true;
      break;
    }

  }

  if ( converged )
  {
    cout << "  The iteration has converged.\n";
  }
  else
  {
    cout << "  The iteration has NOT converged.\n";
  }

  wtime = omp_get_wtime ( ) - wtime;
  cout << "\n";
  cout << "  Elapsed seconds = " << wtime << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "POISSON_OPENMP:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double r8mat_rms ( int nx, int ny, double a[NX][NY] )
{
  int i;
  int j;
  double v;

  v = 0.0;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      v = v + a[i][j] * a[i][j];
    }
  }
  v = sqrt ( v / ( double ) ( nx * ny )  );

  return v;
}
//****************************************************************************80

void rhs ( int nx, int ny, double f[NX][NY] )
{
  double fnorm;
  int i;
  int j;
  double x;
  double y;
//
//  The "boundary" entries of F store the boundary values of the solution.
//  The "interior" entries of F store the right hand sides of the Poisson equation.
//
  for ( j = 0; j < ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        f[i][j] = u_exact ( x, y );
      }
      else
      {
        f[i][j] = - uxxyy_exact ( x, y );
      }
    }
  }

  fnorm = r8mat_rms ( nx, ny, f );

  cout << "  RMS of F = " << fnorm << "\n";

  return;
}
//****************************************************************************80

void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY], 
  int itold, int itnew, double u[NX][NY], double unew[NX][NY] )
{
  int i;
  int it;
  int j;

# pragma omp parallel shared ( dx, dy, f, itnew, itold, nx, ny, u, unew ) private ( i, it, j )

  for ( it = itold + 1; it <= itnew; it++ )
  {
//
//  Save the current estimate.
//
# pragma omp for
    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        u[i][j] = unew[i][j];
      }
    }
//
//  Compute a new estimate.
//
# pragma omp for
    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        if ( i == 0 || j == 0 || i == nx - 1 || j == ny - 1 )
        {
          unew[i][j] = f[i][j];
        }
        else
        { 
          unew[i][j] = 0.25 * ( 
            u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] + f[i][j] * dx * dy );
        }
      }
    }

  }
  return;
}
//****************************************************************************80

void timestamp ( )
{
  double pi = 3.141592653589793;
  double value;

  value = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y );

  return value;
}
# undef NX
# undef NY
