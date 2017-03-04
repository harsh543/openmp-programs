# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
double f ( double x );
double cpu_time ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )
{
  double a = 0.0;
  double b = 10.0;
  double error;
  double exact = 0.49936338107645674464;
  int i;
  int n = 10000000;
  double total;
  double wtime;
  double x;

  timestamp ( );
  cout << "\n";
  cout << "QUAD_OPENMP:\n";
  cout << "  C++ version\n";
  cout << "  Use OpenMP for parallel execution.\n";
  cout << "  Estimate the integral of f(x) from A to B.\n";
  cout << "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n";
  cout << "\n";
  cout << "  A        = " << a << "\n";
  cout << "  B        = " << b << "\n";
  cout << "  N        = " << n << "\n";
  cout << "  Exact    = " << setprecision(16) << setw(24) << exact << "\n";

  wtime = omp_get_wtime ( );

  total = 0.0;

# pragma omp parallel shared ( a, b, n ) private ( i, x )

# pragma omp for reduction ( + : total )

  for ( i = 0; i < n; i++ )
  {
    x = ( ( double ) ( n - i - 1 ) * a 
        + ( double ) (     i     ) * b ) 
        / ( double ) ( n     - 1 );
    total = total + f ( x );
  }

  wtime = omp_get_wtime ( ) - wtime;

  total = ( b - a ) * total / ( double ) n;
  error = fabs ( total - exact );

  cout << "\n";
  cout << "  Estimate = " << setprecision(16) << setw(24) << total << "\n";
  cout << "  Error    = " << error << "\n";
  cout << "  W time   = " << wtime << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "QUAD_OPENMP:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double f ( double x )
{
  double pi = 3.141592653589793;
  double value;

  value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

  return value;
}
//****************************************************************************80

double cpu_time ( )
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

void timestamp ( )

{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
