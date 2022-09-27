#include <R.h>
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <R_ext/Applic.h>  // ST addition
#include <R_ext/Rdynload.h>  // ST addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define SWAP(a,b)   { int t; t=a; a=b; b=t; }  // Macro for swapping

/*
 * Brent Optimization algorithms
 * taken from https://people.sc.fsu.edu/~jburkardt/c_src/brent/brent.html
 * Authors: Richard Brent and John Burkardt (C code)
 */

double zerofinder ( double a, double b, double machep, double t, double f ( double x ) )

/******************************************************************************/
/*
  Purpose:

    ZERO seeks the root of a function F(X) in an interval [A,B].

  Discussion:

    The interval [A,B] must be a change of sign interval for F.
    That is, F(A) and F(B) must be of opposite signs.  Then
    assuming that F is continuous implies the existence of at least
    one value C between A and B for which F(C) = 0.

    The location of the zero is determined to within an accuracy
    of 6 * MACHEPS * abs ( C ) + 2 * T.

    Thanks to Thomas Secretin for pointing out a transcription error in the
    setting of the value of P, 11 February 2013.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2013

  Author:

    Original FORTRAN77 version by Richard Brent.
    C version by John Burkardt.

  Reference:

    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.

  Parameters:

    Input, double A, B, the endpoints of the change of sign interval.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double T, a positive error tolerance.

    Input, double F ( double x ), a user-supplied function whose zero
    is being sought.

    Output, double ZERO, the estimated value of a zero of
    the function F.
*/
{
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double m;
  double p;
  double q;
  double r;
  double s;
  double sa;
  double sb;
  double tol;
/*
  Make local copies of A and B.
*/
  sa = a;
  sb = b;
  fa = f ( sa );
  fb = f ( sb );

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;

  for ( ; ; )
  {
    if ( fabs ( fc ) < fabs ( fb ) )
    {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * machep * fabs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( fabs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if ( sa == c )
      {
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      }

      if ( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
        p < fabs ( 0.5 * s * q ) )
      {
        d = p / q;
      }
      else
      {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if ( tol < fabs ( d ) )
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    fb = f ( sb );

    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}

/******************************************************************************/

/*
 * END Brent Optimization algorithms
 */

/*
 * PELT code, modified and reused from
 * Rebecca Killick<r.killick@lancs.ac.uk
 * https://github.com/rkillick/changepoi
 */

// Cost functions

void meanvar_norm(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){
	double l = *end - *start;
	double x = SS[ *end ]  - SS[ *start ];
	double x2 = SS[ *n + *end ] - SS[ *n + *start ];
	double sigsq=(x2-((x*x)/l))/l;
	if(sigsq<=0){
	  sigsq=0.00000000001;
	}
	if(*MBIC == 0){
	  *cost = l*(log(2*M_PI)+log(sigsq)+1);
	}else{
	  *cost = l*(log(2*M_PI)+log(sigsq)+1)+log(l);
	}
}

void meanvar_poisson(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){
	double l = *end - *start;
	double x = SS[ *end ]  - SS[ *start ];

	if(x==0){
		*cost = 0;
	}else{
	  if(*MBIC == 0){
		*cost = 2*x*(log(l)-log(x));
	  }else{
		*cost = 2*x*(log(l)-log(x))+log(l);
	  }
	}
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 * Extensions to PELT code
 * Author: Michael P Schneider
 * michael.schneider@cruk.cam.ac.uk
 */

void meanvar_negbin_mm(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc) {

	/* number of samples / n */
	double l = *end - *start;
	/* sum_{i=1}^n y_i */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* sum_{i=1}^n y_i^2 */
	double x2 = SS[ *n + *end ]  - SS[ *n + *start ];
	/* estimator of mean */
	double m = x/l;
	/* offset index for raw values */
	int index = *n + *n;

	/* Find y_max */
	double y_max = 0;
	for (int i = (*start+1); i <= *end; i++) {
		if(y_max < SS[index + i]){
			y_max = SS[index + i];
		}
	}

	/* Computation of variance */
    /* note this is not necessarily numerically stable, but reasonably effective */
	double ssq = (x2 / l) - (m*m);
	/* see Anraku & Yanagimoto -> sample variance */
	ssq = ssq * (l/(l-1));
	/* if above line is commented out, we use population variance */
	if(ssq <= 0){
		ssq = 1e-10;
	}

	/* Method of Moments estimate of overdispersion */
	double c = (ssq - m) / (m * m);
	if(c <= -1./y_max){
		c = -1./y_max + 1e-8;
	}

	if(ssq <= m){
		c = -1./y_max + 1e-8;
	}

	/* Cost is minus the maximum log-likelihood */
	/* Compute log likelihood */
	double ll = 0;

	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(c*m));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(c*j);
		}
		ll -= lgamma(SS[index + i] + 1.);
	}
	/* end of likelihood computation */

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
      /* printf("COST-MM2: %f [%f] (%f - %f)/ %d - %f)\n", *cost, c, ssq,m, *start, l); */
}

void meanvar_negbin_bcmm(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* number of samples / n */
	double l = *end - *start;
	/* sum_{i=1}^n y_i */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* sum_{i=1}^n y_i^2 */
	double x2 = SS[ *n + *end ]  - SS[ *n + *start ];
	/* estimator of m */
	double m = x/l;

	/* offset index for raw values */
    int index = *n + *n;

  /* Find y_max */
  double y_max = 0;
  for (int i = (*start+1); i <= *end; i++) {
    if(y_max < SS[index + i]){
      y_max = SS[index + i];
    }
  }

  /* Computation of variance */
  /* note this is not necessarily numerically stable, but reasonably effective */
  double ssq = (x2 / l) - (m*m);
  /* see MM method for details of below line */
  ssq = ssq * (l/(l-1));
  if(ssq <= 0){
      ssq = 1e-10;
  }

	/* Bias-corrected Method of Moments estimate of dispersion parameter */
	double c = (ssq - m) / ((m*m) - (ssq/l));
	if(c <= -1./y_max){
		c = -1./y_max + 1e-8;
	}

	if(ssq <= m){
		c = -1./y_max + 1e-8;
	}
	/* Cost is minus the maximum log-likelihood */

	/* Compute log likelihood */
	double ll = 0;

	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(fmax(-1.+1e-6, c*m)));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(fmax(-1.+1e-6, c*j));
		}
		ll -= lgamma(SS[index + i] + 1.);
	}
	/* end of likelihood computation */

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
	/* printf("COST-BCMM: %f [%f] (%f - %f)/ %d - %f)\n", *cost, c, ssq,m, *start, l); */
}

void meanvar_negbin_ml(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* Maximum likelihood estimator of Negative Binomial */

	/* number of samples / n */
	double l = *end - *start;
	/* \sum_{i=1}^n y_i	 */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* estimator of m */
	double m = x/l;
	/* offset index for raw values */
	int index = *n + *n;

	/* Find y_max */
	double y_max = 0;
	for (int i = (*start+1); i <= *end; i++) {
		if(y_max < SS[index + i]){
			y_max = SS[index + i];
		}
	}

	/* Numerically identify c that maximises log likelihood */
	/* \delta logL / \delta c */
	double F_dll(double c, double m, double l, int index, double *SS, int *start, int *end){
		double dll = 0;
		dll += (l*(1./(c*c)) * log1p(c*m)) - (l * m / c);

		for (int i = (*start+1); i <= *end; i++) {
			for (int j = 0; j < SS[index + i]; j++) {
				dll += j / (1. + (c * j));
			}
		}
		return(dll);
	}

	/* evaluate F with the infered parameter (m) and unknown x (corresponding to value of c) */
	double F_dll_instance(double x){
		return F_dll(x, m, l, index, SS, start, end);
	}

	double c = zerofinder((-1./y_max)+1e-6, 10, 1e-10, 1e-6, F_dll_instance);

	/* Cost is minus the maximum log-likelihood */
	/* Compute log likelihood */
	double ll = 0;

	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(fmax(-1.+1e-6, c*m)));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(fmax(-1.+1e-6, c*j));
		}
		ll -= lgamma(SS[index + i] + 1.);
	}

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
	/* printf("COST-ML: %f [%f] (NA - %f)/ %d - %f)\n", *cost, c, m, *start, l); */
}

void meanvar_negbin_cml(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* Cox and Reid’s (1987) approximate conditional inference */
	/* Cox, David Roxbee, and Nancy Reid. "Parameter orthogonality and approximate conditional inference." Journal of the Royal Statistical Society: Series B (Methodological) 49.1 (1987): 1-18. */

	/* more modern reference
	 * Robinson, Mark D., and Gordon K. Smyth. "Small-sample estimation of negative binomial dispersion, with applications to SAGE data." Biostatistics 9.2 (2008): 321-332.
	 */

	/* number of samples / n */
	double l = *end - *start;
	/* \sum_{i=1}^n y_i	 */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* estimator of m */
	double m = x/l;
	/* offset index for raw values */
	int index = *n + *n;

	/* Find y_max */
	double y_max = 0;
	for (int i = (*start+1); i <= *end; i++) {
		if(y_max < SS[index + i]){
			y_max = SS[index + i];
		}
	}

	double F_dcml(double phi, int index, double n, double lambda, double *SS, int *start, int *end){

		double dll = 0;
		dll += (n * (1./(phi*phi)) * log1p(phi*m)) - (n * m / phi);

		for (int i = (*start+1); i <= *end; i++) {
			for (int j = 0; j < SS[index + i]; j++) {
				dll += j / (1. + (phi * j));
			}
		}

		/* derivative of -0.5 * log * j_lambda,lambda term */
		double dll_cml = dll + (0.5 * (lambda / (1. + (phi * lambda))));
		return(dll_cml);
	}

	double F_dcml_instance(double x){
		return F_dcml(x, index, l, m, SS, start, end);
	}

	/* c_ml > -1./y_max */
	/* we find the maximum by looking at the root of the first derivative - we skip checking for the second derivative */
	double c = zerofinder((-1./y_max)+1e-6, 10, 1e-10, 1e-6, F_dcml_instance);

	/* Cost is minus the maximum log-likelihood */

	/* Compute log likelihood */
	double ll = 0;

	/* int index = *n + *n; */
	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(c*m));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(c*j);
		}
		ll -= lgamma(SS[index + i] + 1.);
	}
	/* subtract Fisher term */
	ll = ll - 0.5 * log(l / (m * (1. + c * m)));

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
	/* printf("COST-CML: %f [%f] (NA - %f)/ %d - %f)\n", *cost, c, m, *start, l); */
}

void meanvar_negbin_alpha(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* Maximum likelihood estimator of Negative Binomial */

	/* number of samples / n */
	double l = *end - *start;
	/* \sum_{i=1}^n y_i	 */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* estimator of m */
	double m = x/l;
	/* offset index for raw values */
	int index = *n + *n;

	double c = *alpha;

	/* Cost is minus the maximum log-likelihood */
	/* Compute log likelihood */
	double ll = 0;

	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(fmax(-1.+1e-6, c*m)));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(fmax(-1.+1e-6, c*j));
		}
		ll -= lgamma(SS[index + i] + 1.);
	}

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
    /* printf("COST-NB-alpha: %f [%f] (NA - %f)/ %d - %f)\n", *cost, c, m, *start, l); */
}

void meanvar_negbin(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* Maximum likelihood estimator of Negative Binomial */

	/* number of samples / n */
	double l = *end - *start;
	/* \sum_{i=1}^n y_i	 */
	double x = SS[ *end ]  - SS[ *start ];
	if(x==0){
		*cost = 0;
		return;
	}
	/* estimator of m */
	double m = x/l;
    double m_old = x/l;

    /* find m_gap closest to copy number state, i.e. min_{m_gap} abs(m - m_gap) */
    double m_gap = 0.5;
    double delta = fabs(m_gap - m);
    for(int i = 1; i <= 50; i++) {
      if(fabs(i * *rpc - m) < delta){
        delta = fabs(i * *rpc - m);
        m_gap = i * *rpc;
      }
    }
    m = m_gap;

	/* offset index for raw values */
	int index = *n + *n;

	double c = *alpha;

	/* Cost is minus the maximum log-likelihood */
	/* Compute log likelihood */
	double ll = 0;

	for (int i = (*start+1); i <= *end; i++) {
		ll += (SS[index + i] * log(m)) - ((SS[index + i] + (1./c)) * log1p(fmax(-1.+1e-6, c*m)));
		for (int j = 0; j < SS[index + i]; j++) {
			ll += log1p(fmax(-1.+1e-6, c*j));
		}
		ll -= lgamma(SS[index + i] + 1.);
	}

    if(*MBIC == 0){
      *cost = -1. * ll;
    }else{
      *cost = -1. * ll + log(l);
    }
    /* printf("COST-NB-full: %f [%f] (%f - %f)/ %d - %f)\n", *cost, c, m_gap, m_old, *start, l); */
}

void state_space(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* State space segmentation */
	double l = *end - *start;

	/* p is the state space truncation level or the maximum state space level */
	double max_ll = -1.0/0.0;
	double ll = 0;

	int index = 0;
	for (int i = 0; i <= *p; i++){
		index = i * *n;
		ll = SS[ *end + index ]  - SS[ *start + index ];
		if(ll > max_ll) {
			max_ll = ll;
		}
	}

	/* Cost is minus the maximum log-likelihood */
    if(*MBIC == 0){
      *cost = -1. * max_ll;
    }else{
      *cost = -1. * max_ll + log(l);
    }
	/* printf("STATESPACE: %f / %d - %f)\n", *cost, *start, l); */
}

void state_space_penalized(double *SS, int *size, int *n, int *p, int *minorder, int *optimalorder, int *maxorder, int *start, int *end, double *cost, double *tol, int *error, double *shape, int *MBIC, double *alpha, double *rpc){

	/* State space segmentation */
	double l = *end - *start;

	/* p is the state space truncation level or the maximum state space level */
	double max_ll = -1.0/0.0;
	double ll = 0;

	int index = 0;
	for (int i = 0; i <= *p; i++){
		index = i * *n;
		ll = SS[ *end + index ]  - SS[ *start + index ];
		if(ll > max_ll) {
			max_ll = ll;
		}
	}

	/* Cost is minus the maximum log-likelihood */
    double gap_penalty = log(1. + 0.00);
    printf("STATESPACE: %f (%f | %f) / %d - %f)\n", *cost, gap_penalty, max_ll, *start, l);
    *cost = -1. * max_ll + gap_penalty;
}

/* Helper functions from changepoint package */

//Find the maximum case
void max_which(double *array, int n, double *maxval, int *maxid){
  //array - values for which to find the maximum
  //n - number of items to search
  //maxval - maximum value
  //maxid - index of the maximum value in the array
  int i;
  *maxid = 0;
  *maxval = array[*maxid];
  for(i = 1; i < n; i++){
    if(array[i] > *maxval){
      *maxid = i;
      *maxval = array[i];
    }
  }
  return;
}

//Find the minimum case
void min_which(double *array, int *n, double *minval, int *minid){
  //array - values for which to find the minimum
  //n - number of items to search
  //minval - minimum value
  //minid - index of the minimum value in the array
  int i;
  *minid = 0;
  *minval = array[*minid];
  for(i = 1; i < *n; i++){
    if(array[i] < *minval){
      *minid = i;
      *minval = array[i];
    }
  }
  return;
}

void order_vec( int a[], int n ){
  int i, j;
  for(i = 0; i < n; i++){  // Make a pass through the array for each element
    for(j = 1; j < (n-i); j++){  		// Go through the array beginning to end
      if(a[j-1] > a[j]){       // If the the first number is greater, swap it
        SWAP(a[j-1],a[j]);
      }
    }
  }
}
