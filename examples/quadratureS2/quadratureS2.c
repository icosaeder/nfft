/* $Id$
 *
 * fastsumS2 - Fast summation of spherical radial functions
 *
 * Copyright (C) 2005 Jens Keiner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/* Include standard C headers. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

/* Include FFTW header. */
#include <fftw3.h>

/* Include NFFT 3 library header. */
#include "nfft3.h"

/* Include NFFT 3 utilities headers. */
#include "util.h"

/** Enumeration for parameter values */
enum boolean {NO = 0, YES = 1};

/** Enumeration for quadrature grid types */
enum gridtype {GRID_GAUSS_LEGENDRE = 0, GRID_CLENSHAW_CURTIS = 1,
  GRID_HEALPIX = 2, GRID_EQUIDISTRIBUTION_7_1_11 = 3};

/**
 * The main program.
 *
 * \param argc The number of arguments
 * \param argv An array containing the arguments as C-strings
 *
 * \return Exit code
 */
int main (int argc, char **argv)
{
  int tc;                      /**< The index variable for testcases          */
  int tc_max;                  /**< The number of testcases                   */

  int *N;                      /**< The array containing the cut-off degrees  *
                                    \f$N\f$                                   */
  int N_max;                   /**< The maximum cut-off degree \f$N\f$ for the*
                                    current testcase                          */
  int next_power_of_two_N_max; /**< Next greater power of two with respect to *
                                    N_max                                     */
  int *p;                      /**< The array containing the grid size
                                    parameters                                */
  int p_max;                   /**< The maximum grid size parameter           */
  int iN;                      /**< Index variable for cut-off degrees        */
  int iN_max;                  /**< The maximum number of cut-off degrees     */

  int use_nfsft;               /**< Whether to use the NFSFT algorithm or not */
  int use_nfft;                /**< Whether to use the NFFT algorithm or not  */
  int use_fpt;                 /**< Whether to use the FPT algorithm or not   */
  int cutoff;                  /**< The current NFFT cut-off parameter        */
  double threshold;            /**< The current NFSFT threshold parameter     */

  int gridtype;                /**< The type of quadrature grid to be used    */
  int repetitions;             /**< The number of repetitions to be performed */

  double *w;                   /**< The quadrature weights                    */
  double *x;                   /**< The quadrature nodes                      */
  complex *f_ref;              /**< The reference function values             */
  complex *f;                  /**< The function values                       */
  complex *f_hat_ref;          /**< The reference spherical Fourier           *
                                    coefficients                              */
  complex *f_hat;              /**< The spherical Fourier coefficients        */

  nfsft_plan plan;             /**< The NFSFT plan                            */
  nfsft_plan plan_gen;         /**< The NFSFT plan                            */

  double t;                    /**< The computation time needed for a single  *
                                    run                                       */
  double t_avg;                /**< The average computation time needed       */
  double err_infty_avg;        /**< The average error \f$E_\infty\f$          */
  double err_2_avg;            /**< The average error \f$E_2\f$               */

  int i;                       /**< A loop variable                           */
  int j;                       /**< A loop variable                           */
  int k;                       /**< A loop variable                           */
  int n;                       /**< A loop variable                           */
  int d;                       /**< A loop variable                           */

  int m_theta;                 /**< The current number of different           *
                                    colatitudinal angles (for grids)          */
  int m_phi;                   /**< The current number of different           *
                                    longitudinal angles (for grids).          */
  int m_total;                 /**< The total number nodes.                   */
  int m_theta_max;             /**< The maximum number of different           *
                                    colatitudinal angles (for grids).         */
  int m_phi_max;               /**< The maximum number of different           *
                                    longitudinal angles (for grids).          */
  int m_total_max;          /**< The total maximum number of grid nodes.   */
  double *theta;               /**< An array for saving the angles theta of a *
                                    grid                                      */
  double *phi;                 /**< An array for saving the angles phi of a   *
                                    grid                                      */
  fftw_plan fplan;             /**< An FFTW plan for computing Clenshaw-Curtis
                                    quadrature weights                        */
  int nside;                   /**< The size parameter for the HEALPix grid   */
  int d2;
  int gamma;
  double thetai;
  double gammai;
  FILE *file;

  /* Read the number of testcases. */
  fscanf(stdin,"testcases=%d\n",&tc_max);
  fprintf(stdout,"%d\n",tc_max);

  /* Process each testcase. */
  for (tc = 0; tc < tc_max; tc++)
  {
    /* Check if the fast transform shall be used. */
    fscanf(stdin,"nfsft=%d\n",&use_nfsft);
    fprintf(stdout,"%d\n",use_nfsft);
    if (use_nfsft != NO)
    {
      /* Check if the NFFT shall be used. */
      fscanf(stdin,"nfft=%d\n",&use_nfft);
      fprintf(stdout,"%d\n",use_nfsft);
      if (use_nfft != NO)
      {
        /* Read the cut-off parameter. */
        fscanf(stdin,"cutoff=%d\n",&cutoff);
        fprintf(stdout,"%d\n",cutoff);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        cutoff = 1;
      }
      /* Check if the fast polynomial transform shall be used. */
      fscanf(stdin,"fpt=%d\n",&use_fpt);
      fprintf(stdout,"%d\n",use_fpt);
      if (use_fpt != NO)
      {
        /* Read the NFSFT threshold parameter. */
        fscanf(stdin,"threshold=%lf\n",&threshold);
        fprintf(stdout,"%lf\n",threshold);
      }
      else
      {
        /* TODO remove this */
        /* Initialize unused variable with dummy value. */
        threshold = 1000.0;
      }
    }
    else
    {
      /* TODO remove this */
      /* Set dummy values. */
      use_nfft = NO;
      use_fpt = NO;
      cutoff = 3;
      threshold = 1000.0;
    }

    /* Read the quadrature grid type. */
    fscanf(stdin,"gridtype=%d\n",&gridtype);
    fprintf(stdout,"%d\n",gridtype);

    /* Read the number of repetitions. */
    fscanf(stdin,"repetitions=%d\n",&repetitions);
    fprintf(stdout,"%d\n",repetitions);

    /* Initialize maximum cut-off degree and grid size parameter. */
    N_max = 0;
    p_max = 0;

    /* Read the number of cut-off degrees. */
    fscanf(stdin,"bandwidths=%d\n",&iN_max);
    fprintf(stdout,"%d\n",iN_max);

    /* Allocate memory for the cut-off degrees and grid size parameters. */
    N = (int*) malloc(iN_max*sizeof(int));
    p = (int*) malloc(iN_max*sizeof(int));

    /* Read the cut-off degrees and grid size parameters. */
    for (iN = 0; iN < iN_max; iN++)
    {
      /* Read cut-off degree and grid size parameter. */
      fscanf(stdin,"%d %d\n",&N[iN],&p[iN]);
      fprintf(stdout,"%d %d\n",N[iN],p[iN]);
      N_max = MAX(N_max,N[iN]);
      p_max = MAX(p_max,p[iN]);
    }

    /* Determine the maximum number of nodes. */
    switch (gridtype)
    {
      case GRID_GAUSS_LEGENDRE:
        m_theta_max = p_max+1;
        m_phi_max = 2*p_max+2;
        m_total_max = m_theta_max * m_phi_max;
        break;
      case GRID_CLENSHAW_CURTIS:
        m_theta_max = 2*p_max+1;
        m_phi_max = 2*p_max+2;
        m_total_max = m_theta_max * m_phi_max;
        break;
      case GRID_HEALPIX:
        m_theta_max = 1;
        m_phi_max = 12*p_max*p_max;
        m_total_max = m_theta_max * m_phi_max;
        break;
      case GRID_EQUIDISTRIBUTION_7_1_11:
        m_theta_max = 1;
        m_phi_max = 2+4*((int)floor((p_max+1)/2.0))*((int)floor(p_max/2.0));
        m_total_max = m_theta_max * m_phi_max;
        break;
    }

    /* Allocate memory for data structures. */
    w = (double*) malloc(m_theta_max*sizeof(double));
    x = (double*) malloc(2*m_total_max*sizeof(double));
    f_ref = (complex*) malloc(m_total_max*sizeof(complex));
    f = (complex*) malloc(m_total_max*sizeof(complex));
    f_hat_ref = (complex*) malloc(NFSFT_F_HAT_SIZE(N_max)*sizeof(complex));
    f_hat = (complex*) malloc(NFSFT_F_HAT_SIZE(N_max)*sizeof(complex));

    /* Do precomputation. */
    nfsft_precompute(N_max,threshold, 0U |
      ((use_nfsft==NO)?(NFSFT_NO_FAST_ALGORITHM):(0U)));

    /* Process all cut-off bandwidths. */
    for (iN = 0; iN < iN_max; iN++)
    {
      switch (gridtype)
      {
        case GRID_GAUSS_LEGENDRE:
          /* Calculate grid dimensions. */
          m_theta = p[iN] + 1;
          m_phi = 2*p[iN] + 2;
          m_total = m_theta*m_phi;

          /* Read quadrature weights. */
          for (k = 0; k < m_theta; k++)
          {
            fscanf(stdin,"%le\n",&w[k]);
            w[k] *= (2.0*PI)/((double)m_phi);
          }

          /* Allocate memory to store the grid's angles. */
          theta = (double*) malloc(m_theta*sizeof(double));
          phi = (double*) malloc(m_phi*sizeof(double));

          /* Read angles theta. */
          for (k = 0; k < m_theta; k++)
          {
            fscanf(stdin,"%le\n",&theta[k]);
          }

          /* Generate the grid angles phi. */
          for (n = 0; n < m_phi; n++)
          {
            phi[n] = n/((double)m_phi);
            phi[n] -= ((phi[n]>=0.5)?(1.0):(0.0));
          }

          /* Generate the grid's nodes. */
          d = 0;
          for (k = 0; k < m_theta; k++)
          {
            for (n = 0; n < m_phi; n++)
            {
              x[2*d] = phi[n];
              x[2*d+1] = theta[k];
              d++;
            }
          }

          /* Free the arrays for the grid's angles. */
          free(theta);
          free(phi);

          break;

        case GRID_CLENSHAW_CURTIS:
          /* Calculate grid dimensions. */
          m_theta = 2*p[iN] + 1;
          m_phi = 2*p[iN] + 2;
          m_total = m_theta*m_phi;

          /* Allocate memory to store the grid's angles. */
          theta = (double*) malloc(m_theta*sizeof(double));
          phi = (double*) malloc(m_phi*sizeof(double));

          /* Generate the grid angles theta. */
          for (k = 0; k < m_theta; k++)
          {
            theta[k] = k/((double)2*(m_theta-1));
          }

          /* Generate the grid angles phi. */
          for (n = 0; n < m_phi; n++)
          {
            phi[n] = n/((double)m_phi);
            phi[n] -= ((phi[n]>=0.5)?(1.0):(0.0));
          }

          /* Generate quadrature weights. */
          fplan = fftw_plan_r2r_1d(p[iN]+1, w, w, FFTW_REDFT00, 0U);
          for (k = 0; k < p[iN]+1; k++)
          {
            w[k] = -2.0/(4*k*k-1);
          }
          fftw_execute(fplan);
          w[0] *= 0.5;

          for (k = 0; k < p[iN]+1; k++)
          {
            w[k] *= (2.0*PI)/((double)(m_theta-1)*m_phi);
            w[m_theta-1-k] = w[k];
          }
          fftw_destroy_plan(fplan);

          /* Generate the grid's nodes. */
          d = 0;
          for (k = 0; k < m_theta; k++)
          {
            for (n = 0; n < m_phi; n++)
            {
              x[2*d] = phi[n];
              x[2*d+1] = theta[k];
              d++;
            }
          }

          /* Free the arrays for the grid's angles. */
          free(theta);
          free(phi);

          break;

        case GRID_HEALPIX:
          /* Calculate grid dimensions. */
          nside = next_power_of_2(ceil((2.0*N[iN])));
          m_theta = 1;
          m_phi = 12*nside*nside;
          m_total = m_theta*m_phi;

          d = 0;
          for (k = 1; k <= nside-1; k++)
          {
            for (n = 0; n <= 4*k-1; n++)
            {
              x[2*d+1] = 1 - (k*k)/((double)(3.0*nside*nside));
              x[2*d] =  ((n+0.5)/(4*k));
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }

          d2 = d-1;

          for (k = nside; k <= 3*nside; k++)
          {
            for (n = 0; n <= 4*nside-1; n++)
            {
              x[2*d+1] = 2.0/(3*nside)*(2*nside-k);
              x[2*d] = (n+((k%2==0)?(0.5):(0.0)))/(4*nside);
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }

          for (k = 1; k <= nside-1; k++)
          {
            for (n = 0; n <= 4*k-1; n++)
            {
              x[2*d+1] = -x[2*d2+1];
              x[2*d] =  x[2*d2];
              d++;
              d2--;
            }
          }

          for (d = 0; d < m_total; d++)
          {
            x[2*d+1] = acos(x[2*d+1])/(2.0*PI);
          }

          w[0] = (4.0*PI)/(m_total);
          break;

        case GRID_EQUIDISTRIBUTION_7_1_11:
          /* Calculate grid dimensions. */
          gamma = 2*N[iN];
          m_theta = 1;
          m_phi = 2+4*((int)floor((gamma+1)/2.0))*((int)floor(gamma/2.0));
          m_total = m_theta*m_phi;
          x[0] = 0.0;
          x[1] = 0.0;
          d = 1;
          for (k = 1; k <= gamma-1; k++)
          {
            thetai = (k*PI)/gamma;
            gammai = ((k<=(gamma/2.0))?(4*k):(4*(gamma-k)));
            for(n = 1; n <= gammai; n++)
            {
              x[2*d+1] = thetai/(2.0*PI);
              x[2*d] = ((n-0.5)*((2.0*PI)/gammai))/(2.0*PI);
              x[2*d] -= (x[2*d]>=0.5)?(1.0):(0.0);
              d++;
            }
          }
          x[2*d+1] = 0.5;
          x[2*d] = 0.0;
          w[0] = (4.0*PI)/(m_total);
          break;

        default:
          break;
      }

      /* Generate random function samples by sampling a bandlimited function. */
      nfsft_init_guru(&plan_gen,N[iN],m_total,
        NFSFT_NORMALIZED | ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
                      ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)), cutoff);
      plan_gen.f_hat = f_hat_ref;
      plan_gen.x = x;
      plan_gen.f = f_ref;
      nfsft_precompute_x(&plan_gen);
      for (k = 0; k < plan_gen.N_total; k++)
      {
        f_hat_ref[k] = 0.0;
      }
      for (k = 0; k < N[iN]; k++)
      {
        for (n = -k; n <= k; n++)
        {
          f_hat_ref[NFSFT_INDEX(k,n,&plan_gen)] =
          drand48()-0.5 + I*(drand48()-0.5);
        }
      }
      if (use_nfsft != NO)
      {
        /* Execute the NFSFT transformation. */
        nfsft_trafo(&plan_gen);
      }
      else
      {
        /* Execute the direct NDSFT transformation. */
        ndsft_trafo(&plan_gen);
      }
      nfsft_finalize(&plan_gen);

      /* Init transform plan. */
      nfsft_init_guru(&plan,N[iN],m_total, NFSFT_NORMALIZED |
        ((use_nfft!=NO)?(0U):(NFSFT_USE_NDFT)) |
        ((use_fpt!=NO)?(0U):(NFSFT_USE_DPT)), cutoff);
      plan.f_hat = f_hat;
      plan.x = x;
      plan.f = f;
      nfsft_precompute_x(&plan);

      /* Initialize cumulative time variable. */
      t_avg = 0.0;
      err_infty_avg = 0.0;
      err_2_avg = 0.0;

      /* Cycle through all runs. */
      for (i = 0; i < repetitions; i++)
      {
        /* Copy exact funtion values to working array. */
        memcpy(f,f_ref,m_total*sizeof(complex));

        /* Initialize time measurement. */
        t = second();

        /* Multiplication with the quadrature weights. */
        /*fprintf(stderr,"\n");*/
        d = 0;
        for (k = 0; k < m_theta; k++)
        {
          for (n = 0; n < m_phi; n++)
          {
            /*fprintf(stderr,"f_ref[%d] = %le + I*%le,\t f[%d] = %le + I*%le,  \t w[%d] = %le\n",
            d,creal(f_ref[d]),cimag(f_ref[d]),d,creal(f[d]),cimag(f[d]),k,
            w[k]);*/
            f[d] *= w[k];
            d++;
          }
        }

        /*fprintf(stderr,"\n");
        d = 0;
        for (d = 0; d < grid_total; d++)
        {
          fprintf(stderr,"f[%d] = %le + I*%le, theta[%d] = %le, phi[%d] = %le\n",
                  d,creal(f[d]),cimag(f[d]),d,x[2*d+1],d,x[2*d]);
        }*/

        /* Check if the fast NFSFT algorithm shall be tested. */
        if (use_nfsft != NO)
        {
          /* Execute the adjoint NFSFT transformation. */
          nfsft_adjoint(&plan);
        }
        else
        {
          /* Execute the adjoint direct NDSFT transformation. */
          ndsft_adjoint(&plan);
        }

        /* Multiplication with the Fourier-Legendre coefficients. */
        /*for (k = 0; k <= m[im]; k++)
        {
          for (n = -k; n <= k; n++)
          {
            fprintf(stderr,"f_hat[%d,%d] = %le\t + I*%le\n",k,n,
                    creal(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]),
                    cimag(f_hat[NFSFT_INDEX(k,n,&plan_adjoint)]));
          }
        }*/

        if (use_nfsft != NO)
        {
          /* Execute the NFSFT transformation. */
          nfsft_trafo(&plan);
        }
        else
        {
          /* Execute the direct NDSFT transformation. */
          ndsft_trafo(&plan);
        }

        t_avg += second() - t;

        err_infty_avg += error_l_infty_complex(f, f_ref, m_total);
        err_2_avg += error_l_2_complex(f, f_ref, m_total);

        for (d = 0; d < m_total; d++)
        {
          fprintf(stderr,"f_ref[%d] = %le + I*%le,\t f[%d] = %le + I*%le\n",
            d,creal(f_ref[d]),cimag(f_ref[d]),d,creal(f[d]),cimag(f[d]));
        }
      }

      /* Calculate average time needed. */
      t_avg = t_avg/((double)repetitions);

      /* Calculate the average error. */
      err_infty_avg = err_infty_avg/((double)repetitions);

      /* Calculate the average error. */
      err_2_avg = err_2_avg/((double)repetitions);

      /* Print out the error measurements. */
      fprintf(stdout,"%+le %+le %+le\n", t_avg, err_infty_avg, err_2_avg);
      fprintf(stderr,"%d: %4d %4d %+le %+le %+le\n", tc, N[iN], p[iN],
        t_avg, err_infty_avg, err_2_avg);

      /* Finalize the NFSFT plans */
      nfsft_finalize(&plan);
    } /* for (im = 0; im < im_max; im++) - Process all cut-off
       * bandwidths.*/
    fprintf(stderr,"\n");

    /* Delete precomputed data. */
    nfsft_forget();

    /* Free data arrays. */
    free(w);
    free(x);
    free(f_ref);
    free(f);
    free(f_hat_ref);
    free(f_hat);

    /* Free memory for cut-off bandwidths and grid size parameters. */
    free(N);
    free(p);
  } /* for (tc = 0; tc < tc_max; tc++) - Process each testcase. */

  /* Return exit code for successful run. */
  return EXIT_SUCCESS;
}