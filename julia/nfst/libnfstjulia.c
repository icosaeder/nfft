#include "config.h"

#include <stdio.h>

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3.h"
#include "infft.h"

nfst_plan* jnfst_alloc(void){
	nfst_plan* p = nfft_malloc(sizeof(nfst_plan));
	return p;
}

void jnfst_init(nfst_plan* p, NFFT_INT d, NFFT_INT* N, NFFT_INT M, NFFT_INT* n, NFFT_INT m, unsigned int f1, unsigned int f2){
    nfst_init_guru(p,d,N,M,n,m,f1,f2);
}

double* jnfst_set_x(nfst_plan* p, double* X){
	NFFT_INT M = p->M_total;
	NFFT_INT d = p->d;
	NFFT_INT r,c;
	for (r = 0; r < M; r++)
		for (c = 0; c < d; c++)
			p->x[d*r+c] = X[d*r+c];
	nfst_precompute_one_psi(p);
	return p->x;
}

// setting Fourier coefficients and returning pointer for access by Julia
double* jnfst_set_fhat(nfst_plan* p,double* f_hat){
	NFFT_INT n = p->N_total;
	NFFT_INT k;
	for (k=0;k<n;k++)
		p->f_hat[k] = f_hat[k];
	return p->f_hat;
}

// setting values and returning pointer for access by Julia
double* jnfst_set_f(nfst_plan* p,double* f){
	NFFT_INT M = p->M_total;
	NFFT_INT j;
	for (j=0;j<M;j++)
		p->f[j] = f[j];
	return p->f;
}

// nfst trafo, return pointer to values for access by Julia if pointer isn't set
double* jnfst_trafo(nfst_plan* p){
	nfst_trafo(p);
	return p->f;
}

// nfst adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double* jnfst_adjoint(nfst_plan* p){
	nfst_adjoint(p);
	return p->f_hat;
}

// nfst trafo, return pointer to values for access by Julia if pointer isn't set
double* jnfst_trafo_direct(nfst_plan* p){
	nfst_trafo_direct(p);
	return p->f;
}

// nfst adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double* jnfst_adjoint_direct(nfst_plan* p){
	nfst_adjoint_direct(p);
	return p->f_hat;
}

// nfst plan finalizer
void jnfst_finalize(nfst_plan* p){
	nfst_finalize(p);
	nfft_free(p);
}
