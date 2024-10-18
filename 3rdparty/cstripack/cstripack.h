/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef NFFT3_H
#define NFFT3_H

/* #include <f2c.h> */

NFFT_INT addnod_(NFFT_INT *nst, NFFT_INT *k, double *x, double *y, double *z__,
  NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew, NFFT_INT *ier);

double areas_(double *v1, double *v2, double *v3);

NFFT_INT bdyadd_(NFFT_INT *kk, NFFT_INT *i1, NFFT_INT *i2, NFFT_INT *list, NFFT_INT *lptr,
  NFFT_INT *lend, NFFT_INT *lnew);

NFFT_INT bnodes_(NFFT_INT *n, NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend,
  NFFT_INT *nodes, NFFT_INT *nb, NFFT_INT *na, NFFT_INT *nt);

NFFT_INT circum_(double *v1, double *v2, double *v3, double *c__, NFFT_INT *ier);

NFFT_INT covsph_(NFFT_INT *kk, NFFT_INT *n0, NFFT_INT *list, NFFT_INT *lptr,
  NFFT_INT *lend, NFFT_INT *lnew);

NFFT_INT crlist_(NFFT_INT *n, NFFT_INT *ncol, double *x, double *y, double *z__,
  NFFT_INT *list, NFFT_INT *lend, NFFT_INT *lptr, NFFT_INT *lnew, NFFT_INT *ltri,
  NFFT_INT *listc, NFFT_INT *nb, double *xc, double *yc, double *zc, double *rc,
  NFFT_INT *ier);

NFFT_INT delarc_(NFFT_INT *n, NFFT_INT *io1, NFFT_INT *io2, NFFT_INT * list,
  NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew, NFFT_INT *ier);

NFFT_INT delnb_(NFFT_INT *n0, NFFT_INT *nb, NFFT_INT *n, NFFT_INT *list,
  NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew, NFFT_INT *lph);

NFFT_INT delnod_(NFFT_INT *k, NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT *list,
 NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew, NFFT_INT *lwk, NFFT_INT *iwk,
 NFFT_INT *ier);

NFFT_INT edge_(NFFT_INT *in1, NFFT_INT *in2, double *x, double *y, double *z__, NFFT_INT *lwk,
  NFFT_INT *iwk, NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *ier);

NFFT_INT getnp_(double *x, double *y, double *z__, NFFT_INT *list, NFFT_INT *lptr,
  NFFT_INT *lend, NFFT_INT *l, NFFT_INT *npts, double *df, NFFT_INT *ier);

NFFT_INT insert_(NFFT_INT *k, NFFT_INT *lp, NFFT_INT *list, NFFT_INT *lptr,
  NFFT_INT *lnew);

NFFT_INT inside_(double *p, NFFT_INT *lv, double *xv, double *yv, double *zv, NFFT_INT *
  nv, NFFT_INT *listv, NFFT_INT *ier);

NFFT_INT intadd_(NFFT_INT *kk, NFFT_INT *i1, NFFT_INT *i2, NFFT_INT *i3, NFFT_INT *list,
  NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew);

NFFT_INT intrsc_(double *p1, double *p2, double *cn, double *p, NFFT_INT *ier);

NFFT_INT jrand_(NFFT_INT *n, NFFT_INT *ix, NFFT_INT *iy, NFFT_INT *iz);

NFFT_INT left_(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
  double *x0, double *y0, double *z0);

NFFT_INT lstptr_(NFFT_INT *lpl, NFFT_INT *nb, NFFT_INT *list, NFFT_INT *lptr);

NFFT_INT nbcnt_(NFFT_INT *lpl, NFFT_INT *lptr);

NFFT_INT nearnd_(double *p, NFFT_INT *ist, NFFT_INT *n, double *x, double *y,
  double *z__, NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, double *al);

NFFT_INT optim_(double *x, double *y, double *z__, NFFT_INT *na, NFFT_INT *list,
  NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *nit, NFFT_INT *iwk, NFFT_INT *ier);

NFFT_INT scoord_(double *px, double *py, double *pz, double *plat, double *plon, double *pnrm);

double store_(double *x);

NFFT_INT swap_(NFFT_INT *in1, NFFT_INT *in2, NFFT_INT *io1, NFFT_INT *	io2,
  NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lp21);

NFFT_INT swptst_(NFFT_INT *n1, NFFT_INT *n2, NFFT_INT *n3, NFFT_INT *n4, double *x,
  double *y, double *z__);

NFFT_INT trans_(NFFT_INT *n, double *rlat, double *rlon, double *x, double *y, double *z__);

NFFT_INT trfind_(NFFT_INT *nst, double *p, NFFT_INT *n, double *x, double *y, double *z__,
  NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, double *b1, double *b2,
  double *b3, NFFT_INT *i1, NFFT_INT *i2, NFFT_INT *i3);

NFFT_INT trlist_(NFFT_INT *n, NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend,
  NFFT_INT *nrow, NFFT_INT *nt, NFFT_INT *ltri, NFFT_INT *ier);

NFFT_INT trlprt_(NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT *iflag,
  NFFT_INT *nrow, NFFT_INT *nt, NFFT_INT *ltri, NFFT_INT *lout);

NFFT_INT trmesh_(NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT	*list,
  NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lnew, NFFT_INT *near__, NFFT_INT *next,
  double *dist, NFFT_INT *ier);

NFFT_INT trplot_(NFFT_INT *lun, double *pltsiz, double *elat, double *elon, double *a,
  NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT *list, NFFT_INT *lptr,
  NFFT_INT *lend, char *title, NFFT_INT *numbr, NFFT_INT *ier, short title_len);

NFFT_INT trprnt_(NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT *iflag,
  NFFT_INT *list, NFFT_INT *lptr, NFFT_INT *lend, NFFT_INT *lout);

NFFT_INT vrplot_(NFFT_INT *lun, double *pltsiz, double *elat, double *elon, double *a,
  NFFT_INT *n, double *x, double *y, double *z__, NFFT_INT *nt, NFFT_INT *listc,
  NFFT_INT *lptr, NFFT_INT *lend, double *xc, double *yc, double *zc, char *title,
  NFFT_INT *numbr, NFFT_INT *ier, short title_len);

#endif
