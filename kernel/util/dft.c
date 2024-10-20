/*
 * Copyright (c) 2024 Konstantin Ryabinin
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

#include "infft.h"

#define MAKE_FFTW_DIMS(hownamy, istride, ostride)                                     \
    fftw_iodim64 *dims = (fftw_iodim64 *)alloca((size_t)rank * sizeof(fftw_iodim64)); \
    NFFT_INT i = rank - 1;                                                            \
    dims[i].n = n[i];                                                                 \
    dims[i].is = dims[i].os = 1;                                                      \
    for (--i; i >= 0; --i)                                                            \
    {                                                                                 \
        dims[i].n = n[i];                                                             \
        dims[i].is = dims[i].os = n[i + 1] * dims[i + 1].is;                          \
    }                                                                                 \
    fftw_iodim64 howmany_dims;                                                        \
    howmany_dims.n = howmany;                                                         \
    howmany_dims.is = istride;                                                        \
    howmany_dims.os = ostride;

fftw_plan Y(plan_dft)(NFFT_INT rank, const NFFT_INT *n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags)
{
    MAKE_FFTW_DIMS(1, 1, 1);
    return FFTW(plan_guru64_dft)((int)rank, dims, 1, &howmany_dims, in, out, sign, flags);
}

fftw_plan Y(plan_dft_1d)(NFFT_INT n0, fftw_complex *in, fftw_complex *out, int sign, unsigned flags)
{
    return Y(plan_dft)(1, &n0, in, out, sign, flags);
}

fftw_plan Y(plan_dft_2d)(NFFT_INT n0, NFFT_INT n1, fftw_complex *in, fftw_complex *out, int sign, unsigned flags)
{
    const NFFT_INT n[2] = { n0, n1 };
    return Y(plan_dft)(2, n, in, out, sign, flags);
}

fftw_plan Y(plan_r2r)(NFFT_INT rank, const NFFT_INT *n, double *in, double *out,
                      const fftw_r2r_kind *kind, unsigned flags)
{
    MAKE_FFTW_DIMS(1, 1, 1);
    return FFTW(plan_guru_r2r)((int)rank, dims, 1, &howmany_dims, in, out, kind, flags);
}

fftw_plan Y(plan_many_r2r)(NFFT_INT rank, const NFFT_INT *n, NFFT_INT howmany,
                           double *in, const NFFT_INT *inembed,
                           NFFT_INT istride, NFFT_INT idist,
                           double *out, const NFFT_INT *onembed,
                           NFFT_INT ostride, NFFT_INT odist,
                           const fftw_r2r_kind *kind, unsigned flags)
{
    // WARNING: it is assumed, that
    // inembed == NULL
    // onembed == NULL
    // idist == 1
    // odist == 1
    // which is the case in NFFT.
    // Guru interface cannot handle these parameters.
    MAKE_FFTW_DIMS(hownamy, istride, ostride);
    return FFTW(plan_guru_r2r)((int)rank, dims, 1, &howmany_dims, in, out, kind, flags);
}
