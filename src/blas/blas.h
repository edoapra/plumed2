/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib BLAS library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference BLAS implementation.

The reference BLAS implementation is available from http://www.netlib.org/blas 

BLAS does not come with a formal named "license", but a general statement that 

"The reference BLAS is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified BLAS files as the original netlib versions, so do what you want with them.
However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_blas_blas_h
#define __PLUMED_blas_blas_h
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Header definitions for the standard BLAS library.
 *
 * This is the subset of BLAS routines used for the
 * linear algebra operations in Gromacs.
 * Do NOT use this for other purposes - we only provide this as a
 * simple fallback/reference implementation when no optimized BLAS
 * is present. If you need an implementation for your own code
 * there are several much faster versions out there.
 *
 * All routines are compatible with the BLAS reference implementation,
 * meaning they assume fortran-style matrix row/column organization.
 *
 * There is plenty of documentation for these routines available
 * at http://www.netlib.org/blas , so there is no point in repeating
 * it here.
 */
#ifndef GMX_BLAS_H
#define GMX_BLAS_H

/*! \cond */


/* These are not required by this file, but by the internal BLAS
 * implementation.  In principle, they could be included in each file
 * that requires them, but this is simpler.  Since the header is internal
 * to the linearyalgebra/ module, the added complexity may not be worth it. */
#include "real.h"

#ifndef __PLUMED_BLAS_RETURNS_FLOAT
#define __PLUMED_BLAS_RETURNS_FLOAT float
#endif
#ifdef __PLUMED_HAS_ILP64
typedef long Integer;
#else
typedef int Integer;
#endif
#if ! defined (__PLUMED_HAS_EXTERNAL_BLAS)
#include "def_internal.h"
namespace PLMD{
namespace blas{
#else
namespace PLMD{
namespace blas{
}
}
#include "def_external.h"
extern "C"{
#endif
#if 0
}
#endif

/* Double precision versions */
double
    PLUMED_BLAS_F77_FUNC(dasum, DASUM) (Integer *n, double *dx, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(daxpy, DAXPY) (Integer *n, double *da, double *dx, Integer *incx, double *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(dcopy, DCOPY) (Integer *n, double *dx, Integer *incx, double *dy, Integer *incy);

double
    PLUMED_BLAS_F77_FUNC(ddot, DDOT) (Integer *n, double *dx, Integer *incx, double *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(dgemm, DGEMM) (const char *transa, const char *transb, Integer *m, Integer *n, Integer *k,
                            double *alpha, double *a, Integer *lda, double *b, Integer *ldb,
                            double *beta, double *c, Integer *ldc);

void
    PLUMED_BLAS_F77_FUNC(dgemv, DGEMV) (const char *trans, Integer *m, Integer *n, double *alpha, double *a, Integer *lda,
                            double *x, Integer *incx, double *beta, double *y, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(dger, DGER) (Integer *m, Integer *n, double *alpha, double *x, Integer *incx,
                          double *y, Integer *incy, double *a, Integer *lda);

double
    PLUMED_BLAS_F77_FUNC(dnrm2, DNRM2) (Integer  *n, double *x, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(drot, DROT) (Integer *n, double *dx, Integer *incx,
                          double *dy, Integer *incy, double *c, double *s);

void
    PLUMED_BLAS_F77_FUNC(dscal, DSCAL) (Integer *n, double *fact, double *dx, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(dswap, DSWAP) (Integer *n, double *dx, Integer *incx, double *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(dsymv, DSYMV) (const char *uplo, Integer *n, double *alpha, double *a, Integer *lda,
                            double *x, Integer *incx, double *beta, double *y, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(dsyr2, DSYR2) (const char *uplo, Integer *n, double *alpha, double *x, Integer *incx,
                            double *y, Integer *incy, double *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(dsyr2k, DSYR2K) (const char *uplo, const char *trans, Integer *n, Integer *k, double *alpha, double *a,
                              Integer *lda, double *b, Integer *ldb, double *beta, double *c, Integer *ldc);

void
    PLUMED_BLAS_F77_FUNC(dtrmm, DTRMM) (const char *side, const char *uplo, const char *transa, const char *diag, Integer *m, Integer *n,
                            double *alpha, double *a, Integer *lda, double *b, Integer *ldb);

void
    PLUMED_BLAS_F77_FUNC(dtrmv, DTRMV) (const char *uplo, const char *trans, const char *diag, Integer *n,
                            double *a, Integer *lda, double *x, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(dtrsm, DTRSM) (const char *side, const char *uplo, const char *transa, const char *diag, Integer *m, Integer *n,
                            double *alpha, double *a, Integer *lda, double *b, Integer *ldb);

Integer
    PLUMED_BLAS_F77_FUNC(idamax, IDAMAX) (Integer *n, double *dx, Integer *incx);



/* Single precision versions */
__PLUMED_BLAS_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(sasum, SASUM) (Integer *n, float *dx, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(saxpy, SAXPY) (Integer *n, float *da, float *dx, Integer *incx, float *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(scopy, SCOPY) (Integer *n, float *dx, Integer *incx, float *dy, Integer *incy);

__PLUMED_BLAS_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(sdot, SDOT) (Integer *n, float *dx, Integer *incx, float *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(sgemm, SGEMM) (const char *transa, const char *transb, Integer *m, Integer *n, Integer *k,
                            float *alpha, float *a, Integer *lda, float *b, Integer *ldb,
                            float *beta, float *c, Integer *ldc);

void
    PLUMED_BLAS_F77_FUNC(sgemv, SGEMV) (const char *trans, Integer *m, Integer *n, float *alpha, float *a, Integer *lda,
                            float *x, Integer *incx, float *beta, float *y, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(sger, SGER) (Integer *m, Integer *n, float *alpha, float *x, Integer *incx,
                          float *y, Integer *incy, float *a, Integer *lda);

__PLUMED_BLAS_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(snrm2, SNRM2) (Integer  *n, float *x, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(srot, SROT) (Integer *n, float *dx, Integer *incx,
                          float *dy, Integer *incy, float *c, float *s);

void
    PLUMED_BLAS_F77_FUNC(sscal, SSCAL) (Integer *n, float *fact, float *dx, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(sswap, SSWAP) (Integer *n, float *dx, Integer *incx, float *dy, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(ssymv, SSYMV) (const char *uplo, Integer *n, float *alpha, float *a, Integer *lda,
                            float *x, Integer *incx, float *beta, float *y, Integer *incy);

void
    PLUMED_BLAS_F77_FUNC(ssyr2, SSYR2) (const char *uplo, Integer *n, float *alpha, float *x, Integer *incx,
                            float *y, Integer *incy, float *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(ssyr2k, SSYR2K) (const char *uplo, const char *trans, Integer *n, Integer *k, float *alpha, float *a,
                              Integer *lda, float *b, Integer *ldb, float *beta, float *c, Integer *ldc);

void
    PLUMED_BLAS_F77_FUNC(strmm, STRMM) (const char *side, const char *uplo, const char *transa, const char *diag, Integer *m, Integer *n,
                            float *alpha, float *a, Integer *lda, float *b, Integer *ldb);

void
    PLUMED_BLAS_F77_FUNC(strmv, STRMV) (const char *uplo, const char *trans, const char *diag, Integer *n,
                            float *a, Integer *lda, float *x, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(strsm, STRSM) (const char *side, const char *uplo, const char *transa, const char *diag, Integer *m, Integer *n,
                            float *alpha, float *a, Integer *lda, float *b, Integer *ldb);

Integer
    PLUMED_BLAS_F77_FUNC(isamax, ISAMAX) (Integer *n, float *dx, Integer *incx);


}
#if ! defined (__PLUMED_HAS_EXTERNAL_BLAS)
}
#endif

/*! \endcond */

#endif /* GMX_BLAS_H */
#endif
