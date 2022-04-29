/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
These files are semi-automatic translations by f2c from the original netlib LAPACK library.
The source has been modified to (mostly) use modern C formatting, and to get rid of
compiler warnings. Any errors in doing this should be blamed on the GROMACS developers, and
not the reference LAPACK implementation.

The reference LAPACK implementation is available from http://www.netlib.org/lapack 

LAPACK does not come with a formal named "license", but a general statement saying:

"The reference LAPACK is a freely-available software package. It is available from netlib
via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software
packages (and has been). We only ask that proper credit be given to the authors."

While the rest of GROMACS is LGPL, we think it's only fair to give you the same rights to
our modified LAPACK files as the original netlib versions, so do what you want with them.

However, be warned that we have only tested that they to the right thing in the cases used
in GROMACS (primarily full & sparse matrix diagonalization), so in most cases it is a much
better idea to use the full reference implementation.

Erik Lindahl, 2008-10-07.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_lapack_lapack_h
#define __PLUMED_lapack_lapack_h
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
 * Header definitions for the standard LAPACK library.
 *
 * This is the subset of LAPACK routines used for the
 * linear algebra operations in Gromacs. Most of the execution time
 * will be spent in the BLAS routines, which you hopefully have an
 * optimized version of. Gromacs includes reference implementations
 * of both BLAS and LAPACK so it compiles everywhere, but you should
 * really try to find a vendor or otherwise optimized version at least
 * of BLAS for better performance.
 *
 * Do NOT use this code for other purposes - we only provide this as a
 * simple fallback/reference implementation when no optimized BLAS
 * is present. If you need an implementation for your own code
 * there are several much faster versions out there.
 *
 * All routines are compatible with the LAPACK/BLAS reference implementations,
 * meaning they assume fortran-style matrix row/column organization.
 *
 * There is plenty of documentation for these routines available
 * at http://www.netlib.org/lapack , so there is no point in repeating
 * it here.
 */
#ifndef GMX_LAPACK_H
#define GMX_LAPACK_H

/*! \cond */


/* These are not required by this file, but by the internal LAPACK
 * implementation.  In principle, they could be included in each file
 * that requires them, but this is simpler.  Since the header is internal
 * to the linearyalgebra/ module, the added complexity may not be worth it. */
#include "real.h"

#ifndef __PLUMED_LAPACK_RETURNS_FLOAT
#define __PLUMED_LAPACK_RETURNS_FLOAT float
#endif
#ifdef __PLUMED_HAS_ILP64
#define Integer long
#else
#define Integer int
#endif
#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)
#include "def_internal.h"
namespace PLMD{
namespace lapack{
#else
#include "def_external.h"
extern "C"{
#endif
#if 0
}
#endif
/* Double precision */

void
    PLUMED_BLAS_F77_FUNC(dbdsdc, DBDSDC) (const char *uplo, const char *compq, Integer *n, double *d, double *e, double *u,
                              Integer *ldu, double *vt, Integer *ldvt, double *q, Integer *iq, double *work,
                              Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgetf2, DGETF2) (Integer *m, Integer *n, double *a, Integer *lda, Integer *ipiv, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlamrg, DLAMRG) (Integer *n1, Integer *n2, double *a, Integer *dtrd1, Integer *dtrd2, Integer *index);

void
    PLUMED_BLAS_F77_FUNC(dlarnv, DLARNV) (Integer *idist, Integer *iseed, Integer *n, double *x);

void
    PLUMED_BLAS_F77_FUNC(dlasd0, DLASD0) (Integer *n, Integer *sqre, double *d, double *e, double *u,
                              Integer *ldu, double *vt, Integer *ldvt, Integer *smlsiz, Integer *iwork,
                              double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasda, DLASDA) (Integer *icompq, Integer *smlsiz, Integer *n, Integer *sqre, double *d, double *e,
                              double *u, Integer *ldu, double *vt, Integer *k, double *difl, double *difr,
                              double *z, double *poles, Integer *givptr, Integer *givcol, Integer *ldgcol,
                              Integer *perm, double *givnum, double *c, double *s,
                              double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasq6, DLASQ6) (Integer *i0, Integer *n0, double *z, Integer *pp, double *dmin, double *dmin1,
                              double *dmin2, double *dn, double *dnm1, double *dnm2);

void
    PLUMED_BLAS_F77_FUNC(dorgl2, DORGL2) (Integer *m, Integer *n, Integer *k, double *a, Integer *lda,
                              double *tau, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dbdsqr, DBDSQR) (const char *uplo, Integer *n, Integer *ncvt, Integer *nru, Integer *ncc, double *d,
                              double *e, double *vt, Integer *ldvt, double *u, Integer *ldu,
                              double *c, Integer *ldc, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgetrf, DGETRF) (Integer *m, Integer *n, double *a, Integer *lda, Integer *ipiv, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgetri, DGETRI) (Integer *n, double *a, Integer *lda, Integer *ipiv, double *work,
                              Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgetrs, DGETRS) (const char *trans, Integer *n, Integer *nrhs,   double *a, Integer *lda, Integer *ipiv,
                              double *b, Integer *ldb, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dtrtri, DTRTRI) (const char *uplo, const char *diag, Integer *n, double *a, Integer *lda, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dtrti2, DTRTI2) (const char *uplo, const char *diag, Integer *n, double *a, Integer *lda, Integer *info);

double
    PLUMED_BLAS_F77_FUNC(dlange, DLANGE) (const char *norm, Integer *m, Integer *n, double *a, Integer *lda, double *work);

void
    PLUMED_BLAS_F77_FUNC(dlarrbx, DLARRBX) (Integer *n, double *d, double *l, double *ld, double *lld, Integer *ifirst,
                                Integer *ilast, double *rtol1, double *rtol2, Integer *offset, double *w,
                                double *wgap, double *werr, double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasd1, DLASD1) (Integer *nl, Integer *nr, Integer *sqre, double *d, double *alpha, double *beta,
                              double *u, Integer *ldu, double *vt, Integer *ldvt, Integer *idxq, Integer *iwork,
                              double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasdq, DLASDQ) (const char *uplo, Integer *sqre, Integer *n, Integer *ncvt, Integer *nru, Integer *ncc,
                              double *d, double *e, double *vt, Integer *ldvt, double *u, Integer *ldu,
                              double *c, Integer *ldc, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasr, DLASR) (const char *side, const char *pivot, const char *direct, Integer *m, Integer *n, double *c,
                            double *s, double *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(dorglq, DORGLQ) (Integer *m, Integer *n, Integer *k, double *a, Integer *lda,
                              double *tau, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dormtr, DORMTR) (const char *side, const char *uplo, const char *trans, Integer *m, Integer *n, double *a,
                              Integer *lda, double *tau, double *c, Integer *ldc,
                              double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgebd2, DGEBD2) (Integer *m, Integer *n, double *a, Integer *lda, double *d, double *e,
                              double *tauq, double *taup, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlabrd, DLABRD) (Integer *m, Integer *n, Integer *nb, double *a, Integer *lda, double *d,
                              double *e, double *tauq, double *taup, double *x,
                              Integer *ldx, double *y, Integer *ldy);

double
    PLUMED_BLAS_F77_FUNC(dlanst, DLANST) (const char *norm, Integer *n, double *d, double *e);

double
    PLUMED_BLAS_F77_FUNC(dlansy, DLANSY) (const char *norm, const char *uplo, Integer *n, double *a, Integer *lda, double *work);

void
    PLUMED_BLAS_F77_FUNC(dlarrex, DLARREX) (const char *range, Integer *n, double *vl, double *vu, Integer *il, Integer *iu,
                                double *d, double *e, double *tol, Integer *nsplit,
                                Integer *isplit, Integer *m, double *w, Integer *iblock, Integer *indexw,
                                double *gersch, double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasd2, DLASD2) (Integer *nl, Integer *nr, Integer *sqre, Integer *k, double *d, double *z,
                              double *alpha, double *beta, double *u, Integer *ldu, double *vt,
                              Integer *ldvt, double *dsigma, double *u2, Integer *ldu2, double *vt2,
                              Integer *ldvt2, Integer *idxp, Integer *idx, Integer *idxc,
                              Integer *idxq, Integer *coltyp, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasdt, DLASDT) (Integer *n, Integer *lvl, Integer *nd, Integer *inode, Integer *ndiml,
                              Integer *ndimr, Integer *msub);

void
    PLUMED_BLAS_F77_FUNC(dlasrt, DLASRT) (const char *id, Integer *n, double *d, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasrt2, DLASRT2) (const char *id, Integer *n, double *d, Integer *key, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ilasrt2, ILASRT2) (const char *id, Integer *n, Integer *d, Integer *key, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dorgqr, DORGQR) (Integer *m, Integer *n, Integer *k, double *a, Integer *lda, double *tau,
                              double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dstebz, DSTEBZ) (const char *range, const char *order, Integer *n, double *vl, double *vu,
                              Integer *il, Integer *iu, double *abstol, double *d, double *e,
                              Integer *m, Integer *nsplit, double *w, Integer *iblock, Integer *isplit,
                              double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dsteqr, DSTEQR) (const char *compz, Integer *n, double *d__, double *e,
                              double *z__,  Integer *ldz, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgebrd, DGEBRD) (Integer *m, Integer *n, double *a, Integer *lda, double *d, double *e,
                              double *tauq, double *taup, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlacpy, DLACPY) (const char *uplo, Integer *m, Integer *n, double *a, Integer *lda, double *b, Integer *ldb);

double
    PLUMED_BLAS_F77_FUNC(dlapy2, DLAPY2) (double * x, double * y);


void
    PLUMED_BLAS_F77_FUNC(dlarrfx, DLARRFX) (Integer *n, double *d, double *l, double *ld, double *lld, Integer *ifirst,
                                Integer *ilast, double *w, double *sigma, double *dplus, double *lplus,
                                double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasd3, DLASD3) (Integer *nl, Integer *nr, Integer *sqre, Integer *k, double *d, double *q, Integer *ldq,
                              double *dsigma, double *u, Integer *ldu, double *u2, Integer *ldu2,
                              double *vt, Integer *ldvt, double *vt2, Integer *ldvt2, Integer *idxc,
                              Integer *ctot, double *z, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlaset, DLASET) (const char *uplo, Integer *m, Integer *n, double *alpha,
                              double *beta, double *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(dlassq, DLASSQ) (Integer *n, double *x, Integer *incx, double *scale, double *sumsq);

void
    PLUMED_BLAS_F77_FUNC(dorm2l, DORM2L) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, double *a, Integer *lda,
                              double *tau, double *c, Integer *ldc, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dstegr, DSTEGR) (const char *jobz, const char *range, Integer *n, double *d, double *e, double *vl,
                              double *vu, Integer *il, Integer *iu, double *abstol, Integer *m, double *w,
                              double *z, Integer *ldz, Integer *isuppz, double *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ssteqr, SSTEQR) (const char *compz, Integer *n, float *d__, float *e,
                              float *z__,  Integer *ldz, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgelq2, DGELQ2) (Integer *m, Integer *n, double *a, Integer *lda, double *tau, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlae2, DLAE2) (double *a, double *b, double *c, double *rt1, double *rt2);

void
    PLUMED_BLAS_F77_FUNC(dlaev2, DLAEV2) (double *a, double *b, double *c, double *rt1, double *rt2,
                              double *cs1, double *cs2);

void
    PLUMED_BLAS_F77_FUNC(dlar1vx, DLAR1VX) (Integer *n, Integer *b1, Integer *bn, double *sigma, double *d, double *l, double *ld,
                                double *lld, double *eval, double *gersch, double *z, double *ztz, double *mingma,
                                Integer *r, Integer *isuppz, double *work);

void
    PLUMED_BLAS_F77_FUNC(dlarrvx, DLARRVX) (Integer *n, double *d, double *l, Integer *isplit, Integer *m, double *w,
                                Integer *iblock, Integer *indexw, double *gersch, double *tol, double *z, Integer *ldz,
                                Integer *isuppz, double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasd4, DLASD4) (Integer *n, Integer *i, double *d, double *z, double *delta,
                              double *rho, double *sigma, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasq1, DLASQ1) (Integer *n, double *d, double *e, double *work, Integer *info);


void
    PLUMED_BLAS_F77_FUNC(dlasv2, DLASV2) (double *f, double *g, double *h, double *ssmin, double *ssmax,
                              double *snr, double *csr, double *snl, double *csl);

void
    PLUMED_BLAS_F77_FUNC(dorm2r, DORM2R) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, double *a,
                              Integer *lda, double *tau, double *c, Integer *ldc, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dstein, DSTEIN) (Integer *n, double *d, double *e, Integer *m, double *w, Integer *iblock, Integer *isplit,
                              double *z, Integer *ldz, double *work, Integer *iwork, Integer *ifail, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgelqf, DGELQF) (Integer *m, Integer *n, double *a, Integer *lda, double *tau,
                              double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlaebz, DLAEBZ) (Integer *ijob, Integer *nitmax, Integer *n, Integer *mmax, Integer *minp, Integer *nbmin,
                              double *abstol, double *reltol, double *pivmin, double *d, double *e,
                              double *e2, Integer *nval, double *ab, double *c, Integer *mout, Integer *nab,
                              double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlarf, DLARF) (const char *side, Integer *m, Integer *n, double *v, Integer *incv, double *tau,
                            double *c, Integer *ldc, double *work);

void
    PLUMED_BLAS_F77_FUNC(dlartg, DLARTG) (double *f, double *g, double *cs, double *sn, double *r);

void
    PLUMED_BLAS_F77_FUNC(dlasd5, DLASD5) (Integer *i, double *d, double *z, double *delta,
                              double *rho, double *dsigma, double *work);

void
    PLUMED_BLAS_F77_FUNC(dlasq2, DLASQ2) (Integer *n, double *z, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasq3, DLASQ3) (Integer *i0, Integer *n0, double *z, Integer *pp, double *dmin,
                              double *sigma, double *desig, double *qmax, Integer *nfail,
                              Integer *iter, Integer *ndiv, Integer *ieee);

void
    PLUMED_BLAS_F77_FUNC(dlaswp, DLASWP) (Integer *n, double *a, Integer *lda, Integer *k1, Integer *k2, Integer *ipiv, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(dormbr, DORMBR) (const char *vect, const char *side, const char *trans, Integer *m, Integer *n, Integer *k,
                              double *a, Integer *lda, double *tau, double *c, Integer *ldc, double *work,
                              Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dsterf, DSTERF) (Integer *n, double *d, double *e, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgeqr2, DGEQR2) (Integer *m, Integer *n, double *a, Integer *lda, double *tau,
                              double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlaed6, DLAED6) (Integer *kniter, Integer *orgati, double *rho, double *d,
                              double *z, double *finit, double *tau, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlarfb, DLARFB) (const char *side, const char *trans, const char *direct, const char *storev, Integer *m, Integer *n,
                              Integer *k, double *v, Integer *ldv, double *t, Integer *ldt, double *c,
                              Integer *ldc, double *work, Integer *ldwork);

void
    PLUMED_BLAS_F77_FUNC(dlaruv, DLARUV) (Integer *iseed, Integer *n, double *x);

void
    PLUMED_BLAS_F77_FUNC(dlasd6, DLASD6) (Integer *icompq, Integer *nl, Integer *nr, Integer *sqre, double *d, double *vf,
                              double *vl, double *alpha, double *beta, Integer *idxq, Integer *perm,
                              Integer *givptr, Integer *givcol, Integer *ldgcol, double *givnum, Integer *ldgnum,
                              double *poles, double *difl, double *difr, double *z, Integer *k,
                              double *c, double *s, double *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlatrd, DLATRD) (const char *uplo, Integer *n, Integer *nb, double *a, Integer *lda, double *e,
                              double * tau, double *w, Integer *ldw);

void
    PLUMED_BLAS_F77_FUNC(dorml2, DORML2) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, double *a,
                              Integer *lda, double *tau, double *c, Integer *ldc, double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dstevr, DSTEVR) (const char *jobz, const char *range, Integer *n, double *d, double *e, double *vl,
                              double *vu, Integer *il, Integer *iu, double *abstol, Integer *m, double *w,
                              double *z, Integer *ldz, Integer *isuppz, double *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dsytrd, DSYTRD) (const char *uplo, Integer *n, double *  a, Integer *lda, double *d,
                              double *e, double *tau, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dsyevr, DSYEVR) (const char *jobz, const char *range, const char *uplo, Integer *n,
                              double *a, Integer *lda, double *vl, double *vu, Integer *
                              il, Integer *iu, double *abstol, Integer *m, double *w,
                              double *z__, Integer *ldz, Integer *isuppz, double *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dormql, DORMQL) (const char *side, const char *trans, Integer *m, Integer *n,
                              Integer *k, double *a, Integer *lda, double *tau, double *
                              c, Integer *ldc, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dormqr, DORMQR) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, double *a,
                              Integer *lda, double *tau, double *c, Integer *ldc,
                              double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dorgbr, DORGBR) (const char *vect, Integer *m, Integer *n, Integer *k, double *a, Integer *lda,
                              double *tau, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasq5, DLASQ5) (Integer *i0, Integer *n0, double *z, Integer *pp, double *tau, double *dmin,
                              double *dmin1, double *dmin2, double *dn, double *dnm1,
                              double *dnm2, Integer *ieee);

void
    PLUMED_BLAS_F77_FUNC(dlasd8, DLASD8) (Integer *icompq, Integer *k, double *d, double *z, double *vf, double *vl,
                              double *difl, double *difr, Integer *lddifr, double *dsigma,
                              double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlascl, DLASCL) (const char *type, Integer *kl, Integer *ku, double *cfrom, double *cto, Integer *m,
                              Integer *n, double *a, Integer *lda, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlarft, DLARFT) (const char *direct, const char *storev, Integer *n, Integer *k, double *v,
                              Integer *ldv, double *tau, double *t, Integer *ldt);

void
    PLUMED_BLAS_F77_FUNC(dlagts, DLAGTS) (Integer *job, Integer *n, double *a, double *b, double *c, double *d,
                              Integer *in, double *y, double *tol, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgesdd, DGESDD) (const char *jobz, Integer *m, Integer *n, double *a, Integer *lda, double *s, double *u,
                              Integer *ldu, double *vt, Integer *ldvt, double *work, Integer *lwork,
                              Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dsytd2, DSYTD2) (const char *uplo, Integer *n, double *a, Integer *lda, double *d,
                              double *e, double *tau, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dormlq, DORMLQ) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, double *a, Integer *lda,
                              double *tau, double *c, Integer *ldc, double *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dorg2r, DORG2R) (Integer *m, Integer *n, Integer *k, double *a, Integer *lda, double *tau,
                              double *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlasq4, DLASQ4) (Integer *i0, Integer *n0, double *z, Integer *pp, Integer *n0in, double *dmin,
                              double *dmin1, double *dmin2, double *dn, double *dn1,
                              double *dn2, double *tau, Integer *ttype);

void
    PLUMED_BLAS_F77_FUNC(dlasd7, DLASD7) (Integer *icompq, Integer *nl, Integer *nr, Integer *sqre, Integer *k, double *d, double *z,
                              double *zw, double *vf, double *vfw, double *vl, double *vlw,
                              double *alpha, double *beta, double *dsigma, Integer *idx, Integer *idxp,
                              Integer *idxq, Integer *perm, Integer *givptr, Integer *givcol, Integer *ldgcol,
                              double *givnum, Integer *ldgnum, double *c, double *s, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dlas2, DLAS2) (double *f, double *g, double *h, double *ssmin, double *ssmax);

void
    PLUMED_BLAS_F77_FUNC(dlarfg, DLARFG) (Integer *n, double *alpha, double *x, Integer *incx, double *tau);

void
    PLUMED_BLAS_F77_FUNC(dlagtf, DLAGTF) (Integer *n, double *a, double *lambda, double *b, double *c,
                              double *tol, double *d, Integer *in, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(dgeqrf, DGEQRF) (Integer *m, Integer *n, double *a, Integer *lda, double *tau,
                              double *work, Integer *lwork, Integer *info);



/* Single precision */

void
    PLUMED_BLAS_F77_FUNC(sbdsdc, SBDSDC) (const char *uplo, const char *compq, Integer *n, float *d, float *e, float *u,
                              Integer *ldu, float *vt, Integer *ldvt, float *q, Integer *iq, float *work,
                              Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgetf2, SGETF2) (Integer *m, Integer *n, float *a, Integer *lda, Integer *ipiv, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slamrg, SLAMRG) (Integer *n1, Integer *n2, float *a, Integer *dtrd1, Integer *dtrd2, Integer *index);

void
    PLUMED_BLAS_F77_FUNC(slarnv, SLARNV) (Integer *idist, Integer *iseed, Integer *n, float *x);

void
    PLUMED_BLAS_F77_FUNC(slasd0, SLASD0) (Integer *n, Integer *sqre, float *d, float *e, float *u,
                              Integer *ldu, float *vt, Integer *ldvt, Integer *smlsiz, Integer *iwork,
                              float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasda, SLASDA) (Integer *icompq, Integer *smlsiz, Integer *n, Integer *sqre, float *d, float *e,
                              float *u, Integer *ldu, float *vt, Integer *k, float *difl, float *difr,
                              float *z, float *poles, Integer *givptr, Integer *givcol, Integer *ldgcol,
                              Integer *perm, float *givnum, float *c, float *s,
                              float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasq6, SLASQ6) (Integer *i0, Integer *n0, float *z, Integer *pp, float *dmin, float *dmin1,
                              float *dmin2, float *dn, float *dnm1, float *dnm2);

void
    PLUMED_BLAS_F77_FUNC(sorgl2, SORGL2) (Integer *m, Integer *n, Integer *k, float *a, Integer *lda,
                              float *tau, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sbdsqr, SBDSQR) (const char *uplo, Integer *n, Integer *ncvt, Integer *nru, Integer *ncc, float *d,
                              float *e, float *vt, Integer *ldvt, float *u, Integer *ldu,
                              float *c, Integer *ldc, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgetrf, SGETRF) (Integer *m, Integer *n, float *a, Integer *lda, Integer *ipiv, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgetri, SGETRI) (Integer *n, float *a, Integer *lda, Integer *ipiv, float *work,
                              Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgetrs, SGETRS) (const char *trans, Integer *n, Integer *nrhs,   float *a, Integer *lda, Integer *ipiv,
                              float *b, Integer *ldb, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(strtri, STRTRI) (const char *uplo, const char *diag, Integer *n, float *a, Integer *lda, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(strti2, STRTI2) (const char *uplo, const char *diag, Integer *n, float *a, Integer *lda, Integer *info);

__PLUMED_LAPACK_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(slange, SLANGE) (const char *norm, Integer *m, Integer *n, float *a, Integer *lda, float *work);

void
    PLUMED_BLAS_F77_FUNC(slarrbx, SLARRBX) (Integer *n, float *d, float *l, float *ld, float *lld, Integer *ifirst,
                                Integer *ilast, float *rtol1, float *rtol2, Integer *offset, float *w,
                                float *wgap, float *werr, float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasd1, SLASD1) (Integer *nl, Integer *nr, Integer *sqre, float *d, float *alpha, float *beta,
                              float *u, Integer *ldu, float *vt, Integer *ldvt, Integer *idxq, Integer *iwork,
                              float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasdq, SLASDQ) (const char *uplo, Integer *sqre, Integer *n, Integer *ncvt, Integer *nru, Integer *ncc,
                              float *d, float *e, float *vt, Integer *ldvt, float *u, Integer *ldu,
                              float *c, Integer *ldc, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasr, SLASR) (const char *side, const char *pivot, const char *direct, Integer *m, Integer *n, float *c,
                            float *s, float *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(sorglq, SORGLQ) (Integer *m, Integer *n, Integer *k, float *a, Integer *lda,
                              float *tau, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sormtr, SORMTR) (const char *side, const char *uplo, const char *trans, Integer *m, Integer *n, float *a,
                              Integer *lda, float *tau, float *c, Integer *ldc,
                              float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgebd2, SGEBD2) (Integer *m, Integer *n, float *a, Integer *lda, float *d, float *e,
                              float *tauq, float *taup, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slabrd, SLABRD) (Integer *m, Integer *n, Integer *nb, float *a, Integer *lda, float *d,
                              float *e, float *tauq, float *taup, float *x,
                              Integer *ldx, float *y, Integer *ldy);

__PLUMED_LAPACK_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(slanst, SLANST) (const char *norm, Integer *n, float *d, float *e);

__PLUMED_LAPACK_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(slansy, SLANSY) (const char *norm, const char *uplo, Integer *n, float *a, Integer *lda, float *work);

void
    PLUMED_BLAS_F77_FUNC(slarrex, SLARREX) (const char *range, Integer *n, float *vl, float *vu, Integer *il, Integer *iu,
                                float *d, float *e, float *tol, Integer *nsplit,
                                Integer *isplit, Integer *m, float *w, Integer *iblock, Integer *indexw,
                                float *gersch, float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasd2, SLASD2) (Integer *nl, Integer *nr, Integer *sqre, Integer *k, float *d, float *z,
                              float *alpha, float *beta, float *u, Integer *ldu, float *vt,
                              Integer *ldvt, float *dsigma, float *u2, Integer *ldu2, float *vt2,
                              Integer *ldvt2, Integer *idxp, Integer *idx, Integer *idxc,
                              Integer *idxq, Integer *coltyp, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasdt, SLASDT) (Integer *n, Integer *lvl, Integer *nd, Integer *inode, Integer *ndiml,
                              Integer *ndimr, Integer *msub);

void
    PLUMED_BLAS_F77_FUNC(slasrt, SLASRT) (const char *id, Integer *n, float *d, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasrt2, SLASRT2) (const char *id, Integer *n, float *d, Integer *key, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sorgqr, SORGQR) (Integer *m, Integer *n, Integer *k, float *a, Integer *lda, float *tau,
                              float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sstebz, SSTEBZ) (const char *range, const char *order, Integer *n, float *vl, float *vu,
                              Integer *il, Integer *iu, float *abstol, float *d, float *e,
                              Integer *m, Integer *nsplit, float *w, Integer *iblock, Integer *isplit,
                              float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgebrd, SGEBRD) (Integer *m, Integer *n, float *a, Integer *lda, float *d, float *e,
                              float *tauq, float *taup, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slacpy, SLACPY) (const char *uplo, Integer *m, Integer *n, float *a, Integer *lda, float *b, Integer *ldb);

__PLUMED_LAPACK_RETURNS_FLOAT
    PLUMED_BLAS_F77_FUNC(slapy2, SLAPY2) (float * x, float * y);

void
    PLUMED_BLAS_F77_FUNC(slarrfx, SLARRFX) (Integer *n, float *d, float *l, float *ld, float *lld, Integer *ifirst,
                                Integer *ilast, float *w, float *sigma, float *dplus, float *lplus,
                                float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasd3, SLASD3) (Integer *nl, Integer *nr, Integer *sqre, Integer *k, float *d, float *q, Integer *ldq,
                              float *dsigma, float *u, Integer *ldu, float *u2, Integer *ldu2,
                              float *vt, Integer *ldvt, float *vt2, Integer *ldvt2, Integer *idxc,
                              Integer *ctot, float *z, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slaset, SLASET) (const char *uplo, Integer *m, Integer *n, float *alpha,
                              float *beta, float *a, Integer *lda);

void
    PLUMED_BLAS_F77_FUNC(slassq, SLASSQ) (Integer *n, float *x, Integer *incx, float *scale, float *sumsq);

void
    PLUMED_BLAS_F77_FUNC(sorm2l, SORM2L) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, float *a, Integer *lda,
                              float *tau, float *c, Integer *ldc, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sstegr, SSTEGR) (const char *jobz, const char *range, Integer *n, float *d, float *e, float *vl,
                              float *vu, Integer *il, Integer *iu, float *abstol, Integer *m, float *w,
                              float *z, Integer *ldz, Integer *isuppz, float *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgelq2, SGELQ2) (Integer *m, Integer *n, float *a, Integer *lda, float *tau, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slae2, SLAE2) (float *a, float *b, float *c, float *rt1, float *rt2);

void
    PLUMED_BLAS_F77_FUNC(slaev2, SLAEV2) (float *a, float *b, float *c, float *rt1, float *rt2,
                              float *cs1, float *cs2);

void
    PLUMED_BLAS_F77_FUNC(slar1vx, SLAR1VX) (Integer *n, Integer *b1, Integer *bn, float *sigma, float *d, float *l, float *ld,
                                float *lld, float *eval, float *gersch, float *z, float *ztz, float *mingma,
                                Integer *r, Integer *isuppz, float *work);

void
    PLUMED_BLAS_F77_FUNC(slarrvx, SLARRVX) (Integer *n, float *d, float *l, Integer *isplit, Integer *m, float *w,
                                Integer *iblock, Integer *indexw, float *gersch, float *tol, float *z, Integer *ldz,
                                Integer *isuppz, float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasd4, SLASD4) (Integer *n, Integer *i, float *d, float *z, float *delta,
                              float *rho, float *sigma, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasq1, SLASQ1) (Integer *n, float *d, float *e, float *work, Integer *info);


void
    PLUMED_BLAS_F77_FUNC(slasv2, SLASV2) (float *f, float *g, float *h, float *ssmin, float *ssmax,
                              float *snr, float *csr, float *snl, float *csl);

void
    PLUMED_BLAS_F77_FUNC(sorm2r, SORM2R) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, float *a,
                              Integer *lda, float *tau, float *c, Integer *ldc, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sstein, SSTEIN) (Integer *n, float *d, float *e, Integer *m, float *w, Integer *iblock, Integer *isplit,
                              float *z, Integer *ldz, float *work, Integer *iwork, Integer *ifail, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgelqf, SGELQF) (Integer *m, Integer *n, float *a, Integer *lda, float *tau,
                              float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slaebz, SLAEBZ) (Integer *ijob, Integer *nitmax, Integer *n, Integer *mmax, Integer *minp, Integer *nbmin,
                              float *abstol, float *reltol, float *pivmin, float *d, float *e,
                              float *e2, Integer *nval, float *ab, float *c, Integer *mout, Integer *nab,
                              float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slarf, SLARF) (const char *side, Integer *m, Integer *n, float *v, Integer *incv, float *tau,
                            float *c, Integer *ldc, float *work);

void
    PLUMED_BLAS_F77_FUNC(slartg, SLARTG) (float *f, float *g, float *cs, float *sn, float *r);

void
    PLUMED_BLAS_F77_FUNC(slasd5, SLASD5) (Integer *i, float *d, float *z, float *delta,
                              float *rho, float *dsigma, float *work);

void
    PLUMED_BLAS_F77_FUNC(slasq2, SLASQ2) (Integer *n, float *z, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasq3, SLASQ3) (Integer *i0, Integer *n0, float *z, Integer *pp, float *dmin,
                              float *sigma, float *desig, float *qmax, Integer *nfail,
                              Integer *iter, Integer *ndiv, Integer *ieee);

void
    PLUMED_BLAS_F77_FUNC(slaswp, SLASWP) (Integer *n, float *a, Integer *lda, Integer *k1, Integer *k2, Integer *ipiv, Integer *incx);

void
    PLUMED_BLAS_F77_FUNC(sormbr, SORMBR) (const char *vect, const char *side, const char *trans, Integer *m, Integer *n, Integer *k,
                              float *a, Integer *lda, float *tau, float *c, Integer *ldc, float *work,
                              Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ssterf, SSTERF) (Integer *n, float *d, float *e, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgeqr2, SGEQR2) (Integer *m, Integer *n, float *a, Integer *lda, float *tau,
                              float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slaed6, SLAED6) (Integer *kniter, Integer *orgati, float *rho, float *d,
                              float *z, float *finit, float *tau, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slarfb, SLARFB) (const char *side, const char *trans, const char *direct, const char *storev, Integer *m, Integer *n,
                              Integer *k, float *v, Integer *ldv, float *t, Integer *ldt, float *c,
                              Integer *ldc, float *work, Integer *ldwork);

void
    PLUMED_BLAS_F77_FUNC(slaruv, SLARUV) (Integer *iseed, Integer *n, float *x);

void
    PLUMED_BLAS_F77_FUNC(slasd6, SLASD6) (Integer *icompq, Integer *nl, Integer *nr, Integer *sqre, float *d, float *vf,
                              float *vl, float *alpha, float *beta, Integer *idxq, Integer *perm,
                              Integer *givptr, Integer *givcol, Integer *ldgcol, float *givnum, Integer *ldgnum,
                              float *poles, float *difl, float *difr, float *z, Integer *k,
                              float *c, float *s, float *work, Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slatrd, SLATRD) (const char *uplo, Integer *n, Integer *nb, float *a, Integer *lda, float *e,
                              float * tau, float *w, Integer *ldw);

void
    PLUMED_BLAS_F77_FUNC(sorml2, SORML2) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, float *a,
                              Integer *lda, float *tau, float *c, Integer *ldc, float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sstevr, SSTEVR) (const char *jobz, const char *range, Integer *n, float *d, float *e, float *vl,
                              float *vu, Integer *il, Integer *iu, float *abstol, Integer *m, float *w,
                              float *z, Integer *ldz, Integer *isuppz, float *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ssytrd, SSYTRD) (const char *uplo, Integer *n, float *  a, Integer *lda, float *d,
                              float *e, float *tau, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ssyevr, SSYEVR) (const char *jobz, const char *range, const char *uplo, Integer *n,
                              float *a, Integer *lda, float *vl, float *vu, Integer *
                              il, Integer *iu, float *abstol, Integer *m, float *w,
                              float *z__, Integer *ldz, Integer *isuppz, float *work,
                              Integer *lwork, Integer *iwork, Integer *liwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sormql, SORMQL) (const char *side, const char *trans, Integer *m, Integer *n,
                              Integer *k, float *a, Integer *lda, float *tau, float *
                              c, Integer *ldc, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sormqr, SORMQR) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, float *a,
                              Integer *lda, float *tau, float *c, Integer *ldc,
                              float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sorgbr, SORGBR) (const char *vect, Integer *m, Integer *n, Integer *k, float *a, Integer *lda,
                              float *tau, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasq5, SLASQ5) (Integer *i0, Integer *n0, float *z, Integer *pp, float *tau, float *dmin,
                              float *dmin1, float *dmin2, float *dn, float *dnm1,
                              float *dnm2, Integer *ieee);

void
    PLUMED_BLAS_F77_FUNC(slasd8, SLASD8) (Integer *icompq, Integer *k, float *d, float *z, float *vf, float *vl,
                              float *difl, float *difr, Integer *lddifr, float *dsigma,
                              float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slascl, SLASCL) (const char *type, Integer *kl, Integer *ku, float *cfrom, float *cto, Integer *m,
                              Integer *n, float *a, Integer *lda, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slarft, SLARFT) (const char *direct, const char *storev, Integer *n, Integer *k, float *v,
                              Integer *ldv, float *tau, float *t, Integer *ldt);

void
    PLUMED_BLAS_F77_FUNC(slagts, SLAGTS) (Integer *job, Integer *n, float *a, float *b, float *c, float *d,
                              Integer *in, float *y, float *tol, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgesdd, SGESDD) (const char *jobz, Integer *m, Integer *n, float *a, Integer *lda, float *s, float *u,
                              Integer *ldu, float *vt, Integer *ldvt, float *work, Integer *lwork,
                              Integer *iwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(ssytd2, SSYTD2) (const char *uplo, Integer *n, float *a, Integer *lda, float *d,
                              float *e, float *tau, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sormlq, SORMLQ) (const char *side, const char *trans, Integer *m, Integer *n, Integer *k, float *a, Integer *lda,
                              float *tau, float *c, Integer *ldc, float *work, Integer *lwork, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sorg2r, SORG2R) (Integer *m, Integer *n, Integer *k, float *a, Integer *lda, float *tau,
                              float *work, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slasq4, SLASQ4) (Integer *i0, Integer *n0, float *z, Integer *pp, Integer *n0in, float *dmin,
                              float *dmin1, float *dmin2, float *dn, float *dn1,
                              float *dn2, float *tau, Integer *ttype);

void
    PLUMED_BLAS_F77_FUNC(slasd7, SLASD7) (Integer *icompq, Integer *nl, Integer *nr, Integer *sqre, Integer *k, float *d, float *z,
                              float *zw, float *vf, float *vfw, float *vl, float *vlw,
                              float *alpha, float *beta, float *dsigma, Integer *idx, Integer *idxp,
                              Integer *idxq, Integer *perm, Integer *givptr, Integer *givcol, Integer *ldgcol,
                              float *givnum, Integer *ldgnum, float *c, float *s, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(slas2, SLAS2) (float *f, float *g, float *h, float *ssmin, float *ssmax);

void
    PLUMED_BLAS_F77_FUNC(slarfg, SLARFG) (Integer *n, float *alpha, float *x, Integer *incx, float *tau);

void
    PLUMED_BLAS_F77_FUNC(slagtf, SLAGTF) (Integer *n, float *a, float *lambda, float *b, float *c,
                              float *tol, float *d, Integer *in, Integer *info);

void
    PLUMED_BLAS_F77_FUNC(sgeqrf, SGEQRF) (Integer *m, Integer *n, float *a, Integer *lda, float *tau,
                              float *work, Integer *lwork, Integer *info);


}
#if ! defined(__PLUMED_HAS_EXTERNAL_LAPACK)
}
#endif

/*! \endcond */

#endif /* GMX_LAPACK_H */
#endif
