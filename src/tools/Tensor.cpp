/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2021 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Tensor.h"
#include "Exception.h"

#include "lapack/lapack.h"

namespace PLMD {

void TensorGenericAux::local_dsyevr(const char *jobz, const char *range, const char *uplo, Integer *n,
                                    double *a, Integer *lda, double *vl, double *vu, Integer *
                                    il, Integer *iu, double *abstol, Integer *m, double *w,
                                    double *z__, Integer *ldz, Integer *isuppz, double *work,
                                    Integer *lwork, Integer *iwork, Integer *liwork, Integer *info) {
  plumed_lapack_dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z__,ldz,isuppz,work,lwork,iwork,liwork,info);
}


}

