/**
*> \brief \b DGBTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGBTRS + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrs.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGBTRS solves a system of linear equations
*>    A * X = B  or  A**T * X = B
*> with a general band matrix A using the LU factorization computed
*> by DGBTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations.
*>          = 'N':  A * X = B  (No transpose)
*>          = 'T':  A**T* X = B  (Transpose)
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The number of subdiagonals within the band of A.  KL >= 0.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The number of superdiagonals within the band of A.  KU >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*>          Details of the LU factorization of the band matrix A, as
*>          computed by DGBTRF.  U is stored as an upper triangular band
*>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*>          the multipliers used during the factorization are stored in
*>          rows KL+KU+2 to 2*KL+KU+1.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices; for 1 <= i <= N, row i of the matrix was
*>          interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleGBcomputational
*
*  =====================================================================
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*  -- Translated from FORTRAN by Shinsuke Ogawa July 2013 --
**/

#include "common.h"

void dgbtrs(
      char TRANS,
      integer N,
      integer KL,
      integer KU,
      integer NRHS,
      const double* AB,
      integer LDAB,
      const integer* IPIV,
      double* B,
      integer LDB,
      integer* INFO)
{
      const double ONE = 1.0;

      logical LNOTI, NOTRAN;
      integer I, J, KD, L, LM;

/**
*     Test the input parameters.
**/
      *INFO = 0;
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         *INFO = -1;
      } else if ( N < 0) {
         *INFO = -2;
      } else if (KL < 0) {
         *INFO = -3;
      } else if (KU < 0) {
         *INFO = -4;
      } else if (NRHS < 0) {
         *INFO = -5;
      } else if (LDAB < 2*KL+KU+1) {
         *INFO = -7;
      } else if (LDB < MAX( 1, N )) {
         *INFO = -10;
      }
      if ( *INFO != 0 ) {
         xerbla( "DGBTRS", -(*INFO) );
         return;
      }
/**
*     Quick return if possible
**/
      if ( N == 0 || NRHS == 0 )
         return;

      KD = KU + KL + 1;
      LNOTI = KL > 0;

      if ( NOTRAN ) {
/**
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
**/
         if ( LNOTI ) {
            for (J = 1; J <= N - 1; J++) {
               LM = MIN( KL, N-J );
               L = IPIV( J );
               if ( L != J )
                  dswap( NRHS, &B( L, 1 ), LDB, &B( J, 1 ), LDB );
               dger( LM, NRHS, -ONE, &AB( KD+1, J ), 1, &B( J, 1 ), LDB, &B( J+1, 1 ), LDB );
            }
         }

         for (I = 1; I <= NRHS; I++) {
/**
*           Solve U*X = B, overwriting B with X.
**/
            dtbsv( 'U', 'N', 'N', N, KL+KU, AB, LDAB, &B( 1, I ), 1 );
         }

      } else {
/**
*        Solve A**T*X = B.
**/
         for (I = 1; I <= NRHS; I++) {
/**
*           Solve U**T*X = B, overwriting B with X.
**/
            dtbsv( 'U', 'T', 'N', N, KL+KU, AB, LDAB, &B( 1, I ), 1 );
         }
/**
*        Solve L**T*X = B, overwriting B with X.
**/
         if ( LNOTI ) {
            for (J = N - 1; J >= 1; J--) {
               LM = MIN( KL, N-J );
               dgemv( 'T', LM, NRHS, -ONE, &B( J+1, 1 ), LDB, &AB( KD+1, J ), 1, ONE, &B( J, 1 ), LDB );
               L = IPIV( J );
               if ( L != J )
                  dswap( NRHS, &B( L, 1 ), LDB, &B( J, 1 ), LDB );
            }
         }
      }
}
