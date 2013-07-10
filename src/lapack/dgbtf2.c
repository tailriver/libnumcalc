/**
*> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGBTF2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, KL, KU, LDAB, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGBTF2 computes an LU factorization of a real m-by-n band matrix A
*> using partial pivoting with row interchanges.
*>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
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
*> \param[in,out] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*>          On entry, the matrix A in band storage, in rows KL+1 to
*>          2*KL+KU+1; rows 1 to KL of the array need not be set.
*>          The j-th column of A is stored in the j-th column of the
*>          array AB as follows:
*>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*>
*>          On exit, details of the factorization: U is stored as an
*>          upper triangular band matrix with KL+KU superdiagonals in
*>          rows 1 to KL+KU+1, and the multipliers used during the
*>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*>          See below for further details.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*>               has been completed, but the factor U is exactly
*>               singular, and division by zero will occur if it is used
*>               to solve a system of equations.
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
*> \date September 2012
*
*> \ingroup doubleGBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The band storage scheme is illustrated by the following example, when
*>  M = N = 6, KL = 2, KU = 1:
*>
*>  On entry:                       On exit:
*>
*>      *    *    *    +    +    +       *    *    *   u14  u25  u36
*>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*>
*>  Array elements marked * are not used by the routine; elements marked
*>  + need not be set on entry, but are required by the routine to store
*>  elements of U, because of fill-in resulting from the row
*>  interchanges.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*  -- Translated from FORTRAN by Shinsuke Ogawa July 2013 --
**/

#include "common.h"

void dgbtf2(
      integer M,
      integer N,
      integer KL,
      integer KU,
      double* AB,
      integer LDAB,
      integer* IPIV,
      integer* INFO)
{
      const double ONE = 1.0;
      const double ZERO = 0.0;

      integer I, J, JP, JU, KM, KV;

/**
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
**/
      KV = KU + KL;

/**
*     Test the input parameters.
**/
      *INFO = 0;
      if ( M < 0 ) {
         *INFO = -1;
      } else if ( N < 0 ) {
         *INFO = -2;
      } else if ( KL < 0 ) {
         *INFO = -3;
      } else if ( KU < 0 ) {
         *INFO = -4;
      } else if ( LDAB < KL+KV+1 ) {
         *INFO = -6;
      }
      if ( *INFO != 0 ) {
         xerbla( "DGBTF2", -(*INFO) );
         return;
      }

/**
*     Quick return if possible
**/
      if ( M == 0 || N == 0 )
         return;

/**
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
**/
      for (J = KU + 2; J <= MIN( KV, N ); J++) {
         for (I = KV - J + 2; I <= KL; I++) {
            AB( I, J ) = ZERO;
         }
      }

/**
*     JU is the index of the last column affected by the current stage
*     of the factorization.
**/
      JU = 1;

      for (J = 1; J <= MIN( M, N ); J++) {
/**
*        Set fill-in elements in column J+KV to zero.
**/
         if ( J+KV <= N ) {
            for (I = 1; I <= KL; I++) {
               AB( I, J+KV ) = ZERO;
            }
         }
/**
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
**/
         KM = MIN( KL, M-J );
         JP = idamax( KM+1, &AB( KV+1, J ), 1 ) + 1;
         IPIV( J ) = JP + J - 1;
         if ( AB( KV+JP, J ) != ZERO ) {
            JU = MAX( JU, MIN( J+KU+JP-1, N ) );
/**
*           Apply interchange to columns J to JU.
**/
            if ( JP != 1 )
               dswap( JU-J+1, &AB( KV+JP, J ), LDAB-1, &AB( KV+1, J ), LDAB-1 );

            if ( KM > 0 ) {
/**
*              Compute multipliers.
**/
               dscal( KM, ONE / AB( KV+1, J ), &AB( KV+2, J ), 1 );
/**
*              Update trailing submatrix within the band.
**/
               if ( JU > J )
                  dger( KM, JU-J, -ONE, &AB( KV+2, J ), 1, &AB( KV, J+1 ), LDAB-1, &AB( KV+1, J+1 ), LDAB-1 );
            }
         } else {
/**
*           If pivot is zero, set INFO to the index of the pivot
*           unless a zero pivot has already been found.
**/
            if ( *INFO == 0 )
               *INFO = J;
         }
      }
}
