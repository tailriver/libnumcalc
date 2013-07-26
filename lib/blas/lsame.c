/**
*> \brief \b LSAME
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION LSAME(CA,CB)
* 
*       .. Scalar Arguments ..
*       CHARACTER CA,CB
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CA
*> \verbatim
*>          CA is CHARACTER*1
*> \endverbatim
*>
*> \param[in] CB
*> \verbatim
*>          CB is CHARACTER*1
*>          CA and CB specify the single characters to be compared.
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
*> \ingroup aux_blas
*
*  =====================================================================
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*  -- Translated from FORTRAN by Shinsuke Ogawa July 2013 --
**/

#include "common.h"


logical lsame(
      char CA,
      char CB)
{
      const char ZCODE = 'Z';

/**
*     Test if the characters are equal
**/
      if (CA == CB) return CA == CB;

/**
*     Now test for equivalence if both characters are alphabetic.
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
**/

      if (ZCODE == 90 || ZCODE == 122) {
/**
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
**/
          if (CA >= 97 && CA <= 122) CA = CA - 32;
          if (CB >= 97 && CB <= 122) CB = CB - 32;
      } else if (ZCODE == 233 || ZCODE == 169) {
/**
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
**/
          if ((CA >= 129 && CA <= 137) || (CA >= 145 && CA <= 153) || (CA >= 162 && CA <= 169))
              CA = CA + 64;
          if ((CB >= 129 && CB <= 137) || (CB >= 145 && CB <= 153) || (CB >= 162 && CB <= 169))
              CB = CB + 64;
      } else if (ZCODE == 218 || ZCODE == 250) {
/**
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
**/
          if (CA >= 225 && CA <= 250) CA = CA - 32;
          if (CB >= 225 && CB <= 250) CB = CB - 32;
      }
      return CA == CB;
}
