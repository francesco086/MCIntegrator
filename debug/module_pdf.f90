MODULE pdf
   IMPLICIT NONE

   REAL (KIND=8), PRIVATE, PARAMETER :: PI=3.141592653589793238462643383279502884197169399375105820974944592d0

   INTEGER, PRIVATE :: ndim

   REAL(KIND=8), PRIVATE :: C_gaussian

CONTAINS

   SUBROUTINE PDF_setNDim(n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      ndim=n
   END SUBROUTINE PDF_setNDim

   SUBROUTINE PDF_setGaussianC(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: C
      C_gaussian=C
   END SUBROUTINE PDF_setGaussianC

   FUNCTION PDF_getGaussianC()
      IMPLICIT NONE
      REAL(KIND=8) :: PDF_getGaussianC
      PDF_getGaussianC=C_gaussian
   END FUNCTION PDF_getGaussianC

   FUNCTION PDF_gaussian(x)
      IMPLICIT NONE
      REAL(KIND=8) :: PDF_gaussian
      REAL(KIND=8), INTENT(IN) :: x(:)

      PDF_gaussian=DEXP(-C_gaussian*DOT_PRODUCT(x,x))/(DSQRT(PI/C_gaussian)**ndim)

   END FUNCTION PDF_gaussian

END MODULE pdf
