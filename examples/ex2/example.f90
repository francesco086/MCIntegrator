MODULE func
   IMPLICIT NONE
   
   CONTAINS
   
      FUNCTION gaussian(x)
         REAL(KIND=8) :: gaussian
         REAL(KIND=8), INTENT(IN) :: x(:)
         gaussian=DEXP(-x(1)*x(1))
      END FUNCTION
      
      FUNCTION decaying_exponential(x)
         REAL(KIND=8) :: decaying_exponential
         REAL(KIND=8), INTENT(IN) :: x(:)
         decaying_exponential=DEXP(-DABS(x(1)))*0.5d0
      END FUNCTION

END MODULE func


PROGRAM ex2
   USE mcintegrator
   USE func
   IMPLICIT NONE
   TYPE (MCI) :: imc1, imc2
   REAL(KIND=8) :: average1, error1, average2, error2
   
   
   PRINT*, "We want to integrate"
   PRINT*, "   f(x) = exp(-x^2)"
   PRINT*, "between -Infinity and +Infinity (result is sqrt(pi)=1.772453850905516)"
   PRINT*, 
   
   PRINT*, 
   PRINT*, "There are two possibilities to obtain this result with MCI:"
   PRINT*, 
   PRINT*, "1) Integrate in a range [-L,L] such that f(-L)=f(L)~=0"
   PRINT*, "   We choose L=10."
   CALL imc1%initialize(1)   !1-dimensional integral
   CALL imc1%setIRange((/-10.d0,10.d0/))   !integral range is [-10,10]
   CALL imc1%setObservable(gaussian)   !integrate the gaussian function
   CALL imc1%integrate(1000000_8,average1,error1)   !integrate by sampling 1000000 points
   PRINT*, "   Result is: ", average1, " +- ", error1
   
   PRINT*, 
   PRINT*, "2) Leave the range as [-Infinity,+Infinity] (as it is by default), and sample from ",&
              "the probability density function (IMPORTANT: it is normalized)"
   PRINT*, "      g(x) = exp(-|x|)/2"
   CALL imc2%initialize(1)   !1-dimensional integral
   CALL imc2%setSamplingFunction(decaying_exponential)   !set the sampling function
   CALL imc2%setObservable(gaussian)   !integrate the gaussian function
   CALL imc2%integrate(1000000_8,average2,error2)   !integrate by sampling 1000000 points
   PRINT*, "   Result is: ", average2, " +- ", error2
   
   CALL imc1%terminate()
   CALL imc2%terminate()

END PROGRAM ex2
