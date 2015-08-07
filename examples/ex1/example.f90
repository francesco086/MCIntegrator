PROGRAM ex1
   USE mcintegrator
   IMPLICIT NONE
   INTERFACE
      FUNCTION x3(x)
         REAL(KIND=8) :: x3
         REAL(KIND=8), INTENT(IN) :: x(:)
      END FUNCTION x3
   END INTERFACE
   TYPE(MCI) :: imc
   REAL(KIND=8) :: average, error

   PRINT *, "We want to numerically integrate the function"
   PRINT *, "  f(x) = x^3"
   PRINT *, "between -1 and 3 (result is 20)."

   CALL imc%initialize(1)   !intialize a MCI object for a 1 dimensional integral
   CALL imc%setIRange((/-1.d0,3.d0/))  !integral domain is [-1,3]
   CALL imc%setObservable(x3)   !set the function to be integrated
   CALL imc%integrate(1000000_8,average,error)   !start the numerical computation of the integral
   
   PRINT *, 
   PRINT *, "Result with MCIntegrator:"
   PRINT *, "  ", average, " +- ", error

   CALL imc%terminate()  !free the memory stored by imc

END PROGRAM ex1



FUNCTION x3(x)
   REAL(KIND=8) :: x3
   REAL(KIND=8), INTENT(IN) :: x(:)
   x3=x(1)*x(1)*x(1)
END FUNCTION x3
