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

PROGRAM ex3
   USE func
   USE mcintegrator
   IMPLICIT NONE
   INTERFACE
      FUNCTION x3(x)
         REAL(KIND=8) :: x3
         REAL(KIND=8), INTENT(IN) :: x(:)
      END FUNCTION x3
   END INTERFACE
   TYPE(MCI) :: imc
   INTEGER :: i
   REAL(KIND=8) :: target_acc_rate, time0, time1, acc_rate, mrt2_step(1:1)
   REAL(KIND=8) :: average(1:9), error(1:9), eff_0
   
   PRINT*, "We want to integrate"
   PRINT*, "   f(x) = exp(-x^2)"
   PRINT*, "between -Infinity and +Infinity (result is sqrt(pi)=1.772453850905516)"
   PRINT*, "In this example we want to show how the target acceptance rate can affect the result."
   
   CALL imc%initialize(1)   !intialize a MCI object for a 1 dimensional integral
   CALL imc%setObservable(gaussian)   !integrate the gaussian function
   CALL imc%setSamplingFunction(decaying_exponential)   !set the sampling function
   PRINT*, 
   
   target_acc_rate=-0.05
   DO i = 1, 9, 1
      target_acc_rate=target_acc_rate+0.1d0
      CALL imc%setTargetAcceptanceRate(target_acc_rate)   !set the acceptance rate
      CALL CPU_TIME(time0)
      CALL imc%integrate(10000000_8,average(i),error(i))   !start the numerical computation of the integral
      CALL CPU_TIME(time1)
      acc_rate=imc%getAcceptanceRate()
      mrt2_step=imc%getMRT2Step()
      PRINT*, "Acceptance Rate = ", acc_rate
      PRINT*, "MRT2 Step       = ", mrt2_step
      PRINT*, "                     ", average(i), " +- ", error(i)
      IF (i==1) eff_0=1.d0/(error(i)*error(i)*(time1-time0))
      PRINT*, "Efficiency =      ", (1.d0/(error(i)*error(i)*(time1-time0)))/eff_0
      PRINT*, 
   END DO

   CALL imc%terminate()
   
END PROGRAM ex3


FUNCTION x3(x)
   REAL(KIND=8) :: x3
   REAL(KIND=8), INTENT(IN) :: x(:)
   x3=x(1)*x(1)*x(1)
END FUNCTION x3
