! Set of functions for estimating the average and error of a set of
! data, correlated or uncorrelated.
! If data are uncorrelated the clocking technique with the plateau
! criterium is used.
! All functions return a 2-dimensional array, the first value is the
! average and the second the standard deviation


MODULE estimators
   IMPLICIT NONE

   !LOGICAL, PRIVATE :: debug_mode=.FALSE.
   

CONTAINS

   FUNCTION correlated_estimator(N,x)
      IMPLICIT NONE
      REAL(KIND=8) :: correlated_estimator(1:2)
      INTEGER, PARAMETER :: MIN_BLOCKS=6, MAX_BLOCKS=50
      INTEGER, PARAMETER :: MAX_PLATEAU_AVERAGE=4
      INTEGER(KIND=8), INTENT(IN) :: N
      REAL(KIND=8), INTENT(IN) :: x(1:N)
      INTEGER :: i1, i2, i_min
      REAL(KIND=8) :: av(MIN_BLOCKS:MAX_BLOCKS), err(MIN_BLOCKS:MAX_BLOCKS)
      REAL(KIND=8) :: foo(1:2), delta, accdelta(MIN_BLOCKS+MAX_PLATEAU_AVERAGE:MAX_BLOCKS-MAX_PLATEAU_AVERAGE)

      !IF (debug_mode) THEN
      !   PRINT *, "correlated_estimator:"
      !   PRINT *, "#blocks     #estimate           #error"
      !END IF
      DO i1 = MIN_BLOCKS, MAX_BLOCKS, 1
         foo=block_estimator(N,x,i1)
         av(i1)=foo(1) ; err(i1)=foo(2)
         !IF (debug_mode) THEN
         !   PRINT *, i1, av(i1), err(i1)
         !END IF
      END DO

      delta=0.d0
      DO i1 = 1, MAX_PLATEAU_AVERAGE, 1
      DO i2 = MIN_BLOCKS+MAX_PLATEAU_AVERAGE, MAX_BLOCKS-MAX_PLATEAU_AVERAGE, 1
         SELECT CASE(i1)
         CASE(1)
            delta=(-0.5d0*err(i2-1)+0.5d0*err(i2+1))
         CASE(2)
            delta=((1.d0/12.d0)*err(i2-2)-(2.d0/3.d0)*err(i2-1)+ &
                  (2.d0/3.d0)*err(i2+1)-(1.d0/12.d0)*err(i2+2))
         CASE(3)
            delta=(-(1.d0/60.d0)*err(i2-3)+(3.d0/20.d0)*err(i2-2)-0.75d0*err(i2-1)+ &
                  0.75d0*err(i2+1)-(3.d0/20.d0)*err(i2+2)+(1.d0/60.d0)*err(i2+3))
         CASE(4)
            delta=((1.d0/280.d0)*err(i2-4)-(4.d0/105.d0)*err(i2-3)+0.2d0*err(i2-2)- &
                  0.8d0*err(i2-1)+0.8d0*err(i2+1)-0.2d0*err(i2+2)+(4.d0/105.d0)*err(i2+3)-(1.d0/280.d0)*err(i2+4))
         END SELECT
         accdelta(i2)=accdelta(i2)+delta
      END DO
      END DO
      
      i_min=MIN_BLOCKS+MAX_PLATEAU_AVERAGE
      DO i2 = MIN_BLOCKS+MAX_PLATEAU_AVERAGE, MAX_BLOCKS-MAX_PLATEAU_AVERAGE, 1
         IF (DABS(accdelta(i2))<DABS(accdelta(i_min))) i_min=i2
      END DO

      correlated_estimator(1)=0.2d0*(av(i_min-2)+av(i_min-1)+av(i_min)+av(i_min+1)+av(i_min+2))
      correlated_estimator(2)=0.2d0*(err(i_min-2)+err(i_min-1)+err(i_min)+err(i_min+1)+err(i_min+2))
      
   END FUNCTION correlated_estimator


   FUNCTION block_estimator(N,x,nblocks)
      IMPLICIT NONE
      REAL(KIND=8) :: block_estimator(1:2)
      INTEGER(KIND=8), INTENT(IN) :: N
      REAL(KIND=8), INTENT(IN) :: x(1:N)
      INTEGER, INTENT(IN) :: nblocks !number of blocks
      INTEGER :: i1
      REAL(KIND=8) :: av(1:nblocks), norm

      !IF (debug_mode) THEN
      !   PRINT *, "block_estimator: number of data elements = ", N
      !   PRINT *, "block_estimator: number of blocks = ", nblocks
      !   PRINT *, "block_estimator: number of data for block = ", N/nblocks
      !END IF
      norm=1.d0/REAL(N/nblocks,8)
      DO i1 = 1, nblocks, 1
         av(i1)=SUM(x((i1-1)*(N/nblocks)+1:i1*(N/nblocks)))*norm
         !IF (debug_mode) THEN
         !   PRINT *, "block_estimator: compute block average from ", (i1-1)*N/nblocks+1, &
         !            " to ", i1*N/nblocks
         !END IF
      END DO
      !IF (debug_mode)  THEN
      !   PRINT *, "block_estimator: block averages"
      !   PRINT *, av
      !END IF
      block_estimator=uncorrelated_estimator(INT(nblocks,8),av)
   
   END FUNCTION block_estimator


   FUNCTION uncorrelated_estimator(N,x)
      IMPLICIT NONE
      REAL(KIND=8) :: uncorrelated_estimator(1:2)
      INTEGER(KIND=8), INTENT(IN) :: N
      REAL(KIND=8), INTENT(IN) :: x(1:N)
      REAL(KIND=8) :: av2, norm

      IF (N<2) STOP "MPI error uncorrelated_estimator() : N must be larger than 1"

      !IF (debug_mode) THEN
      !   PRINT *, "uncorrelated_estimator: number of data elements = ", N
      !END IF
      norm=1.d0/REAL(N,8)
      uncorrelated_estimator(1)=SUM(x)*norm
      av2=DOT_PRODUCT(x,x)*norm
      uncorrelated_estimator(2)=DSQRT((av2-uncorrelated_estimator(1)*uncorrelated_estimator(1))/REAL(N-1,8))

   END FUNCTION uncorrelated_estimator

END MODULE estimators
