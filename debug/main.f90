MODULE tools
   IMPLICIT NONE
   

CONTAINS

   FUNCTION one(x)
      IMPLICIT NONE
      REAL(KIND=8) :: one
      REAL(KIND=8), INTENT(IN) :: x(:)
   
      one=1.d0
   
   END FUNCTION one

   FUNCTION x2(x)
      IMPLICIT NONE
      REAL(KIND=8) :: x2
      REAL(KIND=8), INTENT(IN) :: x(:)
   
      x2=DOT_PRODUCT(x,x)
   
   END FUNCTION x2

END MODULE tools



PROGRAM test
   USE estimators
	USE mcintegrator
   USE pdf
   USE tools
	IMPLICIT NONE
   INTEGER :: i1
   INTEGER(KIND=8), PARAMETER :: ndata1=10000
   REAL(KIND=8) :: data1(1:ndata1)
   REAL(KIND=8) :: estimator(1:2)
   INTEGER(KIND=4), PARAMETER :: nblocks1=13
   TYPE(MCI) :: mcintegral
   INTEGER, PARAMETER :: ndim1=1
   REAL(KIND=8) :: irange1(1:2)
   REAL(KIND=8) :: x1(1:ndim1)
   REAL(KIND=8) :: mrt2step1(1:ndim1)
   REAL(KIND=8) :: accrate1
   INTEGER(KIND=8) :: nmc1=1000000
   REAL(KIND=8) :: est1(1:2)

   PRINT *, "Test protocol"
   PRINT *,
   
   
   PRINT *, " > > > TEST ESTIMATORS"
   PRINT *, 
   PRINT *, "All data are equal to 2.5"
   data1=2.5d0
   estimator=uncorrelated_estimator(ndata1,data1)
   PRINT *, "- uncorrelated_estimator() -> ", estimator
   estimator=block_estimator(ndata1,data1,nblocks1)
   PRINT *, "- block_estimator(",INT(nblocks1,2),") -> ", estimator
   estimator=correlated_estimator(ndata1,data1)
   PRINT *, "- correlated_estimator() -> ", estimator
   PRINT *, 
   
   PRINT *, "Data are sampled randomly from 0.5 to 2."
   CALL RANDOM_NUMBER(data1)
   data1=data1*1.5d0+0.5d0
   estimator=uncorrelated_estimator(ndata1,data1)
   PRINT *, "- uncorrelated_estimator() -> ", estimator
   estimator=block_estimator(ndata1,data1,nblocks1)
   PRINT *, "- block_estimator(",INT(nblocks1,2),") -> ", estimator
   estimator=correlated_estimator(ndata1,data1)
   PRINT *, "- correlated_estimator() -> ", estimator
   PRINT *, 

   PRINT *, 
   PRINT *, " > > > TEST MC_INTEGRAL"
   PRINT *, 
   CALL mcintegral%initialize(ndim1)
   PRINT *, "- mc_integral%initialize(",INT(ndim1,2),") -> number of dimensions = ", mcintegral%getNDim()
   PRINT *, "- integral range = ", mcintegral%getIRange()
   PRINT *, "- x = ", mcintegral%getX()
   PRINT *, "- stepmrt2 = ", mcintegral%getMRT2step()
   PRINT *, "- acceptance rate = ", mcintegral%getTargetAcceptanceRate()

   PRINT *, 
   irange1=(/-3.d0,3.d0/)
   CALL mcintegral%setIRange(irange1)
   PRINT *, "- mc_integral%setIRange() -> integral range = ", mcintegral%getIRange()

   PRINT *, 
   PRINT *, "- mc_integral%setX() -> "
   x1=1*irange1(1)-2.d0
   DO i1 = 1, 10, 1
      x1=x1+1.d0; CALL mcintegral%setX(x1)
      PRINT *, x1, " -> ", mcintegral%getX()
   END DO

   PRINT *, 
   mrt2step1=DABS(irange1(2)-irange1(1))*0.25d0
   CALL mcintegral%setMRT2Step(mrt2step1)
   PRINT *, "- mc_integral%setMRT2Step() -> ", mcintegral%getMRT2step()

   PRINT *, 
   accrate1=0.65d0
   CALL mcintegral%setTargetAcceptanceRate(accrate1)
   PRINT *, "- mc_integral%setTargetAcceptanceRate() -> ", mcintegral%getTargetAcceptanceRate()

   !PRINT *, 
   !CALL PDF_setNDim(ndim1)
   !CALL PDF_setGaussianC(0.1d0)
   !CALL mcintegral%setSamplingFunction()
   !PRINT *, "- mc_integral%setSamplingFunction() -> "
   !x1=-3.d0
   !DO i1 = 1, 11, 1
   !   PRINT *, x1, " -> ", mcintegral%pdf(x1)
   !   x1=x1+0.6d0
   !END DO

   PRINT *, 
   CALL mcintegral%setObservable(x2)
   PRINT *, "- mc_integral%setObservable() -> "
   x1=-3.d0
   DO i1 = 1, 11, 1
      PRINT *, x1, " -> ", mcintegral%obs(x1)
      x1=x1+0.6d0
   END DO

   PRINT *, 
   CALL mcintegral%findMRT2Step()
   PRINT *, "- mc_integral%findMRT2Step() -> ", mcintegral%getMRT2step()

   PRINT *, 
   CALL mcintegral%initialDecorrelation()
   PRINT *, "- mc_integral%initalDecorrelation()"

   PRINT *, 
   CALL mcintegral%integrate(nmc1,est1(1),est1(2))
   PRINT *, "- mc_integral%integrate() -> ", est1(1), " +- ", est1(2)

   CALL mcintegral%terminate()


END PROGRAM test


