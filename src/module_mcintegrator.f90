MODULE mcintegrator
   USE estimators
	IMPLICIT NONE

   ! Initial MRT2 step, in case it is not set
   REAL(KIND=8), PRIVATE, PARAMETER :: INITIAL_STEP=0.1d0

   ! Abstract interface for the sampling function
   ABSTRACT INTERFACE
      FUNCTION mcfunction(x)
         IMPLICIT NONE
         REAL(KIND=8) :: mcfunction
         REAL(KIND=8), INTENT(IN) :: x(:)
      END FUNCTION mcfunction
   END INTERFACE

   TYPE MCI
      LOGICAL, PRIVATE :: flaginit=.FALSE.   !flag that tells if the class has been initialized or not
      INTEGER, PRIVATE :: ndim   !number of dimensions
      REAL(KIND=8), PRIVATE, ALLOCATABLE :: irange(:,:)   !integration range
      REAL(KIND=8), PRIVATE, ALLOCATABLE :: x(:)   !walker position
      REAL(KIND=8), PRIVATE, ALLOCATABLE :: mrt2step(:)    !random step
      REAL(KIND=8), PRIVATE :: targetaccrate=0.5d0   !target acceptance rate
      INTEGER(KIND=8), PRIVATE :: acc, rej   !accepted and rejected moves
      PROCEDURE(mcfunction), POINTER, NOPASS :: pdf   !the MC integration will be done sampling from this pdf
      LOGICAL, PRIVATE :: flagnopdf
      REAL(KIND=8), PRIVATE :: pdfx   !value of the pdf function in x
      PROCEDURE(mcfunction), POINTER, NOPASS :: obs   !the MC integration will be done accumulating this observable
      REAL(KIND=8), PRIVATE :: obsx   !value of the obs function in x
      INTEGER(KIND=8), PRIVATE :: ridx  !running index, which keeps track of the index of datax
      REAL(KIND=8), PRIVATE, ALLOCATABLE :: datax(:)  !array that will contain all the measured observable
      REAL(KIND=8), PRIVATE :: vol  !integrated volume 

   CONTAINS

      PROCEDURE :: initialize, terminate

      PROCEDURE :: setIRange, setX, setMRT2Step, setTargetAcceptanceRate
      PROCEDURE :: setSamplingFunction, setObservable

      PROCEDURE :: getNDim
      PROCEDURE :: getIRange, getX, getMRT2Step, getTargetAcceptanceRate

      PROCEDURE :: integrate
      PROCEDURE :: getAcceptanceRate
      PROCEDURE :: findMRT2Step
      PROCEDURE :: initialDecorrelation
      PROCEDURE, PRIVATE :: doStepMRT2
      PROCEDURE, PRIVATE :: applyPBC
      PROCEDURE, PRIVATE :: computeNewObservable, confirmOldObservable
      PROCEDURE, PRIVATE :: resetAccRejCounters
      PROCEDURE, PRIVATE :: sample

   END TYPE


CONTAINS

   SUBROUTINE integrate(this,nmc,average,error)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER(KIND=8), INTENT(IN) :: nmc
      REAL(KIND=8), INTENT(OUT) :: average, error
      INTEGER :: i1
      REAL(KIND=8) :: foo(1:2)

      ! Find the optimal MRT2 step
      CALL this%findMRT2Step()

      ! Initial decorrelation
      CALL this%initialDecorrelation()

      ! Allocation of the data array
      ALLOCATE(this%datax(1:nmc))

      ! Monte Carlo sampling
      CALL this%sample(nmc,.TRUE.)

      ! Estimation of the integral
      foo(1:2)=correlated_estimator(nmc,this%datax(1:nmc))
      IF (this%flagnopdf) foo=foo*this%vol
      average=foo(1) ; error=foo(2)

      DEALLOCATE(this%datax)

   END SUBROUTINE integrate


   SUBROUTINE initialDecorrelation(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER(KIND=8), PARAMETER :: MIN_NMC=100
      LOGICAL :: flag_loop
      REAL(KIND=8) :: newestimate(1:2), oldestimate(1:2)

      ALLOCATE(this%datax(1:MIN_NMC))

      CALL this%sample(MIN_NMC,.TRUE.)
      oldestimate=correlated_estimator(MIN_NMC,this%datax(1:MIN_NMC))
      flag_loop=.TRUE.
      DO WHILE(flag_loop)
         CALL this%sample(MIN_NMC,.TRUE.)
         newestimate=correlated_estimator(MIN_NMC,this%datax(1:MIN_NMC))
         IF (DABS(newestimate(1)-oldestimate(1))<=newestimate(2)+oldestimate(2)) THEN
            ! Same estimators, the markov chain is now sampling correctly
            flag_loop=.FALSE.
         ELSE
            ! Estimators too different, the markov chain is not ready yet
            flag_loop=.TRUE.
         END IF
         oldestimate=newestimate
      END DO

      DEALLOCATE(this%datax)
      
   END SUBROUTINE initialDecorrelation


   SUBROUTINE sample(this,npoints,flag_obs)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER(KIND=8), INTENT(IN) :: npoints !number of points that should be sampled
      LOGICAL, INTENT(IN) :: flag_obs !store the observables in datax or not?
      LOGICAL :: flag_acc
      INTEGER(KIND=8) :: i1

      this%ridx=1
      IF (this%flagnopdf) THEN
         this%pdfx=this%vol
      ELSE
         this%pdfx=this%pdf(this%x)   
      END IF
      CALL this%resetAccRejCounters()
      DO i1 = 1, npoints, 1
         CALL this%doStepMRT2(flag_acc)
         IF (flag_obs) THEN
            IF (flag_acc) THEN
               CALL this%computeNewObservable()
            ELSE
               CALL this%confirmOldObservable()
            END IF
         END IF
      END DO
      
   END SUBROUTINE sample


   SUBROUTINE findMRT2Step(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER(KIND=8), PARAMETER :: MIN_STAT=100
      INTEGER(KIND=8), PARAMETER :: MIN_CONS=4
      REAL(KIND=8), PARAMETER :: TOLERANCE=0.05d0
      INTEGER(KIND=8), PARAMETER :: MAX_NUM_ATTEMPTS=10
      INTEGER(KIND=8) :: i1, i2, counter
      REAL(KIND=8) :: rate
      
      CALL this%resetAccRejCounters()
      IF (this%flagnopdf) THEN
         this%pdfx=this%vol
      ELSE
         this%pdfx=this%pdf(this%x)   
      END IF
      i1=0
      counter=0
      DO WHILE (i1<MIN_CONS)
         counter=counter+1
         DO i2 = 1, MIN_STAT, 1
            CALL this%doStepMRT2()
         END DO
         CALL this%getAcceptanceRate(rate)
         IF (rate>this%targetaccrate+TOLERANCE) THEN
            ! acceptance too high, step too small
            i1=0
            this%mrt2step=this%mrt2step*MIN(2.d0,rate/this%targetaccrate)
         ELSE IF (rate<this%targetaccrate-TOLERANCE) THEN
            ! acceptance too low, step too large
            i1=0
            this%mrt2step=this%mrt2step*MAX(0.5d0,rate/this%targetaccrate)
         ELSE
            i1=i1+1
         END IF
         CALL this%resetAccRejCounters()
         ! The while loop is running since too long
         IF (counter>MAX_NUM_ATTEMPTS) i1=MIN_CONS
         ! The mrt2step is +Infinity
         IF (MAXVAL(this%mrt2step)>MAXVAL(this%irange(2,:)-this%irange(1,:))) THEN
            this%mrt2step=MAXVAL(this%irange(2,:)-this%irange(1,:))
            i1=MIN_CONS
         END IF
         ! The mrt2step is 0.
         IF (MINVAL(this%mrt2step)<TINY(0.d0)) THEN
            this%mrt2step=TINY(0.e0)
            i1=MIN_CONS
         END IF
      END DO
      
   END SUBROUTINE findMRT2Step


   SUBROUTINE doStepMRT2(this,flag_acc)
      IMPLICIT NONE
      CLASS(MCI) :: this
      LOGICAL, INTENT(INOUT), OPTIONAL :: flag_acc
      REAL(KIND=8) :: propx(1:this%ndim), eta(1:this%ndim), proppdf
      
      CALL RANDOM_NUMBER(eta)
      eta=eta-0.5d0
      propx=this%x+this%mrt2step*eta
      CALL this%applyPBC(propx)
      IF (this%flagnopdf) THEN
         proppdf=this%vol
      ELSE
         proppdf=this%pdf(propx)
      END IF
      CALL RANDOM_NUMBER(eta(1))
      IF (proppdf/this%pdfx>eta(1)) THEN
         !accepted
         this%acc=this%acc+1
         this%x=propx
         this%pdfx=proppdf
         IF (PRESENT(flag_acc)) flag_acc=.TRUE.
      ELSE
         !rejected
         this%rej=this%rej+1
         IF (PRESENT(flag_acc)) flag_acc=.FALSE.
      END IF
      
   END SUBROUTINE doStepMRT2
   

   SUBROUTINE applyPBC(this,v)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), OPTIONAL, INTENT(INOUT) :: v(1:this%ndim)
      INTEGER :: i1

      IF (PRESENT(v)) THEN
         DO i1 = 1, this%ndim, 1
            DO WHILE (v(i1)<this%irange(1,i1))
               v(i1)=v(i1)+(this%irange(2,i1)-this%irange(1,i1))
            END DO
            DO WHILE (v(i1)>this%irange(2,i1))
               v(i1)=v(i1)-(this%irange(2,i1)-this%irange(1,i1))
               END DO
         END DO
      ELSE
         DO i1 = 1, this%ndim, 1
            DO WHILE (this%x(i1)<this%irange(1,i1))
               this%x(i1)=this%x(i1)+(this%irange(2,i1)-this%irange(1,i1))
            END DO
            DO WHILE (this%x(i1)>this%irange(2,i1))
               this%x(i1)=this%x(i1)-(this%irange(2,i1)-this%irange(1,i1))
               END DO
         END DO
      END IF
      
   END SUBROUTINE applyPBC


   SUBROUTINE getAcceptanceRate(this,rate)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), INTENT(OUT) :: rate

      rate=REAL(this%acc,8)/REAL(this%acc+this%rej,8)
      
   END SUBROUTINE getAcceptanceRate


   SUBROUTINE resetAccRejCounters(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      
      this%acc=0
      this%rej=0

   END SUBROUTINE resetAccRejCounters


   SUBROUTINE computeNewObservable(this)
      IMPLICIT NONE
      CLASS(MCI) :: this

      IF (this%flagnopdf) THEN
         this%obsx=this%obs(this%x)
      ELSE
         this%obsx=this%obs(this%x)/this%pdfx
      END IF
      this%datax(this%ridx)=this%obsx
      this%ridx=this%ridx+1
      
   END SUBROUTINE computeNewObservable


   SUBROUTINE confirmOldObservable(this)
      IMPLICIT NONE
      CLASS(MCI) :: this

      this%datax(this%ridx)=this%obsx
      this%ridx=this%ridx+1
      
   END SUBROUTINE confirmOldObservable


   !!!!!!!!  SET PROPERTIES OF THE MC INTEGRATION  !!!!!!!!


   SUBROUTINE setObservable(this,observable)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTERFACE
         FUNCTION observable(x)
            IMPLICIT NONE
            REAL(KIND=8) :: observable
            REAL(KIND=8), INTENT(IN) :: x(:)
         END FUNCTION observable
      END INTERFACE

      this%obs=>observable
      
   END SUBROUTINE setObservable


   SUBROUTINE setSamplingFunction(this,sampling_function)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTERFACE
         FUNCTION sampling_function(x)
            IMPLICIT NONE
            REAL(KIND=8) :: sampling_function
            REAL(KIND=8), INTENT(IN) :: x(:)
         END FUNCTION sampling_function
      END INTERFACE

      this%pdf=>sampling_function
      this%flagnopdf=.FALSE.
      
   END SUBROUTINE setSamplingFunction


   SUBROUTINE setIRange(this,irange)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), INTENT(IN) :: irange(1:2,1:this%ndim)
      INTEGER :: i1

      this%irange=irange
      DO i1 = 1, this%ndim, 1
         this%x(i1)=(irange(2,i1)+irange(1,i1))*0.5d0
      END DO

      this%vol=1.d0
      DO i1 = 1, this%ndim, 1
         this%vol=this%vol*(this%irange(2,i1)-this%irange(1,i1))
      END DO
      
   END SUBROUTINE setIRange


   SUBROUTINE setX(this,x0)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), INTENT(IN) :: x0(1:this%ndim)

      this%x=x0
      CALL this%applyPBC()
      
   END SUBROUTINE setX


   SUBROUTINE setMRT2Step(this,step)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), INTENT(IN) :: step(1:this%ndim)

      this%mrt2step=step
      
   END SUBROUTINE setMRT2Step


   SUBROUTINE setTargetAcceptanceRate(this,targetaccrate)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8), INTENT(IN) :: targetaccrate

      this%targetaccrate=targetaccrate
      
   END SUBROUTINE setTargetAcceptanceRate


   FUNCTION getTargetAcceptanceRate(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8) :: getTargetAcceptanceRate

      getTargetAcceptanceRate=this%targetaccrate
   
   END FUNCTION getTargetAcceptanceRate


   FUNCTION getMRT2Step(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8) :: getMRT2Step(1:this%ndim)

      getMRT2Step=this%mrt2step
   
   END FUNCTION getMRT2Step


   FUNCTION getX(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8) :: getX(1:this%ndim)

      getX=this%x
   
   END FUNCTION getX


   FUNCTION getNDim(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER :: getNDim
   
      getNDim=this%ndim
   
   END FUNCTION getNDim


   FUNCTION getIRange(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      REAL(KIND=8) :: getIRange(1:2,1:this%ndim)

      getIRange=this%irange
   
   END FUNCTION getIRange


   SUBROUTINE initialize(this,ndim)
      IMPLICIT NONE
      CLASS(MCI) :: this
      INTEGER, INTENT(IN) :: ndim
      INTEGER :: i1

      IF (this%flaginit) STOP "MCI error: Cannot initialize already initialized object."

      ! ndim
      this%ndim=ndim
      ! irange
      ALLOCATE(this%irange(1:2,1:this%ndim))
      DO i1 = 1, this%ndim, 1
         this%irange(1,i1)=-HUGE(0.d0) ; this%irange(2,i1)=HUGE(0.d0)
      END DO
      ! vol
      this%vol=1.d0
      DO i1 = 1, this%ndim, 1
         this%vol=this%vol*(this%irange(2,i1)-this%irange(1,i1))
      END DO
      ! x
      ALLOCATE(this%x(1:this%ndim))
      this%x=0.d0
      ! mrt2step
      ALLOCATE(this%mrt2step(1:this%ndim))
      this%mrt2step=INITIAL_STEP
      ! pdf
      this%flagnopdf=.TRUE.

      this%flaginit=.TRUE.
      
   END SUBROUTINE initialize


   SUBROUTINE terminate(this)
      IMPLICIT NONE
      CLASS(MCI) :: this
      
      IF (.NOT. this%flaginit) STOP "MCI error: Cannot terminate a not initialized object."
      DEALLOCATE(this%irange)
      DEALLOCATE(this%x)
      DEALLOCATE(this%mrt2step)
      this%flaginit=.FALSE.
      IF (ALLOCATED(this%datax)) DEALLOCATE(this%datax)
      
   END SUBROUTINE terminate

	
END MODULE mcintegrator


