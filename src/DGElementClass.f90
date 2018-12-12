MODULE DGElementClass
  USE NodalDGStorageClass
  USE FluxRoutines
  IMPLICIT NONE

  TYPE DGElement
     REAL(KIND=RP)                             :: delta_x,xL,xR,lambdamax
     INTEGER                                   :: nEqn
     REAL(KIND=RP),ALLOCATABLE, DIMENSION(:)   :: QR,QL,FstarR,FstarL
     REAL(KIND=RP),ALLOCATABLE, DIMENSION(:,:) :: Q,Q_dot
  END type DGElement

CONTAINS

  SUBROUTINE ConstructDGElement(this,nEqn,xL,xR,N)
    IMPLICIT NONE
    INTEGER        ,INTENT(IN)    :: N
    TYPE(DGElement),INTENT(INOUT) :: this
    REAL(KIND=RP)  ,INTENT(IN)    :: xL,xR
    INTEGER        ,INTENT(IN)    :: nEqn

    ALLOCATE(this%QR(nEqn),this%QL(nEqn),this%FstarR(nEqn),this&
         &%FstarL(nEqn))
    ALLOCATE(this%Q(0:N,nEqn),this%Q_dot(0:N,nEqn))
    this%nEqn     = nEqn
    this%xL       = xL
    this%xR       = xR
    this%delta_x  = xR - xL
    this%QR       = 0.0_RP
    this%QL       = 0.0_RP
    this%FstarR   = 0.0_RP
    this%FstarL   = 0.0_RP
    this%Q        = 0.0_RP
    this%Q_dot    = 0.0_RP
    
  END SUBROUTINE ConstructDGElement

  !///////////////////////////////////////!

  SUBROUTINE LocalTimeDerivative(this,DG)
    IMPLICIT NONE
    TYPE(DGElement)      ,INTENT(INOUT) :: this
    TYPE(NodalDGStorage) ,INTENT(IN)    :: DG
    !
    INTEGER                                   :: j,N,nEqn
    REAL(KIND=RP),DIMENSION(0:DG%N,this%nEqn) :: F,F_prime
    !
    N    = DG%N
    nEqn = this%nEqn
    DO j = 0,N
       CALL EulerAnalyticFlux(this%Q,F,N,nEqn)
    END DO

    CALL SystemDGDerivative(this%FstarR,this%FstarL,F,F_prime,DG&
         &%D_hat,DG%weights,DG%l_at_one,DG%l_at_minus_one,nEqn,N)
    this%Q_dot = (-2.0_RP/this%delta_x)*F_prime

  END SUBROUTINE LocalTimeDerivative

  !//////////////////////////////////////////!

  SUBROUTINE SystemDGDerivative(FR,FL,F,Fprime,D_hat,weights,l_one&
       &,l_minus_one,nEqn,N)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(nEqn)     ,INTENT(IN)  :: FR,FL
    REAL(KIND=RP),DIMENSION(0:N,nEqn) ,INTENT(IN)  :: F
    REAL(KIND=RP),DIMENSION(0:N,nEqn) ,INTENT(OUT) :: Fprime
    REAL(KIND=RP),DIMENSION(0:N,0:N)  ,INTENT(IN)  :: D_hat
    REAL(KIND=RP),DIMENSION(0:N)      ,INTENT(IN)  :: weights,l_one&
         &,l_minus_one
    INTEGER                           ,INTENT(IN)  :: nEqn,N
    !local
    INTEGER :: i,j

    Fprime = 0.0_RP
    
    !calculate volume terms
    DO i = 1,nEqn
       CALL CalculateVolume(Fprime(:,i),F(:,i),D_hat,N)
    END DO

    !Surface terms

    DO j = 0,N
       DO i = 1,nEqn
          Fprime(j,i) = Fprime(j,i) + (FR(i)*l_one(j) + FL(i)&
               &*l_minus_one(j))/weights(j)
       END DO
    END DO

  END SUBROUTINE SystemDGDerivative

  !///////////////////////////////////////////!

  SUBROUTINE DestructElement(this)
    IMPLICIT NONE
    TYPE(DGElement),INTENT(INOUT) :: this

    this%nEqn     = 0
    this%xL       = 0.0_RP
    this%xR       = 0.0_RP
    this%delta_x  = 0.0_RP
    DEALLOCATE(this%QR,this%QL,this%FstarR,this%FstarL,this%Q,this&
         &%Q_dot)
  END SUBROUTINE DestructElement
END MODULE DGElementClass

          
    
    
  
