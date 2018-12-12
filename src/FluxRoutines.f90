MODULE FluxRoutines
  USE Diverses
  IMPLICIT NONE
  !Contains many different flux functions and Riemann solvers
  REAL(KIND=RP)         :: mu,g=9.812_RP,R=287.058_RP,gamma,Pr

CONTAINS

  SUBROUTINE setPhysPar(muVal,gammaVal,PrVal)
    IMPLICIT NONE
    REAL(KIND=RP),INTENT(IN) ::muVal,gammaVal,PrVal

    mu=muVal
    gamma=gammaVal
    Pr=PrVal
  END SUBROUTINE setPhysPar

  !//////////////////////////////////////////////////////
  
  SUBROUTINE calculateVolume(Fout,Fin,Dhat,N)
    IMPLICIT NONE

    REAL(KIND=RP), DIMENSION(0:N)    , INTENT(OUT)    :: Fout
    REAL(KIND=RP), DIMENSION(0:N)    , INTENT(IN)     :: Fin
    REAL(KIND=RP), DIMENSION(0:N,0:N), INTENT(IN)     :: Dhat
    INTEGER                          , INTENT(IN)     :: N
    
    !local variables
    INTEGER :: i

    DO i=0,N
       Fout(i)=dot_product(Dhat(i,:),Fin)
    END DO
    
    
  END SUBROUTINE calculateVolume
  
 !////////////////////////////////////////////////////////////////////

  SUBROUTINE RiemannSolver(QR,QL,Fout,nEqn,N,lambdamax)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(nEqn),INTENT(IN)  :: QR,QL
    REAL(KIND=RP),DIMENSION(nEqn),INTENT(OUT) :: Fout
    INTEGER                      ,INTENT(IN)  :: nEqn,N
    REAL(KIND=RP)                ,INTENT(IN)  :: lambdamax

    !local variables

    REAL(KIND=RP),DIMENSION(nEqn)   :: Fsharp,FL,FR
    INTEGER                         :: i,j

    CALL EulerAnalyticFluxPoint(QR,FL,N,nEqn)
    CALL EulerAnalyticFluxPoint(QL,FR,N,nEqn)
    Fsharp=0.5_RP*(FL+FR)
    Fout=Fsharp-0.5_RP*lambdamax*(QL-QR)
    
  END SUBROUTINE RiemannSolver

  !////////////////////////////////////////////////////////////////////

  SUBROUTINE EulerAnalyticFlux(Q,Fout,N,nEqn)
    IMPLICIT NONE
    INTEGER                          ,INTENT(IN) :: N, nEqn
    REAL(KIND=RP),DIMENSION(0:N,nEqn),INTENT(IN) :: Q
    REAL(KIND=RP),DIMENSION(0:N,nEqn),INTENT(OUT):: Fout

    !local variables
    REAL(KIND=RP),DIMENSION(0:n) :: p
    p=(gamma-1.0_RP)*(Q(:,3)-0.5_RP*(Q(:,2)**2)/Q(:,1))
    
    Fout(:,1)=Q(:,2)
    Fout(:,2)=Q(:,2)**2/Q(:,1)+p
    Fout(:,3)=Q(:,2)/Q(:,1)*(Q(:,3)+p)
    
  END SUBROUTINE EulerAnalyticFlux

  SUBROUTINE EulerAnalyticFluxPoint(Q,Fout,N,nEqn)
    IMPLICIT NONE
    INTEGER                      ,INTENT(IN) :: N, nEqn
    REAL(KIND=RP),DIMENSION(nEqn),INTENT(IN) :: Q
    REAL(KIND=RP),DIMENSION(nEqn),INTENT(OUT):: Fout

    !local variables
    REAL(KIND=RP) :: p
    p=(gamma-1.0_RP)*(Q(3)-0.5_RP*(Q(2)**2)/Q(1))
    
    Fout(1)=Q(2)
    Fout(2)=Q(2)**2/Q(1)+p
    Fout(3)=Q(2)/Q(1)*(Q(3)+p)
    
  END SUBROUTINE EulerAnalyticFluxPoint
  
END MODULE FluxRoutines
