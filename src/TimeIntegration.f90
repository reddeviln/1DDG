MODULE Timeintegration
  USE DGMeshClass
  IMPLICIT NONE
  REAL(KIND=RP)                             :: CFL,lambdamax
  REAL(KIND=RP), DIMENSION(:,:) , ALLOCATABLE :: points
  REAL(KIND=RP), DIMENSION(5)               ::  a=(/0.0_rp, -567301805773.0_rp/1357537059087.0_rp,&
       &-2404267990393.0_rp/2016746695238.0_rp, -3550918686646.0_rp/2091501179385.0_rp,&
       &-1275806237668.0_rp/842570457699.0_rp/), b=(/0.0_rp, 1432997174477.0_rp/9575080441755.0_rp,&
       &2526269341429.0_rp/6820363962896.0_rp, 2006345519317.0_rp/3224310063776.0_rp,&
       &2802321613138.0_rp/2924317926251.0_rp /), c=(/1432997174477.0_rp/9575080441755.0_rp, &
       &5161836677717.0_rp/13612068292357.0_rp,1720146321549.0_rp/2090206949498.0_rp, 3134564353537.0_rp/4481467310338.0_rp,&
       &2277821191437.0_rp/14882151754819.0_rp /)
CONTAINS

  SUBROUTINE ConstructSimulation(this, NQ, N, nEqn, xmin, xmax, CFLval)
    IMPLICIT NONE
    TYPE(DGMesh) , INTENT(INOUT) :: this
    INTEGER      , INTENT(IN)    :: NQ, N, nEqn
    REAL(KIND=RP), INTENT(IN)    :: xmin, xmax, CFLval

    !local variables

    REAL(KIND=RP)                   :: dx
    INTEGER                         :: i, K
    REAL(KIND=RP) , DIMENSION(0:NQ) :: nodes

    K = NQ
    dx=(xmax-xmin)/real(NQ,kind=RP)
    DO i=0,NQ
       nodes(i) = i*dx+xmin
    END DO
    CFL = CFLval
    write(*, fmt="(1x,a)", advance="no") "Constructing Mesh ..."
    CALL ConstructMesh1D(this,nodes,NQ,N,nEqn)
    print*, "Done"
    write(*, fmt="(1x,a)", advance="no")"Imposing Initial Condition ..."
    CALL InitialCondition(this)
    print*, "Done"
    
  END SUBROUTINE ConstructSimulation

  !//////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE InitialCondition(this)
    IMPLICIT NONE
    TYPE(DGMesh) , INTENT(INOUT) :: this

    !local variables

    REAL(KIND=RP) , DIMENSION(0:this%DG%N)              :: x,w
    REAL(KIND=RP) , DIMENSION(0:this%DG%N*this%K)       :: Qplot, pplot
    INTEGER                                             :: i
    ALLOCATE(points(1:this%K, 0:this%DG%N))
    Qplot=-1.0_RP
    pplot=0
    CALL LegendreGaussLobattoNodesandWeights(this%DG%N,x,w)
    DO i=1,this%K
       points(i,:)=(this%e(i)%xR-this%e(i)%xL)*0.5_RP*x+(this%e(i)%xR&
            &+this%e(i)%xL)*0.5_RP
       this%e(i)%Q(:,1)=2.0_RP+0.5_RP*sin(2.0_RP*pi*points(i,:))
       this%e(i)%Q(:,2)=this%e(i)%Q(:,1)
       this%e(i)%Q(:,3)=this%e(i)%Q(:,1)*(1+1/((gamma-1.0_RP)*this&
            &%e(i)%Q(:,1)))
       Qplot((i-1)*this%DG%N:i*this%DG%N)=this%e(i)%Q(:,1)
       pplot((i-1)*this%DG%N:i*this%DG%N)=points(i,:)
    END DO
    OPEN(file='../Plots/initial.dat', unit=15)
    CALL ExportToTecplot_1D(pplot,qplot,this%DG%N,15,"rho")
    CLOSE(15)
  END SUBROUTINE InitialCondition
    
  
  
END MODULE Timeintegration


