MODULE DGMeshClass
  USE DGElementClass
  IMPLICIT NONE

  TYPE NodePointers
     INTEGER :: eLeft,eRight
  END type NodePointers

  TYPE DGMesh
     INTEGER                                     :: K
     TYPE(NodalDGStorage)                        :: DG
     TYPE(DGElement)   ,ALLOCATABLE,DIMENSION(:) :: e
     TYPE(NodePointers),ALLOCATABLE,DIMENSION(:) :: p
     REAL(KIND=RP)                               :: lambdamax
  END type DGMesh

CONTAINS

  SUBROUTINE ConstructMesh1D(this,x_nodes,K,N,nEqn)
    IMPLICIT NONE
    TYPE(DGMesh)                ,INTENT(INOUT) :: this
    REAL(KIND=RP),DIMENSION(0:K),INTENT(IN)    :: x_nodes
    INTEGER                     ,INTENT(IN)    :: K,N,nEqn
    !local
    INTEGER :: i

    this%K = K
    ALLOCATE(this%e(K),this%p(K))
    CALL ConstructNodalDGStorage(this%DG,N)
    DO i = 1,K
       CALL ConstructDGElement(this%e(i),nEqn,x_nodes(i-1),x_nodes(i),N)
    END DO
    !Set the node pointers
    DO i = 1,K-1
       this%p(i)%eLeft  = i
       this%p(i)%eRight = i+1
    END DO
    this%p(K)%eLeft  = K
    this%p(K)%eRight = 1

  END SUBROUTINE ConstructMesh1D

  !/////////////////////////////////////////!

  SUBROUTINE GlobalTimeDerivative(this,t)
    IMPLICIT NONE
    TYPE(DGMesh) ,INTENT(INOUT) :: this
    REAL(KIND=RP),INTENT(IN)    :: t
    !Local
    REAL(KIND=RP),DIMENSION(this%e(1)%nEqn) :: QL_ext,QR_ext,F
    REAL(KIND=RP)                           :: x_temp
    INTEGER                                 :: i,j,idL,idR,nEqn,N

    N    = this%DG%N
    nEqn = this%e(1)%nEqn
    
    !boundaryconditions here if needed
    DO i=1,this%K
       this%e(i)%QL=this%e(i)%Q(0,:)
       this%e(i)%QR=this%e(i)%Q(N,:)
    END DO
    
    !Riemannproblem
    
    DO j = 1,this%K
       idL = this%p(j)%eLeft
       idR = this%p(j)%eRight
       CALL RiemannSolver(this%e(idL)%QR,this%e(idR)%QL,F,this%e(idR)%nEqn,N,max(this%e(idL)%lambdamax&
            ,this%e(idR)%lambdamax))
       this%e(idR)%FstarL = -F
       this%e(idL)%FstarR = F
    END DO

    !Time derivative on each element

    DO j = 1,this%K
       CALL LocalTimeDerivative(this%e(j),this%DG)
    END DO
    
  END SUBROUTINE GlobalTimeDerivative
  
  !////////////////////////////////////////////////!

  SUBROUTINE DestructDGMesh(this)
    IMPLICIT NONE
    TYPE(DGMesh),INTENT(INOUT) :: this

    INTEGER :: i

    DO i = 1,this%K
       CALL DestructElement(this%e(i))
    END DO
    CALL DestructNodalDGStorage(this%DG)
    DEALLOCATE(this%e,this%p)

  END SUBROUTINE DestructDGMesh
END MODULE DGMeshClass

    
