PROGRAM Driver
  USE TimeIntegration
  IMPLICIT NONE
  REAL(KIND=RP) :: muVal=0.2_RP,gammaVal=1.4_RP,PrVal=0.7_RP,CFLVal=0.3_RP,t=0.0_RP,tend=1.0_RP
  REAL(KIND=RP) :: xmin=-1.0_RP,xmax=1.0_RP,ymin=-1.0_RP,ymax=1.0_RP&
       &,zmin=-1.0_RP,zmax=1.0_RP
  TYPE(DGMesh)  :: Simulation
  INTEGER,PARAMETER       :: NQ=8,N=3,nEqn=5
  INTEGER       :: i,m=0,j=0
  CHARACTER(len=3) :: numChar
  CHARACTER(len=24) :: fname ='../Plots/Movies/UXXX.tec'
  REAL(KIND=RP),DIMENSION(1:nq**3,0:N,0:N,0:N) :: Qplot
  print*,"THIS IS A DGCODE BY NILS DORNBUSCH inspired by 'FORTRAN NOTES'&
       &"," by ANDREW WINTERS"
  print*,"Starting everything up..."
  CALL setPhysPar(muVal,gammaVal,PrVal)
  print*,"physical constants set successfully"
  CALL ConstructSimulation(Simulation,NQ,N,nEqn,xmin,xmax,CFLVal)
END PROGRAM Driver

