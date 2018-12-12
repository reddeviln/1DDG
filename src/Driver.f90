PROGRAM Driver
  USE TimeIntegration
  IMPLICIT NONE
  REAL(KIND=RP) :: muVal=0.2_RP,gammaVal=1.4_RP,PrVal=0.7_RP,CFLVal=0.9_RP,t=0.0_RP,tend=1.0_RP
  REAL(KIND=RP) :: xmin=0.0_RP,xmax=1.0_RP, dt
  TYPE(DGMesh)  :: Simulation
  INTEGER,PARAMETER       :: NQ=8,N=3,nEqn=3
  INTEGER       :: i,m=0,j=0
  CHARACTER(len=3) :: numChar
  CHARACTER(len=26) :: fname ='../Plots/Movies/UXXX.curve'
  REAL(KIND=RP),DIMENSION(1:nq**3,0:N,0:N,0:N) :: Qplot
  print*,"THIS IS A DGCODE BY NILS DORNBUSCH inspired by 'FORTRAN NOTES'&
       &"," by ANDREW WINTERS"
  print*,"Starting everything up..."
  CALL setPhysPar(muVal,gammaVal,PrVal)
  print*,"physical constants set successfully"
  CALL ConstructSimulation(Simulation,NQ,N,nEqn,xmin,xmax,CFLVal)
  print*,"Finished setting up. Starting the time loop."

  DO WHILE(tend-t>epsilon(t))
     CALL getLambdaMaxGlobally(Simulation)
     print*,"t:",t
     print*, "lambda:",Simulation%lambdamax
     print*,"sum Energy 1st Element:", sum(Simulation%e(1)%Q(:,3))
     dt=CFL/Simulation%lambdamax*((xmax-xmin)/real(NQ,kind=RP))/real(N+1,kind=RP)
     print*,"dt:",dt
     IF(tend<(t+dt)) THEN
        dt=tend-t
        print*,"dt adjusted to", dt
     END IF
     CALL getRungeKuttaStep(Simulation,t,dt)
     t=t+dt
     IF(MODULO(j,5).EQ.0) THEN
        m=m+1
        WRITE(numChar,'(i3)')m
        IF (m.GE.100) THEN
           fName(18:20) = numChar
        ELSEif(m.GE.10) then
           fName(18:18)    = "0"
           fName(19:20)  = numChar(2:3)
        ELSE
           fName(18:19) = "00"
           fName(20:20)=numChar(3:3)
        END IF
        CALL preparePlot(Simulation,fName,1)
     END IF
     j=j+1
  END DO
  CALL preparePlot(Simulation,"../Plots/endTime.curve",1)
print*,"Done without a problem! See you next time."
  
END PROGRAM Driver

