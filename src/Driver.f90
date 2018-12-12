PROGRAM Driver
  USE TimeIntegration
  IMPLICIT NONE
  REAL(KIND=RP) :: muVal=0.2_RP,gammaVal=1.4_RP,PrVal=0.7_RP,CFLVal=0.9_RP,t=0.0_RP,tend=1.0_RP
  REAL(KIND=RP) :: xmin=0.0_RP,xmax=1.0_RP, dt
  TYPE(DGMesh)  :: Simulation
  INTEGER       :: NQ,N,nEqn=3
  INTEGER       :: i,l,m=0,j=0, numSteps=3
  CHARACTER(len=3) :: numChar
  CHARACTER(len=26) :: fname ='../Plots/Movies/UXXX.curve'
  REAL(KIND=RP),DIMENSION(:,:),allocatable :: Solution
  REAL(KIND=RP),DIMENSION(:,:), allocatable :: Fehler
  REAL(KIND=RP),DIMENSION(:), allocatable :: EOC,pSolution
  print*,"-------------------------------------------------------------"
  print*,"------------THIS IS A DGCODE WRITTEN BY NILS DORNBUSCH-------"
  print*,"------It is inspired by 'FORTRAN NOTES' by ANDREW WINTERS----"
  print*,"-------------------------------------------------------------"
  print*,"Starting everything up..."
  CALL setPhysPar(muVal,gammaVal,PrVal)
  print*,"physical constants set successfully"
  Print*,"-------------------------------------------------------------"
  print*,"Performing ",numSteps," step(s)."
  ALLOCATE(Fehler(nEqn,numSteps))
  ALLOCATE(EOC(numSteps))
  DO i=1,numSteps
     print*,"Step ",i,"of ",numSteps
     N=3
     NQ=2*i
     ALLOCATE(Solution(0:N*NQ,nEqn))
     ALLOCATE(pSolution(0:N*NQ))
     CALL ConstructSimulation(Simulation,NQ,N,nEqn,xmin,xmax,CFLVal)
     print*,"Finished setting up. Starting the time loop."
     
     DO WHILE(tend-t>epsilon(t))
        CALL getLambdaMaxGlobally(Simulation)
        print*,"t:",t
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
     print*,"Step ",i," finished successfully."
     DO l=1,Simulation%K
        Solution((l-1)*N:l*N,:)=Simulation%e(l)%Q(:,:)
        pSolution((l-1)*N:l*N)=points(l,:)
     END DO
     Fehler(1,i)=maxval(abs(Solution(:,1)-(2+0.5*sin(2*pi*(pSolution-tend)))))
     Fehler(2,i)=maxval(abs(Solution(:,2)-(2+0.5*sin(2*pi*(pSolution-tend)))))
     Fehler(3,i)=maxval(abs(Solution(:,3)-(1+0.25*sin(2*pi*(pSolution-tend))+1/(gammaVal-1))))
     DEALLOCATE(Solution,pSolution)
     CALL DestructDGMesh(Simulation)
     DEALLOCATE(points)
     t=0
     print*,"-----------------------------------------------------"
  END DO
  print*,"Fehler:",maxval(Fehler,dim=1)
  EOC=0
  DO i=2,numSteps
     EOC(i)=log(maxval(Fehler(:,i))/maxval(Fehler(:,i-1)))/log((i-1.0_RP)/i)
  END DO
  print*, "EOC:", EOC
  
END PROGRAM Driver
   
