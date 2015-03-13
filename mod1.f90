!-------------------------------------------------------------------------------------
Subroutine Rand_Field(initial)
  Use global
  Implicit None

  Integer                :: initial   ! Defines the type of initial cond
  Integer                :: i,j,dummy ! dummy is a grbage index
  Integer                :: IOstatus1
  Real*8                 :: rand, dum_rad      ! used for random number

6002 format(I3,7E20.10,I10,I3)

  ! This is for reading the initial omega fields from a files.
  ! The M and N should match the size of the data file
  If (initial==1) then
     ! Data for bilayer
     If (bilayer==1) then
        If (iter==1) then 
           Open(4,file='./results/shapes/bilayer_M50_N50.dat',status='old')
           Do i=1,N
              Do j=1,M
                 Read(4,*) dummy,wA(i,j),wB(i,j),wC(i,j)
                 Call random_number(rand)
                 wC(i,j)=wC(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wA(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wB(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wE(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wD(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0
                 
              End Do
           End Do
           Close(4)
        ElseIf (iter>1) then
           Open(4,file='./results/omega.bilayer',status='old')
           Do i=1,N
              Do j=1,M
                 Read(4,*) dummy,wA(i,j),wB(i,j),wC(i,j),wE(i,j),wD(i,j)
                 Call random_number(rand)
                 wC(i,j)=wC(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wA(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wB(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wE(i,j)=wE(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wD(i,j)=wD(i,j)+(2.0*rand-1.0)*1.0d0
              End Do
           End Do
           Close(4)
        end If
     End If
     ! Data for disk
     If (disk==1) then
        If (iter==1) then
           Open(4,file='./results/shapes/disk_M50_N50.dat',status='old')
           Do i=1,N
              Do j=1,M
                 Read(4,*) dummy,wA(i,j),wB(i,j),wC(i,j)
                 Call random_number(rand)
                 wC(i,j)=wC(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wA(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wB(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wE(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wD(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0

              End Do
           End Do
           Close(4)
        ElseIf (iter>1) then
           Open(4,file='./results/omega.disk',status='old')
           Do i=1,N
              Do j=1,M
                 Read(4,*) dummy,wA(i,j),wB(i,j),wC(i,j),wE(i,j),wD(i,j)
                 Call random_number(rand)
                 wC(i,j)=wC(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wA(i,j)=wA(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wB(i,j)=wB(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wE(i,j)=wE(i,j)+(2.0*rand-1.0)*1.0d0
                 Call random_number(rand)
                 wD(i,j)=wD(i,j)+(2.0*rand-1.0)*1.0d0
              End Do
           End Do
           Close(4)
        end If
     end If
  ElseIf (initial==0) then
     ! This is for a random initial condition
     Do i=1,N
        Do j=1,M
           Call random_number(rand)
           wA(i,j)=0.5*rand
           Call random_number(rand)
           wB(i,j)=0.5*rand
           Call random_number(rand)
           wC(i,j)=0.5*rand
           Call random_number(rand)
           wE(i,j)=0.5*rand
           Call random_number(rand)
           wD(i,j)=0.5*rand
        End Do
     End Do
  ElseIf (initial==2) then
     ! This is to do a specific type of initial condition
     ! like a sin function of a line or what ever
     Do i=1,N
        Do j=1,18
            wC(i,j) = 300.0
            wA(i,j) = 0.0
            wB(i,j) = 0.0
            wD(i,j) = 0.0
            wE(i,j) = 0.0
        End do
        Do j=19,22
            wC(i,j) = 0.0
            wA(i,j) = 300.0
            wB(i,j) = 0.0
            wD(i,j) = 0.0
            wE(i,j) = 0.0
        End Do
        Do j=23,26
            wC(i,j) = 0.0
            wA(i,j) = 0.0
            wB(i,j) = 300.0
            wD(i,j) = 0.0
            wE(i,j) = 0.0
        End Do
        Do j=27,30
            wC(i,j) = 0.0
            wA(i,j) = 0.0
            wB(i,j) = 0.0
            wD(i,j) = 300.0
            wE(i,j) = 0.0
        End Do
        Do j=31,34
            wC(i,j) = 0.0
            wA(i,j) = 0.0
            wB(i,j) = 0.0
            wD(i,j) = 0.0
            wE(i,j) = 300.0
        End Do
        Do j=35,M
            wC(i,j) = 300.0
            wA(i,j) = 0.0
            wB(i,j) = 0.0
            wD(i,j) = 0.0
            wE(i,j) = 0.0
        End Do
    End do

     !Do i=1,N
      !  Do j=1,3
       !    wB(i,j)=-1.0
       !    wC(i,j)=1.0
       ! End Do
     !End Do
     !wD=wB


  end If
  Return

End Subroutine Rand_Field
!------------------------------- Calculating Phi-----------------------------------------------
Subroutine phi(dummy)
  Use global
  Implicit None
  
  ! The integration here is the Simpsons rule
  Integer                :: dummy
  Integer                :: i,j,s
 
  pA=0.0d0
  pB=0.0d0
  pC=0.0d0
  pE=0.0d0
  pD=0.0d0
 
  ! Calculating PhiA
  !----------------------------------------------------
  Do i=1,N
     Do j=1,M
        Do s=0,NA
           If ((s==0).or.(s==NA)) then
              pA(i,j)=pA(i,j)+(qA(i,j,s)*qAD(i,j,(NA-s)))
           ElseIf ((Mod(s,2).ne.0).and.(s.ne.NA)) then
              pA(i,j)=pA(i,j)+4.0*(qA(i,j,s)*qAD(i,j,(NA-s)))
           ElseIf ((Mod(s,2)==0).and.(s.ne.0).and.(s.ne.NA)) then
              pA(i,j)=pA(i,j)+2.0*(qA(i,j,s)*qAD(i,j,(NA-s)))
           Else
              Print*,"Something gone wrong in pA calculation"
           End If
        End Do
        pA(i,j)=pA(i,j)*delt/3.0d0
     End Do
  End Do
  !----------------------------------------------------
  ! Calculating PhiB
  !----------------------------------------------------
  Do i=1,N
     Do j=1,M
        Do s=0,NB
           If ((s==0).or.(s==NB)) then
              pB(i,j)=pB(i,j)+(qB(i,j,s)*qBD(i,j,(NB-s)))
           ElseIf ((Mod(s,2).ne.0).and.(s.ne.NB)) then
              pB(i,j)=pB(i,j)+4.0*(qB(i,j,s)*qBD(i,j,(NB-s)))
           ElseIf ((Mod(s,2)==0).and.(s.ne.0).and.(s.ne.NB)) then
              pB(i,j)=pB(i,j)+2.0*(qB(i,j,s)*qBD(i,j,(NB-s)))
           Else
              Print*,"Something gone wrong in pB calculation"
           End If
        End Do
        pB(i,j)=pB(i,j)*delt/3.0d0
     End Do
  End Do
  !----------------------------------------------------
  ! Calculating PhiC
  !----------------------------------------------------
  Do i=1,N
     Do j=1,M
        Do s=0,NC
           If ((s==0).or.(s==NC)) then
              pC(i,j)=pC(i,j)+(qC(i,j,s)*qC(i,j,(NC-s)))
           ElseIf ((Mod(s,2).ne.0).and.(s.ne.NC)) then
              pC(i,j)=pC(i,j)+4.0*(qC(i,j,s)*qC(i,j,(NC-s)))
           ElseIf ((Mod(s,2)==0).and.(s.ne.0).and.(s.ne.NC)) then
              pC(i,j)=pC(i,j)+2.0*(qC(i,j,s)*qC(i,j,(NC-s)))
           Else
              Print*,"Something gone wrong in pC calculation"
           End If
        End Do
        pC(i,j)=pC(i,j)*delt/3.0d0
     End Do
  End Do
  !----------------------------------------------------
  ! Calculating PhiE
  !----------------------------------------------------
  Do i=1,N
     Do j=1,M
        Do s=0,NE
           If ((s==0).or.(s==NE)) then
              pE(i,j)=pE(i,j)+(qE(i,j,s)*qED(i,j,(NE-s)))
           ElseIf ((Mod(s,2).ne.0).and.(s.ne.NE)) then
              pE(i,j)=pE(i,j)+4.0*(qE(i,j,s)*qED(i,j,(NE-s)))
           ElseIf ((Mod(s,2)==0).and.(s.ne.0).and.(s.ne.NE)) then
              pE(i,j)=pE(i,j)+2.0*(qE(i,j,s)*qED(i,j,(NE-s)))
           Else
              Print*,"Something gone wrong in pE calculation"
           End If
        End Do
        pE(i,j)=pE(i,j)*delt/3.0d0
     End Do
  End Do 
  !----------------------------------------------------
  ! Calculating PhiD
  !----------------------------------------------------
  Do i=1,N
     Do j=1,M
        Do s=0,ND
           If ((s==0).or.(s==ND)) then
              pD(i,j)=pD(i,j)+(qD(i,j,s)*qDD(i,j,(ND-s)))
           ElseIf ((Mod(s,2).ne.0).and.(s.ne.ND)) then
              pD(i,j)=pD(i,j)+4.0*(qD(i,j,s)*qDD(i,j,(ND-s)))
           ElseIf ((Mod(s,2)==0).and.(s.ne.0).and.(s.ne.ND)) then
              pD(i,j)=pD(i,j)+2.0*(qD(i,j,s)*qDD(i,j,(ND-s)))
           Else
              Print*,"Something gone wrong in pD calculation"
           End If
        End Do
        pD(i,j)=pD(i,j)*delt/3.0d0
     End Do
  End Do 

  !Do i=1,20
     !do j=1,M
        !pA(i,j)=0.0
        !pB(i,j)=0.0
        !pE(i,j)=0.0
        !pD(i,j)=0.0
     !End do
  !End Do

  pA=exp(muABDE)*(pA)
  pB=exp(muABDE)*(pB)
  pC=exp(kappaC*muC)*(pC)/kappaC
  pE=exp(muABDE)*(pE)
  pD=exp(muABDE)*(pD)
 
  Return

End Subroutine phi
!---------------------------- Calculating eta -------------------------------------------------
Subroutine Pressure( )
  Use global
  Implicit None

  Integer              :: i,j
 
  Do i=1,N
     Do j=1,M
        eta(i,j)=eta(i,j)-(1.0-(pA(i,j)+pB(i,j)+pC(i,j)+pE(i,j)+pD(i,j)))
     End Do
  End Do
  
  Return

End Subroutine Pressure
!---------------------------- Calculating eta2 -------------------------------------------------
! eta2 is used for pinning the system in a configuration
Subroutine Pressure2( )
  Use global
  Implicit None

  Integer              :: i,j
  Integer              :: dum_rad

  print*,"Hello"
 
  ! Ntip and Mtip are defined in the main.f90
  Do i=Ntip,Ntip
     Do j=Mtip,Mtip
        eta2(i,j)=eta2(i,j)-10.0*(pB(i,j)-pD(i,j))
     End Do
  End Do
 
  Return
End Subroutine Pressure2
!------------------- Calculating Partition Function -------------------------------------------
Subroutine Q_partition( )
  Use global
  Implicit None

  Integer            :: i,j

  !++++++++++++++++ 2-D trapz rule++++++++++++++++++++++++++++
 
  Q_ABDE=0.0d0
  Q_C=0.0d0
  !Q_ED=0.0d0

  Do i=2,N-1
     Do j=2,M-1
        Q_ABDE=Q_ABDE+qE(i,j,NE)*((real(i-1)*delr)+(r_0))*delr*delz
        Q_C=Q_C+qC(i,j,NC)*((real(i-1)*delr)+(r_0))*delr*delz
     End Do
  End Do

  i=1
  Do j=2,M-1
        Q_ABDE=Q_ABDE+(0.5d0*qE(i,j,NE)*((real(i-1)*delr)+(r_0))*delr*delz)
        Q_C=Q_C+(0.5d0*qC(i,j,NC)*((real(i-1)*delr)+(r_0))*delr*delz)
  End Do

  i=N
  Do j=2,M-1
     Q_ABDE=Q_ABDE+(0.5d0*qE(i,j,NE)*((real(i-1)*delr)+(r_0))*delr*delz)
     Q_C=Q_C+(0.5d0*qC(i,j,NC)*((real(i-1)*delr)+(r_0))*delr*delz)
  End Do

  j=1
  Do i=2,N-1
     Q_ABDE=Q_ABDE+(0.5d0*qE(i,j,NE)*((real(i-1)*delr)+(r_0))*delr*delz)
     Q_C=Q_C+(0.5d0*qC(i,j,NC)*((real(i-1)*delr)+(r_0))*delr*delz)
  End Do

  j=M
  Do i=2,N-1
        Q_ABDE=Q_ABDE+(0.5d0*qE(i,j,NE)*((real(i-1)*delr)+(r_0))*delr*delz)
        Q_C=Q_C+(0.5d0*qC(i,j,NC)*((real(i-1)*delr)+(r_0))*delr*delz)
  End Do


  Q_ABDE=Q_ABDE+0.25d0*delr*delz*(qE(1,1,NE)*(real(1-1)*delr+r_0)+qE(1,M,NE)*(real(1-1)*delr+r_0) &
        +qE(N,1,NE)*(real(N-1)*delr+(r_0))+qE(N,M,NE)*(real(N-1)*delr+(r_0)))

  !Q_ABDE=Q_ABDE+0.25d0*delr*delz*(qAD(1,1,NA)*(real(1-1)*delr+r_0)+qAD(1,M,NA)*(real(1-1)*delr+r_0) &
  !      +qAD(N,1,NA)*(real(N-1)*delr+(r_0))+qAD(N,M,NA)*(real(N-1)*delr+(r_0)))

  Q_C=Q_C+0.25d0*delr*delz*(qC(1,1,NC)*(real(1-1)*delr+r_0)+qC(1,M,NC)*(real(1-1)*delr+r_0) &
       +qC(N,1,NC)*(real(N-1)*delr+(r_0))+qC(N,M,NC)*(real(N-1)*delr+(r_0)))


 
  
  Q_ABDE=Q_ABDE*2.0*pi
  Q_C=Q_C*2.0*pi
  Q_ABDE=Q_ABDE/vol
  Q_C=Q_C/vol

  Return
End Subroutine Q_partition
!-------------------------------------------------------------------------------------
Subroutine New_Fields( )
  Use global
  Implicit None

  Integer                       :: i,j

  ! These are used for conv condition
  conv_w=0.0d0
  conv_p=0.0d0
  Do i=1,N
     Do j=1,M

        dwA(i,j)=xAB*pB(i,j)+xAC*pC(i,j)+xAE*pE(i,j)+xAD*pD(i,j)-eta2(i,j)+eta(i,j)-wA(i,j)
        dwB(i,j)=xAB*pA(i,j)+xBC*pC(i,j)+xBE*pE(i,j)+xBD*pD(i,j)+eta2(i,j)+eta(i,j)-wB(i,j)
        dwC(i,j)=xAC*pA(i,j)+xBC*pB(i,j)+xCE*pE(i,j)+xCD*pD(i,j)+eta(i,j)-wC(i,j)
        dwE(i,j)=xAE*pA(i,j)+xBE*pB(i,j)+xCE*pC(i,j)+xED*pD(i,j)-eta2(i,j)+eta(i,j)-wE(i,j)
        dwD(i,j)=xAD*pA(i,j)+xBD*pB(i,j)+xCD*pC(i,j)+xED*pE(i,j)+eta2(i,j)+eta(i,j)-wD(i,j)

        dpp(i,j)=1.0d0-(pA(i,j)+pB(i,j)+pC(i,j)+pE(i,j)+pD(i,j))


        conv_w=conv_w+abs(dwA(i,j)+dwB(i,j)+dwC(i,j)+dwE(i,j)+dwD(i,j))
        conv_p=conv_p+abs(dpp(i,j))
        ! updating the omega condition
        wA(i,j)=wA(i,j)+sig*dwA(i,j)-sig2*dpp(i,j)
        wB(i,j)=wB(i,j)+sig*dwB(i,j)-sig2*dpp(i,j)
        wC(i,j)=wC(i,j)+sig*dwC(i,j)-sig2*dpp(i,j)
        wE(i,j)=wE(i,j)+sig*dwE(i,j)-sig2*dpp(i,j)
        wD(i,j)=wD(i,j)+sig*dwD(i,j)-sig2*dpp(i,j)
     End Do
  End Do
 
  conv_w=conv_w/(N*M)
  conv_p=conv_p/(N*M)

  Return  
End Subroutine New_Fields
!********************* Writting the data ************************
Subroutine profile(msg)
  Use global
  Implicit None

  Integer              :: i,j,ii,msg
  Real                 :: cr,cz,ct,epsct
  Real                 :: px,py,pz
  Character(len=50)    :: outfile,orderpar

  If (msg==1) then
     epsct=(2.0*pi)/N
     open(2,file='./results/profile1.dat',status='replace')
     Do i=1,N
        Do j=1,M
           write(2,*) (i*delr),(j*delz),pA(i,j),pB(i,j),pC(i,j),pE(i,j),pD(i,j)
           write(2,*)
        End Do
     End Do
     close(2)
     
  ElseIf (msg==2) then
     
     write(orderpar,'(f10.4)') OP

     outfile='./results/outphi/profile/'//TRIM(orderpar)//'.dat'
  
     epsct=(2.0*pi)/N
     open(2,file=outfile)
     Do i=1,N
        Do j=1,M
           write(2,*) i,j,pA(i,j),pB(i,j),pC(i,j),pE(i,j),pD(i,j),muC,muABDE
           write(2,*)
        End Do
     End Do
     close(2)

     open(4,file='./results/xyz.dat',status='replace')
     ct=0.0
     Do i=1,N
     px=(((i-1)*delr+r_0))*cos(ct)
     py=(((i-1)*delr+r_0))*sin(ct)
     pz=(i-1)*delz
     write(4,*) ct,(i-1)*delr,(i-1)*delz
     ct=ct+epsct
     end Do
     close(4)
     open(8,file='./results/ABCDE.dat',status='replace')
     Do ii=1,N ! theta
     Do i=1,N ! r
     Do j=1,N ! z
     write(8,*) pA(i,j),pB(i,j),pC(i,j),pE(i,j),pD(i,j)
     End Do
     End Do
     End Do
     close(8)
  End If

  Return

End Subroutine profile
!******************** The Interaction Matrix  *************************
Subroutine XMatrix( )
  Use global
  Implicit None

  Integer           :: errorflag

  !Here I am putting the interaction parameters into a matrix
  XM(1,1)=0.0d0
  XM(2,1)=xAB
  XM(3,1)=xAC
  XM(4,1)=xAE
  XM(5,1)=xAD

  XM(1,2)=xAB
  XM(2,2)=0.0d0
  XM(3,2)=xBC
  XM(4,2)=xBE
  XM(5,2)=xBD


  XM(1,3)=xAC
  XM(2,3)=xBC
  XM(3,3)=0.0d0
  XM(4,3)=xCE
  XM(5,3)=xCD

  XM(1,4)=xAE
  XM(2,4)=xBE
  XM(3,4)=xCE
  XM(4,4)=0.0d0
  XM(5,4)=xED

  XM(1,5)=xAD
  XM(2,5)=xBD
  XM(3,5)=xCD
  XM(4,5)=xED
  XM(5,5)=0.0d0


  Return

End Subroutine XMatrix
!*********************************************
Subroutine cleanme( )
  Implicit None
  
  ! Emptying data files ***********************
  open(8,file='./results/profile1.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/xyz.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/ABCDE.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/profile4.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/phi1Dr.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/phi1Dz.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/phi.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/omega.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/b.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  open(8,file='./results/a.dat',status='old')
  write(8,*) 
  write(8,*)
  close(8)
  !*********************************************


End Subroutine cleanme
