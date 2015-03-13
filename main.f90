!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!This is the main program which will calculate the free-energy for a given
! structure. I am using the Grand-Canonical ensemble. This is the code for the
! AB/C system in the 2D cylindrical geometry, with the (r) and (z) axis.
!                                                 Written By: Ashkan Dehghan
!                                                          McMaster University
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Program Main
  Use global
  use omp_lib
  Implicit None

  ! Declaring the variables
  Integer             :: once,msg1,msg2
  Integer             :: s,s2,i,j,initial,jj        ! Some counting indices
  Integer             :: muED_up,muED_down          ! For steping up and down.
  Real*8              :: ptot,phiA,phiB,phiC,phiE,phiD,frac
  Real*8              :: fEnow,fEold                ! Used for calculating the dfE
  Real*8              :: D_r,D_z                    ! Box size in the r and z direction is Rg^2
  Real*8              :: Area,Tip_R
  Real*8              :: F_int_0, AB_0

  ! Type in here the delF and omega if your
  ! using the difference method
  ! Bilayer delfE and excess copolymer
  F_int_0=0.0
  AB_0=0.0
  
  ! Reading in the parameters
  open(1,file='./data/data.dat',status='old')
  read(1,*) D_r,D_z
  read(1,*) N,M
  read(1,*) NA,NB,NC
  read(1,*) NE,ND
  read(1,*) sig,sig2
  read(1,*) xAB,xAC,xAE,xAD
  read(1,*) xBC,xBE,xBD
  read(1,*) xCE,xCD
  read(1,*) xED
  read(1,*) muABDE,muC
  close(1)

  ! The step sizes are defined here
  delt=1.0/(NA+NB)
  delr=D_r/real(N-1)
  delz=D_z/real(M-1)
  
  ! Distance from the center of cylinder 
  r_0=5.0

  ! Initial Set Up 0=Random 1=Read from file 2=make your own condition
  initial=1
  ! If you have chosen 1 for Initial, then you can pick what initial configuration
  ! you can have, 1=on 0=off
  !()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
  ! If you turn on ten_find then it will find the mu for the tensionless bilayer

  ten_find=0

  ! up + down - 
  muED_up=0
  muED_down=0

  ! If ten_find is on, turn off bilayer
  bilayer=0
  once=1
  disk=1

  !()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
  

  ! Cleaning the data files
  If (disk==1) then
     open(8,file='./results/fE_disk.dat',status='old')
     write(8,*) 
     write(8,*)
     close(8)
  ElseIf (bilayer==1) then
     open(8,file='./results/fE_bilayer.dat',status='old')
     write(8,*) 
     write(8,*)
     close(8)
  End If

  ! This will make a interaction parameter matrix, its located in mod1.f90
  Call XMatrix( )

  
  ! This is the pinning condition position
  If (disk==1) then
     Ntip=20
     Mtip=1
     ! In Rg
     Tip_r=(Ntip-1)*delr
  End If

  ! Setting kappa and fA
  kappaC=real(NC)/real(NA+NB)
  fracA=real(NA)/real(NA+NB)
  fracB=real(NB)/real(NA+NB)
  fracE=real(NE)/real(NE+ND)
  fracD=real(ND)/real(NE+ND)
  ! This will allocate the proper memory
  Call Allocate_mod(1)
  iter=1
  Do

     ! This where I will put the secant method
     If (ten_find==1) Call secant( )

     ! For initial configuration of the omega fields
     Call Rand_Field(initial)
     iter=iter+1

     ! setting the incomp and pinning to zero initially
     eta2=0.0
     eta=0.0
     ! calculates the homogenous free energy
     Call fE_homo( )
     
     ! Here the volume and Area are defined, depending on the configuration
     vol=2.0*pi*(0.5*(M-1)*delz*(((N-1)*delr+r_0)**2-(r_0)**2))
     If (bilayer==1) Area=pi*(((N-1)*delr+r_0)**2-(r_0)**2)
     If (disk==1) Area=pi*((Tip_r+r_0)**2-(r_0)**2)
   
     ! Cleans all the data fiels
     Call cleanme( )
     
     s2=0
     fE_old=0.0d0
    
     ! do stop either in scf-conv is reached or s=100000
     in:Do s=1,100000

        !open(4,file='./results/a.dat',status='old',access='append')
        ! Here I do the solving of the diff-eq They are in modA.f90 and so on.
        Call qA_forward( )  !First A block
        Call qB_forward( )  !First B block
        Call qD_forward( )  !Second B block
        Call qE_forward( )  !Second A block

        Call qC_forward( )  !Solvent

        Call qED_forward( ) !Second A block
        Call qDD_forward( ) !Second B block
        Call qBD_forward( ) !First B block
        Call qAD_forward( ) !First A block


        ! Calcularing the partition function
        Call Q_partition( )
        ! Calculating the phi
        Call phi( )
        Call Pressure( )
        ! pressure 2 is the pinning condition, only needed for the disk
        If (disk==1) Call Pressure2( )
        Call FreeEnergy( )
        Call totalphi(ptot,phiA,phiB,phiC,phiE,phiD)
        Call New_Fields( )
        
        !****************************************************************
        !**********************Print During the Run**********************
       print*,"Dfe: ",((fE-fE_hom)*vol)/Area,"phiA: ",phiA,"phiB: ",phiB,"phiC: ",phiC
       print*,"fE:", fE, "fE_hom:", fE_hom
  
        
        write(4,*) real(s),fE,dfffE
        
        ! Here I write some temp data fiels for plotting 
        If (s>s2) then
6002       format(I3,7E20.10,I10,I3)
           Call profile(1)
           open(8,file='./results/phi.dat',status='replace')
           open(2,file='./results/phi1Dr.dat',status='replace')
           open(7,file='./results/phi1Dz.dat',status='replace')
           If (disk==1) open(4,file='./results/omega.disk',status='replace')
           If (bilayer==1) open(4,file='./results/omega.bilayer',status='replace')
           Do i=1,N
              Do j=1,M
                 write(2,*) (i*delr),(j*delz),pA(i,M/2),pB(i,M/2),pC(i,M/2),pE(i,M/2),pD(i,M/2)
                 write(7,*) (i*delr),(j*delz),pA(N/2,j),pB(N/2,j),pC(N/2,j),pE(N/2,j),pD(N/2,j)
                 write(4,6002) 1,wA(i,j),wB(i,j),wC(i,j),wE(i,j),wD(i,j)
                 write(8,*) (i*delr),(j*delz),pA(i,j),pB(i,j),pC(i,j),pE(i,j),pD(i,j)
                 write(8,*)
              End Do
           End Do
           s2=s2+10
           close(4)
           close(8)
           close(2)
           close(7)
        End If
    
        ! Conv conditions
        If ((conv_p<1.0d-4).and.(conv_w<1.0d-4).and.(dfffE<1.0d-4)) exit in
        
     End Do in
     close(4)
  
     OP=((phiA+phiB)-(phiE+phiD))/((phiA+phiB)+(phiE+phiD))

     ! Prints to screen:
     print*,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
     If (bilayer==1) then
        print*,"Shape chosen is: Bilayer"
        print*,"order_parameter",OP
        print*,"muABDE=",muABDE,"  ","muC=",muC
        print*,"fE",((fE-fE_hom)*vol)/Area
     ElseIf (disk==1) then
        print*,"Shape chosen is: Disk"
        print*,"order_parameter",OP
        print*,"muABDE=",muABDE,"  ","muC=",muC
        print*,"fE",((fE-fE_hom)*vol)/Area
     End If
     print*,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    



     !*************************** Writting to the Data Fields *****************************************
     If (disk==1) then
        open(8,file='./results/fE_disk.dat',status='old',access='append')
        write(8,*) ND,((fE-fE_hom)*vol)/Area,((phiA+phiB+phiE+phiD)-phiAB_hom)*vol/Area,(r_0+Tip_r) &
             ,(phiA+phiB),(phiE+phiD),Area,muC,muABDE
        close(8)
     end If

     If (bilayer==1) then
        open(8,file='./results/fE_bilayer.dat',status='old',access='append')
        write(8,*)  ((fE-fE_hom)*vol)/Area,((phiA+phiB+phiE+phiD)-phiAB_hom)*vol/Area,(r_0+Tip_r) &
             ,(phiA+phiB),(phiE+phiD),Area,muC
        close(8)
     end If

     Call profile(2)
     !************************************************************************************************

     ! Scanning over mu or radius depending on the choice of initial condition
     If ((bilayer==1).and.(once==1)) stop

     If ((disk==1).and.(once==1)) stop

     !If (disk==1) then
      
      !  If (muED_up==1) then
       !    muED=muED+0.1
       !    If (OP<(-0.99d0)) stop
       ! ElseIf (muED_down==1) then
       !    muED=muED-0.1
       !    If (OP>0.99d0) stop
       ! End If
       !   If ((muED_up==0).and.(muED_down==0)) then
       !      stop
       ! End If

     !End If


     If (bilayer==1) then
        If (muC<-5.0) stop
        muC=muC-0.1
     End If

  End Do

  Call Allocate_mod(0)

End Program Main
