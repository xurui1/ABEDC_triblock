Subroutine secant()
  Use global
  Implicit None

  Integer             :: s,s2,i,j,ii,msg,msg1
  Real*8              :: mu1,mu2,mu3
  Real*8              :: fE1,fE2,fE3


  mu1=muC
  mu2=muC+0.01


  bilayer=1
  msg1=1
  If (disk==1) then
     disk=0
  else
     msg1=0
  end If
  
  Call Rand_Field(1)

  ii=1
  msg=0
  scnt:Do
  
     eta2=0.0
     eta=0.0
     Call fE_homo( )
     vol=2.0*pi*(0.5*(M-1)*delz*(((N-1)*delr+r_0)**2-(r_0)**2))
     
     
     s2=0
     fE_old=0.0d0
     Call cleanme()
     !*********************************************************************************************************************
     ! do stop either in scf-conv is reached or s=100000
     in:Do s=1,100000
        ! Here I do the solving of the diff-eq They are in modA.f90 and so on.
        Call qA_forward( )   
        Call qB_forward( ) 
        Call qD_forward( )
        Call qE_forward( )


        Call qC_forward( )


        Call qED_forward( )
        Call qDD_forward( )
        Call qBD_forward( )
        Call qAD_forward( )  



        Call Q_partition( )
        Call phi( )
        Call Pressure( )
        Call FreeEnergy( )
        Call New_Fields( )   
        !****************************************************************
        !**********************Print During the Run**********************
        print*,Conv_w,fE-fE_hom,muC,ii
        !****************************************************************
        !****************************************************************
        ! Here I write some temp data fiels for plotting 
        If (s>s2) then
6002       format(I3,7E20.10,I10,I3)
           Call profile( )
           open(8,file='./results/phi.dat',status='replace')
           open(2,file='./results/phi1Dr.dat',status='replace')
           open(7,file='./results/phi1Dz.dat',status='replace')
           If (disk==1) open(4,file='./results/omega.disk',status='replace')
           If (bilayer==1) open(4,file='./results/omega.bilayer',status='replace')
           Do i=1,N
              Do j=1,M
                 write(2,*) (i*delr),(j*delz),pA(i,1),pB(i,1),pC(i,1),pE(i,1),pD(i,1)
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
        If ((conv_p<1.0d-4).and.(conv_w<1.0d-4).and.(dfffE<1.0d-4)) exit in
     End Do in
     !*********************************************************************************************************************
    
     If (ii==1) then
        fE1=fE-fE_hom
        muC=mu2
     End If
     If (ii==2) then
        fE2=fE-fE_hom
        mu3=mu2-(fE2*((mu2-mu1)/(fE2-fE1)))
        muC=mu3
     End If
     If (ii==3) then
        fE3=fE-fE_hom
        If (abs(fE3)<1.0d-5) then
           exit scnt
        Else
           mu1=mu2
           mu2=mu3
           muC=mu1
           ii=1
           msg=1
        End If
     End If

     print*,muC

     If (msg==0) ii=ii+1
     msg=0

     If (ii>3) print*,"something is wrong in the secant method mod"
     
  End Do scnt
  muC=mu3
 
  If ((msg1==1).and.(disk==0)) then
     disk=1
     bilayer=0
  End If
  
  !print*,"==============================================================="
  !print*,"The muC for a tenssionless bilayer is:"
  !print*,"muC=",muC,"deltafE=",fE-fE_hom
  !print*,"==============================================================="

  Return
  
End Subroutine secant
