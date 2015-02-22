! Here I am calculating the free energy, integration is done using the Trapoz-Rule
Subroutine FreeEnergy(dummy)
  Use global
  Implicit None

  Integer                   :: dummy
  Integer                   :: i,j,ii,jj
  Real*8                    :: F1,F2,F3,F4,F5,F6,F7,F8,F9
  Real*8                    :: FF1,FF2,FF3,FF4,FF5,FF6,FF7,FF8,FF9
  Real*8                    :: fEhom
  Real*8                    :: p_vect(5),w_vect(5)

  fE=0.0d0

  !*******************Bulk************************
  w_vect=0.0d0
  p_vect=0.0d0

  F1=0.0d0
  FF1=0.0d0

  Do i=2,N-1
     Do j=2,M-1
        p_vect(1)=pA(i,j)
        p_vect(2)=pB(i,j)
        p_vect(3)=pC(i,j)
        p_vect(4)=pE(i,j)
        p_vect(5)=pD(i,j)

        w_vect(1)=wA(i,j)
        w_vect(2)=wB(i,j)
        w_vect(3)=wC(i,j)
        w_vect(4)=wE(i,j)
        w_vect(5)=wD(i,j)
        
        Do ii=1,5
           Do jj=ii,5
              F1=F1+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz  
           End Do
           FF1=FF1+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz  
        End Do
     End Do
  End Do
  !*******************Side 1************************
  w_vect=0.0d0
  p_vect=0.0d0
  F2=0.0d0
  FF2=0.0d0
  Do i=1,1
     Do j=2,M-1
        p_vect(1)=pA(i,j)
        p_vect(2)=pB(i,j)
        p_vect(3)=pC(i,j)
        p_vect(4)=pE(i,j)
        p_vect(5)=pD(i,j)
 
        w_vect(1)=wA(i,j)
        w_vect(2)=wB(i,j)
        w_vect(3)=wC(i,j)
        w_vect(4)=wE(i,j)
        w_vect(5)=wD(i,j)

        Do ii=1,5
           Do jj=ii,5
              F2=F2+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
           End Do
           FF2=FF2+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
        End Do
     End Do
  End Do
  !*******************Side 2************************
  w_vect=0.0d0
  p_vect=0.0d0
  F3=0.0d0
  FF3=0.0d0
  Do i=N,N
     Do j=2,M-1
        p_vect(1)=pA(i,j)
        p_vect(2)=pB(i,j)
        p_vect(3)=pC(i,j)
        p_vect(4)=pE(i,j)
	p_vect(5)=pD(i,j)

        w_vect(1)=wA(i,j)
        w_vect(2)=wB(i,j)
        w_vect(3)=wC(i,j)
	w_vect(4)=wE(i,j)
	w_vect(5)=wD(i,j)
        
        Do ii=1,5
           Do jj=ii,5
              F3=F3+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
           End Do
           FF3=FF3+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
        End Do
     End Do
  End Do
  !*******************Side 3************************
  w_vect=0.0d0
  p_vect=0.0d0
  F4=0.0d0
  FF4=0.0d0
  Do i=2,N-1
     Do j=1,1
        p_vect(1)=pA(i,j)
        p_vect(2)=pB(i,j)
        p_vect(3)=pC(i,j)
        p_vect(4)=pE(i,j)
	p_vect(5)=pD(i,j)

        w_vect(1)=wA(i,j)
        w_vect(2)=wB(i,j)
        w_vect(3)=wC(i,j)
	w_vect(4)=wE(i,j)
	w_vect(5)=wD(i,j)
        
        Do ii=1,5
           Do jj=ii,5
              F4=F4+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
           End Do
           FF4=FF4+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
        End Do
     End Do
  End Do
  !*******************Side 4************************
  w_vect=0.0d0
  p_vect=0.0d0
  F5=0.0d0
  FF5=0.0d0
  Do i=2,N-1
     Do j=M,M
        p_vect(1)=pA(i,j)
        p_vect(2)=pB(i,j)
        p_vect(3)=pC(i,j)
        p_vect(4)=pE(i,j)
	p_vect(5)=pD(i,j)

        w_vect(1)=wA(i,j)
        w_vect(2)=wB(i,j)
        w_vect(3)=wC(i,j)
	w_vect(4)=wE(i,j)
	w_vect(5)=wD(i,j)
        
        Do ii=1,5
           Do jj=ii,5
              F5=F5+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
           End Do
           FF5=FF5+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.5d0  
        End Do
     End Do
  End Do
  !*******************Corrner 1************************
  i=1
  j=1
  w_vect=0.0d0
  p_vect=0.0d0
  F6=0.0d0
  FF6=0.0d0
  p_vect(1)=pA(i,j)
  p_vect(2)=pB(i,j)
  p_vect(3)=pC(i,j)
  p_vect(4)=pE(i,j)
  p_vect(5)=pD(i,j)
  
  w_vect(1)=wA(i,j)
  w_vect(2)=wB(i,j)
  w_vect(3)=wC(i,j)
  w_vect(4)=wE(i,j)
  w_vect(5)=wD(i,j)

  Do ii=1,5
     Do jj=ii,5
        F6=F6+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
     End Do
     FF6=FF6+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
  End Do
   !*******************Corrner 2************************
  i=N
  j=1
  w_vect=0.0d0
  p_vect=0.0d0
  F7=0.0d0
  FF7=0.0d0
  p_vect(1)=pA(i,j)
  p_vect(2)=pB(i,j)
  p_vect(3)=pC(i,j)
  p_vect(4)=pE(i,j)
  p_vect(5)=pD(i,j)
  
  w_vect(1)=wA(i,j)
  w_vect(2)=wB(i,j)
  w_vect(3)=wC(i,j)
  w_vect(4)=wE(i,j)
  w_vect(5)=wD(i,j)
  
  Do ii=1,5
     Do jj=ii,5
        F7=F7+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
     End Do
     FF7=FF7+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
  End Do
  !*******************Corrner 3************************
  i=1
  j=M
  w_vect=0.0d0
  p_vect=0.0d0
  F8=0.0d0
  FF8=0.0d0
  p_vect(1)=pA(i,j)
  p_vect(2)=pB(i,j)
  p_vect(3)=pC(i,j)
  p_vect(4)=pE(i,j)
  p_vect(5)=pD(i,j)
  
  w_vect(1)=wA(i,j)
  w_vect(2)=wB(i,j)
  w_vect(3)=wC(i,j)
  w_vect(4)=wE(i,j)
  w_vect(5)=wD(i,j)

  Do ii=1,5
     Do jj=ii,5
        F8=F8+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
     End Do
     FF8=FF8+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
  End Do
  !*******************Corner 4************************
  i=N
  j=M
  w_vect=0.0d0
  p_vect=0.0d0
  F9=0.0d0
  FF9=0.0d0
  p_vect(1)=pA(i,j)
  p_vect(2)=pB(i,j)
  p_vect(3)=pC(i,j)
  p_vect(4)=pE(i,j)
  p_vect(5)=pD(i,j)

  w_vect(1)=wA(i,j)
  w_vect(2)=wB(i,j)
  w_vect(3)=wC(i,j)
  w_vect(4)=wE(i,j)
  w_vect(5)=wD(i,j)

  Do ii=1,5
     Do jj=ii,5
        F9=F9+p_vect(ii)*p_vect(jj)*XM(ii,jj)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
     End Do
     FF9=FF9+p_vect(ii)*w_vect(ii)*delr*((real(i-1)*delr)+(r_0))*delz*0.25d0  
  End Do

  ! Calculating the Free Energy

  fE=(F1+F2+F3+F4+F5+F6+F7+F8+F9)-(FF1+FF2+FF3+FF4+FF5+FF6+FF7+FF8+FF9)
  fE=fE*2.0*pi
  fE=(fE/vol)-(exp(muABDE)*Q_ABDE)-(exp(muC*kappaC)*Q_C/kappaC)

   !-(exp(muED*kappaED)*Q_ED/kappaED)

  ! printing individual parts of fE
  !print*,-(exp(muAB)*Q_AB)-(exp(muC*kappaC)*Q_C/kappaC)-(exp(muED*kappaED)*Q_ED/kappaED), &
       !(F1+F2+F3+F4+F5+F6+F7+F8+F9)*2.0*pi/vol, &
       !(FF1+FF2+FF3+FF4+FF5+FF6+FF7+FF8+FF9)*2.0*pi/vol
  
  dfffE=0.0d0
  dfffE=abs(fE-fE_old)
  fE_old=fE


  Return
End Subroutine FreeEnergy
