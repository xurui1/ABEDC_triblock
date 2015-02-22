!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This is where I am solving the diff equation using the EDI method
! I am not a good programmer, so this is messy
! Since I wrote this code a 2D xy planer system, some of the subroutines are
! still called x and y and I dont wanna change things around.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! This subroutine makes the lower, middle and upper diagonal elements of the matrix
! when scanning over the r-direction (or index i)
Subroutine E_Matrix_y(ii)
  Use global
  Implicit None
  
  Integer                    :: i,j,ii
  Real*8                     :: alphaE,betaEL,betaEU

  ! Here U means upper and L means lower the other one is the middle element
  DiagEx=0.0d0
  DiagEUx=0.0d0
  DiagELx=0.0d0
  Do i=1,N
     DiagEx(i)=1.0d0+(delt/delr**2)+((delt/2.0d0)*wE(i,ii))
  End Do

  Do i=2,N-1
     DiagEUx(i)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(i-1))+r_0)*4.0d0*delr))
  End Do

  Do i=1,N-2
     DiagELx(i)=-(delt/(2.0d0*delr**2))+(delt/(((delr*real(i+1-1))+r_0)*4.0d0*delr))
  End Do

  ! The corners are different because we are using a 0-der boundary condition
  DiagEUx(1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(1-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(1-1))+r_0)*4.0d0*delr))

  DiagELx(N-1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(N-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(N-1))+r_0)*4.0d0*delr))

    
  Return

End Subroutine E_Matrix_y
!--------------------------------------------------------------------------------
! Same thing here, I am calculating the matrix, this is when we are scanning 
! over the z-direction or (j)
Subroutine E_Matrix_x(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaE,betaE

  DiagEy=0.0d0
  DiagEUy=0.0d0
  DiagELy=0.0d0

  Do i=1,M
     DiagEy(i)=1.0d0+delt/delz**2
  End Do
  Do i=1,M-1
     DiagEUy(i)=-delt/(2.0d0*delz**2)
  End Do
  Do i=1,M-1
     DiagELy(i)=-delt/(2.0d0*delz**2)
  End Do
  ! Samething here, cornenrs are special, due to bc 
  DiagEUy(1)=2.0d0*DiagEUy(1)
  DiagELy(M-1)=2.0d0*DiagELy(M-1)

  Return

End Subroutine E_Matrix_x
!-------------------------------------qE------------------------------------------------
! This is the forwad propagator. the qE_0 is a intermidiate step, you can ignor it
Subroutine qE_forward( )
  Use global
  Implicit None

  Integer                  :: errorflag
  Integer                  :: s,i,j
  Integer                  :: errorflagE
  Real*8                   :: gamma,bettaU,bettaL,betta

  qE_0=0.0d0
  qE=0.0d0
  ! Initializing the qs
  Do i=1,N
     Do j=1,M
        qE_0(i,j)=qD(i,j,ND)
        qE(i,j,0)=qD(i,j,ND)
     End Do
  End Do
 
  Do s=1,NE

     bEy=0.0d0
     !**********************Scan Over z***************
     Do i=1,N
        
        Call E_Matrix_x(i)  
        ! Depending on the bc, you have to populate the right hand side accordingly.
        If (i==1) then
           Do j=1,M 

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bEy(j)=gamma*qE_0(i,j)+bettaU*qE_0(i+1,j)+bettaL*qE_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bEy(j)=gamma*qE_0(i,j)+bettaU*qE_0(i-1,j)+bettaL*qE_0(i-1,j)

           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bEy(j)=gamma*qE_0(i,j)+bettaU*qE_0(i+1,j)+bettaL*qE_0(i-1,j)

           End Do
        End If
     
        DiagEdy=DiagEy
        DiagEUdy=DiagEUy
        DiagELdy=DiagELy
        ! This is an external library for solving the matrix, using LEPECK
        !Call DGTSV(M,1,DiagELdy,DiagEdy,DiagEUdy,bEy,M,errorflagE)
        Call TDMA(M,DiagELdy,DiagEdy,DiagEUdy,bEy)
      
        Do j=1,M
           qE(i,j,s)=bEy(j)
        End Do

     End Do
     
     bEy=0.0d0
     Do i=1,N
        Do j=1,M
           qE_0(i,j)=qE(i,j,s)
        End Do
     End Do

     bEx=0.0d0
     !**********************Scan Over r***************

     Do j=1,M

        Call E_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bEx(i)=gamma*qE_0(i,j)+betta*qE_0(i,j+1)+betta*qE_0(i,j+1)
        
           End do
        ElseIf (j==M) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bEx(i)=gamma*qE_0(i,j)+betta*qE_0(i,j-1)+betta*qE_0(i,j-1)
      
           End do
        Else
           Do i=1,N   
     
              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bEx(i)=gamma*qE_0(i,j)+betta*qE_0(i,j+1)+betta*qE_0(i,j-1)

           End Do
        End If
     
        
        DiagEdx=DiagEx
        DiagEUdx=DiagEUx
        DiagELdx=DiagELx
        !Call DGTSV(N,1,DiagELdx,DiagEdx,DiagEUdx,bEx,N,errorflagE)
        Call TDMA(N,DiagELdx,DiagEdx,DiagEUdx,bEx)

        Do i=1,N
           qE(i,j,s)=bEx(i)
        End Do
     End Do
    
     Do i=1,N
        Do j=1,M
           qE_0(i,j)=qE(i,j,s)
        End Do
     End Do
     
  End Do

  Return

End Subroutine qE_forward
!-------------------------------------qED------------------------------------------------
! This is for the complementary prpagator, qED(agor)
Subroutine qED_forward( )
  Use global
  Implicit None

  Integer                  :: s,i,j
  Integer                  :: errorflagE
  Real*8                   :: gamma,bettaU,bettaL,betta

  ! Same as before, but the initial condition is from qD
  qED_0=0.0d0
  qED=0.0d0
  
  Do i=1,N
     Do j=1,M
        qED(i,j,0)=1.0
        qED_0(i,j)=1.0
     End Do
  End Do
    
  Do s=1,NE
     
     bEy=0.0d0
     !**********************Scan Over z***************
     Do i=1,N
        Call E_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bEy(j)=gamma*qED_0(i,j)+bettaU*qED_0(i+1,j)+bettaL*qED_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))  
    
              bEy(j)=gamma*qED_0(i,j)+bettaU*qED_0(i-1,j)+bettaL*qED_0(i-1,j)
          
           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wE(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bEy(j)=gamma*qED_0(i,j)+bettaU*qED_0(i+1,j)+bettaL*qED_0(i-1,j)

           End Do
        End If
        
        DiagEdy=DiagEy
        DiagEUdy=DiagEUy
        DiagELdy=DiagELy
        !Call DGTSV(M,1,DiagELdy,DiagEdy,DiagEUdy,bEy,M,errorflagE)
        Call TDMA(M,DiagELdy,DiagEdy,DiagEUdy,bEy)

        Do j=1,M
           qED(i,j,s)=bEy(j)
        End Do
     End Do

     bEy=0.0d0
     Do i=1,N
        Do j=1,M
           qED_0(i,j)=qED(i,j,s)
        End Do
     End Do

     bEx=0.0d0
     !**********************Scan Over r***************
     Do j=1,M
        Call E_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bEx(i)=gamma*qED_0(i,j)+betta*qED_0(i,j+1)+betta*qED_0(i,j+1)

           End Do
        ElseIf (j==M) then
           Do i=1,N   

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bEx(i)=gamma*qED_0(i,j)+betta*qED_0(i,j-1)+betta*qED_0(i,j-1)

           End Do
        Else
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bEx(i)=gamma*qED_0(i,j)+betta*qED_0(i,j+1)+betta*qED_0(i,j-1)

           End Do
        End If
        

        DiagEdx=DiagEx
        DiagEUdx=DiagEUx
        DiagELdx=DiagELx
        !Call DGTSV(N,1,DiagELdx,DiagEdx,DiagEUdx,bEx,N,errorflagE)     
        Call TDMA(N,DiagELdx,DiagEdx,DiagEUdx,bEx)     
    
        Do i=1,N
           qED(i,j,s)=bEx(i)
        End Do

     End Do
     
     Do i=1,N
        Do j=1,M
           qED_0(i,j)=qED(i,j,s)
        End Do
     End Do
     

  End Do
  
  Return

End Subroutine qED_forward
