!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This is where I am solving the diff equation using the ADI method
! I am not a good programmer, so this is messy
! Since I wrote this code a 2D xy planer system, some of the subroutines are
! still called x and y and I dont wanna change things around.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! This subroutine makes the lower, middle and upper diagonal elements of the matrix
! when scanning over the r-direction (or index i)
Subroutine A_Matrix_y(ii)
  Use global
  Implicit None
  
  Integer                    :: i,j,ii
  Real*8                     :: alphaA,betaAL,betaAU

  ! Here U means upper and L means lower the other one is the middle element
  DiagAx=0.0d0
  DiagAUx=0.0d0
  DiagALx=0.0d0
  Do i=1,N
     DiagAx(i)=1.0d0+(delt/delr**2)+((delt/2.0d0)*wA(i,ii))
  End Do

  Do i=2,N-1
     DiagAUx(i)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(i-1))+r_0)*4.0d0*delr))
  End Do

  Do i=1,N-2
     DiagALx(i)=-(delt/(2.0d0*delr**2))+(delt/(((delr*real(i+1-1))+r_0)*4.0d0*delr))
  End Do

  ! The corners are different because we are using a 0-der boundary condition
  DiagAUx(1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(1-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(1-1))+r_0)*4.0d0*delr))

  DiagALx(N-1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(N-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(N-1))+r_0)*4.0d0*delr))

    
  Return

End Subroutine A_Matrix_y
!--------------------------------------------------------------------------------
! Same thing here, I am calculating the matrix, this is when we are scanning 
! over the z-direction or (j)
Subroutine A_Matrix_x(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaA,betaA

  DiagAy=0.0d0
  DiagAUy=0.0d0
  DiagALy=0.0d0

  Do i=1,M
     DiagAy(i)=1.0d0+delt/delz**2
  End Do
  Do i=1,M-1
     DiagAUy(i)=-delt/(2.0d0*delz**2)
  End Do
  Do i=1,M-1
     DiagALy(i)=-delt/(2.0d0*delz**2)
  End Do
  ! Samething here, corners are special, due to bc
  DiagAUy(1)=2.0d0*DiagAUy(1)
  DiagALy(M-1)=2.0d0*DiagALy(M-1)

  Return

End Subroutine A_Matrix_x
!-------------------------------------qA------------------------------------------------
! This is the forwad propagator. the qA_0 is a intermidiate step, you can ignore it
Subroutine qA_forward( )
  Use global
  Implicit None

  Integer                  :: errorflag
  Integer                  :: s,i,j
  Integer                  :: errorflagA
  Real*8                   :: gamma,bettaU,bettaL,betta

  qA_0=0.0d0
  qA=0.0d0
  ! Initializing the qs
  Do i=1,N
     Do j=1,M
        qA_0(i,j)=1.0d0
        qA(i,j,0)=1.0d0
     End Do
  End Do
 
  Do s=1,NA

     bAy=0.0d0
     !**********************Scan Over z***************
     Do i=1,N
        
        Call A_Matrix_x(i)  
        ! Depending on the bc, you have to populate the right hand side accordingly.
        If (i==1) then
           Do j=1,M 

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bAy(j)=gamma*qA_0(i,j)+bettaU*qA_0(i+1,j)+bettaL*qA_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bAy(j)=gamma*qA_0(i,j)+bettaU*qA_0(i-1,j)+bettaL*qA_0(i-1,j)

           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bAy(j)=gamma*qA_0(i,j)+bettaU*qA_0(i+1,j)+bettaL*qA_0(i-1,j)

           End Do
        End If
     
        DiagAdy=DiagAy
        DiagAUdy=DiagAUy
        DiagALdy=DiagALy
        ! This is an external library for solving the matrix, using LAPACK
        !Call DGTSV(M,1,DiagALdy,DiagAdy,DiagAUdy,bAy,M,errorflagA)
        ! this is the house code
        Call TDMA(M,DiagALdy,DiagAdy,DiagAUdy,bAy)
             
        Do j=1,M
           qA(i,j,s)=bAy(j)
        End Do

     End Do
     
     bAy=0.0d0
     Do i=1,N
        Do j=1,M
           qA_0(i,j)=qA(i,j,s)
        End Do
     End Do

     bAx=0.0d0
     !**********************Scan Over r***************

     Do j=1,M

        Call A_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bAx(i)=gamma*qA_0(i,j)+betta*qA_0(i,j+1)+betta*qA_0(i,j+1)
        
           End do
        ElseIf (j==M) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bAx(i)=gamma*qA_0(i,j)+betta*qA_0(i,j-1)+betta*qA_0(i,j-1)
      
           End do
        Else
           Do i=1,N   
     
              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bAx(i)=gamma*qA_0(i,j)+betta*qA_0(i,j+1)+betta*qA_0(i,j-1)

           End Do
        End If
     
        
        DiagAdx=DiagAx
        DiagAUdx=DiagAUx
        DiagALdx=DiagALx
        !Call DGTSV(N,1,DiagALdx,DiagAdx,DiagAUdx,bAx,N,errorflagA)
        Call TDMA(N,DiagALdx,DiagAdx,DiagAUdx,bAx)

        Do i=1,N
           qA(i,j,s)=bAx(i)
        End Do
     End Do
    
     Do i=1,N
        Do j=1,M
           qA_0(i,j)=qA(i,j,s)
        End Do
     End Do
     
  End Do
!Do i = 1,N
 !  Do j=1,M
  !    Do s=1,NA
   !      print*,i,' ',j,' ',s,' qA:',qA(i,j,s)
    !  End Do
  ! End Do
!End Do


  Return

End Subroutine qA_forward
!-------------------------------------qAD------------------------------------------------
! This is for the complementary prpagator, qAD(agor)
Subroutine qAD_forward( )
  Use global
  Implicit None

  Integer                  :: s,i,j
  Integer                  :: errorflagA
  Real*8                   :: gamma,bettaU,bettaL,betta

  ! Same as before, but the initial condition is from qB
  qAD_0=0.0d0
  qAD=0.0d0
  
  Do i=1,N
     Do j=1,M
        qAD(i,j,0)=qBD(i,j,NB)
        qAD_0(i,j)=qBD(i,j,NB)
     End Do
  End Do
    
  Do s=1,NA
     
     bAy=0.0d0
     !**********************Scan Over z***************
     Do i=1,N
        Call A_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bAy(j)=gamma*qAD_0(i,j)+bettaU*qAD_0(i+1,j)+bettaL*qAD_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))  
    
              bAy(j)=gamma*qAD_0(i,j)+bettaU*qAD_0(i-1,j)+bettaL*qAD_0(i-1,j)
          
           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wA(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bAy(j)=gamma*qAD_0(i,j)+bettaU*qAD_0(i+1,j)+bettaL*qAD_0(i-1,j)

           End Do
        End If
        
        DiagAdy=DiagAy
        DiagAUdy=DiagAUy
        DiagALdy=DiagALy
        !Call DGTSV(M,1,DiagALdy,DiagAdy,DiagAUdy,bAy,M,errorflagA)
        Call TDMA(M,DiagALdy,DiagAdy,DiagAUdy,bAy)

        Do j=1,M
           qAD(i,j,s)=bAy(j)
        End Do
     End Do

     bAy=0.0d0
     Do i=1,N
        Do j=1,M
           qAD_0(i,j)=qAD(i,j,s)
        End Do
     End Do

     bAx=0.0d0
     !**********************Scan Over r***************
     Do j=1,M
        Call A_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bAx(i)=gamma*qAD_0(i,j)+betta*qAD_0(i,j+1)+betta*qAD_0(i,j+1)

           End Do
        ElseIf (j==M) then
           Do i=1,N   

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bAx(i)=gamma*qAD_0(i,j)+betta*qAD_0(i,j-1)+betta*qAD_0(i,j-1)

           End Do
        Else
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bAx(i)=gamma*qAD_0(i,j)+betta*qAD_0(i,j+1)+betta*qAD_0(i,j-1)

           End Do
        End If
        

        DiagAdx=DiagAx
        DiagAUdx=DiagAUx
        DiagALdx=DiagALx
        !Call DGTSV(N,1,DiagALdx,DiagAdx,DiagAUdx,bAx,N,errorflagA)      
        Call TDMA(N,DiagALdx,DiagAdx,DiagAUdx,bAx)

        Do i=1,N
           qAD(i,j,s)=bAx(i)
        End Do

     End Do
     
     Do i=1,N
        Do j=1,M
           qAD_0(i,j)=qAD(i,j,s)
        End Do
     End Do
     

  End Do
!Do i=1,N
 !   Do j=1,M
  !      Do s=1,NA
   !         print*,i,j,s,qAD(i,j,s)
    !    End Do
   ! End Do
!End Do

  Return

End Subroutine qAD_forward
