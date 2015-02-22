Subroutine B_Matrix_y(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaB,betaBL,betaBU

  
  DiagBx=0.0d0
  DiagBUx=0.0d0
  DiagBLx=0.0d0

  Do i=1,N
     DiagBx(i)=1.0d0+(delt/delr**2)+((delt/2.0d0)*wB(i,ii))
  End Do

  Do i=2,N-1
     DiagBUx(i)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(i-1))+r_0)*4.0d0*delr))
  End Do

  Do i=1,N-2
     DiagBLx(i)=-(delt/(2.0d0*delr**2))+(delt/(((delr*real(i+1-1))+r_0)*4.0d0*delr))
  End Do

  DiagBUx(1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(1-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(1-1))+r_0)*4.0d0*delr))

  DiagBLx(N-1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(N-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(N-1))+r_0)*4.0d0*delr))


  Return

End Subroutine B_Matrix_y
!--------------------------------------------------------------------------------
Subroutine B_Matrix_x(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaB,betaB

  DiagBy=0.0d0
  DiagBUy=0.0d0
  DiagBLy=0.0d0

  Do i=1,M
     DiagBy(i)=1.0d0+delt/delz**2
  End Do
  Do i=1,M-1
     DiagBUy(i)=-delt/(2.0d0*delz**2)
  End Do
  Do i=1,M-1
     DiagBLy(i)=-delt/(2.0d0*delz**2)
  End Do

  DiagBUy(1)=2.0d0*DiagBUy(1)
  DiagBLy(M-1)=2.0d0*DiagBLy(M-1)

  Return

End Subroutine B_Matrix_x
!-------------------------------------qB------------------------------------------------
Subroutine qB_forward( )
  Use global
  Implicit None

  Integer                  :: errorflag
  Integer                  :: s,i,j
  Integer                  :: errorflagB
  Real*8                   :: gamma,bettaU,bettaL,betta

  qB_0=0.0d0
  qB=0.0d0

  Do i=1,N
     Do j=1,M
        qB_0(i,j)=qA(i,j,NA)
        qB(i,j,0)=qA(i,j,NA)
     End Do
  End Do
 
  Do s=1,NB

     bBy=0.0d0
     !**********************Scan Over x***************
     Do i=1,N
        
        Call B_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M 

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bBy(j)=gamma*qB_0(i,j)+bettaU*qB_0(i+1,j)+bettaL*qB_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bBy(j)=gamma*qB_0(i,j)+bettaU*qB_0(i-1,j)+bettaL*qB_0(i-1,j)

           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bBy(j)=gamma*qB_0(i,j)+bettaU*qB_0(i+1,j)+bettaL*qB_0(i-1,j)

           End Do
        End If
     
        DiagBdy=DiagBy
        DiagBUdy=DiagBUy
        DiagBLdy=DiagBLy
        !Call DGTSV(M,1,DiagBLdy,DiagBdy,DiagBUdy,bBy,M,errorflagB)
        Call TDMA(M,DiagBLdy,DiagBdy,DiagBUdy,bBy)

        Do j=1,M
           qB(i,j,s)=bBy(j)
        End Do

     End Do
     
     bBy=0.0d0
     Do i=1,N
        Do j=1,M
           qB_0(i,j)=qB(i,j,s)
        End Do
     End Do

     bBx=0.0d0
     !**********************Scan Over y***************

     Do j=1,M

        Call B_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bBx(i)=gamma*qB_0(i,j)+betta*qB_0(i,j+1)+betta*qB_0(i,j+1)
        
           End do
        ElseIf (j==M) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bBx(i)=gamma*qB_0(i,j)+betta*qB_0(i,j-1)+betta*qB_0(i,j-1)
      
           End do
        Else
           Do i=1,N   
     
              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bBx(i)=gamma*qB_0(i,j)+betta*qB_0(i,j+1)+betta*qB_0(i,j-1)

           End Do
        End If
     
        
        DiagBdx=DiagBx
        DiagBUdx=DiagBUx
        DiagBLdx=DiagBLx
        !Call DGTSV(N,1,DiagBLdx,DiagBdx,DiagBUdx,bBx,N,errorflagB)
        Call TDMA(N,DiagBLdx,DiagBdx,DiagBUdx,bBx)

        Do i=1,N
           qB(i,j,s)=bBx(i)
        End Do
     End Do
    
     Do i=1,N
        Do j=1,M
           qB_0(i,j)=qB(i,j,s)
        End Do
     End Do
     
  End Do
  
  Return

End Subroutine qB_forward
!-------------------------------------qBD------------------------------------------------
Subroutine qBD_forward( )
  Use global
  Implicit None

  Integer                  :: s,i,j
  Integer                  :: errorflagB
  Real*8                   :: gamma,bettaU,bettaL,betta



  qBD_0=0.0d0
  qBD=0.0d0
  
  Do i=1,N
     Do j=1,M
        qBD(i,j,0)=qDD(i,j,ND)
        qBD_0(i,j)=qDD(i,j,ND)
     End Do
  End Do
    
  Do s=1,NB
     
     bBy=0.0d0
     !**********************Scan Over x***************
     Do i=1,N
        Call B_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bBy(j)=gamma*qBD_0(i,j)+bettaU*qBD_0(i+1,j)+bettaL*qBD_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))  
    
              bBy(j)=gamma*qBD_0(i,j)+bettaU*qBD_0(i-1,j)+bettaL*qBD_0(i-1,j)
          
           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wB(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bBy(j)=gamma*qBD_0(i,j)+bettaU*qBD_0(i+1,j)+bettaL*qBD_0(i-1,j)

           End Do
        End If
        
        DiagBdy=DiagBy
        DiagBUdy=DiagBUy
        DiagBLdy=DiagBLy
        !Call DGTSV(M,1,DiagBLdy,DiagBdy,DiagBUdy,bBy,M,errorflagB)
        Call TDMA(M,DiagBLdy,DiagBdy,DiagBUdy,bBy)

        Do j=1,M
           qBD(i,j,s)=bBy(j)
        End Do
     End Do

     bBy=0.0d0
     Do i=1,N
        Do j=1,M
           qBD_0(i,j)=qBD(i,j,s)
        End Do
     End Do

     bBx=0.0d0
     !**********************Scan Over y***************
     Do j=1,M
        Call B_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bBx(i)=gamma*qBD_0(i,j)+betta*qBD_0(i,j+1)+betta*qBD_0(i,j+1)

           End Do
        ElseIf (j==M) then
           Do i=1,N   

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bBx(i)=gamma*qBD_0(i,j)+betta*qBD_0(i,j-1)+betta*qBD_0(i,j-1)

           End Do
        Else
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bBx(i)=gamma*qBD_0(i,j)+betta*qBD_0(i,j+1)+betta*qBD_0(i,j-1)

           End Do
        End If
        

        DiagBdx=DiagBx
        DiagBUdx=DiagBUx
        DiagBLdx=DiagBLx
        !Call DGTSV(N,1,DiagBLdx,DiagBdx,DiagBUdx,bBx,N,errorflagB)     
        Call TDMA(N,DiagBLdx,DiagBdx,DiagBUdx,bBx)     
    
        Do i=1,N
           qBD(i,j,s)=bBx(i)
        End Do

     End Do
     
     Do i=1,N
        Do j=1,M
           qBD_0(i,j)=qBD(i,j,s)
        End Do
     End Do
     

  End Do
  
  Return

End Subroutine qBD_forward
