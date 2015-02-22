Subroutine D_Matrix_y(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaD,betaDL,betaDU

  
  DiagDx=0.0d0
  DiagDUx=0.0d0
  DiagDLx=0.0d0

  Do i=1,N
     DiagDx(i)=1.0d0+(delt/delr**2)+((delt/2.0d0)*wD(i,ii))
  End Do

  Do i=2,N-1
     DiagDUx(i)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(i-1))+r_0)*4.0d0*delr))
  End Do

  Do i=1,N-2
     DiagDLx(i)=-(delt/(2.0d0*delr**2))+(delt/(((delr*real(i+1-1))+r_0)*4.0d0*delr))
  End Do

  DiagDUx(1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(1-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(1-1))+r_0)*4.0d0*delr))

  DiagDLx(N-1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(N-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(N-1))+r_0)*4.0d0*delr))


  Return

End Subroutine D_Matrix_y
!--------------------------------------------------------------------------------
Subroutine D_Matrix_x(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaD,betaD

  DiagDy=0.0d0
  DiagDUy=0.0d0
  DiagDLy=0.0d0

  Do i=1,M
     DiagDy(i)=1.0d0+delt/delz**2
  End Do
  Do i=1,M-1
     DiagDUy(i)=-delt/(2.0d0*delz**2)
  End Do
  Do i=1,M-1
     DiagDLy(i)=-delt/(2.0d0*delz**2)
  End Do

  DiagDUy(1)=2.0d0*DiagDUy(1)
  DiagDLy(M-1)=2.0d0*DiagDLy(M-1)

  Return

End Subroutine D_Matrix_x
!-------------------------------------qD------------------------------------------------
Subroutine qD_forward( )
  Use global
  Implicit None

  Integer                  :: errorflag
  Integer                  :: s,i,j
  Integer                  :: errorflagD
  Real*8                   :: gamma,bettaU,bettaL,betta

  qD_0=0.0d0
  qD=0.0d0

  Do i=1,N
     Do j=1,M
        qD_0(i,j)=qB(i,j,NB)
        qD(i,j,0)=qB(i,j,NB)
     End Do
  End Do
 
  Do s=1,ND

     bDy=0.0d0
     !**********************Scan Over x***************
     Do i=1,N
        
        Call D_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M 

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bDy(j)=gamma*qD_0(i,j)+bettaU*qD_0(i+1,j)+bettaL*qD_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bDy(j)=gamma*qD_0(i,j)+bettaU*qD_0(i-1,j)+bettaL*qD_0(i-1,j)

           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bDy(j)=gamma*qD_0(i,j)+bettaU*qD_0(i+1,j)+bettaL*qD_0(i-1,j)

           End Do
        End If
     
        DiagDdy=DiagDy
        DiagDUdy=DiagDUy
        DiagDLdy=DiagDLy
        !Call DGTSV(M,1,DiagDLdy,DiagDdy,DiagDUdy,bDy,M,errorflagD)
        Call TDMA(M,DiagDLdy,DiagDdy,DiagDUdy,bDy)
      
        Do j=1,M
           qD(i,j,s)=bDy(j)
        End Do

     End Do
     
     bDy=0.0d0
     Do i=1,N
        Do j=1,M
           qD_0(i,j)=qD(i,j,s)
        End Do
     End Do

     bDx=0.0d0
     !**********************Scan Over y***************

     Do j=1,M

        Call D_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bDx(i)=gamma*qD_0(i,j)+betta*qD_0(i,j+1)+betta*qD_0(i,j+1)
        
           End do
        ElseIf (j==M) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bDx(i)=gamma*qD_0(i,j)+betta*qD_0(i,j-1)+betta*qD_0(i,j-1)
      
           End do
        Else
           Do i=1,N   
     
              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bDx(i)=gamma*qD_0(i,j)+betta*qD_0(i,j+1)+betta*qD_0(i,j-1)

           End Do
        End If
     
        
        DiagDdx=DiagDx
        DiagDUdx=DiagDUx
        DiagDLdx=DiagDLx
        !Call DGTSV(N,1,DiagDLdx,DiagDdx,DiagDUdx,bDx,N,errorflagD)
        Call TDMA(N,DiagDLdx,DiagDdx,DiagDUdx,bDx)

        Do i=1,N
           qD(i,j,s)=bDx(i)
        End Do
     End Do
    
     Do i=1,N
        Do j=1,M
           qD_0(i,j)=qD(i,j,s)
        End Do
     End Do
     
  End Do
  
  Return

End Subroutine qD_forward
!-------------------------------------qDD------------------------------------------------
Subroutine qDD_forward( )
  Use global
  Implicit None

  Integer                  :: s,i,j
  Integer                  :: errorflagD
  Real*8                   :: gamma,bettaU,bettaL,betta



  qDD_0=0.0d0
  qDD=0.0d0
  
  Do i=1,N
     Do j=1,M
        qDD(i,j,0)=qED(i,j,NE)
        qDD_0(i,j)=qED(i,j,NE)
     End Do
  End Do
    
  Do s=1,ND
     
     bDy=0.0d0
     !**********************Scan Over x***************
     Do i=1,N
        Call D_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bDy(j)=gamma*qDD_0(i,j)+bettaU*qDD_0(i+1,j)+bettaL*qDD_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))  
    
              bDy(j)=gamma*qDD_0(i,j)+bettaU*qDD_0(i-1,j)+bettaL*qDD_0(i-1,j)
          
           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wD(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))

              bDy(j)=gamma*qDD_0(i,j)+bettaU*qDD_0(i+1,j)+bettaL*qDD_0(i-1,j)

           End Do
        End If
        
        DiagDdy=DiagDy
        DiagDUdy=DiagDUy
        DiagDLdy=DiagDLy
        !Call DGTSV(M,1,DiagDLdy,DiagDdy,DiagDUdy,bDy,M,errorflagD)
        Call TDMA(M,DiagDLdy,DiagDdy,DiagDUdy,bDy)

        Do j=1,M
           qDD(i,j,s)=bDy(j)
        End Do
     End Do

     bDy=0.0d0
     Do i=1,N
        Do j=1,M
           qDD_0(i,j)=qDD(i,j,s)
        End Do
     End Do

     bDx=0.0d0
     !**********************Scan Over y***************
     Do j=1,M
        Call D_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bDx(i)=gamma*qDD_0(i,j)+betta*qDD_0(i,j+1)+betta*qDD_0(i,j+1)

           End Do
        ElseIf (j==M) then
           Do i=1,N   

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bDx(i)=gamma*qDD_0(i,j)+betta*qDD_0(i,j-1)+betta*qDD_0(i,j-1)

           End Do
        Else
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)

              bDx(i)=gamma*qDD_0(i,j)+betta*qDD_0(i,j+1)+betta*qDD_0(i,j-1)

           End Do
        End If
        

        DiagDdx=DiagDx
        DiagDUdx=DiagDUx
        DiagDLdx=DiagDLx
        !Call DGTSV(N,1,DiagDLdx,DiagDdx,DiagDUdx,bDx,N,errorflagD)     
        Call TDMA(N,DiagDLdx,DiagDdx,DiagDUdx,bDx)     
    
        Do i=1,N
           qDD(i,j,s)=bDx(i)
        End Do

     End Do
     
     Do i=1,N
        Do j=1,M
           qDD_0(i,j)=qDD(i,j,s)
        End Do
     End Do
     

  End Do
  
  Return

End Subroutine qDD_forward
