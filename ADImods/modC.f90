Subroutine C_Matrix_y(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaC,betaCL,betaCU

  
  DiagCx=0.0d0
  DiagCUx=0.0d0
  DiagCLx=0.0d0

  Do i=1,N
     DiagCx(i)=1.0d0+(delt/delr**2)+((delt/2.0d0)*wC(i,ii))
  End Do

  Do i=2,N-1
     DiagCUx(i)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(i-1))+r_0)*4.0d0*delr))
  End Do

  Do i=1,N-2
     DiagCLx(i)=-(delt/(2.0d0*delr**2))+(delt/(((delr*real(i+1-1))+r_0)*4.0d0*delr))
  End Do

  DiagCUx(1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(1-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(1-1))+r_0)*4.0d0*delr))

  DiagCLx(N-1)=-(delt/(2.0d0*delr**2))-(delt/(((delr*real(N-1))+r_0)*4.0d0*delr)) &
       -(delt/(2.0d0*delr**2))+(delt/(((delr*real(N-1))+r_0)*4.0d0*delr))


  Return

End Subroutine C_Matrix_y
!--------------------------------------------------------------------------------
Subroutine C_Matrix_x(ii)
  Use global
  Implicit None

  Integer                    :: i,j,ii
  Real*8                     :: alphaC,betaC

  DiagCy=0.0d0
  DiagCUy=0.0d0
  DiagCLy=0.0d0

  Do i=1,M
     DiagCy(i)=1.0d0+delt/delz**2
  End Do
  Do i=1,M-1
     DiagCUy(i)=-delt/(2.0d0*delz**2)
  End Do
  Do i=1,M-1
     DiagCLy(i)=-delt/(2.0d0*delz**2)
  End Do

  DiagCUy(1)=2.0d0*DiagCUy(1)
  DiagCLy(M-1)=2.0d0*DiagCLy(M-1)

  Return

End Subroutine C_Matrix_x
!-------------------------------------qC------------------------------------------------
Subroutine qC_forward( )
  Use global
  Implicit None

  Integer                  :: errorflag
  Integer                  :: s,i,j
  Integer                  :: errorflagC
  Real*8                   :: gamma,bettaU,bettaL,betta

  qC_0=0.0d0
  qC=0.0d0

  Do i=1,N
     Do j=1,M
        qC_0(i,j)=1.0
        qC(i,j,0)=1.0
     End Do
  End Do
 
  Do s=1,NC

     bCy=0.0d0
     !**********************Scan Over x***************
     Do i=1,N
        
        Call C_Matrix_x(i)  
        
        If (i==1) then
           Do j=1,M 

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wC(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bCy(j)=gamma*qC_0(i,j)+bettaU*qC_0(i+1,j)+bettaL*qC_0(i+1,j)

           End Do
        ElseIf (i==N) then
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wC(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bCy(j)=gamma*qC_0(i,j)+bettaU*qC_0(i-1,j)+bettaL*qC_0(i-1,j)

           End Do
        Else
           Do j=1,M

              gamma=1.0d0-(delt/delr**2)-((delt/2.0d0)*wC(i,j))
              bettaL=(delt/(2.0d0*delr**2))-(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              bettaU=(delt/(2.0d0*delr**2))+(delt/(((real(i-1)*delr)+(r_0))*4.0d0*delr))
              
              bCy(j)=gamma*qC_0(i,j)+bettaU*qC_0(i+1,j)+bettaL*qC_0(i-1,j)

           End Do
        End If
     
        DiagCdy=DiagCy
        DiagCUdy=DiagCUy
        DiagCLdy=DiagCLy
        !Call DGTSV(M,1,DiagCLdy,DiagCdy,DiagCUdy,bCy,M,errorflagC)
        Call TDMA(M,DiagCLdy,DiagCdy,DiagCUdy,bCy)
      
        Do j=1,M
           qC(i,j,s)=bCy(j)
        End Do

     End Do
     
     bCy=0.0d0
     Do i=1,N
        Do j=1,M
           qC_0(i,j)=qC(i,j,s)
        End Do
     End Do

     bCx=0.0d0
     !**********************Scan Over y***************

     Do j=1,M

        Call C_Matrix_y(j)  
        
        If (j==1) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bCx(i)=gamma*qC_0(i,j)+betta*qC_0(i,j+1)+betta*qC_0(i,j+1)
        
           End do
        ElseIf (j==M) then
           Do i=1,N

              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bCx(i)=gamma*qC_0(i,j)+betta*qC_0(i,j-1)+betta*qC_0(i,j-1)
      
           End do
        Else
           Do i=1,N   
     
              gamma=1.0d0-delt/delz**2
              betta=delt/(2.0d0*delz**2)
              bCx(i)=gamma*qC_0(i,j)+betta*qC_0(i,j+1)+betta*qC_0(i,j-1)

           End Do
        End If
     
        
        DiagCdx=DiagCx
        DiagCUdx=DiagCUx
        DiagCLdx=DiagCLx
        !Call DGTSV(N,1,DiagCLdx,DiagCdx,DiagCUdx,bCx,N,errorflagC)
        Call TDMA(N,DiagCLdx,DiagCdx,DiagCUdx,bCx)

        Do i=1,N
           qC(i,j,s)=bCx(i)
        End Do
     End Do
    
     Do i=1,N
        Do j=1,M
           qC_0(i,j)=qC(i,j,s)
        End Do
     End Do
     
  End Do
  
  Return

End Subroutine qC_forward
