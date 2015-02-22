!-------------------------------------------------------------------------------------
Subroutine totalphi(ptot_in,phiA,phiB,phiC,phiE,phiD)
  Use global
  Implicit None

  Integer            :: dummy
  Integer            :: i,j
  Real*8             :: ptot_in
  Real*8             :: phiA,phiB,phiC,phiE,phiD

  phiC=0.0d0
  phiA=0.0d0
  phiB=0.0d0 
  phiE=0.0d0
  phiD=	0.0d0
  ptot_in=0.0
  !-----------------------------------------------------------------
  
  Do i=2,N-1
     Do j=2,N-1
        phiA=phiA+pA(i,j)*delr*(real(i-1)*delr+r_0)*delz
        phiB=phiB+pB(i,j)*delr*(real(i-1)*delr+r_0)*delz
        phiC=phiC+pC(i,j)*delr*(real(i-1)*delr+r_0)*delz
        phiE=phiE+pE(i,j)*delr*(real(i-1)*delr+r_0)*delz
        phiD=phiD+pD(i,j)*delr*(real(i-1)*delr+r_0)*delz
     End Do
  End Do


  i=1
  Do j=2,N-1
     phiA=phiA+0.5d0*pA(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiB=phiB+0.5d0*pB(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiC=phiC+0.5d0*pC(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiE=phiE+0.5d0*pE(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiD=phiD+0.5d0*pD(i,j)*delr*(real(i-1)*delr+r_0)*delz
  End Do
  
  i=N
  Do j=2,N-1
     phiA=phiA+0.5d0*pA(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiB=phiB+0.5d0*pB(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiC=phiC+0.5d0*pC(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiE=phiE+0.5d0*pE(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiD=phiD+0.5d0*pD(i,j)*delr*(real(i-1)*delr+r_0)*delz
  End Do
  
  j=1
  Do i=2,N-1
     phiA=phiA+0.5d0*pA(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiB=phiB+0.5d0*pB(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiC=phiC+0.5d0*pC(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiE=phiE+0.5d0*pE(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiD=phiD+0.5d0*pD(i,j)*delr*(real(i-1)*delr+r_0)*delz
  End Do
 
  j=N
  Do i=2,N-1
     phiA=phiA+0.5d0*pA(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiB=phiB+0.5d0*pB(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiC=phiC+0.5d0*pC(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiE=phiE+0.5d0*pE(i,j)*delr*(real(i-1)*delr+r_0)*delz
     phiD=phiD+0.5d0*pD(i,j)*delr*(real(i-1)*delr+r_0)*delz
  End Do

  phiA=phiA+0.25d0*delr*delz*((pA(1,1)*(real(1-1)*delr+r_0))+(pA(1,N)*(real(1-1)*delr+r_0)))
  phiA=phiA+0.25d0*delr*delz*((pA(N,1)*(real(N-1)*delr+r_0))+(pA(N,N)*(real(N-1)*delr+r_0)))

  phiB=phiB+0.25d0*delr*delz*((pB(1,1)*(real(1-1)*delr+r_0))+(pB(1,N)*(real(1-1)*delr+r_0)) &
       +(pB(N,1)*(real(N-1)*delr+r_0))+(pB(N,N)*(real(N-1)*delr+r_0)))

  phiC=phiC+0.25d0*delr*delz*((pC(1,1)*(real(1-1)*delr+r_0))+(pC(1,N)*(real(1-1)*delr+r_0)) &
       +(pC(N,1)*(real(N-1)*delr+r_0))+(pC(N,N)*(real(N-1)*delr+r_0)))


  phiE=phiE+0.25d0*delr*delz*((pE(1,1)*(real(1-1)*delr+r_0))+(pE(1,N)*(real(1-1)*delr+r_0)) &
       +(pE(N,1)*(real(N-1)*delr+r_0))+(pE(N,N)*(real(N-1)*delr+r_0)))

  phiD=phiD+0.25d0*delr*delz*((pD(1,1)*(real(1-1)*delr+r_0))+(pD(1,N)*(real(1-1)*delr+r_0)) &
       +(pD(N,1)*(real(N-1)*delr+r_0))+(pD(N,N)*(real(N-1)*delr+r_0)))
  
  !------------------------------------------------------------------
  phiA=2.0*pi*phiA/vol
  phiB=2.0*pi*phiB/vol
  phiC=2.0*pi*phiC/vol
  phiE=2.0*pi*phiE/vol
  phiD=2.0*pi*phiD/vol

  ptot_in=phiA+phiB+phiC+phiE+phiD


  Return

End Subroutine totalphi
