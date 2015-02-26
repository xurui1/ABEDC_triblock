Subroutine fE_homo()
  Use global
  Implicit None

  Integer               :: i,j
  Real*8                :: pA_ave,pB_ave,pC_ave,pE_ave,pD_ave
  Real*8                :: wA_ave,wB_ave,wC_ave,wE_ave,wD_ave
  Real*8                :: dwA_ave,dWB_ave,dwC_ave,dwE_ave,dwD_ave,dpp_ave
  Real*8                :: eta_ave
  Real*8                :: f_int, f_omeg
  Real*8                :: p_vect(5),w_vect(5)

  f_int=0.0
  f_omeg=0.0

  dwA_ave=0.0
  dwB_ave=0.0
  dwC_ave=0.0
  dwE_ave=0.0
  dwD_ave=0.0

  eta_ave=0.0

  pA_ave=0.002
  pB_ave=pA_ave
  pE_ave=pA_ave
  pD_ave=pA_ave
  pC_ave=1.0-(pA_ave+pB_ave+pE_ave+pD_ave)

  wA_ave=xAB*pB_ave+xAC*pC_ave+xAE*pE_ave+xAD*pD_ave+eta_ave
  wB_ave=xAB*pA_ave+xBC*pC_ave+xBE*pE_ave+xBD*pD_ave+eta_ave
  wC_ave=xAC*pA_ave+xBC*pB_ave+xCE*pE_ave+xCD*pD_ave+eta_ave
  wE_ave=xAE*pA_ave+xBE*pB_ave+xCE*pC_ave+xED*pD_ave+eta_ave
  wD_ave=xAD*pA_ave+xBD*pB_ave+xCD*pC_ave+xED*pE_ave+eta_ave

  Do i=1,10000000
     
     eta_ave=eta_ave-0.05*(1.0-(pA_ave+pB_ave+pC_ave+pE_ave+pD_ave))

     pA_ave=exp(muABDE-wA_ave*fracA-wB_ave*fracB-wE_ave*fracE &
                -wD_ave*fracD)*(fracA/2.0)
     pB_ave=exp(muABDE-wA_ave*fracA-wB_ave*fracB-wE_ave*fracE &
                -wD_ave*fracD)*(fracB/2.0)
     pC_ave=exp(kappaC*(muC-wC_ave))
     pE_ave=exp(muABDE-wA_ave*fracA-wB_ave*fracB-wE_ave*fracE &
                -wD_ave*fracD)*(fracE/2.0)
     pD_ave=exp(muABDE-wA_ave*fracA-wB_ave*fracB-wE_ave*fracE &
                -wD_ave*fracD)*(fracD/2.0)

     dwA_ave=(xAB*pB_ave+xAC*pC_ave+xAE*pE_ave+xAD*pD_ave+eta_ave)-wA_ave
     dwB_ave=(xAB*pA_ave+xBC*pC_ave+xBE*pE_ave+xBD*pD_ave+eta_ave)-wB_ave
     dwC_ave=(xAC*pA_ave+xBC*pB_ave+xCE*pE_ave+xCD*pD_ave+eta_ave)-wC_ave
     dwE_ave=(xAE*pA_ave+xBE*pB_ave+xCE*pC_ave+xED*pD_ave+eta_ave)-wE_ave
     dwD_ave=(xAD*pA_ave+xBD*pB_ave+xCD*pC_ave+xED*pE_ave+eta_ave)-wD_ave

     dpp_ave=1.0-(pA_ave+pB_ave+pC_ave+pE_ave+pD_ave)
     
     wA_ave=wA_ave+0.005*dwA_ave
     wB_ave=wB_ave+0.005*dwB_ave
     wC_ave=wC_ave+0.005*dwC_ave
     wE_ave=wE_ave+0.005*dwE_ave
     wD_ave=wD_ave+0.005*dwD_ave

  End Do

  phiAB_hom=pA_ave+pB_ave+pE_ave+pD_ave
  phiC_hom=pC_ave

 
  p_vect(1)=pA_ave
  p_vect(2)=pB_ave
  p_vect(3)=pC_ave
  p_vect(4)=pE_ave
  p_vect(5)=pD_ave

  w_vect(1)=wA_ave
  w_vect(2)=wB_ave
  w_vect(3)=wC_ave
  w_vect(4)=wE_ave
  w_vect(5)=wD_ave

  Do i=1,5
     Do j=i,5
        f_int=f_int+p_vect(i)*p_vect(j)*XM(i,j)
     End Do
     f_omeg=f_omeg+p_vect(i)*w_vect(i)
  End Do
  

  fE_hom=f_int-f_omeg &
       -(exp(muABDE-wA_ave*fracA-wB_ave*fracB-wE_ave*fracE-wD_ave*fracD)) &
       -(exp(kappaC*(muC-wC_ave))/kappaC)


  Return

End Subroutine fE_homo
