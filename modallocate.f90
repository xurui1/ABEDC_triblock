Subroutine Allocate_mod(order)
  Use global
  Implicit None

  Integer        :: order


  If (order==1) then
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !_________________________________Allocation_A________________________________________________________
     Allocate(qA(1:N,1:M,0:NA),qA_0(0:N+1,0:M+1),qAD(1:N,1:M,0:NA),qAD_0(0:N+1,0:M+1))
     Allocate(wA(1:N,1:M),dwA(1:N,1:M),pA(1:N,1:M),bAx(1:N),bAy(1:M))
     Allocate(DiagAx(1:N),DiagAUx(1:N-1),DiagALx(1:N-1))
     Allocate(DiagAdx(1:N),DiagAUdx(1:N-1),DiagALdx(1:N-1))
     Allocate(DiagAy(1:M),DiagAUy(1:M-1),DiagALy(1:M-1))
     Allocate(DiagAdy(1:M),DiagAUdy(1:M-1),DiagALdy(1:M-1))
     !_________________________________Allocation_B_________________________________________________________
     Allocate(qB(1:N,1:M,0:NB),qB_0(0:N+1,0:M+1),qBD(1:N,1:M,0:NB),qBD_0(0:N+1,0:M+1))
     Allocate(wB(1:N,1:M),dwB(1:N,1:M),pB(1:N,1:M),bBx(1:N),bBy(1:M))
     Allocate(DiagBx(1:N),DiagBUx(1:N-1),DiagBLx(1:N-1))
     Allocate(DiagBdx(1:N),DiagBUdx(1:N-1),DiagBLdx(1:N-1))
     Allocate(DiagBy(1:M),DiagBUy(1:M-1),DiagBLy(1:M-1))
     Allocate(DiagBdy(1:M),DiagBUdy(1:M-1),DiagBLdy(1:M-1))
     !_________________________________Allocation_C_________________________________________________________
     Allocate(qC(1:N,1:M,0:NC),qC_0(0:N+1,0:M+1))
     Allocate(wC(1:N,1:M),dwC(1:N,1:M),pC(1:N,1:M),bCx(1:N),bCy(1:M))
     Allocate(DiagCx(1:N),DiagCUx(1:N-1),DiagCLx(1:N-1))
     Allocate(DiagCdx(1:N),DiagCUdx(1:N-1),DiagCLdx(1:N-1))
     Allocate(DiagCy(1:M),DiagCUy(1:M-1),DiagCLy(1:M-1))
     Allocate(DiagCdy(1:M),DiagCUdy(1:M-1),DiagCLdy(1:M-1))
     !_________________________________Allocation_E________________________________________________________
     Allocate(qE(1:N,1:M,0:NE),qE_0(0:N+1,0:M+1),qED(1:N,1:M,0:NE),qED_0(0:N+1,0:M+1))
     Allocate(wE(1:N,1:M),dwE(1:N,1:M),pE(1:N,1:M),bEx(1:N),bEy(1:M))
     Allocate(DiagEx(1:N),DiagEUx(1:N-1),DiagELx(1:N-1))
     Allocate(DiagEdx(1:N),DiagEUdx(1:N-1),DiagELdx(1:N-1))
     Allocate(DiagEy(1:M),DiagEUy(1:M-1),DiagELy(1:M-1))
     Allocate(DiagEdy(1:M),DiagEUdy(1:M-1),DiagELdy(1:M-1))
     !_________________________________Allocation_D_________________________________________________________
     Allocate(qD(1:N,1:M,0:ND),qD_0(0:N+1,0:M+1),qDD(1:N,1:M,0:ND),qDD_0(0:N+1,0:M+1))
     Allocate(wD(1:N,1:M),dwD(1:N,1:M),pD(1:N,1:M),bDx(1:N),bDy(1:M))
     Allocate(DiagDx(1:N),DiagDUx(1:N-1),DiagDLx(1:N-1))
     Allocate(DiagDdx(1:N),DiagDUdx(1:N-1),DiagDLdx(1:N-1))
     Allocate(DiagDy(1:M),DiagDUy(1:M-1),DiagDLy(1:M-1))
     Allocate(DiagDdy(1:M),DiagDUdy(1:M-1),DiagDLdy(1:M-1))



     !_________________________________Allocation_Rest_______________________________________________________
     Allocate(dpp(1:N,1:M),eta(1:N,1:M),eta2(1:N,1:M))
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ElseIf (order==0) then
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DeAllocation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !_________________________________DeAllocation_A______________________________________________________
     DeAllocate(qA,qA_0,qAD,qAD_0)
     DeAllocate(wA,dwA,pA,bAx,bAy)
     DeAllocate(DiagAx,DiagAUx,DiagALx)
     DeAllocate(DiagAdx,DiagAUdx,DiagALdx)
     DeAllocate(DiagAy,DiagAUy,DiagALy)
     DeAllocate(DiagAdy,DiagAUdy,DiagALdy)
     !_________________________________DeAllocation_B_______________________________________________________
     DeAllocate(qB,qB_0,qBD,qBD_0)
     DeAllocate(wB,dwB,pB,bBx,bBy)
     DeAllocate(DiagBx,DiagBUx,DiagBLx)
     DeAllocate(DiagBdx,DiagBUdx,DiagBLdx)
     DeAllocate(DiagBy,DiagBUy,DiagBLy)
     DeAllocate(DiagBdy,DiagBUdy,DiagBLdy)
     !_________________________________DeAllocation_C_______________________________________________________
     DeAllocate(qC,qC_0)
     DeAllocate(wC,dwC,pC,bCx,bCy)
     DeAllocate(DiagCx,DiagCUx,DiagCLx)
     DeAllocate(DiagCdx,DiagCUdx,DiagCLdx)
     DeAllocate(DiagCy,DiagCUy,DiagCLy)
     DeAllocate(DiagCdy,DiagCUdy,DiagCLdy)
     !_________________________________DeAllocation_E______________________________________________________
     DeAllocate(qE,qE_0,qED,qED_0)
     DeAllocate(wE,dwE,pE,bEx,bEy)
     DeAllocate(DiagEx,DiagEUx,DiagELx)
     DeAllocate(DiagEdx,DiagEUdx,DiagELdx)
     DeAllocate(DiagEy,DiagEUy,DiagELy)
     DeAllocate(DiagEdy,DiagEUdy,DiagELdy)
     !_________________________________DeAllocation_D_______________________________________________________
     DeAllocate(qD,qD_0,qDD,qDD_0)
     DeAllocate(wD,dwD,pD,bDx,bDy)
     DeAllocate(DiagDx,DiagDUx,DiagDLx)
     DeAllocate(DiagDdx,DiagDUdx,DiagDLdx)
     DeAllocate(DiagDy,DiagDUy,DiagDLy)
     DeAllocate(DiagDdy,DiagDUdy,DiagDLdy)
     !_________________________________DeAllocation_Rest____________________________________________________
     DeAllocate(dpp,eta,eta2)
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Else
     Print*,"Something has gone wrong in the allocation module"
  End If

  Return
End Subroutine Allocate_mod
