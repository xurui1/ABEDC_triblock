Module global
  Implicit None
  Save

  Integer                    :: N,M
  Real*8,parameter           :: pi=3.14159265d0
  Real*8                     :: delr,delz,delt,sig,sig2
  Real*8                     :: Q_ABDE,Q_C,fE,fE_old
  Real*8                     :: xAB,xAC,xBC,xAE,xAD,xBE,xBD,xCE,xCD,xED

  Real*8                     :: Vol
  Real*8                     :: Conv_w,Conv_p,dfffE,fE_hom
  Real*8                     :: kappaC,fracA,fracE,fracB,fracD,muABDE,muC,r_0
  Real*8                     :: XM(5,5)

  Real*8                     :: phiAB_hom,phiC_hom

  Integer                    :: bilayer,disk,ten_find
  Integer                    :: Mtip,Ntip,iter

  Real*8                     :: OP ! Order parameter

  Real*8,Allocatable         :: eta(:,:),dpp(:,:),eta2(:,:)
  

  !*************** A-Block ************************
  Integer                    :: NA
  Real*8,Allocatable         :: bAx(:),bAy(:),wA(:,:),dwA(:,:),pA(:,:)
  Real*8,Allocatable         :: qA(:,:,:),qA_0(:,:),qAD(:,:,:),qAD_0(:,:)
  Real*8,Allocatable         :: DiagAx(:),DiagAUx(:),DiagALx(:)
  Real*8,Allocatable         :: DiagAdx(:),DiagAUdx(:),DiagALdx(:)
  Real*8,Allocatable         :: DiagAy(:),DiagAUy(:),DiagALy(:)
  Real*8,Allocatable         :: DiagAdy(:),DiagAUdy(:),DiagALdy(:)
  !*************** B-Block ************************
  Integer                    :: NB
  Real*8,Allocatable         :: bBx(:),bBy(:),wB(:,:),dwB(:,:),pB(:,:)
  Real*8,Allocatable         :: qB(:,:,:),qB_0(:,:),qBD(:,:,:),qBD_0(:,:)
  Real*8,Allocatable         :: DiagBx(:),DiagBUx(:),DiagBLx(:)
  Real*8,Allocatable         :: DiagBdx(:),DiagBUdx(:),DiagBLdx(:)
  Real*8,Allocatable         :: DiagBy(:),DiagBUy(:),DiagBLy(:)
  Real*8,Allocatable         :: DiagBdy(:),DiagBUdy(:),DiagBLdy(:)
 !*************** C-Chain ************************
  Integer                    :: NC
  Real*8,Allocatable         :: bCx(:),bCy(:),wC(:,:),pC(:,:),dwC(:,:)
  Real*8,Allocatable         :: qC(:,:,:),qC_0(:,:)
  Real*8,Allocatable         :: DiagCx(:),DiagCUx(:),DiagCLx(:)
  Real*8,Allocatable         :: DiagCdx(:),DiagCUdx(:),DiagCLdx(:)
  Real*8,Allocatable         :: DiagCy(:),DiagCUy(:),DiagCLy(:)
  Real*8,Allocatable         :: DiagCdy(:),DiagCUdy(:),DiagCLdy(:)
  !*************** E-Block ************************
  Integer                    :: NE
  Real*8,Allocatable         :: bEx(:),bEy(:),wE(:,:),dwE(:,:),pE(:,:)
  Real*8,Allocatable         :: qE(:,:,:),qE_0(:,:),qED(:,:,:),qED_0(:,:)
  Real*8,Allocatable         :: DiagEx(:),DiagEUx(:),DiagELx(:)
  Real*8,Allocatable         :: DiagEdx(:),DiagEUdx(:),DiagELdx(:)
  Real*8,Allocatable         :: DiagEy(:),DiagEUy(:),DiagELy(:)
  Real*8,Allocatable         :: DiagEdy(:),DiagEUdy(:),DiagELdy(:)
  !*************** D-Block ************************
  Integer                    :: ND
  Real*8,Allocatable         :: bDx(:),bDy(:),wD(:,:),dwD(:,:),pD(:,:)
  Real*8,Allocatable         :: qD(:,:,:),qD_0(:,:),qDD(:,:,:),qDD_0(:,:)
  Real*8,Allocatable         :: DiagDx(:),DiagDUx(:),DiagDLx(:)
  Real*8,Allocatable         :: DiagDdx(:),DiagDUdx(:),DiagDLdx(:)
  Real*8,Allocatable         :: DiagDy(:),DiagDUy(:),DiagDLy(:)
  Real*8,Allocatable         :: DiagDdy(:),DiagDUdy(:),DiagDLdy(:)

end Module global




