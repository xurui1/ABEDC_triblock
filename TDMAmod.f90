Subroutine TDMA(N_in, DiagL_in, Diag_in, DiagU_in, b_in)
  Use global
  Implicit None

  integer                :: N_in,i
  Real*8                   :: DiagL_in(N_in-1),DiagU_in(N_in-1),Diag_in(N_in),b_in(N_in)
  Real*8,Allocatable       :: d(:),c(:)

  Allocate(d(N_in),c(N_in))

  c(1)=DiagU_in(1)/Diag_in(1)
  d(1)=b_in(1)/Diag_in(1)

  Do i=2,(N_in-1)
     c(i)=DiagU_in(i)/(Diag_in(i)-(c(i-1)*DiagL_in(i-1)))
  End Do

  Do i=2,N_in
     d(i)=(b_in(i)-(d(i-1)*DiagL_in(i-1)))/(Diag_in(i)-(c(i-1)*DiagL_in(i-1)))
  End Do

  b_in(N_in)=d(N_in)
  Do i=(N_in-1),1,-1
     b_in(i)=d(i)-(c(i)*b_in(i+1))
  End Do

  DeAllocate(d,c)

End Subroutine TDMA
