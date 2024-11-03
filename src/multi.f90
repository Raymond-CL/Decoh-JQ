module multi

  use const
  use qcd, only : Nc,CA,CF,TR
  implicit none
  private
  integer , parameter :: Nf = 3d0
  real(wp), parameter :: Lqcd = 0.247
  real(wp), parameter :: b = (11d0*CA-2d0*Nf)/(12d0*PI)
  real(wp), parameter :: a = 1d0/4d0+5d0*Nf/(54d0*PI*b)
  real(wp), parameter :: AA = sqrt(144d0/(33d0-2d0*nf))
  real(wp), parameter :: BB = (33d0+2d0*nf/9d0)/(33d0-2d0*nf)
  public :: getpt,nMLL2

contains

  function getpt(ptfinal,R,e,Qm) result(ptinit)
  real(wp), intent(in) :: ptfinal,R,e,Qm
  real(wp) :: ptinit
  real(wp) :: ptin,ptout
  ptin = ptfinal
  ptout = ptfinal + e * nMLL2(ptfinal*R,Qm)
  do while(abs(ptin-ptout).gt.1d0)
    ptin = ptout
    ptout = ptfinal + e * nMLL2(ptin*R,Qm)
  enddo
  ptinit = ptout
  return
  end function getpt

  function nLL1(Q) result(n)
  real(wp), intent(in) :: Q
  real(wp) :: n
  real(wp) :: c,temp
  c = 0.017d0
  temp = sqrt( 2d0*CA/PI/b * log(Q*Q/Lqcd/Lqcd) )
  n = c * exp(temp)
  return
  end function nLL1

  function nMLL1(Q) result(n)
  real(wp), intent(in) :: Q
  real(wp) :: n
  real(wp) :: c,temp1,temp2
  c = 0.04d0
  temp1 = sqrt( 2d0*CA/PI/b * log(Q*Q/Lqcd/Lqcd) )
  temp2 = a * log(1d0/b/log(Q*Q/Lqcd/Lqcd))
  n = c * exp(temp1 + temp2)
  return
  end function nMLL1  

  function nLL2(Q,Qm) result(n)
  real(wp), intent(in) :: Q,Qm
  real(wp) :: n
  real(wp) :: c,x1,x2
  real(wp) :: BI0,BI1,BK0,BK1
  real(wp), external :: besi0,besi1,besk0,besk1
  c = 0.0023d0
  x1 = AA*sqrt(log(Q/Lqcd))
  x2 = AA*sqrt(log(Qm/Lqcd))
  BI0 = besi0 (x2)
  BI1 = besi1 (x1)
  BK0 = besk0 (x2)
  BK1 = besk1 (x1)
  n = c * x1 * (BI1*BK0 + BK1*BI0)
  return
  end function nLL2

  function nMLL2(Q,Qm) result(n)
  real(wp), intent(in) :: Q,Qm
  real(wp) :: n
  real(wp) :: c,x1,x2
  real(wp) :: BI0,BI1,BK0,BK1
  real(wp) :: val(3)
  integer :: stat
  c = 0.23d0
  x1 = AA*sqrt(log(Q/Lqcd))
  x2 = AA*sqrt(log(Qm/Lqcd))
  call ribesl(x1,BB-1d0,3,1,val,stat)
  BI1 = val(3)
  call rkbesl(x2,BB-1d0,2,1,val,stat)
  BK0 = val(2)
  call rkbesl(x1,BB-1d0,3,1,val,stat)
  BK1 = val(3)
  call ribesl(x2,BB-1d0,2,1,val,stat)
  BI0 = val(2)
  n = c * x1 * ((x2/x1)**BB) * (BI1*BK0 + BK1*BI0)
  return
  end function nMLL2

end module multi