module kinematics

  use const
  implicit none
  public
  real(wp) :: CME
  real(wp) :: ptTmin,ptTmax
  real(wp) :: yTmin ,yTmax
  real(wp) :: ptAmin,ptAmax
  real(wp) :: yAmin ,yAmax
  real(wp) :: ptTrig,ptAsso
  real(wp) :: yTrig,yAsso
  real(wp) :: pT,pTjet,pTq,pTg
  real(wp) :: x1,x2
  real(wp) :: x,y,th
  real(wp) :: mans,mant,manu
  real(wp) :: mufac,muren
  real(wp) :: as
  logical :: fillnow

contains

  pure function ampA(n1,n2,d) result(res)
  real(wp), intent(in) :: n1,n2,d
  real(wp) :: res
  res = 4d0/9d0*(n1*n1+n2*n2)/d/d
  return
  end function ampA

  pure function ampB(n,d1,d2) result(res)
  real(wp), intent(in) :: n,d1,d2
  real(wp) :: res
  res = ampA(n,d1,d2) + ampA(n,d2,d1) - 8d0/27d0*n*n/d1/d2
  return
  end function ampB

  pure function ampC(n1,n2,d) result(res)
  real(wp), intent(in) :: n1,n2,d
  real(wp) :: res
  res = 16d0/81d0*(n1*n1+n2*n2)/n1/n2 - ampA(n1,n2,d)
  return
  end function ampC

  pure function ampD(n1,n2,n3) result(res)
  real(wp), intent(in) :: n1,n2,n3
  real(wp) :: res
  res = 3d0-n1*n2/n3**2-n1*n3/n2**2-n2*n3/n1**2
  return
  end function ampD

end module kinematics