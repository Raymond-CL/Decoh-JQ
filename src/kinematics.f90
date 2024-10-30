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

  pure function ampA(s,t,u) result(res)
  real(wp), intent(in) :: s,t,u
  real(wp) :: res
  res = 4d0/9d0 * (s*s + t*t)/u/u
  return
  end function ampA

  pure function ampB(s,t,u) result(res)
  real(wp), intent(in) :: s,t,u
  real(wp) :: res
  res = 4d0/9d0 * (s*s + t*t)/u/u
  return
  end function ampB

end module kinematics