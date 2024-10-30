module kinematics
  use, intrinsic :: iso_fortran_env, only : wp => real64
  implicit none
  private :: wp
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
end module kinematics