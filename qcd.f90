module qcd

  use, intrinsic :: iso_fortran_env, only : wp => real64
  implicit none
  public
  real(wp), parameter :: PI = 4d0*atan(1d0)
  real(wp), parameter :: PI2 = PI*PI
  real(wp), parameter :: PI3 = PI*PI*PI
  integer, parameter :: Nc = 3
  real(wp), parameter :: Ca = dble(Nc)
  real(wp), parameter :: Cf = (dble(Nc)**2-1d0)/2d0/dble(Nc)
  real(wp), parameter :: Tr = 1d0/2d0
  real(wp), parameter :: AsmZ = 0.1181d0
  integer, public :: Nf,qcdloop
  real(wp), protected :: b0,b1,b2
  real(wp), public :: Lqcd,Lqcd2

contains

  subroutine setqcd(x,y)
  integer, intent(in) :: x,y
  Nf = x
  qcdloop = y
  b0 = (33d0-2d0*dble(Nf))/(12d0*PI)
  b1 = (153d0-19d0*dble(Nf))/(24d0*PI2)
  b2 = (2857d0-5033d0/9d0*dble(Nf)+325d0/27d0*dble(Nf)**2) &
     /(128d0*PI3)
  Lqcd = 0d0
  if(qcdloop .eq. 1) then
    if(Nf.eq.3) Lqcd = 0.247d0
    if(Nf.eq.4) Lqcd = 0.155d0
    if(Nf.eq.5) Lqcd = 0.0884d0
    if(Nf.eq.6) Lqcd = 0.0456d0
  elseif(qcdloop .eq. 2) then
    if(Nf.eq.3) Lqcd = 0.741d0
    if(Nf.eq.4) Lqcd = 0.438d0
    if(Nf.eq.5) Lqcd = 0.228d0
    if(Nf.eq.6) Lqcd = 0.0987d0
  elseif(qcdloop .eq. 3) then
    if(Nf.eq.3) Lqcd = 0.641d0
    if(Nf.eq.4) Lqcd = 0.389d0
    if(Nf.eq.5) Lqcd = 0.210d0
    if(Nf.eq.6) Lqcd = 0.0953d0
  endif
  Lqcd2 = Lqcd**2
  end subroutine setqcd

  function alphas(q) result(as)
  real(wp), intent(in) :: q
  real(wp) :: as
  real(wp) :: q2,l2,t,l
  q2 = q**2
  l2 = lqcd**2
  t = log(q2/l2)
  l = log(t)
  as = 1d0/(b0*t)
  if(qcdloop.ge.2) then
    as = as * -b1*l/(b0**2*t)
  endif
  if(qcdloop.ge.3) then
    as = as * +(b1**2*(l**2-l-1d0)+b0*b2)/(b0**4*t**2)
  endif
  return
  end function alphas

  function Pqq(z) result(res)
  real(wp), intent(in) :: z
  real(wp) :: res
  res = CF*(1d0+z*z)/(1d0-z)
  if(z.eq.1d0) res = res + CF*3d0/2d0
  return
  end function Pqq

  function Pgq(z) result(res)
  real(wp), intent(in) :: z
  real(wp) :: res
  res = CF*(1d0+(1d0-z)**2)/z
  return
  end function Pgq

  function Pqg(z) result(res)
  real(wp), intent(in) :: z
  real(wp) :: res
  res = Tr*(z*z+(1d0-z)**2)
  return
  end function Pqg

  function Pgg(z) result(res)
  real(wp), intent(in) :: z
  real(wp) :: res
  real(wp) :: beta0
  beta0 = 11d0/3d0*CA-4d0/3d0*Tr*Nf
  res = 2d0*CA*(z/(1d0-z)+(1d0-z)/z+z*(1d0-z))
  if(z.eq.1d0) res = res + beta0/2d0
  return
  end function Pgg

end module qcd
