module bdmps

  use const
  use qcd, only : CA,CF
  implicit none
  real(wp) :: wcq,wcg

contains

  function De(wc,eps,i) result(res)
  real(wp), intent(in) :: wc,eps
  integer, intent(in) :: i
  real(wp) :: res
  real(wp) :: as,al,temp
  as = 0.2d0
  if(i.eq.0) then
    al = 2d0*CA/PI * as
  else
    al = 2d0*CF/PI * as
  endif
  temp = al**2 * wc / 2d0 / eps
  res = sqrt(temp) * exp( -PI * temp ) / eps
  end function De

end module bdmps