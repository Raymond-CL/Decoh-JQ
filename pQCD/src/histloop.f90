module histloop

  use const
  implicit none
  private
  real(wp) :: dmin,dmax
  logical :: islin
  integer, public :: id,nd
  real(wp), public :: dL,dM,dR,bin
  public :: setloop,setbin

contains

  subroutine setloop(a,b,n,scale)
  real(wp), intent(in) :: a,b
  integer, intent(in) :: n
  logical :: scale
  dmin = a
  dmax = b
  nd = n
  islin = scale
  if(islin) then
    bin = (dmax-dmin)/nd
  else
    if(dmin.le.0d0) dmin = 1e-10
    bin = log(dmax/dmin)/nd
  endif
  end subroutine setloop

  subroutine setbin
  if(islin) then
    dL = dmin + (id-1)*bin
    dR = dmin + id*bin
  else
    dL = dmin*exp((id-1)*bin)
    dR = dmin*exp(id*bin)
  endif
  dM = (dL+dR)/2d0
  end subroutine setbin

end module histloop