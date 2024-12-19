module eventcounter
  use const
  implicit none
  real(wp), public :: eff
  integer(ip), public :: accevent,totevent
  real(wp), public :: tot_xsec
  logical, public :: docount
end module eventcounter