module const
  use, intrinsic :: iso_fortran_env, only : sp => real32
  use, intrinsic :: iso_fortran_env, only : dp => real64
  use, intrinsic :: iso_fortran_env, only : ip => int64
  use, intrinsic :: iso_fortran_env, only : stdin => input_unit
  use, intrinsic :: iso_fortran_env, only : stdout => output_unit
  implicit none
  public
  integer, parameter :: wp = dp     ! working precision
  real(wp), parameter :: PI = 4d0*atan(1d0)
  real(wp), parameter :: twoPI = 2d0*PI
  real(wp), parameter :: PI2 = PI*PI
end module const



module phyconst
  use const
  implicit none
  ! conversion units
  public
  real(wp), parameter :: planck = 6.62607015d-34   ! Planck constant (J.s)
  real(wp), parameter :: planckbar = planck/twoPI  ! reduced Planck (J.s)
  real(wp), parameter :: charge = 1.602176634d-19  ! chrg mag (J/eV)
  real(wp), parameter :: hbar = planckbar/charge   ! hbar (eV.s)
  real(wp), parameter :: light = 299792458d0      ! light speed (m/s)
  real(wp), parameter :: hbarc = hbar*light      ! conversion unit (eV.m)
  real(wp), parameter :: giga = 1d+9
  real(wp), parameter :: mili = 1d-3
  real(wp), parameter :: micro = 1d-6
  real(wp), parameter :: nano = 1d-9
  real(wp), parameter :: pico = 1d-12
  real(wp), parameter :: femto = 1d-15
  real(wp), parameter :: fm2barn = 1d-2         ! fm^2->barns (barn/fm^2)
  real(wp), parameter :: gevfm = hbarc/giga/femto  ! common conversion (GeV.fm)
  real(wp), parameter :: gev2barn = gevfm*gevfm*fm2barn ! (GeV^2.barn)
end module phyconst