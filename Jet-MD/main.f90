module md

  ! ------------------------------------
  ! Jet multiplicity distribution module
  ! ------------------------------------
  ! module to generate table of jet multiplicity distribution 
  ! Pn(y), where n is number of partons and y is scale difference
  ! Requires:
  !   - dqag() subroutine from quadpack_double.f90
  !     https://people.math.sc.edu/Burkardt/f_src/quadpack_double/quadpack_double.html
  ! Calculates:
  !   - P_n(y) = sum_{i=1}^{n-1} i/(n-1) * P_{n-i}(y) * S_{i}(y)
  !   - S_i(y) = \int_0^y dy' (y-y') * \gamma(y') P_i(y')
  !   - using an iterative algorithm, we can generate 
  !     a table of Pn(y) up to n = Nn
  !   - Pny_table(:,0) stores y values
  !   - Sny_table(:,n-1) is used to stores integration results
  ! Parameters:
  !   - Nf, Lqcd: for Nf=3, Lqcd~0.247
  !       we can use other values of Nf and Lqcd
  !   - Q0: stopping scale
  !       usually set to Q0=Lqcd, but can be a bit larger
  !       to prevent singularity
  !       or set to heavy quark mass
  !       or QGP medium scale
  !   - ymax: maximum scale range
  !       ymax=8 is about Q=745, which is about a
  !       1800 GeV jet with R=0.4, more than enough
  !       for current energies
  !   - Nn: number of partons
  !       determines how large the Pny table will be
  !       e.g. set Nn=20 to see the probability of a jet
  !       that can spawn 20 partons from Q to Q0
  !   - Ny: number of y values
  !       also determines how large the Pny table will be
  !       the larger the Ny, the more accurate the
  !       interpolation result will be
  ! Note:
  !   - if set working precision wp=>real32 or wp=>real128, 
  !     the single/quadriple precision version of dqag()
  !     is needed, and can be found in same link.

  use iso_fortran_env, only: wp=>real64
  implicit none
  ! color factors
  real(wp), parameter :: Nc=3d0
  integer,  parameter :: Nf=3
  real(wp), parameter :: b=11d0/3d0*Nc-2d0/3d0*Nf
  ! QCD scales
  real(wp), parameter :: Lqcd=0.25d0;
  real(wp), parameter :: Q0=0.5d0
  real(wp), parameter :: lambda=log(Q0/Lqcd)
  ! scale range y=log(Q/Q0)
  real(wp), parameter :: ymin=0d0
  real(wp), parameter :: ymax=8d0
  ! table dimensions
  integer,  parameter :: Nn=2
  integer,  parameter :: Ny=16
  ! data table
  real(wp) :: Pny_tab(0:Ny,0:Nn)
  real(wp) :: Sny_tab(0:Ny,0:Nn-1)
  ! tmp var to keep track of current n and y
  integer :: n_now
  real(wp) :: y_now

contains

  ! initialize tables
  subroutine initialize
    integer :: n,y
    do y = 0,Ny
      Pny_tab(y,0) = y * ymax / Ny
    end do
    do n = 1,Nn
      do y = 0,Ny
        Pny_tab(y,n) = 0d0
        Sny_tab(y,n-1) = 0d0
      end do
    end do
  end subroutine initialize

  ! print Pny
  subroutine print_Pny(u)
    integer,intent(in) :: u
    integer :: n,y
    do y = 0,Ny
      do n = 0,Nn
        write(u,'(es12.3)',advance='no') Pny_tab(y,n)
      enddo
      write(u,*)
    enddo
  end subroutine print_Pny

  ! print Sny
  subroutine print_Sny(u)
    integer,intent(in) :: u
    integer :: n,y
    do y = 0,Ny
      do n = 0,Nn-1
        write(u,'(es12.3)',advance='no') Sny_tab(y,n)
      enddo
      write(u,*)
    enddo
  end subroutine print_Sny

  pure function gam(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    res = 4d0 * Nc / b / (yp + lambda)
    return
  end function gam

  function Sny_int(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    res = (y_now - yp) * gam(yp)
    return
  end function Sny_int

  subroutine set_Sny(n)
    integer, intent(in) :: n
    integer :: iy
    real(wp) :: ypmin,ypmax
    real(wp) :: epsabs,epsrel
    integer :: key
    real(wp) :: result,abserr
    integer :: neval,ier
    integer, parameter :: limit = 1000
    integer, parameter :: lenw = 10 * limit
    integer :: last,iwork(limit),work(lenw)
    n_now = n
    ypmin = ymin
    epsabs = 0d0
    epsrel = 1d-10
    key = 6
    do iy = 0,Ny
      ypmax = Pny_tab(iy,0)
      y_now = ypmax
      call dqag ( Sny_int, ypmin, ypmax, epsabs, epsrel, &
        key, result, abserr, neval, ier, &
        limit, lenw, last, iwork, work )
      Sny_tab(iy,n-1) = result
    end do
  end subroutine set_Sny

  subroutine set_Pny(n)
    integer, intent(in) :: n
    integer :: iy
    do iy = 0,Ny
      Pny_tab(iy,n) = exp( - Sny_tab(iy,n-1) )
    end do
  end subroutine set_Pny

end module md

program main
  use md
  use iso_fortran_env, only: stdout=>output_unit
  implicit none
  call timestamp
  call initialize
  call set_Sny(1)
  call set_Pny(1)
  call print_Sny(stdout)
  call print_Pny(stdout)
  call timestamp
end program main