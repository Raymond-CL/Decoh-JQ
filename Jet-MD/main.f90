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
  !   - Pny_table(:,0) stores y values (the 0th column)
  !   - Pny_table(:,1:n) stores P_n values (from 1st column on-wards)
  !   - Sny_table(:,n-1) is used to stores temporary integration results
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
  ! QCD parameters
  real(wp), parameter :: Nc=3d0
  real(wp), parameter :: CA=Nc
  real(wp), parameter :: CF=(Nc**2-1d0)/(2d0*Nc)
  integer,  parameter :: Nf=3
  real(wp), parameter :: b=11d0/3d0*Nc-2d0/3d0*Nf
  real(wp), parameter :: Lqcd=0.245748d0
  ! hadronization (perturbative) scales
  real(wp), parameter :: Q0=1.45d0 !0.5d0
  real(wp), parameter :: lambda=log(Q0/Lqcd)
  ! scale range y=log(Q/Q0)
  real(wp), parameter :: ymin=0d0
  real(wp), parameter :: ymax=8d0
  ! table dimensions
  integer,  parameter :: Nn=1000
  integer,  parameter :: Ny=80
  real(wp), parameter :: ybin=(ymax-ymin)/Ny
  ! data table
  real(wp) :: Pg(0:Ny,0:Nn)
  real(wp) :: Pq(0:Ny,0:Nn)
  real(wp) :: Sg(0:Ny,1:Nn-1)
  real(wp) :: Sq(0:Ny,1:Nn-1)
  real(wp) :: nbarg(0:Ny)
  real(wp) :: nbarq(0:Ny)
  ! integration accuracy
  integer :: ikey = 1 ! 1-6, 6 for highest accuracy
  ! tmp var to keep track of current n and y
  integer :: n_now
  real(wp) :: y_now

contains

  ! initialize tables
  subroutine initialize
    integer :: in,iy
    do iy = 0,Ny
      y_now = ymin + iy * ybin
      Pg(iy,0) = y_now
      Pq(iy,0) = y_now
    end do
    do in = 1,Nn
      do iy = 0,Ny
        Pg(iy,in) = 0d0
        Pq(iy,in) = 0d0
      end do
    end do
    do in = 1,Nn-1
      do iy = 0,Ny
        Sg(iy,in) = 0d0
        Sq(iy,in) = 0d0
      end do
    end do
  end subroutine initialize

  subroutine print_Pg(u)
    integer, intent(in) :: u
    integer :: iy,in
    ! write(u,*) "Pg(y,n) table"
    do iy = 0,Ny
      do in = 0,Nn
        write(u,'(es12.4)',advance='no') Pg(iy,in)
      enddo
      write(u,*)
    enddo
    write(u,*)
  end subroutine print_Pg

  ! print Pny cumulative
  ! subroutine print_Pny_2(u)
  !   integer, intent(in) :: u
  !   integer :: in,iy
  !   real(wp) :: tmp
  !   write(u,*) "P_ny cumulative"
  !   do iy = 0,Ny
  !     write(u,'(es12.4)',advance='no') Pny_tab(iy,0)
  !     tmp = 0d0
  !     do in = 1,Nn
  !       tmp = tmp + Pny_tab(iy,in)
  !       write(u,'(es12.4)',advance='no') tmp
  !     enddo
  !     write(u,*)
  !   enddo
  !   write(u,*)
  ! end subroutine print_Pny_2

  subroutine print_Sg(u)
    integer,intent(in) :: u
    integer :: iy,in
    write(u,*) "Sg(y,n) table"
    do iy = 0,Ny
      write(u,'(12x)',advance='no')
      do in = 1,Nn-1
        write(u,'(es12.4)',advance='no') Sg(iy,in)
      enddo
      write(u,*)
    enddo
    write(u,*)
  end subroutine print_Sg

  subroutine print_nbarg(u)
    integer,intent(in) :: u
    integer :: iy
    write(u,*) "<ng>(y) table"
    do iy = 0,Ny
      write(u,'(2es12.4)') Pg(iy,0), nbarg(iy)
    enddo
    write(u,*)
  end subroutine print_nbarg

  ! \gamma_0^2(yp) function (pure may optimize)
  pure function gam2g(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    res = 4d0 * CA / b / (yp + lambda)
    return
  end function gam2g

  ! pure function gam2q(yp) result(res)
  !   real(wp), intent(in) :: yp
  !   real(wp) :: res
  !   res = 4d0 * CF / b / (yp + lambda)
  !   return
  ! end function gam2q

  subroutine set_Pg(n)
    integer, intent(in) :: n
    integer :: iy,in
    real(wp) :: y_cap,total
    if(n.eq.1) then
      do iy = 0,Ny
        y_now = Pg(iy,0)
        y_cap = y_now + lambda
        Pg(iy,n) = exp( gam2g(y_now)*y_cap*(y_now + y_cap*log(lambda/y_cap)) )
      end do
    elseif(n.gt.1) then
      call set_Sg(n-1)
      do iy = 0,Ny
        total = 0d0
        do in = 1,n-1
          total = total + real(in,wp)/(n-1) * Pg(iy,n-in) * Sg(iy,in)
        enddo
        Pg(iy,n) = total
      end do
    endif
  end subroutine set_Pg

  ! interpolation function
  function get_Pg(y,n) result(res)
    real(wp), intent(in) :: y
    integer, intent(in) :: n
    real(wp) :: res
    integer :: i
    real(wp) :: yL, yH, PL, PH
    if(y.lt.ymin .or. y.gt.ymax) res = 0d0
    i = floor((y-ymin)/ybin)
    yL = Pg(i,0)
    yH = Pg(i+1,0)
    PL = Pg(i,n)
    PH = Pg(i+1,n)
    res = PL + (PH - PL) * (y - yL) / (yH - yL)
    return
  end function get_Pg

  ! integrand function
  function Sg_int(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    res = 0d0
    if(n_now.le.1) return
    res = (y_now - yp) * gam2g(yp) * get_Pg(yp,n_now-1)
    return
  end function Sg_int

  ! calculate Sny with QAG integration
  subroutine set_Sg(n)
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
    n_now = n+1
    ypmin = ymin
    epsabs = 0d0
    epsrel = 1d-10
    key = ikey
    do iy = 0,Ny
      y_now = Pg(iy,0)
      ypmax = y_now
      call dqag ( Sg_int, ypmin, ypmax, epsabs, epsrel, &
        key, result, abserr, neval, ier, &
        limit, lenw, last, iwork, work )
      Sg(iy,n) = result
    end do
  end subroutine set_Sg

  subroutine set_nbarg
    integer :: iy,in
    real(wp) :: avg
    do iy = 0,Ny
      do in = 0,Nn
        avg = avg + in * Pg(iy,in)
      end do
      nbarg(iy) = avg
    end do
  end subroutine set_nbarg


end module md

! main program (will be in separate file)
program main
  use md
  use iso_fortran_env, only: stdout=>output_unit
  implicit none
  integer :: u
  integer :: i
  call timestamp

  call initialize
  do i=1,Nn
    call set_Pg(i)
  end do
  call set_nbarg

  ! u = 66
  ! open(unit=u,file='Pyn.dat')
  ! call print_Pg(u)
  ! close(u)
  
  u = stdout
  call print_nbarg(u)

  call timestamp
end program main
