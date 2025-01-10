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
  !   - Pny_table(:,1:n) stores P_n values
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
  ! color factors
  real(wp), parameter :: Nc=3d0
  integer,  parameter :: Nf=3
  real(wp), parameter :: b=11d0/3d0*Nc-2d0/3d0*Nf
  ! QCD scales
  real(wp), parameter :: Lqcd=0.245748d0
  real(wp), parameter :: Q0=39.99d0 !0.5d0
  real(wp), parameter :: lambda=log(Q0/Lqcd)
  ! scale range y=log(Q/Q0)
  real(wp), parameter :: ymin=0d0
  real(wp), parameter :: ymax=8d0
  ! table dimensions
  integer,  parameter :: Nn=3
  integer,  parameter :: Ny=800
  real(wp), parameter :: ybin=(ymax-ymin)/Ny
  ! data table
  real(wp) :: Pny_tab(0:Ny,0:Nn)
  real(wp) :: Sny_tab(0:Ny,0:Nn-1)
  ! tmp var to keep track of current n and y
  integer :: n_now
  real(wp) :: y_now

contains

  ! initialize tables
  subroutine initialize
    integer :: in,iy
    do iy = 0,Ny
      Pny_tab(iy,0) = ymin + iy * ybin
    end do
    do in = 1,Nn
      do iy = 0,Ny
        Pny_tab(iy,in) = 0d0
        Sny_tab(iy,in-1) = 0d0
      end do
    end do
  end subroutine initialize

  ! print Pny
  subroutine print_Pny(u)
    integer, intent(in) :: u
    integer :: in,iy
    write(u,*) "P_ny"
    do iy = 0,Ny
      do in = 0,Nn
        write(u,'(es12.4)',advance='no') Pny_tab(iy,in)
      enddo
      write(u,*)
    enddo
    write(u,*)
  end subroutine print_Pny

  ! print Pny cumulative
  subroutine print_Pny_2(u)
    integer, intent(in) :: u
    integer :: in,iy
    real(wp) :: tmp
    write(u,*) "P_ny cumulative"
    do iy = 0,Ny
      write(u,'(es12.4)',advance='no') Pny_tab(iy,0)
      tmp = 0d0
      do in = 1,Nn
        tmp = tmp + Pny_tab(iy,in)
        write(u,'(es12.4)',advance='no') tmp
      enddo
      write(u,*)
    enddo
    write(u,*)
  end subroutine print_Pny_2

  ! print Sny
  subroutine print_Sny(u)
    integer,intent(in) :: u
    integer :: in,iy
    write(u,*) "S_ny"
    do iy = 0,Ny
      do in = 0,Nn-1
        write(u,'(es12.4)',advance='no') Sny_tab(iy,in)
      enddo
      write(u,*)
    enddo
    write(u,*)
  end subroutine print_Sny

  pure function gam2(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    res = 4d0 * Nc / b / (yp + lambda)
    return
  end function gam2

  function get_Pny(y,n) result(res)
    real(wp), intent(in) :: y
    integer, intent(in) :: n
    real(wp) :: res
    integer :: i
    real(wp) :: yL, yH, PL, PH
    if(y.lt.ymin .or. y.gt.ymax) res = 0d0
    i = floor((y-ymin)/ybin)
    yL = Pny_tab(i,0)
    yH = Pny_tab(i+1,0)
    PL = Pny_tab(i,n)
    PH = Pny_tab(i+1,n)
    res = PL + (PH - PL) * (y - yL) / (yH - yL)
    return
  end function get_Pny

  function Sny_int(yp) result(res)
    real(wp), intent(in) :: yp
    real(wp) :: res
    if(n_now.eq.1) then
      res = (y_now - yp) * gam2(yp)
    elseif(n_now.gt.1) then
      res = (y_now - yp) * gam2(yp) * get_Pny(yp,n_now-1)
    else
      res = 0d0
    endif
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
    key = 1
    do iy = 0,Ny
      y_now = Pny_tab(iy,0)
      ypmax = y_now
      call dqag ( Sny_int, ypmin, ypmax, epsabs, epsrel, &
        key, result, abserr, neval, ier, &
        limit, lenw, last, iwork, work )
      Sny_tab(iy,n-1) = result
    end do
  end subroutine set_Sny

  subroutine set_Pny(n)
    integer, intent(in) :: n
    integer :: iy,in
    if(n.eq.1) then
      do iy = 0,Ny
        Pny_tab(iy,n) = exp( - Sny_tab(iy,n-1) )
      end do
    elseif(n.gt.1) then
      do iy = 0,Ny
        do in = 1,n-1
          Pny_tab(iy,n) = real(in,wp)/(n-1) * Pny_tab(iy,n-in) * Sny_tab(iy,in)
        enddo
      end do
    else
      do iy = 0,Ny
        Pny_tab(iy,n) = 0d0
      end do
    endif
  end subroutine set_Pny

end module md

program main
  use md
  use iso_fortran_env, only: stdout=>output_unit
  implicit none
  integer :: n,u

  integer :: i,in
  double precision :: pt,ptm,ytmp
  ! double precision :: tmp,tot
  ! double precision :: y1,y2,y3
  ! call timestamp
  call initialize
  do n = 1,Nn
    call set_Sny(n)
    call set_Pny(n)
  enddo
  u = stdout
  ! u = 66
  ! open(unit=u,file='Pny.dat')
  ! do n = 1,Nn
  !   y1 = log(60d0/Q0)
  !   y2 = log(300d0/Q0)
  !   y3 = log(1000d0/Q0)
  !   write(u,*) n,get_Pny(y1,n),get_Pny(y2,n),get_Pny(y3,n)
  ! enddo
  ! call print_Sny(u)
  ! call print_Pny(u)

  do i=1,100
    pt = i*10d0 + 40d0
    ptm = pt * 1d0!0.4d0
    ytmp = log(ptm/Q0)
    write(u,'(es12.4)',advance='no') ptm
    write(u,'(es12.4)',advance='no') ytmp
    ! do in = 1,6
      write(u,'(es12.4)',advance='no') get_Pny(ytmp,1)
    ! enddo
    write(u,*)
  enddo
  ! tot = 0d0
  ! do iy = 1,Ny-2
  !   tmp = Pny_tab(iy,1)
  !   tot = tot + 0.5d0 * (7d0-tmp) * gam2(tmp) * tmp
  ! enddo
  ! write(*,*) tot
  ! call print_Pny_2(u)
  ! close(u)
  ! call timestamp
end program main
