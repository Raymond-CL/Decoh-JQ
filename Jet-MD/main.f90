module md
  use iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), parameter :: Nc=3d0
  integer, parameter :: Nf=3
  real(wp), parameter :: b=11d0/3d0*Nc-2d0/3d0*Nf
  real(wp), parameter :: Lqcd=0.25d0;
  real(wp), parameter :: Q0=0.5d0
  real(wp), parameter :: lambda=log(Q0/Lqcd)
  real(wp), parameter :: ymin=0d0
  real(wp), parameter :: ymax=8d0
  integer, parameter :: Nn=4
  integer, parameter :: Ny=16
  real(wp) :: Pny_tab(0:Ny,0:Nn)
  real(wp) :: Sny_tab(0:Ny,1:Nn)

contains

  subroutine initialize
    integer :: n,y
    do y=0,Ny
      Pny_tab(y,0) = y*ymax/Ny
    end do
    do n=1,Nn
      do y=0,Ny
        Pny_tab(y,n) = 0d0
        Sny_tab(y,n) = 0d0
      end do
    end do
  end subroutine initialize

  subroutine print_Pny(u)
    integer,intent(in) :: u
    integer :: n,y
    do y=0,Ny
      do n=0,Nn
        write(u,'(es12.3)',advance='no') Pny_tab(y,n)
      enddo
      write(u,*)
    enddo
  end subroutine print_Pny

  subroutine print_Sny(u)
    integer,intent(in) :: u
    integer :: n,y
    do y=0,Ny
      do n=1,Nn
        write(u,'(es12.3)',advance='no') Sny_tab(y,n)
      enddo
      write(u,*)
    enddo
  end subroutine print_Sny

  ! function get_Pny(n,y) result(res)
  !   integer,intent(in) :: n
  !   real(wp), intent(in) :: y
  !   real(wp) :: res
  !   integer :: iy
  !   iy = int(y/ymax*Ny)
  !   if(iy.eq.0) then
  !     res=Pny_tab(0,n)
  !   elseif(iy.eq.Ny) then
  !     res = Pny_tab(Ny,n)
  !   else
  !   endif
  !   return
  ! end function get_Pny

  ! function fxn ( x )
  !   real(wp) :: fxn
  !   real(wp), intent(in) :: x
  !   ! fxn = exp(-x**2)
  !   fxn = bessel_j0(x)*bessel_j1(x)
  !   return
  ! end function fxn
end module md

program main
  use md
  use iso_fortran_env, only: stdout=>output_unit
  implicit none
  call timestamp
  call initialize
  ! call print_Pny(stdout)
  call timestamp
end program main

! subroutine dqag_test
!   use wrap
!   use iso_fortran_env, only: wp=>real64
!   implicit none
!   real(wp) :: a,b
!   real(wp) :: epsabs,epsrel
!   integer :: key
!   real(wp) :: result,abserr
!   integer :: neval,ier
!   integer, parameter :: limit = 1000
!   integer, parameter :: lenw = 10 * limit
!   integer :: last,iwork(limit),work(lenw)
!   real(wp) :: actual
!   a=0d0
!   b=100d0
!   epsabs=0d0
!   epsrel=1d-10
!   key=6
!   call dqag ( fxn, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
!     limit, lenw, last, iwork, work )
!   ! actual = ( erf(b) - erf(a) ) * sqrt(atan(1d0))
!   actual = 0.5d0 * ( bessel_j0(a)**2 - bessel_j0(b)**2 )
!   write(*,*) 'from ',a,' to ',b
!   write(*,*) 'integrated result =',result
!   write(*,*) 'actual result     =',actual
!   write(*,*) 'estimated error   =',abserr
!   write(*,*) 'actual error      =',abs(actual-result)
!   write(*,*) 'num. of eval.     =',neval
!   write(*,*) 'error code        =',ier
!   return
! end subroutine dqag_test
