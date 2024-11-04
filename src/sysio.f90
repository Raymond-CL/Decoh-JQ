module sysio

  use const
  
contains

  function getu() result(u)
  implicit none
  integer :: u
  logical :: ex
  do u = 10, 100
    inquire(unit=u,opened=ex)
    if(.not. ex) return
  enddo
  end function getu

  subroutine errmsg(code,str)
  integer, intent(in) :: code
  character(len=*), intent(in) :: str
  write(*,*) 'error:',str
  error stop code
  end subroutine errmsg

  subroutine readinput
    use vint
    use optflags
    use kinematics
    implicit none
    integer :: stat,u
    u = getu()
    open(u,iostat=stat,file='input.dat',status='old')
    if(stat.ne.0) call exit(11)
    read(u,*) CME
    read(u,*) ptTmin,ptTmax
    read(u,*) yTmin ,yTmax
    read(u,*) ptAmin,ptAmax
    read(u,*) yAmin ,yAmax
    read(u,*) npt1,itn1
    read(u,*) npt2,itn2
    read(u,*) quench_opt
    read(u,*) decoherent
    close(u)
    ! define some options and flags
    if(quench_opt.eq.0)  ndimn = 3
    if(quench_opt.eq.1)  ndimn = 3
    if(quench_opt.eq.2)  ndimn = 4
    if(quench_opt.eq.3)  ndimn = 7
    ! some fixed variable values need to be set
    prnt = -1
  end subroutine readinput
  
end module sysio