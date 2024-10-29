module const
  use, intrinsic :: iso_fortran_env, only : wp => real64
  implicit none
  real(wp), parameter :: PI = 4d0*atan(1d0)
  real(wp), parameter :: twoPI = 2d0*PI
  real(wp), parameter :: PI2 = PI*PI
end module const



module phyconst
  use, intrinsic :: iso_fortran_env, only : wp => real64
  use const
  implicit none
  ! conversion units
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



module optflags
  implicit none
  ! options and flags
  integer :: quench_opt
end module optflags



module eventcounter
  use, intrinsic :: iso_fortran_env, only : wp => real64
  use, intrinsic :: iso_fortran_env, only : ip => int64
  implicit none
  integer(ip) :: accevent,totevent
  real(wp) :: eff
end module eventcounter



module kinematics
  use, intrinsic :: iso_fortran_env, only : wp => real64
  implicit none
  real(wp) :: CME
  real(wp) :: ptTmin,ptTmax
  real(wp) :: yTmin ,yTmax
  real(wp) :: ptAmin,ptAmax
  real(wp) :: yAmin ,yAmax
  real(wp) :: ptTrig,ptAsso
  real(wp) :: yTrig,yAsso
  real(wp) :: pT,pTjet,pTq,pTg
  real(wp) :: x1,x2
  real(wp) :: x,y,th
  real(wp) :: mans,mant,manu
  real(wp) :: mufac,muren
  real(wp) :: as
  logical :: fillnow
end module kinematics



function getu() result(u)
  implicit none
  integer :: u
  logical :: ex
  do u = 10, 100
    inquire(unit=u,opened=ex)
    if(.not. ex) return
  enddo
end function getu



subroutine readinput
  use vint
  use optflags
  use kinematics
  implicit none
  integer :: stat,u
  interface
    function getu() result(unit)
      implicit none
      integer :: unit
    end function getu
  end interface
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
  close(u)

  quench_opt = 0
  ! 0: no quenching (pp)
  ! 1: constant Delta-E
  ! 2: BDMPS D(e)
  ! 3: BDMPS+hydro geometry qhat
  if(quench_opt.eq.0)  ndimn = 3
  if(quench_opt.eq.1)  ndimn = 3
  if(quench_opt.eq.2)  ndimn = 4
  if(quench_opt.eq.3)  ndimn = 7

  prnt = -1
end subroutine readinput



program main
  use, intrinsic :: iso_fortran_env, only : stdout => output_unit
  use time
  use vint
  use eventcounter
  use optflags
  use kinematics
  use qcd, only : setqcd
  implicit none
  character(len=40), parameter :: pdfgrid = './grids/i2Tn2.00.pds'
  real(wp) :: dmin,dmax,bin,dL,dM,dR
  integer :: id,nd

  call time_start

  call readinput

  call setCT18(pdfgrid)

  call setqcd(5,1)

  ! linear scale
  nd = 48
  dmin = ptAmin;  dmax = ptAmax
  bin = (dmax-dmin)/nd
  ! log scale
  ! nd = 100
  ! dmin = 1d-3;  dmax = 1d0
  ! bin = log(dmax/dmin)/nd

  do id = 1,nd
    ! linear scale
    dL = dmin + (id-1)*bin
    dR = dmin + id*bin
    ! log scale
    ! dL = dmin*exp((id-1)*bin)
    ! dR = dmin*exp(id*bin)
    dM = (dL+dR)/2d0

    limits(1) = yTmin;          limits(ndimn+1) = yTmax
    limits(2) = yAmin;          limits(ndimn+2) = yAmax
    limits(3) = dL;             limits(ndimn+3) = dR
    
    initial = -1
    call vegas(limits(1:2*ndimn),fxn,initial,npt1,itn1,prnt,intres,stddev,chisq)
    initial = +1
    accevent = 0;  totevent = 0
    call vegas(limits(1:2*ndimn),fxn,initial,npt2,itn2,prnt,intres,stddev,chisq)
    eff = dble(accevent)/dble(totevent)*100d0
    write(*,'(5(es11.3),f10.2)') dL,dM,dR,intres/bin,stddev,eff
    !write(*,'(a,2(es15.4),a)') 'vegas result:',intres,stddev/intres,new_line('a')
  enddo

  call time_stop
  call print_time
end program main



function fxn(dx,wgt)
  use const
  use phyconst
  use eventcounter
  use kinematics
  use qcd, only: alphas,nf
  implicit none
  double precision, dimension(:), intent(in) :: dx
  double precision, intent(in) :: wgt
  double precision :: fxn
  double precision, dimension(-6:+6) :: pdf1,pdf2,jff3,jff4
  integer :: i,j
  real(wp) :: dis1t,dis1u,dis2,dis3,dis4t,dis4u,dis5,&
  dis6,dis7tq,dis7tg,dis7uq,dis7ug,dis8!,dis7t,dis7u
  real(wp) :: sig1t,sig1u,sig2,sig3,sig4t,sig4u,sig5,&
  sig6,sig7tq,sig7tg,sig7uq,sig7ug,sig8!,sig7t,sig7u
  !double precision :: dis1,dis2,dis3,dis4,dis5,dis6,dis7,dis8
  !double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8
  real(wp) :: tempfunc
  interface
    function CT18Pdf(iparton,x,Q)
      implicit double precision (a-h,o-z)
    end function CT18Pdf
    function CT18Alphas (Q)
      implicit double precision (a-h,o-z)
    end function CT18Alphas
  end interface

  fxn = 0d0
  totevent = totevent + 1

  yTrig = dx(1)
  yAsso = dx(2)
  pTjet = dx(3)

  pT = pTjet
  pTq = pTjet
  pTg = pTjet

  x1 = pT / CME * (exp(+yTrig) + exp(+yAsso))
  x2 = pT / CME * (exp(-yTrig) + exp(-yAsso))
  if(x1.le.0d0 .or. x1.gt.1d0) return
  if(x2.le.0d0 .or. x2.gt.1d0) return

  mans = +x1*x2*CME*CME
  mant = -x1*CME*pT*exp(-yTrig)
  manu = -x2*CME*pT*exp(+yTrig)

  mufac = pT
  muren = mufac
  !as = CT18Alphas(muren)
  as = alphas(muren)

  pdf1 = 1d0;         pdf2 = 1d0
  jff3 = 1d0;         jff4 = 1d0

  do i = -nf, +nf
    pdf1(i) = CT18Pdf(i,x1,mufac)
    pdf2(i) = CT18Pdf(i,x2,mufac)
  enddo

  ! set distributions
  ! mass ordered parton flavours
  ! i: -6  -5  -4  -3  -2  -1  0  +1  +2  +3  +4  +5  +6
  !    tb  bb  cb  sb  db  ub  g   u   d   s   c   b   t

  ! initialize
  pdf1=1d0; pdf2=1d0; jff3=1d0; jff4=1d0
  
  do i = -nf,+nf
    pdf1(i) = CT18pdf(i,dble(x1),dble(mufac))
    pdf2(i) = CT18pdf(i,dble(x2),dble(mufac))
  enddo
  
  ! q(qb) + q'(qb') -> q(qb) + q'(qb') : only t-channel
  dis1t = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis1t = dis1t + (pdf1(+i)*pdf2(+j)*jff3(+i)*jff4(+j) &
                    +pdf1(+i)*pdf2(-j)*jff3(+i)*jff4(-j) &
                    +pdf1(-i)*pdf2(+j)*jff3(-i)*jff4(+j) &
                    +pdf1(-i)*pdf2(-j)*jff3(-i)*jff4(-j))
  endif
  enddo
  enddo
  sig1t = as*as * 4d0/9d0 * (mans**2 + manu**2)/mant**2

  ! q(qb) + q'(qb') -> q'(qb') + q(qb) : only u-channel
  dis1u = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis1u = dis1u + (pdf1(+i)*pdf2(+j)*jff3(+j)*jff4(+i) &
                    +pdf1(+i)*pdf2(-j)*jff3(-j)*jff4(+i) &
                    +pdf1(-i)*pdf2(+j)*jff3(+j)*jff4(-i) &
                    +pdf1(-i)*pdf2(-j)*jff3(-j)*jff4(-i))
  endif
  enddo
  enddo
  sig1u = as*as * 4d0/9d0 * (mans**2 + mant**2)/manu**2

  ! q(qb) + qb(q) -> q'(qb') + qb'(q') : only s-channel
  dis2 = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis2 = dis2 + (pdf1(+i)*pdf2(-i)*jff3(+j)*jff4(-j) &
                  +pdf1(+i)*pdf2(-i)*jff3(-j)*jff4(+j) &
                  +pdf1(-i)*pdf2(+i)*jff3(+j)*jff4(-j) &
                  +pdf1(-i)*pdf2(+i)*jff3(-j)*jff4(+j))
  endif
  enddo
  enddo
  sig2 = as*as * 4d0/9d0 * (mant**2 + manu**2)/mans**2

  ! q(qb) + q(qb) -> q(qb) + q(qb) : t- and u-channel
  dis3 = 0d0
  do i = 1,nf
    dis3 = dis3 + (pdf1(+i)*pdf2(+i)*jff3(+i)*jff4(+i) &
                  +pdf1(-i)*pdf2(-i)*jff3(-i)*jff4(-i))
  enddo
  sig3 = as*as * 4d0/9d0 * ((mans**2+manu**2)/mant**2 &
                           +(mans**2+mant**2)/manu**2 &
                           -2d0/3d0*mans**2/manu/mant)

  ! q(qb) + qb(q) -> q(qb) + qb(q) : t- and s-channel
  dis4t = 0d0
  do i = 1,nf
    dis4t = dis4t + (pdf1(+i)*pdf2(-i)*jff3(+i)*jff4(-i) &
                    +pdf1(-i)*pdf2(+i)*jff3(-i)*jff4(+i))
  enddo
  sig4t =  as*as * 4d0/9d0 * ((mans**2+manu**2)/mant**2 &
                             +(mant**2+manu**2)/mans**2 &
                             -2d0/3d0*manu**2/mans/mant)

  ! q(qb) + qb(q) -> qb(q) + q(qb) : u- and s-channel
  dis4u = 0d0
  do i = 1,nf
    dis4u = dis4u + (pdf1(+i)*pdf2(-i)*jff3(-i)*jff4(+i) &
                    +pdf1(-i)*pdf2(+i)*jff3(+i)*jff4(-i))
  enddo
  sig4u =  as*as * 4d0/9d0 * ((mans**2+mant**2)/manu**2 &
                             +(mant**2+manu**2)/mans**2 &
                             -2d0/3d0*mant**2/mans/manu)

  ! q(qb) + qb(q) -> g + g : scatter(t-,u-channel) + annihlate(s-channel)
  dis5 = 0d0
  do i = 1,nf
    dis5 = dis5 + (pdf1(+i)*pdf2(-i)*jff3(0)*jff4(0) &
                  +pdf1(-i)*pdf2(+i)*jff3(0)*jff4(0))
  enddo
  sig5 = as*as * (32d0/27d0 * (mant**2+manu**2)/mant/manu &
                 -8d0 / 3d0 * (mant**2+manu**2)/mans**2)

  ! g + g -> q(qb) + qb(q) : scatter(t-,u-channel) + creation(s-channel)
  dis6 = 0d0
  do j=1,nf
    dis6 = dis6 + (pdf1(0)*pdf2(0)*jff3(+j)*jff4(-j) &
                  +pdf1(0)*pdf2(0)*jff3(-j)*jff4(+j))
  enddo
  sig6 = as*as * (1d0/6d0 * (mant**2+manu**2)/mant/manu &
                 -3d0/8d0 * (mant**2+manu**2)/mans**2)

  ! g + q(qb) -> g + q(qb) : triple gluon at t-channel
  dis7tq = 0d0
  do i=1,nf
    dis7tq = dis7tq + (pdf1(0)*pdf2(+i)*jff3(0)*jff4(+i) &
                      +pdf1(0)*pdf2(-i)*jff3(0)*jff4(-i))
  enddo
  sig7tq = as*as * ((manu**2+mans**2)/mant**2 - 4d0/9d0*(mans**2+manu**2)/mans/manu)
  
  ! q(qb) + g -> q(qb) + g : triple gluon at t-channel
  dis7tg = 0d0
  do i=1,nf
    dis7tg = dis7tg + (pdf1(+i)*pdf2(0)*jff3(+i)*jff4(0) &
                      +pdf1(-i)*pdf2(0)*jff3(-i)*jff4(0))
  enddo
  sig7tg = as*as * ((manu**2+mans**2)/mant**2 - 4d0/9d0*(mans**2+manu**2)/mans/manu)

  ! g + q(qb) -> q(qb) + g : triple gluon at u-channel
  dis7ug = 0d0
  do i=1,nf
    dis7ug = dis7ug + (pdf1(0)*pdf2(+i)*jff3(+i)*jff4(0) &
                      +pdf1(0)*pdf2(-i)*jff3(-i)*jff4(0))
  enddo
  sig7ug = as*as * ((mant**2+mans**2)/manu**2 - 4d0/9d0*(mans**2+mant**2)/mans/mant)

  ! q(qb) + g -> g + q(qb) : triple gluon at u-channel
  dis7uq = 0d0
  do i=1,nf
    dis7uq = dis7uq + (pdf1(+i)*pdf2(0)*jff3(0)*jff4(+i) &
                      +pdf1(-i)*pdf2(0)*jff3(0)*jff4(-i))
  enddo
  sig7uq = as*as * ((mant**2+mans**2)/manu**2 - 4d0/9d0*(mans**2+mant**2)/mans/mant)

  ! g + g -> g + g :
  dis8 = pdf1(0)*pdf2(0)*jff3(0)*jff4(0)
  sig8 = as*as * 9d0/2d0*(3d0-mant*manu/mans**2-mans*manu/mant**2-mans*mant/manu**2)

  fxn = fxn + dis1t * sig1t &
            + dis1u * sig1u &
            + dis2  * sig2 &
            + dis3  * sig3 &
            + dis4t * sig4t &
            + dis4u * sig4u &
!           + dis5  * sig5 &
            + dis6  * sig6 &
!           + dis7tq * sig7tq &
            + dis7tg * sig7tg &
!           + dis7uq * sig7uq &
            + dis7ug * sig7ug !&
!           + dis8  * sig8

  fxn = fxn * twoPI * ptTrig * x1*x2 / mans**2
  fxn = fxn / (yAmax-yAmin) * gev2barn / nano

  tempfunc = fxn

  ! separating quark and gluon contributions

  ! assume jet(3) is gluon

  ptTrig = ptAsso !+ 50d0

  x1 = ptTrig / CME * (exp(+yTrig) + exp(+yAsso))
  x2 = ptTrig / CME * (exp(-yTrig) + exp(-yAsso))
  if(x1.le.0d0 .or. x1.gt.1d0) return
  if(x2.le.0d0 .or. x2.gt.1d0) return

  mans = +x1*x2*CME*CME
  mant = -x1*CME*ptTrig*exp(-yTrig) 
  manu = -x2*CME*ptTrig*exp(+yTrig)

  mufac = ptTrig !/ 2d0 ! scale sensitivity multiplier
  muren = mufac
  as = alphas(muren)

  ! set distributions
  ! mass ordered parton flavours
  ! i: -6  -5  -4  -3  -2  -1  0  +1  +2  +3  +4  +5  +6
  !    tb  bb  cb  sb  db  ub  g   u   d   s   c   b   t

  ! initialize
  pdf1=1d0; pdf2=1d0; jff3=1d0; jff4=1d0
  
  do i = -nf,+nf
    pdf1(i) = CT18pdf(i,dble(x1),dble(mufac))
    pdf2(i) = CT18pdf(i,dble(x2),dble(mufac))
  enddo
  
  ! q(qb) + q'(qb') -> q(qb) + q'(qb') : only t-channel
  dis1t = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis1t = dis1t + (pdf1(+i)*pdf2(+j)*jff3(+i)*jff4(+j) &
                    +pdf1(+i)*pdf2(-j)*jff3(+i)*jff4(-j) &
                    +pdf1(-i)*pdf2(+j)*jff3(-i)*jff4(+j) &
                    +pdf1(-i)*pdf2(-j)*jff3(-i)*jff4(-j))
  endif
  enddo
  enddo
  sig1t = as*as * 4d0/9d0 * (mans**2 + manu**2)/mant**2

  ! q(qb) + q'(qb') -> q'(qb') + q(qb) : only u-channel
  dis1u = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis1u = dis1u + (pdf1(+i)*pdf2(+j)*jff3(+j)*jff4(+i) &
                    +pdf1(+i)*pdf2(-j)*jff3(-j)*jff4(+i) &
                    +pdf1(-i)*pdf2(+j)*jff3(+j)*jff4(-i) &
                    +pdf1(-i)*pdf2(-j)*jff3(-j)*jff4(-i))
  endif
  enddo
  enddo
  sig1u = as*as * 4d0/9d0 * (mans**2 + mant**2)/manu**2

  ! q(qb) + qb(q) -> q'(qb') + qb'(q') : only s-channel
  dis2 = 0d0
  do i = 1,nf
  do j = 1,nf
  if(i.ne.j) then
    dis2 = dis2 + (pdf1(+i)*pdf2(-i)*jff3(+j)*jff4(-j) &
                  +pdf1(+i)*pdf2(-i)*jff3(-j)*jff4(+j) &
                  +pdf1(-i)*pdf2(+i)*jff3(+j)*jff4(-j) &
                  +pdf1(-i)*pdf2(+i)*jff3(-j)*jff4(+j))
  endif
  enddo
  enddo
  sig2 = as*as * 4d0/9d0 * (mant**2 + manu**2)/mans**2

  ! q(qb) + q(qb) -> q(qb) + q(qb) : t- and u-channel
  dis3 = 0d0
  do i = 1,nf
    dis3 = dis3 + (pdf1(+i)*pdf2(+i)*jff3(+i)*jff4(+i) &
                  +pdf1(-i)*pdf2(-i)*jff3(-i)*jff4(-i))
  enddo
  sig3 = as*as * 4d0/9d0 * ((mans**2+manu**2)/mant**2 &
                           +(mans**2+mant**2)/manu**2 &
                           -2d0/3d0*mans**2/manu/mant)

  ! q(qb) + qb(q) -> q(qb) + qb(q) : t- and s-channel
  dis4t = 0d0
  do i = 1,nf
    dis4t = dis4t + (pdf1(+i)*pdf2(-i)*jff3(+i)*jff4(-i) &
                    +pdf1(-i)*pdf2(+i)*jff3(-i)*jff4(+i))
  enddo
  sig4t =  as*as * 4d0/9d0 * ((mans**2+manu**2)/mant**2 &
                             +(mant**2+manu**2)/mans**2 &
                             -2d0/3d0*manu**2/mans/mant)

  ! q(qb) + qb(q) -> qb(q) + q(qb) : u- and s-channel
  dis4u = 0d0
  do i = 1,nf
    dis4u = dis4u + (pdf1(+i)*pdf2(-i)*jff3(-i)*jff4(+i) &
                    +pdf1(-i)*pdf2(+i)*jff3(+i)*jff4(-i))
  enddo
  sig4u =  as*as * 4d0/9d0 * ((mans**2+mant**2)/manu**2 &
                             +(mant**2+manu**2)/mans**2 &
                             -2d0/3d0*mant**2/mans/manu)

  ! q(qb) + qb(q) -> g + g : scatter(t-,u-channel) + annihlate(s-channel)
  dis5 = 0d0
  do i = 1,nf
    dis5 = dis5 + (pdf1(+i)*pdf2(-i)*jff3(0)*jff4(0) &
                  +pdf1(-i)*pdf2(+i)*jff3(0)*jff4(0))
  enddo
  sig5 = as*as * (32d0/27d0 * (mant**2+manu**2)/mant/manu &
                 -8d0 / 3d0 * (mant**2+manu**2)/mans**2)

  ! g + g -> q(qb) + qb(q) : scatter(t-,u-channel) + creation(s-channel)
  dis6 = 0d0
  do j=1,nf
    dis6 = dis6 + (pdf1(0)*pdf2(0)*jff3(+j)*jff4(-j) &
                  +pdf1(0)*pdf2(0)*jff3(-j)*jff4(+j))
  enddo
  sig6 = as*as * (1d0/6d0 * (mant**2+manu**2)/mant/manu &
                 -3d0/8d0 * (mant**2+manu**2)/mans**2)

  ! g + q(qb) -> g + q(qb) : triple gluon at t-channel
  dis7tq = 0d0
  do i=1,nf
    dis7tq = dis7tq + (pdf1(0)*pdf2(+i)*jff3(0)*jff4(+i) &
                      +pdf1(0)*pdf2(-i)*jff3(0)*jff4(-i))
  enddo
  sig7tq = as*as * ((manu**2+mans**2)/mant**2 - 4d0/9d0*(mans**2+manu**2)/mans/manu)
  
  ! q(qb) + g -> q(qb) + g : triple gluon at t-channel
  dis7tg = 0d0
  do i=1,nf
    dis7tg = dis7tg + (pdf1(+i)*pdf2(0)*jff3(+i)*jff4(0) &
                      +pdf1(-i)*pdf2(0)*jff3(-i)*jff4(0))
  enddo
  sig7tg = as*as * ((manu**2+mans**2)/mant**2 - 4d0/9d0*(mans**2+manu**2)/mans/manu)

  ! g + q(qb) -> q(qb) + g : triple gluon at u-channel
  dis7ug = 0d0
  do i=1,nf
    dis7ug = dis7ug + (pdf1(0)*pdf2(+i)*jff3(+i)*jff4(0) &
                      +pdf1(0)*pdf2(-i)*jff3(-i)*jff4(0))
  enddo
  sig7ug = as*as * ((mant**2+mans**2)/manu**2 - 4d0/9d0*(mans**2+mant**2)/mans/mant)

  ! q(qb) + g -> g + q(qb) : triple gluon at u-channel
  dis7uq = 0d0
  do i=1,nf
    dis7uq = dis7uq + (pdf1(+i)*pdf2(0)*jff3(0)*jff4(+i) &
                      +pdf1(-i)*pdf2(0)*jff3(0)*jff4(-i))
  enddo
  sig7uq = as*as * ((mant**2+mans**2)/manu**2 - 4d0/9d0*(mans**2+mant**2)/mans/mant)

  ! g + g -> g + g :
  dis8 = pdf1(0)*pdf2(0)*jff3(0)*jff4(0)
  sig8 = as*as * 9d0/2d0*(3d0-mant*manu/mans**2-mans*manu/mant**2-mans*mant/manu**2)

  fxn = fxn + dis5  * sig5 &
            + dis7tq * sig7tq &
            + dis7uq * sig7uq &
            + dis8  * sig8

  fxn = fxn * twoPI * ptTrig * x1*x2 / mans**2
  fxn = fxn / (yAmax-yAmin) * gev2barn / nano

  tempfunc = tempfunc + fxn

  fxn = tempfunc

  accevent = accevent + 1
  return

contains

  pure function ampA(s,t,u) result(res)
  real(wp), intent(in) :: s,t,u
  real(wp) :: res
  res = 4d0/9d0 * (s*s + t*t)/u/u
  return
  end function ampA

end function fxn
