program main
  use time
  use sysio
  use histloop
  use vint
  use eventcounter
  use optflags
  use kinematics
  use qcd, only : setqcd
  implicit none
  character(len=40), parameter :: pdfgrid = './grids/i2Tn2.00.pds'

  ! ************************************
  ! initialization procedures
  ! ************************************
  call time_start
  call readinput
  call setCT18(pdfgrid)
  call setqcd(5,1)
  call setloop(ptAmin,ptAmax,48,.true.)
  
  ! ************************************
  ! looping histogram
  ! ************************************
  tot_xsec = 0d0;
  write(*,*) "|---L---|---M---|---R---|----result----|----error-----|--eff--|"
  do id = 1,nd
    call setbin

    limits(1) = yTmin;          limits(ndimn+1) = yTmax
    limits(2) = yAmin;          limits(ndimn+2) = yAmax
    limits(3) = dL;             limits(ndimn+3) = dR
    if(quench_opt.eq.2) then
      limits(4) = 0d0;             limits(ndimn+4) = 50d0
      ! limits(5) = 0d0;             limits(ndimn+5) = 50d0   ! testing
      ! limits(6) = 0d0;             limits(ndimn+6) = 50d0   ! testing
    endif
    
    ! ************************************
    ! vegas integration
    ! ************************************
    ! warm-up run
    initial = -1
    docount = .false.
    call vegas(limits(1:2*ndimn),fxn,initial,npt1,itn1,prnt,intres,stddev,chisq)
    ! actual run
    initial = +1
    docount = .true.
    accevent = 0;  totevent = 0
    call vegas(limits(1:2*ndimn),fxn,initial,npt2,itn2,prnt,intres,stddev,chisq)
    eff = real(accevent,wp)/real(totevent,wp)*100d0
    write(*,'(3(f8.2),2(es15.5),f8.2)') dL,dM,dR,intres/bin,stddev,eff
  enddo
  write(*,'(a,es15.5)') "total cross-section: ",tot_xsec/itn2

  ! ************************************
  ! setting energy loss options
  ! ************************************
  call time_stop
  call print_time
end program main



function fxn(dx,wgt)
  use phyconst
  use eventcounter
  use kinematics
  use optflags
  use qcd, only: alphas,nf
  use bdmps
  use multi
  implicit none
  double precision, dimension(:), intent(in) :: dx
  double precision, intent(in) :: wgt
  double precision :: fxn
  double precision, dimension(-6:+6) :: pdf1,pdf2
  integer :: i,j
  real(wp) :: dis1,dis2,dis3,dis4,dis5,dis6,dis7tq,dis7tg,dis7uq,dis7ug,dis8
  real(wp) :: sig1,sig2,sig3,sig4,sig5,sig6,sig7tq,sig7tg,sig7uq,sig7ug,sig8
  real(wp) :: fxnq,fxng
  real(wp) :: elossQ,elossG
  real(wp) :: eps
  real(wp) :: Rcone
  real(wp) :: Qmed
  real(wp) :: eps1,eps2,eps3
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
  if(quench_opt.eq.2) then
    eps = dx(4)  !testing
    ! eps1 = dx(4)  ! testing
    ! eps2 = dx(5)  ! testing
    ! eps3 = dx(6)  ! testing
  endif

  ! ************************************
  ! setting energy loss options
  ! ************************************
  elossQ = 0d0
  elossG = 0d0
  if(quench_opt.eq.1) then
    ! coh:10d0; decoh:4.3d0;
    elossQ = 10d0
    elossG = CA/CF * elossQ
  elseif(quench_opt.eq.2) then
    ! coh:15d0; decoh:3d0;
    wcq = 3d0
    wcg = CA/CF * wcq
    elossQ = eps !eps1+eps2+eps3 !eps
    elossG = eps !eps1+eps2+eps3 !eps
  endif
  Rcone = 0.4d0
  Qmed = 0.5d0

  ! ************************************
  ! quark jet differential cross-section
  ! ************************************
  ! get parent pt for coh or decoh case
  if(decoherent .and. quench_opt.ne.0) then
    pTq = getpt(pTjet,Rcone,Qmed,elossG,1)     ! elossG not elossQ
    !ptq = pTjet + eps1+eps2+eps3
  else
    pTq = pTjet + elossQ
  endif

  ! this is a test with avgn cut, not yet tested
  ! avgn = nMLL2(pTq*Rcone,Qmed)
  ! if(avgn.le.1.5d0) then
  !   fxnq = fxnq * De(wcq,eps,1)
  ! else
  !   fxnq = fxnq * De(wcg,eps,0)
  ! endif

  ! momentum fraction, Mandelstam and scales
  x1 = pTq / CME * (exp(+yTrig) + exp(+yAsso))
  x2 = pTq / CME * (exp(-yTrig) + exp(-yAsso))
  if(x1.le.0d0 .or. x1.gt.1d0) return
  if(x2.le.0d0 .or. x2.gt.1d0) return
  mans = +x1*x2*CME*CME
  mant = -x1*CME*pTq*exp(-yTrig)
  manu = -x2*CME*pTq*exp(+yTrig)
  mufac = pTq
  muren = mufac
  as = alphas(muren)  !as = CT18Alphas(muren)

  ! PDFs
  pdf1 = 1d0;         pdf2 = 1d0
  do i = -nf, +nf
    pdf1(i) = CT18Pdf(i,x1,mufac)
    pdf2(i) = CT18Pdf(i,x2,mufac)
  enddo

  ! q + q'    -> q + q'
  ! q + qb'   -> q + qb'
  ! qb + q'   -> qb + q'
  ! qb + qb'  -> qb + qb' : t-channel
  ! q + q'    -> q' + q
  ! q + qb'   -> qb' + q
  ! qb + q'   -> q' + qb
  ! qb + qb'  -> qb' + qb : u-channel
  dis1 = 0d0
  do i = 1,nf
    do j = 1,nf
      if(i.ne.j) then
        dis1 = dis1 + (pdf1(+i)*pdf2(+j) &
                      +pdf1(+i)*pdf2(-j) &
                      +pdf1(-i)*pdf2(+j) &
                      +pdf1(-i)*pdf2(-j))
      endif
    enddo
  enddo
  sig1 = ampA(mans,manu,mant) + ampA(mans,mant,manu)

  ! q + qb    -> q' + qb'
  ! qb + q    -> q' + qb' : s-channel
  ! q + qb    -> qb' + q'
  ! qb + q    -> qb' + q' : s-channel
  dis2 = 0d0
  do i = 1,nf
    do j = 1,nf
      if(i.ne.j) then
        dis2 = dis2 + (pdf1(+i)*pdf2(-i) &
                      +pdf1(-i)*pdf2(+i) )
      endif
    enddo
  enddo
  sig2 = ampA(mant,manu,mans) + ampA(manu,mant,mans)

  ! q + q     -> q + q
  ! qb + qb   -> qb + qb : identical final state
  dis3 = 0d0
  do i = 1,nf
    dis3 = dis3 + (pdf1(+i)*pdf2(+i) &
                  +pdf1(-i)*pdf2(-i))
  enddo
  sig3 = ampB(mans,mant,manu)

  ! q + qb    -> q + qb
  ! qb + q    -> qb + q : t- and s-channel
  ! q + qb    -> qb + q
  ! qb + q    -> q + qb : u- and s-channel
  dis4 = 0d0
  do i = 1,nf
    dis4 = dis4 + (pdf1(+i)*pdf2(-i) &
                  +pdf1(-i)*pdf2(+i))
  enddo
  sig4 = ampB(manu,mans,mant) + ampB(mant,mans,manu) 

  ! g + g     -> q + qb : t-,u-,s-channel
  ! g + g     -> qb + q : t-,u-,s-channel
  dis6 = 0d0
  do j=1,nf
    dis6 = dis6 + (pdf1(0)*pdf2(0))
  enddo
  sig6 = 27d0/32d0 * (ampC(mant,manu,mans) + ampC(manu,mant,mans))
  
  ! q + g     -> q + g
  ! qb + g    -> qb + g : triple gluon at t-channel
  dis7tq = 0d0
  do i=1,nf
    dis7tq = dis7tq + (pdf1(+i)*pdf2(0) &
                      +pdf1(-i)*pdf2(0))
  enddo
  sig7tq = -9d0/4d0 * ampC(mans,manu,mant)

  ! g + q     -> q + g
  ! g + qb    -> qb + g : triple gluon at u-channel
  dis7uq = 0d0
  do i=1,nf
    dis7uq = dis7uq + (pdf1(0)*pdf2(+i) &
                      +pdf1(0)*pdf2(-i))
  enddo
  sig7uq = -9d0/4d0 * ampC(mans,mant,manu)

  ! sum channels
  fxnq = fxnq + dis1  * sig1 &
              + dis2  * sig2 &
              + dis3  * sig3 &
              + dis4  * sig4 &
              + dis6  * sig6 &
              + dis7tq * sig7tq &
              + dis7uq * sig7uq

  ! common factor and unit conversion
  fxnq = fxnq * as*as * twoPI * ptq * x1*x2 / mans**2
  fxnq = fxnq / (yAmax-yAmin) * gev2barn / nano

  ! for BDMPS D(e)
  if(quench_opt.eq.2 .and. .not.decoherent) then
    fxnq = fxnq * De(wcq,eps,1)    !testing
    ! fxnq = fxnq * De(wcq,eps1,1) * De(wcq,eps2,1) * De(wcq,eps3,1)
  elseif(quench_opt.eq.2 .and. decoherent) then
    fxnq = fxnq * De(wcg,eps,0)
    ! fxnq = fxnq * De(wcg,eps1,0) * De(wcg,eps2,0) * De(wcg,eps3,0)
  endif

  ! for average dE or n
  if(xsec_fac.eq.0) then
    fxnq = fxnq
  elseif(xsec_fac.eq.1) then
    fxnq = fxnq * (ptq-pTjet)
  elseif(xsec_fac.eq.2) then
    fxnq = fxnq * nMLL2(ptq*Rcone,Qmed)
  endif

  ! ************************************
  ! gluon jet differential cross-section
  ! ************************************
  ! get parent pt for coh or decoh case
  if(decoherent .and. quench_opt.ne.0) then
    pTg = getpt(pTjet,Rcone,Qmed,elossG,0)
    ! pTg = pTjet + eps1+eps2+eps3
  else
    pTg = pTjet + elossG
  endif

  ! momentum fraction, Mandelstam and scales
  x1 = ptg / CME * (exp(+yTrig) + exp(+yAsso))
  x2 = ptg / CME * (exp(-yTrig) + exp(-yAsso))
  if(x1.le.0d0 .or. x1.gt.1d0) return
  if(x2.le.0d0 .or. x2.gt.1d0) return
  mans = +x1*x2*CME*CME
  mant = -x1*CME*ptg*exp(-yTrig) 
  manu = -x2*CME*ptg*exp(+yTrig)
  mufac = ptg 
  muren = mufac
  as = alphas(muren)

  ! PDFs
  pdf1=1d0;     pdf2=1d0
  do i = -nf, +nf
    pdf1(i) = CT18Pdf(i,x1,mufac)
    pdf2(i) = CT18Pdf(i,x2,mufac)
  enddo
  
  ! q + qb    -> g + g
  ! qb + q    -> g + g : identical final state
  dis5 = 0d0
  do i = 1,nf
    dis5 = dis5 + (pdf1(+i)*pdf2(-i) &
                  +pdf1(-i)*pdf2(+i))
  enddo
  sig5 = 6d0 * ampC(mant,manu,mans)

  ! g + q     -> g + q
  ! g + qb    -> g + qb : triple gluon at t-channel
  dis7tg = 0d0
  do i=1,nf
    dis7tg = dis7tg + (pdf1(0)*pdf2(+i) &
                      +pdf1(0)*pdf2(-i))
  enddo
  sig7tg = -9d0/4d0 * ampC(mans,manu,mant)

  ! q + g     -> g + q
  ! qb + g    -> g + qb : triple gluon at u-channel
  dis7ug = 0d0
  do i=1,nf
    dis7ug = dis7ug + (pdf1(+i)*pdf2(0) &
                      +pdf1(-i)*pdf2(0))
  enddo
  sig7ug = -9d0/4d0 * ampC(mans,mant,manu)

  ! g + g -> g + g : identical final state
  dis8 = pdf1(0)*pdf2(0)
  sig8 = 9d0/2d0 * ampD(mans,mant,manu)

  ! sum channels
  fxng = fxng + dis5  * sig5 &
              + dis7tg * sig7tg &
              + dis7ug * sig7ug &
              + dis8  * sig8

  ! common factor and unit conversion
  fxng = fxng * as*as * twoPI * ptg * x1*x2 / mans**2
  fxng = fxng / (yAmax-yAmin) * gev2barn / nano

  ! for BDMPS D(e)
  if(quench_opt.eq.2) then
    fxng = fxng * De(wcg,eps,0)    !testing
    ! fxng = fxng * De(wcg,eps1,0) * De(wcg,eps2,0) * De(wcg,eps3,0)
    ! fxng = fxng * De(wcg,eps1,0) * De(wcg,eps2,0) * De(wcg,eps3,0)
  endif

  ! for average dE or n
  if(xsec_fac.eq.0) then
    fxng = fxng
  elseif(xsec_fac.eq.1) then
    fxng = fxng * (ptg-pTjet)
  elseif(xsec_fac.eq.2) then
    fxng = fxng * nMLL2(ptg*Rcone,Qmed)
  endif

  ! ************************************
  ! sum quark + gluon jet cross-section
  ! ************************************
  fxn = fxnq + fxng

  ! ************************************
  ! event number and tot_xsec counter
  ! ************************************
  accevent = accevent + 1
  if(docount)  tot_xsec = tot_xsec + fxn*wgt
  return

end function fxn
