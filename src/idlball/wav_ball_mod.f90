
module wav_ball_mod

  use alm_tools
  use healpix_types
  use common_module
  use sincm
  USE fitstools

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  integer :: recfac_count=1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: recfac

  REAL(DP),parameter :: pifac=2d0
  REAL(DP),parameter :: argpi=pifac*3.1415926535897931159980
   REAL(DP),parameter ::  OVFLOW=exp(log(1.0d30))
   REAL(DP),parameter ::  UNFLOW=exp(log(1.0d-30))
   REAL(DP),parameter ::  logOVFLOW=log(1.0d30)


  INTEGER(I4B) :: nf2,feedback=2,intflag=2,nbasis
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f2



contains

  subroutine calc_f2

    IMPLICIT NONE

    INTEGER(I4B) :: i
    REAL(DP) :: x

    nf2=10000000
    ALLOCATE(f2(0:nf2-1))

    x=-1.d0
    !f2(0)=exp(-1.d0/(1.d0-x**2))
    f2(0)=0.d0
    do i=1,nf2-1
       x=2.d0/REAL(nf2,DP)*i-1.d0
       f2(i)=f2(i-1)+exp(-1.d0/(1.d0-x**2))
    enddo

    f2=f2/f2(nf2-1) ! this is phi(u) = integral[f(t)dt,-1,u]/integral[f(t)dt,-1,1]

  end subroutine calc_f2

!====================
  real(DP) function snrv(nn,v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn
    INTEGER(I4B), optional ::nvj
    REAL(DP) :: bb,rvj,v,dr,y,norm

    y = nn*argpi
    norm = sqrt(2d0)

    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)

!print*, 'nvj,dr: ',nvj,dr
   !  if (nn < bb**real(j,DP)) then 
        rvj=v*dr
        snrv=sqrt(dr)*norm*sin(y*rvj) !sqrt(dr)*  !/3d0 !sqrt(ljkv(v,j,bb))*
   !  else
   !    snrv=0d0
   ! endif


    return

  end function snrv

!===============================

  real(DP) function snr(nn,v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn
    INTEGER(I4B), optional ::nvj
    REAL(DP) :: bb,rvj,v,dr

    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)

    rvj=v*dr
    snr=sqrt(2d0)*argpi*nn*sinc(argpi*nn*rvj) 


    return

  end function snr


!===============================
  real(DP) function ljkv(v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B), optional ::nvj
    INTEGER(I4B) :: j,i,nj,nn
    REAL(DP) :: bb,rvj,v,dr

    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)
    rvj=v*dr
    ljkv = dr*rvj**2

    return

  end function ljkv



!=============================
  real(DP) function snrv_sincos(itype,nn,v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn,itype
    INTEGER(I4B), optional ::nvj
    REAL(DP) :: bb,rvj,v,dr,y,norm

    y = nn*argpi
    norm = sqrt(2d0)


    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)

    rvj=v*dr
    if (itype==2) snrv_sincos=sqrt(dr)*norm*sin(y*rvj)

    if (nn==0) norm=1d0
    if (itype==3) snrv_sincos=sqrt(dr)*norm*cos(y*rvj)

    return

  end function snrv_sincos
!================

  real(DP) function snrv_exp(itype,nn,v,j,bb,nvj)
    !this function is to do INVERSE FFT
    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn,itype
    INTEGER(I4B) ::nvj
    REAL(DP) :: bb,rvj,v,dr,y,norm

    y = nn*2d0*pi
    dr=1d0/dble(nvj)
    rvj=v*dr
    norm = sqrt(2d0)


    if (nn==0) norm=1d0
    if (itype==1) snrv_exp=dble(nvj)*norm*cos(y*rvj)
    if (itype==2) snrv_exp=dble(nvj)*norm*sin(y*rvj)

    if (itype<1 .or. itype>2) stop 'snrv_exp - itype must be 1 or 2'

    return

  end function snrv_exp


!====================
  real(DP) function snrv_sin(nn,v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn
    INTEGER(I4B), optional ::nvj
    REAL(DP) :: bb,rvj,v,dr,y,norm

    y = nn*argpi
    norm = sqrt(2d0)

    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)

    rvj=v*dr
    snrv_sin=sqrt(dr)*norm*sin(y*rvj) 

    return

  end function snrv_sin
!====================

  real(DP) function snrv_cos(nn,v,j,bb,nvj)

    IMPLICIT NONE

    INTEGER(I4B) :: j,i,nj,nn
    INTEGER(I4B), optional ::nvj
    REAL(DP) :: bb,rvj,v,dr,y,norm

    y = nn*argpi
    norm = sqrt(2d0)
    if (nn==0) norm=1d0

    dr=1d0/bb**real(j,DP)
    if (present(nvj)) dr=1d0/dble(nvj)

    rvj=v*dr
    snrv_cos=sqrt(dr)*norm*cos(y*rvj) 
    

    return

  end function snrv_cos

!===============================
  ! real(DP) function snr(nn,nshell,ir)

  !   IMPLICIT NONE
  !   INTEGER(I4B) :: nn,nshell,ir
  !   REAL(DP) :: r


  !   call r2dr_vec_ir(nshell , 0d0, 1d0, ir, r)
  !   snr=PI*nn*sinc(PI*nn*r)

  !   return

  ! end function snr


  !==============================================

  subroutine calc_phiphi(phiphi,x,aa)

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: phiphi
    REAL(DP), INTENT(IN) :: x,aa
    INTEGER(I4B) :: ii
    phiphi=0.d0

    if ((x.gt.-1.d0).and.(x.lt.-1d0/aa)) then
       !ii=NINT(REAL(n,DP)/2.d0*(1.d0+3.d0+4.d0*x))
       stop 'STOP, x lt -1/aa, need formula for this case!'
       if (ii.gt.(nf2-1)) ii=nf2-1
       phiphi=f2(ii)
    endif

    if ((x.ge.-1d0/aa).and.(x.le.1d0/aa)) phiphi=1.d0

    if ((x.le.1.d0).and.(x.gt.1d0/aa)) then
       !ii=NINT(REAL(n,DP)/2.d0*(1.d0+4.d0*x-3.d0))
       ii=NINT(REAL(nf2,DP)*(1d0-(aa/(aa-1d0)*(x-1d0/aa))))  !!!the figure 2 in here has disappeared because we try to find an index which runs from 0 to n-1 instead of a number between -1 and 1.
       if (ii.gt.(nf2-1)) ii=nf2-1
       phiphi=f2(ii)
    endif

  end subroutine calc_phiphi

  !==============================================

  subroutine calc_phi(phi,x,aa)

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: phi
    REAL(DP), INTENT(IN) :: x,aa
    REAL(DP) :: xx,phi1,phi2

    xx=x/aa
    call calc_phiphi(phi1,xx,aa)
    xx=x
    call calc_phiphi(phi2,xx,aa)

    !print*,'xx',xx,phi1, phi2

    phi=phi1-phi2
    !print*,phi1,phi2
  end subroutine calc_phi

  !==============================================


  !==============================================

  subroutine calc_gln_stand(j0,nj,lmax,nngrid,gln,aa,gln4fft)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,nngrid
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,0:nngrid,j0:j0+nj-1), INTENT(OUT) :: gln
    INTEGER(I4B) :: j,eln0,eln1,l,n,isfft,nmax
    REAL(DP) :: xx,phi,eln,elnmax,norm
    INTEGER(I4B),optional::gln4fft


    isfft=0
    if (present(gln4fft)) isfft=1

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    nmax=nngrid
    if (isfft==1) nmax=(nngrid+1)/2

    gln=0.d0

    elnmax = real(nmax**2,DP)+REAL(lmax,DP)*(REAL(lmax,DP)+1)

    do j=j0,j0+nj-1
       if (feedback>2) print*,'gln j=',j

       eln0=NINT(aa**(j-1))
       eln1=NINT(aa**(j+1))

       if (eln0.gt.elnmax) print*, 'WARNNING: first eln for j greater than elnmax!!!, eln0, elnmax: ', eln0, elnmax
       if (eln1.gt.elnmax) print*, 'WARNNING: last eln for j greater than elnmax!!!, eln1m elnmax: ', eln0, elnmax

       do n=0,nmax

          do l=0,lmax

             eln = real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1)


             xx=sqrt(real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1))/aa**REAL(j,DP)

             call calc_phi(phi,xx,aa)
             gln(l,n,j)=phi

             !print*,'j, n, l, aa^j, sqrt(n2+l(l+1)),xx:',j,n,l,aa**REAL(j,DP),xx*aa**REAL(j,DP)
          enddo
       enddo
    enddo

    do n=0,nmax
       do l=0,lmax
          norm=sum(gln(l,n,:))
          if (norm.gt.0d0) gln(l,n,:)=gln(l,n,:)/norm
       enddo
    enddo
    gln=sqrt(gln)

    !if isfft then we should store gln such that the radial frequency ordering is 
    !0:nngrid = 1,2 .. N/2, -N/2-1,..-2, -1
    if (isfft==1) then
       do n=nmax+1,nngrid
          gln(0:lmax,n,j0:j0+nj-1) = gln(0:lmax,(nngrid+1)-n,j0:j0+nj-1)
       enddo
    endif

    if ((aa.lt.2d0).and.(j0.eq.0)) gln(0,0:nngrid,0)=1d0

  end subroutine calc_gln_stand

  !==============================================

  subroutine calc_gln_mex(j0,nj,lmax,nmax,gln,aa,p)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p,nmax
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,0:nmax,j0:j0+nj-1), INTENT(OUT) :: gln
    INTEGER(I4B) :: j,l0,l1,l,n


    gln=0.d0

    do j=j0,j0+nj-1
       do n=0,nmax
          do l=0,lmax
             gln(l,n,j)=(sqrt(real(n**2,DP)+REAL(l,DP)*(REAL(l,DP)+1))/aa**j)**(2d0*p)*exp(-REAL(l,DP)**2/aa**(2d0*j))
          enddo
       enddo
    enddo

  end subroutine calc_gln_mex

  !==============================================

  subroutine calc_gln(j0,nj,lmax,nmax,gln,aa,wavtyp,p,nside,smhw,gln4fft)

    INTEGER(I4B), INTENT(IN) :: j0,nj,lmax,p,nmax
    INTEGER(I4B), INTENT(IN), OPTIONAL :: nside,gln4fft
    REAL(DP), INTENT(IN) :: aa
    REAL(DP),  DIMENSION(0:lmax,0:nmax,j0:j0+nj-1), INTENT(OUT) :: gln
    REAL(DP),  DIMENSION(j0:j0+nj-1), INTENT(IN), OPTIONAL :: smhw
    CHARACTER(LEN=128), INTENT(IN) :: wavtyp

    if (present(gln4fft)) then
       if (trim(adjustl(wavtyp)).eq.'standard') call calc_gln_stand(j0,nj,lmax,nmax,gln,aa,gln4fft)
    else
       if (trim(adjustl(wavtyp)).eq.'standard') call calc_gln_stand(j0,nj,lmax,nmax,gln,aa)
    endif

    if (trim(adjustl(wavtyp)).eq.'mexican')  call calc_gln_mex(j0,nj,lmax,nmax,gln,aa,p)
    if ((trim(adjustl(wavtyp)).eq.'smhw').and.(present(nside)).and.(present(smhw)))  stop 'calc_swhm not implemented for needlet in a ball' !call calc_smhw(nside*4,lmax,nj,smhw,gl,.false.,.false.,wavtyp)
    if ((trim(adjustl(wavtyp)).eq.'smhw').and.((.not.present(nside)).or.(.not.present(smhw))))  stop 'SMHW, but not all necessary input to calc_gl!'
    if ((trim(adjustl(wavtyp)).ne.'standard').and.(trim(adjustl(wavtyp)).ne.'mexican').and.(trim(adjustl(wavtyp)).ne.'smhw')) stop 'wavelet type not specified in calc_gl'

  end subroutine calc_gln

  ! subroutine wav_convolve_ball(nside,lmax,nside2,lmax2,nscales,g_real_tot,map,wavmap,me)

  !   use alm_tools

  !   IMPLICIT NONE

  !   INTEGER(I4B), INTENT(IN) :: nside,lmax,nscales,nside2,lmax2
  !   INTEGER(I4B), INTENT(IN), OPTIONAL :: me
  !   REAL(SP), DIMENSION(0:nside**2*12-1), INTENT(IN) :: map
  !   REAL(SP), DIMENSION(0:nside2**2*12-1,0:nscales-1), INTENT(out) :: wavmap
  !   COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  !   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
  !   REAL(DP), DIMENSION(0:lmax2,0:nscales-1), INTENT(IN) :: g_real_tot
  !   REAL(DP), DIMENSION(1:2) :: z = 0.0d0

  !   if (present(me)) then
  !      if (me.eq.0) print*,'allocating'
  !   endif

  !   ALLOCATE(alm(1:1,0:lmax2,0:lmax2),w8ring(1:2*nside,1:1))
  !   w8ring=1d0

  !   if (present(me)) then
  !      if (me.eq.0) print*,'entering anamap'
  !   endif
  !   CALL map2alm(nside, lmax2, lmax2, map, alm, z, w8ring)  

  !   if (present(me)) then
  !      if (me.eq.0) print*,'entering synmap'
  !   endif

  !   CALL synmap_multi_ball(alm(1,:,:),wavmap,nside2,lmax2,g_real_tot,nscales)

  !   if (present(me)) then
  !      if (me.eq.0) print*,'done'
  !   endif

  !   DEALLOCATE(alm,w8ring)

  !   if (present(me)) then
  !      if (me.eq.0) print*,'deallocated'
  !   endif


  ! end subroutine wav_convolve_ball

!============================================

  subroutine synmap_multi_ball(almn,map,nside,nmax, lmax,bb,j0, gln,njmaps,nvmaps)

    use healpix_types
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: nside,lmax,njmaps, j0,nmax,nvmaps
    REAL(DP), DIMENSION(0:nside**2*12-1,0:nvmaps, 0:njmaps-1), INTENT(OUT) :: map
    REAL(DP), DIMENSION(0:lmax,0:nmax, j0:j0+njmaps-1), INTENT(IN) :: gln
    COMPLEX(DPC), DIMENSION(0:lmax,0:lmax,0:nmax), INTENT(IN) :: almn
    real(dp)::bb

    print*, 'calling alm2map_multiball ..'

    CALL alm2map_multi_ball(nside,nmax, lmax,lmax,almn,bb,j0,map,gln,njmaps,nvmaps)

  end subroutine synmap_multi_ball

!============================================

  subroutine alm2map_multi_ball(nsmax,nnmax,  nlmax, nmmax, almn,bb,j0,map, gln ,nmaps,nvmaps, der)

    use healpix_types
    IMPLICIT none
    INTEGER(I4B), INTENT(IN) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax, nnmax,j0,nvmaps
    INTEGER(I4B), INTENT(IN) :: nmmax,nmaps
    INTEGER(I4B), INTENT(IN), OPTIONAL :: der
    COMPLEX(DPC), INTENT(IN),  DIMENSION(0:nlmax,0:nmmax,0:nnmax) :: almn
    REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax, 0:nmaps-1) :: gln
    REAL(DP),     INTENT(OUT), DIMENSION(0:12*nsmax**2-1,0:nvmaps, 0:nmaps-1) :: map

    INTEGER(I4B) :: l, m, ith, indl, scalem, scalel, i         ! alm related
    INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
    integer(i4b):: iv, vv,in,jj,nv,nvj,nnmaxj,nlmaxj,nnj

    REAL(DP) :: cth, sth, dth1, dth2, dst1, OVFLOW, UNFLOW, logOVFLOW,bb
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) :: f2m, fm2, fl2, corfac,vval
    COMPLEX(DPC) :: factor

    CHARACTER(LEN=7), PARAMETER :: code = 'ALMN2MAP_BALL'
    COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almjq
    COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: b_north, b_south
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: b_n, b_s
    INTEGER(I4B) :: status

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: recfac
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    !=======================================================================

    nvj=nvmaps

    !     --- allocates space for arrays ---
    ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)    
    if (status /= 0) stop 'recfac'

    ALLOCATE(b_north(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_north'

    ALLOCATE(b_n(0:nvmaps-1,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_n'

    ALLOCATE(b_south(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_south'

    ALLOCATE(b_s(0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_s'

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) stop 'ring'

    ALLOCATE(almjq(0:nvmaps,0:nmaps-1, 0:nlmax, 0:nmmax),stat = status) 
    if (status /= 0) stop 'almjq'

    if (feedback>2) print*, 'l320: total(gln): ',sum(gln)
    if (feedback>2) print*, 'l321: total(almn): ',sum(almn)
    if (feedback>1) print*,'l322: alm2map_multi_ball l321 A '

    !     ------------ initiate arrays ----------------

    do m = 0, nmmax
       fm2 = DBLE(m) **2
       do l = m, nlmax
          fl2 = DBLE(l+1) **2
          recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       enddo
    enddo
    !     --------------------------------------------

    !     ----- set the whole map to zero ------
    map = 0.0
    !     --------------------------------------

    istart_north = 0
    istart_south = 12*nsmax**2
    OVFLOW=exp(log(1.0d30))
    UNFLOW=exp(log(1.0d-30))
    logOVFLOW=log(1.0d30)

    dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
    dth2 = 2.d0 / (3.d0*DBLE(nsmax))
    dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

    almjq=0


    do l=0,nlmax
       do jj=0,nmaps-1
          do iv=0,nvj-1

             vval= dble(iv)
             do in=0,nnmax                
                almjq(iv, jj,l,0:l) = almjq(iv, jj, l,0:l)+almn(l,0:l,in)*snrv(in, vval,j0+jj,bb,nvj)*gln(l,in, jj)                
             enddo !end in

          enddo  !end iv       
       enddo !end jj
    enddo !end l


    do ith = 1,2*nsmax      ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.d0  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
       endif
       lam_mm = sq4pi_inv ! lambda_00
       scalem=1

       !print*,'alm2map_multi_ball l325 B m loop start '

       b_north=0
       b_south=0

       b_n=0
       b_s=0

       if (mod(ith,50)==0) print*, 'ring number',ith

       do m = 0, nmmax


          f2m = 2.d0 * m

          !           ---------- l = m ----------
          par_lm = 1.d0  ! = (-1)^(l+m)
          if (m .ge. 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac=exp(scalem*logOVFLOW) !Actual lam_mm
          lam_lm = lam_mm*corfac/OVFLOW


          b_n(0:nvmaps-1,0:nmaps-1) = lam_lm *almjq(0:nvmaps-1,0:nmaps-1,m,m) 
          b_s(0:nvmaps-1,0:nmaps-1) = lam_lm *almjq(0:nvmaps-1,0:nmaps-1,m,m)
             
          ! do jj=0,nmaps-1
          !    !if (feedback>10) print*, 'm,jj,nvj: ',m,jj+j0,nvj
          !    do iv=0,nvj-1
          !       vval= dble(iv)
          !       !if (feedback>3 .and. jj==nmaps-1 .and. m==0 .and. ith==1) print*, 'r reconst for j=',jj,':',vval/dble(nvj)

          !       do in=0,nnmax
          !          if (abs(gln(m,in,jj))>1e-16) then
          !             factor = lam_lm * almn(m,m,in)
          !             !if (present(der)) factor = lam_lm * almn(m,m,in)/sth**der
                      
          !             b_n(iv, jj) = b_n(iv,jj) +  factor*snrv(in, vval,j0+jj,bb,nvj)*gln(m,in, jj)
          !             b_s(iv, jj) = b_s(iv,jj) +  factor*snrv(in, vval,j0+jj,bb,nvj)*gln(m,in, jj)
          !          endif
          !       enddo !end in

          !          !b_n(iv, jj) = b_n(iv,jj) +  factor*snr(in,nvmaps,iv)*gln(m,in,j0+jj)**2  !snrv(in, vval,j0+jj,bb)
          !          !b_s(iv, jj) = b_s(iv,jj) +  factor*snr(in,nvmaps,iv)*gln(m,in, j0+jj)**2 !snrv(in, vval,j0+jj,bb)
          !          !print*, 'bb, j0, in, jj, iv, vval, lam_lm, abs(almn), snrv, gln: ',bb, j0,in, jj, iv, vval, lam_lm,abs(almn(in,m,m)),snrv(in, vval,j0+jj,bb),gln(m,in, j0+jj) 
          !          !print*, 'in,iv,j0+jj,nvj,vval',iv,j0+jj,nvj,vval
          !    enddo  !end iv


          ! enddo !end jj

          !if (feedback>4) print*, 'l395: sum(b_n) ',sum(b_n)
          !if (feedback>4) print*, 'l396: sum(b_s) ',sum(b_s)

          !           ---------- l > m ----------
          lam_0 = 0.d0
          lam_1 = 1.d0 
          scalel=0
          a_rec = recfac(m,m)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 

             b_n(0:nvmaps-1,0:nmaps-1) = b_n(0:nvmaps-1,0:nmaps-1) + lam_lm *almjq(0:nvmaps-1,0:nmaps-1,l,m) 
             b_s(0:nvmaps-1,0:nmaps-1) = b_s(0:nvmaps-1,0:nmaps-1) + par_lm * lam_lm *almjq(0:nvmaps-1,0:nmaps-1,l,m)


             ! !!working version
             ! do jj=0,nmaps-1
             !    do iv=0,nvj-1
             !       vval= dble(iv)
             !       do in=0,nnmax
             !          if (abs(gln(l,in,jj))>1e-16) then
             !             factor = lam_lm * almn(l,m,in)
             !             !if (present(der)) factor = lam_lm * almn(l,m,in)/sth**der
                         
             !             b_n(iv, jj) = b_n(iv,jj) +          factor*snrv(in, vval,j0+jj,bb,nvj)*gln(l,in,jj)
             !             b_s(iv, jj) = b_s(iv,jj) + par_lm * factor*snrv(in, vval,j0+jj,bb,nvj)*gln(l,in,jj)
             !          endif
             !          ! b_n(iv, jj) = b_n(iv,jj) +          factor*snr(in,nvmaps,iv)*gln(l,in,j0+jj)**2 !snrv(in, vval,j0+jj,bb)
             !          ! b_s(iv, jj) = b_s(iv,jj) + par_lm * factor*snr(in,nvmaps,iv)*gln(l,in,j0+jj)**2 !snrv(in, vval,j0+jj,bb)
             !       enddo !end in
                   
             !    enddo !end iv
             ! enddo  !end jj
             
             !if (feedback>4) print*, 'l418: sum(b_n) ',sum(b_n)
             !if (feedback>4) print*, 'l419: sum(b_s) ',sum(b_s)

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = recfac(l,m)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) .gt. OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac=exp((scalel+scalem)*logOVFLOW)
             else
                if (abs(lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac=exp((scalel+scalem)*logOVFLOW)
                endif
             endif

          enddo !! end l

          do jj=0,nmaps-1
             !nvj=ceiling(bb**(j0+jj))
             do iv=0,nvj-1
                b_north(m,iv,jj) = b_n(iv,jj) 
                b_south(m,iv,jj) = b_s(iv,jj)
             enddo
          enddo


       enddo !end loop over m

       !if (feedback>2) print*, 'l448: sum(b_north) ',sum(b_north)
       !if (feedback>2) print*, 'l449: sum(b_south) ',sum(b_south)

       !print*,'alm2map_multi_ball l325 C m loop done! '
       !        ---------------------------------------------------------------

       !        obtains the negative m of b(m,theta) (= complex conjugate)


       !        ---------------------------------------------------------------
       !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta)
       !        ---------------------------------------------------------------

       do i=0,nmaps-1
          do vv=0,nvj-1
             call ring_synthesis_ball(nsmax,nlmax,nmmax,b_north(:,vv,i),nph,ring(:),kphi0)   ! north hemisph. + equator
             map(istart_north:istart_north+nph-1,vv,i) = ring(0:nph-1)
             !if (feedback>3) print*, 'l464: j, v, sum(b_north) ',j0+i, vv, sum(b_north(:,vv,i))
             !if (feedback>3) print*, 'l465: j, v, sum(b_south) ',j0+i, vv, sum(b_south(:,vv,i))
          enddo
       enddo
       istart_north = istart_north + nph

       !if (feedback>2) print*, 'l469: ith_ring north, sum(map(istart_north:istart_north+nph-1,:,:)) ',ith, sum(map(istart_north:istart_north+nph-1,:,:))


       !print*,'alm2map_multi_ball l325 north done, ith= ',ith

       istart_south = istart_south - nph
       if (ith .lt. 2*nsmax) then
          do i=0,nmaps-1
             do vv=0,nvj-1
                call ring_synthesis_ball(nsmax,nlmax,nmmax,b_south(:,vv,i),nph,ring(:),kphi0) ! south hemisph. w/o equat
                map(istart_south:istart_south+nph-1,vv,i) = ring(0:nph-1)
             enddo !end vv
          enddo
          !if (feedback>2) print*, 'l482: ith_ring south, sum(map(istart_south:istart_south+nph-1,:,:)) ',ith,sum(map(istart_south:istart_south+nph-1,:,:))

          !print*,'alm2map_multi_ball l325 south done, ith= ',ith
       endif

    enddo    ! loop on cos(theta)

    !if (feedback>1)  print*,'alm2map_multi_ball l325 E, multiball done! sum(map) = ',sum(map)
    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(recfac)
    DEALLOCATE(almjq)
    DEALLOCATE(b_north)
    DEALLOCATE(b_south)
    DEALLOCATE(b_n)
    DEALLOCATE(b_s)
    DEALLOCATE(ring)
    return
  end subroutine alm2map_multi_ball

!=======================================

  subroutine alm2map_sincos_ball(itype,nsmax,nnmax,  nlmax, nmmax, almn,bb,j0,map, gln ,nmaps,nvmaps, der)

    use healpix_types
    IMPLICIT none
    INTEGER(I4B), INTENT(IN) :: nsmax,itype
    INTEGER(I4B), INTENT(IN) :: nlmax, nnmax,j0,nvmaps
    INTEGER(I4B), INTENT(IN) :: nmmax,nmaps
    INTEGER(I4B), INTENT(IN), OPTIONAL :: der
    COMPLEX(DPC), INTENT(IN),  DIMENSION(0:nlmax,0:nmmax,0:nnmax) :: almn
    REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax, 0:nmaps-1) :: gln
    REAL(DP),     INTENT(OUT), DIMENSION(0:12*nsmax**2-1,0:nvmaps, 0:nmaps-1) :: map

    INTEGER(I4B) :: l, m, ith, indl, scalem, scalel, i         ! alm related
    INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
    integer(i4b):: iv, vv,in,jj,nv,nvj,nnmaxj,nlmaxj,nnj

    REAL(DP) :: cth, sth, dth1, dth2, dst1, OVFLOW, UNFLOW, logOVFLOW,bb
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) :: f2m, fm2, fl2, corfac,vval
    COMPLEX(DPC) :: factor

    CHARACTER(LEN=7), PARAMETER :: code = 'ALMN2MAP_BALL'
    COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almjq
    COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: b_north, b_south
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: b_n, b_s
    INTEGER(I4B) :: status

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: recfac
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    !=======================================================================

    nvj=nvmaps

    !     --- allocates space for arrays ---
    ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)    
    if (status /= 0) stop 'recfac'

    ALLOCATE(b_north(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_north'

    ALLOCATE(b_n(0:nvmaps-1,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_n'

    ALLOCATE(b_south(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_south'

    ALLOCATE(b_s(0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_s'

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) stop 'ring'

    ALLOCATE(almjq(0:nvmaps,0:nmaps-1, 0:nlmax, 0:nmmax),stat = status) 
    if (status /= 0) stop 'almjq'

    if (feedback>2) print*, 'l320: total(gln): ',sum(gln)
    if (feedback>2) print*, 'l321: total(almn): ',sum(almn)
    if (feedback>1) print*,'l322: alm2map_multi_ball l321 A '

    !     ------------ initiate arrays ----------------

    do m = 0, nmmax
       fm2 = DBLE(m) **2
       do l = m, nlmax
          fl2 = DBLE(l+1) **2
          recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       enddo
    enddo
    !     --------------------------------------------

    !     ----- set the whole map to zero ------
    map = 0.0
    !     --------------------------------------

    istart_north = 0
    istart_south = 12*nsmax**2
    OVFLOW=exp(log(1.0d30))
    UNFLOW=exp(log(1.0d-30))
    logOVFLOW=log(1.0d30)

    dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
    dth2 = 2.d0 / (3.d0*DBLE(nsmax))
    dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

    almjq=0

    if (itype==2) then
       do l=0,nlmax
          do jj=0,nmaps-1
             do iv=0,nvj-1
                
                vval= dble(iv)
                do in=0,nnmax                
                   almjq(iv, jj,l,0:l) = almjq(iv, jj, l,0:l) + almn(l,0:l,in)*snrv_sin(in, vval,j0+jj,bb,nvj)*gln(l,in, jj)   !sin part
                enddo !end in
                
             enddo  !end iv       
          enddo !end jj
       enddo !end l
    else
       do l=0,nlmax
          do jj=0,nmaps-1
             do iv=0,nvj-1
                
                vval= dble(iv)
                do in=0,nnmax                
                   almjq(iv, jj,l,0:l) = almjq(iv, jj, l,0:l) + almn(l,0:l,in)*snrv_cos(in, vval,j0+jj,bb,nvj)*gln(l,in, jj)   !cos part
                enddo !end in
                
             enddo  !end iv       
          enddo !end jj
       enddo !end l
    endif
    

    do ith = 1,2*nsmax      ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.d0  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
       endif
       lam_mm = sq4pi_inv ! lambda_00
       scalem=1

       !print*,'alm2map_multi_ball l325 B m loop start '

       b_north=0
       b_south=0

       b_n=0
       b_s=0

       if (mod(ith,50)==0) print*, 'ring number',ith


       do m = 0, nmmax


          f2m = 2.d0 * m

          !           ---------- l = m ----------
          par_lm = 1.d0  ! = (-1)^(l+m)
          if (m .ge. 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac=exp(scalem*logOVFLOW) !Actual lam_mm
          lam_lm = lam_mm*corfac/OVFLOW


          b_n(0:nvmaps-1,0:nmaps-1) = lam_lm *almjq(0:nvmaps-1,0:nmaps-1,m,m) 
          b_s(0:nvmaps-1,0:nmaps-1) = lam_lm *almjq(0:nvmaps-1,0:nmaps-1,m,m)

          !           ---------- l > m ----------
          lam_0 = 0.d0
          lam_1 = 1.d0 
          scalel=0
          a_rec = recfac(m,m)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 

             b_n(0:nvmaps-1,0:nmaps-1) = b_n(0:nvmaps-1,0:nmaps-1) + lam_lm *almjq(0:nvmaps-1,0:nmaps-1,l,m) 
             b_s(0:nvmaps-1,0:nmaps-1) = b_s(0:nvmaps-1,0:nmaps-1) + par_lm * lam_lm *almjq(0:nvmaps-1,0:nmaps-1,l,m)


             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = recfac(l,m)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) .gt. OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac=exp((scalel+scalem)*logOVFLOW)
             else
                if (abs(lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac=exp((scalel+scalem)*logOVFLOW)
                endif
             endif

          enddo !! end l

          do jj=0,nmaps-1
             !nvj=ceiling(bb**(j0+jj))
             do iv=0,nvj-1
                b_north(m,iv,jj) = b_n(iv,jj) 
                b_south(m,iv,jj) = b_s(iv,jj)
             enddo
          enddo


       enddo !end loop over m

       do i=0,nmaps-1
          do vv=0,nvj-1
             call ring_synthesis_ball(nsmax,nlmax,nmmax,b_north(:,vv,i),nph,ring(:),kphi0)   ! north hemisph. + equator
             map(istart_north:istart_north+nph-1,vv,i) = ring(0:nph-1)
             !if (feedback>3) print*, 'l464: j, v, sum(b_north) ',j0+i, vv, sum(b_north(:,vv,i))
             !if (feedback>3) print*, 'l465: j, v, sum(b_south) ',j0+i, vv, sum(b_south(:,vv,i))
          enddo
       enddo
       istart_north = istart_north + nph


       istart_south = istart_south - nph
       if (ith .lt. 2*nsmax) then
          do i=0,nmaps-1
             do vv=0,nvj-1
                call ring_synthesis_ball(nsmax,nlmax,nmmax,b_south(:,vv,i),nph,ring(:),kphi0) ! south hemisph. w/o equat
                map(istart_south:istart_south+nph-1,vv,i) = ring(0:nph-1)
             enddo !end vv
          enddo
          !if (feedback>2) print*, 'l482: ith_ring south, sum(map(istart_south:istart_south+nph-1,:,:)) ',ith,sum(map(istart_south:istart_south+nph-1,:,:))

          !print*,'alm2map_multi_ball l325 south done, ith= ',ith
       endif

    enddo    ! loop on cos(theta)

    !if (feedback>1)  print*,'alm2map_multi_ball l325 E, multiball done! sum(map) = ',sum(map)
    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(recfac)
    DEALLOCATE(almjq)
    DEALLOCATE(b_north)
    DEALLOCATE(b_south)
    DEALLOCATE(b_n)
    DEALLOCATE(b_s)
    DEALLOCATE(ring)
    return
  end subroutine alm2map_sincos_ball


!=======================================

  subroutine alm2map_fftw_ball(nsmax,nvmax,nlmax,nmmax,almn,bb,j0,map,gln,nmaps,der)

    use healpix_types
    IMPLICIT none
    INTEGER(I4B), INTENT(IN) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax,j0,nvmax
    INTEGER(I4B), INTENT(IN) :: nmmax,nmaps
    INTEGER(I4B), INTENT(IN), OPTIONAL :: der
    COMPLEX(DPC), INTENT(IN),  DIMENSION(0:nvmax,0:nlmax,0:nmmax) :: almn
    REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nvmax, 0:nmaps-1) :: gln
    REAL(DP),     INTENT(OUT), DIMENSION(0:12*nsmax**2-1,0:nvmax, 0:nmaps-1) :: map

    INTEGER(I4B) :: l, m, ith, indl, scalem, scalel, i         ! alm related
    INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
    integer(i4b):: iv, vv,in,jj,nv,nvj,nvmaxj,nlmaxj,nnj

    REAL(DP) :: cth, sth, dth1, dth2, dst1, OVFLOW, UNFLOW, logOVFLOW,bb
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) :: f2m, fm2, fl2, corfac,vval
    COMPLEX(DPC) :: factor

    CHARACTER(LEN=7), PARAMETER :: code = 'ALMN2MAP_BALL'
    COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almjq
    COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: b_north, b_south
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: b_n, b_s
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fn
    INTEGER(I4B) :: status

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: recfac
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    !=======================================================================


    !     --- allocates space for arrays ---
    ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)    
    if (status /= 0) stop 'recfac'

    ALLOCATE(b_north(0:nmmax,0:nvmax,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_north'

    ALLOCATE(b_n(0:nvmax,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_n'

    ALLOCATE(b_south(0:nmmax,0:nvmax,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_south'

    ALLOCATE(b_s(0:nvmax,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_s'

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) stop 'ring'

    ALLOCATE(almjq(0:nvmax,0:nmaps-1,0:nlmax,0:nmmax),stat = status) 
    if (status /= 0) stop 'almjq'

    if (feedback>2) print*, 'l320: total(gln): ',sum(gln)
    if (feedback>2) print*, 'l321: total(almn): ',sum(almn)
    if (feedback>1) print*,'l322: alm2map_multi_ball l321 A '

    !     ------------ initiate arrays ----------------

    do m = 0, nmmax
       fm2 = DBLE(m) **2
       do l = m, nlmax
          fl2 = DBLE(l+1) **2
          recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       enddo
    enddo
    !     --------------------------------------------

    !     ----- set the whole map to zero ------
    map = 0.0
    !     --------------------------------------

    istart_north = 0
    istart_south = 12*nsmax**2
    OVFLOW=exp(log(1.0d30))
    UNFLOW=exp(log(1.0d-30))
    logOVFLOW=log(1.0d30)

    dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
    dth2 = 2.d0 / (3.d0*DBLE(nsmax))
    dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

    almjq=0

    !test almn and gln
    ! do iv=0,nvmax
    !    print*, 'iv, sum(glm), sum(almn)',sum(gln(:,iv, :)), sum(almn(iv,:,:))
    ! enddo

    !=========================
    !---- this is how it should have been if gln was not reordered 
    !---- such that gln(l, n=0..N/2,-N/2-1, -1,j)
    !====================
    ! do iv=0,nvmax
    !    in=iv  !from n=0, n=N/2
    !    if (iv>(nvmax+1)/2) in=(nvmax+1-iv) !from n=-N/2,-1
    !    fn(iv)=almn(iv,l,m)*gln(l,in, jj) !nvmax=nvmax
    ! enddo
    !===============

    allocate(fn(0:nvmax))
    do m=0,nmmax 
       do l=0,nlmax
          do jj=0,nmaps-1
             fn(0:nvmax)=almn(0:nvmax,l,m)*gln(l,0:nvmax, jj) !nvmax=nvmax
             call complex_fft(fn,backward=.true.) 
             almjq(0:nvmax, jj,l,m)=fn(0:nvmax)
             !if (l==3 .and. m==2) print*, 'l=2,m=2: j, fn = ',jj,real(fn)
          enddo !end j
       enddo ! end l
    enddo !end m
    deallocate(fn)


  
    open(12,file='temp/almjq_from_ball.unf',status='unknown',form='unformatted')
    write(12) almjq
    close(12)

    ! do jj=0,nvmax
    !    print*,'ir, sum(almj1(jj,:,:,:)) = ',sum(almjq(jj,:,:,:))
    ! enddo

    do ith = 1,2*nsmax      ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.d0  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
       endif
       lam_mm = sq4pi_inv ! lambda_00
       scalem=1

       !print*,'alm2map_multi_ball l325 B m loop start '

       b_north=0
       b_south=0

       b_n=0
       b_s=0

       if (mod(ith,50)==0) print*, 'ring number',ith


       do m = 0, nmmax


          f2m = 2.d0 * m

          !           ---------- l = m ----------
          par_lm = 1.d0  ! = (-1)^(l+m)
          if (m .ge. 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac=exp(scalem*logOVFLOW) !Actual lam_mm
          lam_lm = lam_mm*corfac/OVFLOW




          b_n(0:nvmax,0:nmaps-1) = lam_lm *almjq(0:nvmax,0:nmaps-1,m,m) 
          b_s(0:nvmax,0:nmaps-1) = lam_lm *almjq(0:nvmax,0:nmaps-1,m,m)

          !           ---------- l > m ----------
          lam_0 = 0.d0
          lam_1 = 1.d0 
          scalel=0
          a_rec = recfac(m,m)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 

             b_n(0:nvmax,0:nmaps-1) = b_n(0:nvmax,0:nmaps-1) + lam_lm *almjq(0:nvmax,0:nmaps-1,l,m) 
             b_s(0:nvmax,0:nmaps-1) = b_s(0:nvmax,0:nmaps-1) + par_lm * lam_lm *almjq(0:nvmax,0:nmaps-1,l,m)


             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = recfac(l,m)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) .gt. OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac=exp((scalel+scalem)*logOVFLOW)
             else
                if (abs(lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac=exp((scalel+scalem)*logOVFLOW)
                endif
             endif

          enddo !! end l

          do jj=0,nmaps-1
             do iv=0,nvmax
                b_north(m,iv,jj) = b_n(iv,jj) 
                b_south(m,iv,jj) = b_s(iv,jj)
             enddo
          enddo


       enddo !end loop over m

       do i=0,nmaps-1
          do vv=0,nvmax
             call ring_synthesis_ball(nsmax,nlmax,nmmax,b_north(:,vv,i),nph,ring(:),kphi0)   ! north hemisph. + equator
             map(istart_north:istart_north+nph-1,vv,i) = ring(0:nph-1)
             !ap(istart_north:istart_north+nph-1,vv,i,2) = CONJG(ring(0:nph-1)
          enddo
          !if(ith==3) print*, 'i,vv,sum(map(:,:,j))',i,vv,sum(map(:,:,i))
       enddo
       istart_north = istart_north + nph


       istart_south = istart_south - nph
       if (ith .lt. 2*nsmax) then
          do i=0,nmaps-1
             do vv=0,nvmax
                call ring_synthesis_ball(nsmax,nlmax,nmmax,b_south(:,vv,i),nph,ring(:),kphi0) ! south hemisph. w/o equat
                map(istart_south:istart_south+nph-1,vv,i) = ring(0:nph-1)
                !map(istart_south:istart_south+nph-1,vv,i,2) = CONJG(ring(0:nph-1))
             enddo !end vv
          enddo
       endif

    enddo    ! loop on cos(theta)

    ! do vv=0,nvmax
    !    print*, 'j, sum(beta_ballmap)',sum(b_north(:,vv,:)),sum(b_south(:,vv,:))
    ! enddo

    !if (feedback>0) print*,'alm2map_multi_ball l325 E, multiball done! sum(map) = ',sum(map)
    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(recfac)
    DEALLOCATE(almjq)
    DEALLOCATE(b_north)
    DEALLOCATE(b_south)
    DEALLOCATE(b_n)
    DEALLOCATE(b_s)
    DEALLOCATE(ring)
    return
  end subroutine alm2map_fftw_ball


!!================================================================
   subroutine map2alm_ball_r(nsmax, nlmax, nmmax, nnmax, j0,nj,nvmaps,map,bb,gln, nshell, rvec,alm_out, cos_theta_cut, w8ring,imin,imax)
     !=======================================================================
     !     computes the a(l,m) from a map for the HEALPIX pixelisation
     !      for the Temperature field
     !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
     !            = int dtheta lambda_lm(theta) 
     !                  * int dphi T(theta,phi) e^(-i*m*phi)
     !
     !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
     !
     !     * the recurrence of Ylm is the standard one (cf Num Rec)
     !     * the integral over phi is done by FFT
     !
     !     cos_theta_cut (>0) is the cosine of the 
     !     symmetric cut applied to the sky
     !     if it is <= 0 no cut is applied
     !
     !     NB : these a(l,m) have to be multiplied by the pixel solid angle
     !      to give the correct coefficients
     !
     !             -------------------------------
     !         precomputes the Lambda_lm recurrence factor
     !      and is therefore ~50% faster than previous versions
     !     the multiplication by omega_pix is done in the routine
     !             -------------------------------
     !
     !=======================================================================
     IMPLICIT none
     INTEGER(I4B), INTENT(IN) :: nsmax
     INTEGER(I4B), INTENT(IN) :: nlmax
     INTEGER(I4B), INTENT(IN) :: nmmax
     INTEGER(I4B), INTENT(IN) :: nnmax,nvmaps,j0,nshell,nj
     REAL(DP), INTENT(IN) :: bb
     REAL(DP), INTENT(IN), dimension(0:nshell) :: rvec
     REAL(DP),     INTENT(IN), DIMENSION(0:12*nsmax**2-1,0:nvmaps,0:nj-1) :: map
     REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax,0:nj-1) :: gln
     COMPLEX(DPC), INTENT(OUT), DIMENSION(0:nshell,1:1,0:nlmax,0:nlmax) :: alm_out
     COMPLEX(DPC), DIMENSION(0:nnmax,0:nj-1,0:nlmax,0:nlmax) :: alm
     REAL(DP),     INTENT(IN) :: cos_theta_cut
     REAL(DP),     INTENT(IN), DIMENSION(1:2*nsmax,1) :: w8ring
     INTEGER(I4B), OPTIONAL, INTENT(IN) :: imin,imax

     INTEGER(I4B) :: l, m, ith, im_max, scalem, scalel       ! alm related
     INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
     INTEGER(I4B) :: i, j, irmax, irmin,iv,in,ij,nvj,nmaps,jj,nnj

     REAL(DP) :: omega_pix,vval
     REAL(DP) :: cth, sth, dth1, dth2, dst1
     REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
     REAL(DP) :: f2m, fm2, fl2, corfac,factor


     CHARACTER(LEN=7), PARAMETER :: code = 'MAP2ALM'
     COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: phas_n
     COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: phas_s
     COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_x
     INTEGER(I4B) :: status


     REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
     LOGICAL   :: keep_it,keep_n,keep_s

!!!      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ylmx    !xxx
     !=======================================================================


     nmaps=nj
     nvj=nvmaps


     irmax=4*nsmax
     irmin=1

     if (present(imin)) then
        irmin=imin
        irmax=imax
     endif


     ! --- allocates space for arrays ---
     if (allocated(recfac) .eqv. .false.)  then 
        ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)
        if (status /= 0) stop 'recfac'
     endif

     ALLOCATE(phas_n(0:nmmax,0:nnmax,0:nj-1),stat = status) 
     if (status /= 0) stop 'phas_n'

     ALLOCATE(phas_s(0:nmmax,0:nnmax,0:nj-1),stat = status) 
     if (status /= 0) stop 'phas_s'

     ALLOCATE(phas_x(0:nmmax),stat = status)
     if (status /= 0) stop 'phas_x'

     ALLOCATE(ring(0:4*nsmax-1),stat = status) 
     if (status /= 0) stop 'ring'

!!!      ALLOCATE(ylmx(0:nlmax,0:nmmax),stat=status)      !xxx
!!!      if (status /= 0) call die_alloc(code,'ylmx')

     !     ------------ initiate arrays ----------------

     if (recfac_count<2) then
        do m = 0, nmmax
           fm2 = DBLE(m) **2
           do l = m, nlmax
              fl2 = DBLE(l+1) **2
              recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
           enddo
        enddo
     endif
     recfac_count = recfac_count+1

     !     --------------------------------------------

     !     ----- set the whole alm array to zero ------
     alm = 0.0
!!!      ylmx=0.d0     !xxx
     !     -------------------------------------------

     omega_pix = pi / (3.d0 * nsmax * nsmax)
     istart_north = 0
     istart_south = 12*nsmax**2


     dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
     dth2 = 2.d0 / (3.d0*DBLE(nsmax))
     dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

     !     --------------------------------------------

     !-----------------------------------------------------------------------
     !           computes the integral in phi : phas_m(theta)
     !           for each parallele from north to south pole
     !-----------------------------------------------------------------------

     print*,'alm2map_multi_rmap: nnmax,lmax:',nnmax,nlmax

     do ith = 1, 2*nsmax

        if (feedback>1) print*, 'map2alm_r ith ring: ',ith

        phas_x(0:nmmax) = CMPLX(0.d0, 0.d0)   ! arbitrary
        phas_n(0:nmmax,0:nnmax,0:nj-1) = CMPLX(0.d0, 0.d0)   ! North    m >= 0
        phas_s(0:nmmax,0:nnmax,0:nj-1) = CMPLX(0.d0, 0.d0)   ! South    m >= 0

        if (ith .le. nsmax-1) then      ! north polar cap
           nph = 4*ith
           kphi0 = 1 
           cth = 1.d0  - DBLE(ith)**2 * dth1
           sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
        else                            ! tropical band + equat.
           nph = 4*nsmax
           kphi0 = MOD(ith+1-nsmax,2)
           cth = DBLE(2*nsmax-ith) * dth2
           sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
        endif

        keep_it = (ABS(cth).gt.cos_theta_cut) ! part of the sky out of the symmetric cut
        keep_n = (ith.ge.irmin).and.(ith.le.irmax)
        keep_s = (((4*nsmax-ith).ge.irmin).and.((4*nsmax-ith).le.irmax))
        if (.not.(keep_n.or.keep_s)) keep_it=keep_n

        !         if (.not.keep_n) print*,'Not N'
        !         if (.not.keep_s) print*,'Not S'

        if (keep_it.and.keep_n) then

           if (feedback>4)  print*,'north:',ith

           do jj=0,nmaps-1
              do in=0,nnmax

                 ring=0.
                 do iv=0,nvj-1
                    vval= dble(iv)
                    ring(0:nph-1) = ring(0:nph-1) + map(istart_north:istart_north+nph-1,iv,jj)*snrv(in, vval,j0+jj,bb,nvj) * w8ring(ith,1)
                 enddo

                 call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
                 phas_n(:,in,jj)=phas_x

              enddo !end in

           enddo !end jj

           if (feedback>4)  print*,'ring analysis north done:'
        endif

        istart_north = istart_north + nph
        istart_south = istart_south - nph

        if (ith .lt. 2*nsmax .and. (keep_it.and.keep_s)) then
           if (feedback>4)  print*,'south:',ith

           do jj=0,nmaps-1
              do in=0,nnmax

                 ring=0.
                 do iv=0,nvj-1
                    vval= dble(iv)
                    ring(0:nph-1) = ring(0:nph-1) + map(istart_south:istart_south+nph-1,iv,jj)*snrv(in, vval,j0+jj,bb,nvj) * w8ring(ith,1)
                 enddo
                 call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
                 phas_s(:,in,jj)=phas_x

              enddo !end in
           enddo !end jj
           if (feedback>3)  print*,'ring analysis south done:'
        endif

        !-----------------------------------------------------------------------
        !              computes the a_lm by integrating over theta
        !                  lambda_lm(theta) * phas_m(theta)
        !                         for each m and l
        !-----------------------------------------------------------------------

        if (keep_it) then

           if (feedback>3) print*,'recursion'

           lam_mm = sq4pi_inv ! lambda_00
           scalem=1
           im_max = nmmax !16/05/1999 B. D. Wandelt -> proper delta function sampling
           do m = 0, im_max
              f2m = 2.d0 * m
              

              !           ---------- l = m ----------
              par_lm = 1.d0  ! = (-1)^(l+m)
              if (m .ge. 1) then ! lambda_0_0 for m>0
                 lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
              endif

              if (abs(lam_mm).lt.UNFLOW) then
                 lam_mm=lam_mm*OVFLOW
                 scalem=scalem-1
              endif
              corfac = exp(scalem*logOVFLOW)   ! Actual lambda_mm 
              lam_lm = corfac*lam_mm/OVFLOW

              if (feedback>4) print*,'summing over n and v ..m=',m 

              do jj=0,nmaps-1
                 !nnj=min(nnmax,nvj)
                 do in=0,nnmax             
                       factor = (phas_n(m,in,jj) + phas_s(m,in,jj)) 
                       alm(in,jj, m, m) = alm(in, jj, m, m) +  lam_lm * factor 

                 enddo !! end in
              enddo ! jj end

              if (feedback>4) print*,'summing over n and v done! ..m=',m 
!!!            if (ith.eq.2*nsmax) then
!!!               ylmx(m,m)=lam_lm   !xxx
!!!            endif
              !           ---------- l > m ----------
              lam_0 = 0.d0
              lam_1 = 1.d0
              scalel=0
              a_rec = recfac(m,m)
              lam_2 = cth * lam_1 * a_rec
              do l = m+1, nlmax
                 par_lm = - par_lm  ! = (-1)^(l+m)
                 
                 lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors
                 
                 if (feedback>4) print*,'summing over n and v ..l=',l 
                 
                 do jj=0,nmaps-1
                    do in=0,nnmax
                       factor = (phas_n(m,in,jj) + par_lm*phas_s(m,in,jj)) 
                       alm(in,jj, l, m) = alm(in,jj, l, m) + factor * lam_lm
                    enddo !end in
                 enddo !!end jj

                 if (feedback>4) print*,'summing over n and v done! ..l=',l 

!!!               if (ith.eq.2*nsmax) then
!!!                  ylmx(l,m)=lam_lm   !xxx
!!!               endif


                 lam_0 = lam_1 / a_rec
                 lam_1 = lam_2
                 a_rec = recfac(l,m)
                 lam_2 = (cth * lam_1 - lam_0) * a_rec

                 if (abs(lam_2) .gt. OVFLOW) then
                    lam_0=lam_0/OVFLOW
                    lam_1=lam_1/OVFLOW
                    lam_2 = (cth * lam_1 - lam_0) * a_rec
                    scalel=scalel+1
                    corfac=exp((scalem+scalel)*logOVFLOW)
                 else
                    if (abs(lam_2) .lt. UNFLOW) then
                       lam_0=lam_0*OVFLOW
                       lam_1=lam_1*OVFLOW
                       lam_2 = (cth * lam_1 - lam_0) * a_rec 
                       scalel=scalel-1
                       corfac=exp((scalem+scalel)*logOVFLOW)
                    endif
                 endif
              enddo
           enddo

        endif

     enddo


     !     --------------------------------------------------------------
     !     normalisation of the alm (multiplication by pixel solid angle)
     !     --------------------------------------------------------------
     alm = alm * omega_pix



!!!open(12,file='ylmx.unf',status='unknown',form='unformatted')   !xxx
!!!write(12) ylmx(0:nlmax,0:nmmax)
!!!close(12)

     !     --------------------
     !     free memory and exit
     !     --------------------
     DEALLOCATE(phas_n)
     DEALLOCATE(phas_s)
     DEALLOCATE(ring)
     !DEALLOCATE(recfac)
     DEALLOCATE(phas_x)


        ! print*, 'saving final map: ',trim(adjustl(dir))//'arlm_ball_beta.unf'
        ! open(12,file=trim(adjustl(dir))//'arlm_ball_beta.unf',status='unknown',form='unformatted')
        ! write(12) arlm
        ! close(12)



     alm_out=0.0

     do jj=0,nmaps-1
        do in=0,nnmax
           do l = 0, nlmax
              !if ( ((in*in+l*(l+1)) <= bb**(2d0*(j0+jj+1))) .and. ((in*in+l*(l+1)) >= bb**(2d0*(j0+jj-1))) ) then
                 do m = 0, l
                    alm_out(0:nshell-1,1, l, m) = alm_out(0:nshell-1, 1, l, m) + alm(in, jj, l, m)*gln(l,in,jj)*in*argpi*sinc(in*argpi*rvec(0:nshell-1))
                 enddo !end m
              !endif
     enddo !end m
  enddo
enddo


     RETURN
   END subroutine map2alm_ball_r


!=====================================
   subroutine map2alm_ball_rall(nsmax, nlmax, nmmax, nnmax, j0,nj,nvmaps,map,bb,gln, nshell, rvec,alm, cos_theta_cut, w8ring,imin,imax)
     !=======================================================================
     !     computes the a(l,m) from a map for the HEALPIX pixelisation
     !      for the Temperature field
     !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
     !            = int dtheta lambda_lm(theta) 
     !                  * int dphi T(theta,phi) e^(-i*m*phi)
     !
     !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
     !
     !     * the recurrence of Ylm is the standard one (cf Num Rec)
     !     * the integral over phi is done by FFT
     !
     !     cos_theta_cut (>0) is the cosine of the 
     !     symmetric cut applied to the sky
     !     if it is <= 0 no cut is applied
     !
     !     NB : these a(l,m) have to be multiplied by the pixel solid angle
     !      to give the correct coefficients
     !
     !             -------------------------------
     !         precomputes the Lambda_lm recurrence factor
     !      and is therefore ~50% faster than previous versions
     !     the multiplication by omega_pix is done in the routine
     !             -------------------------------
     !
     !=======================================================================
     IMPLICIT none
     INTEGER(I4B), INTENT(IN) :: nsmax
     INTEGER(I4B), INTENT(IN) :: nlmax
     INTEGER(I4B), INTENT(IN) :: nmmax
     INTEGER(I4B), INTENT(IN) :: nnmax,nvmaps,j0,nshell,nj
     REAL(DP), INTENT(IN) :: bb
     REAL(DP), INTENT(IN), dimension(0:nshell) :: rvec
     REAL(SP),     INTENT(IN), DIMENSION(0:12*nsmax**2-1,0:nvmaps,0:nj-1) :: map
     REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax,0:nj-1) :: gln
     COMPLEX(SPC), INTENT(OUT), DIMENSION(0:nshell,1:1,0:nlmax,0:nlmax) :: alm
     REAL(DP),     INTENT(IN) :: cos_theta_cut
     REAL(DP),     INTENT(IN), DIMENSION(1:2*nsmax,1) :: w8ring
     INTEGER(I4B), OPTIONAL, INTENT(IN) :: imin,imax

     INTEGER(I4B) :: l, m, ith, im_max, scalem, scalel       ! alm related
     INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
     INTEGER(I4B) :: i, j, irmax, irmin,iv,in,ij,nvj,nmaps,jj,nnj

     REAL(DP) :: omega_pix,vval
     REAL(DP) :: cth, sth, dth1, dth2, dst1
     REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
     REAL(DP) :: f2m, fm2, fl2, corfac,factor


     CHARACTER(LEN=7), PARAMETER :: code = 'MAP2ALM'
     COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: phas_n
     COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: phas_s
     COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_x
     INTEGER(I4B) :: status


     REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
     LOGICAL   :: keep_it,keep_n,keep_s

!!!      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ylmx    !xxx
     !=======================================================================

     nmaps=nj

     irmax=4*nsmax
     irmin=1

     if (present(imin)) then
        irmin=imin
        irmax=imax
     endif


     ! --- allocates space for arrays ---
     if (allocated(recfac) .eqv. .false.)  then 
        ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)
        if (status /= 0) stop 'recfac'
     endif

     ALLOCATE(phas_n(0:nmmax,0:nnmax,0:nj-1),stat = status) 
     if (status /= 0) stop 'phas_n'

     ALLOCATE(phas_s(0:nmmax,0:nnmax,0:nj-1),stat = status) 
     if (status /= 0) stop 'phas_s'

     ALLOCATE(phas_x(0:nmmax),stat = status)
     if (status /= 0) stop 'phas_x'

     ALLOCATE(ring(0:4*nsmax-1),stat = status) 
     if (status /= 0) stop 'ring'

!!!      ALLOCATE(ylmx(0:nlmax,0:nmmax),stat=status)      !xxx
!!!      if (status /= 0) call die_alloc(code,'ylmx')

     !     ------------ initiate arrays ----------------

     if (recfac_count<2) then
        do m = 0, nmmax
           fm2 = DBLE(m) **2
           do l = m, nlmax
              fl2 = DBLE(l+1) **2
              recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
           enddo
        enddo
     endif
     recfac_count = recfac_count+1

     !     --------------------------------------------

     !     ----- set the whole alm array to zero ------
     alm = 0.0
!!!      ylmx=0.d0     !xxx
     !     -------------------------------------------

     omega_pix = pi / (3.d0 * nsmax * nsmax)
     istart_north = 0
     istart_south = 12*nsmax**2


     dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
     dth2 = 2.d0 / (3.d0*DBLE(nsmax))
     dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

     !     --------------------------------------------

     !-----------------------------------------------------------------------
     !           computes the integral in phi : phas_m(theta)
     !           for each parallele from north to south pole
     !-----------------------------------------------------------------------

     do ith = 1, 2*nsmax

        if (feedback>1) print*, 'map2alm_r ith ring: ',ith

        phas_x(0:nmmax) = CMPLX(0.d0, 0.d0)   ! arbitrary
        phas_n(0:nmmax,0:nnmax,0:nj-1) = CMPLX(0.d0, 0.d0)   ! North    m >= 0
        phas_s(0:nmmax,0:nnmax,0:nj-1) = CMPLX(0.d0, 0.d0)   ! South    m >= 0

        if (ith .le. nsmax-1) then      ! north polar cap
           nph = 4*ith
           kphi0 = 1 
           cth = 1.d0  - DBLE(ith)**2 * dth1
           sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
        else                            ! tropical band + equat.
           nph = 4*nsmax
           kphi0 = MOD(ith+1-nsmax,2)
           cth = DBLE(2*nsmax-ith) * dth2
           sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
        endif

        keep_it = (ABS(cth).gt.cos_theta_cut) ! part of the sky out of the symmetric cut
        keep_n = (ith.ge.irmin).and.(ith.le.irmax)
        keep_s = (((4*nsmax-ith).ge.irmin).and.((4*nsmax-ith).le.irmax))
        if (.not.(keep_n.or.keep_s)) keep_it=keep_n

        !         if (.not.keep_n) print*,'Not N'
        !         if (.not.keep_s) print*,'Not S'

        if (keep_it.and.keep_n) then

           if (feedback>4)  print*,'north:',ith

           do jj=0,nmaps-1
              nvj=ceiling(bb**(j0+jj))
              nnj=  nnmax !min(nnmax,ceiling(bb**(j0+jj+1)))

              do in=1,nnj                

                 ring=0.
                 do iv=1,nvj-1
                    vval= dble(iv)
                    ring(0:nph-1) = ring(0:nph-1) + map(istart_north:istart_north+nph-1,iv,jj)*snrv(in, vval,j0+jj,bb) * w8ring(ith,1)
                 enddo

                 call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
                 phas_n(:,in,jj)=phas_x

              enddo !end in

           enddo !end jj

           if (feedback>4)  print*,'ring analysis north done:'
        endif

        istart_north = istart_north + nph
        istart_south = istart_south - nph

        if (ith .lt. 2*nsmax .and. (keep_it.and.keep_s)) then
           if (feedback>4)  print*,'south:',ith

           do jj=0,nmaps-1
              nvj=ceiling(bb**(j0+jj))
              nnj= nnmax !min(nnmax,ceiling(bb**(j0+jj+1)))

              do in=1,nnj

                 ring=0.
                 do iv=1,nvj-1
                    vval= dble(iv)
                    ring(0:nph-1) = ring(0:nph-1) + map(istart_south:istart_south+nph-1,iv,jj)*snrv(in, vval,j0+jj,bb) * w8ring(ith,1)
                 enddo
                 call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
                 phas_s(:,in,jj)=phas_x

              enddo !end in
           enddo !end jj
           if (feedback>3)  print*,'ring analysis south done:'
        endif

        !-----------------------------------------------------------------------
        !              computes the a_lm by integrating over theta
        !                  lambda_lm(theta) * phas_m(theta)
        !                         for each m and l
        !-----------------------------------------------------------------------

        if (keep_it) then

           if (feedback>3) print*,'recursion'

           lam_mm = sq4pi_inv ! lambda_00
           scalem=1
           im_max = nmmax !16/05/1999 B. D. Wandelt -> proper delta function sampling
           do m = 0, im_max
              f2m = 2.d0 * m
              

              !           ---------- l = m ----------
              par_lm = 1.d0  ! = (-1)^(l+m)
              if (m .ge. 1) then ! lambda_0_0 for m>0
                 lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
              endif

              if (abs(lam_mm).lt.UNFLOW) then
                 lam_mm=lam_mm*OVFLOW
                 scalem=scalem-1
              endif
              corfac = exp(scalem*logOVFLOW)   ! Actual lambda_mm 
              lam_lm = corfac*lam_mm/OVFLOW

              if (feedback>4) print*,'summing over n and v ..m=',m 

              do jj=0,nmaps-1
                 nvj=ceiling(bb**(j0+jj))
                 nnj= nnmax !min(nnmax,ceiling(bb**(j0+jj+1)))

                 do in=1,nnj
                    if ((in*in+m*(m+1)) <= bb**(2d0*(j0+jj+1))) then
                       factor = (phas_n(m,in,jj) + phas_s(m,in,jj))*gln(m,in,jj)
                       alm(1:nshell-1,1, m, m) = alm(1:nshell-1, 1, m, m) +  lam_lm * factor *in*argpi*sinc(in*argpi*rvec(1:nshell-1))
                    endif
                 enddo !! end in

              enddo ! jj end

              if (feedback>4) print*,'summing over n and v done! ..m=',m 
!!!            if (ith.eq.2*nsmax) then
!!!               ylmx(m,m)=lam_lm   !xxx
!!!            endif
              !           ---------- l > m ----------
              lam_0 = 0.d0
              lam_1 = 1.d0
              scalel=0
              a_rec = recfac(m,m)
              lam_2 = cth * lam_1 * a_rec
              do l = m+1, nlmax
                 par_lm = - par_lm  ! = (-1)^(l+m)
                 
                 lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors
                 
                 if (feedback>4) print*,'summing over n and v ..l=',l 
                 
                 do jj=0,nmaps-1
                    nvj=ceiling(bb**(j0+jj))
                    nnj=min(nnmax,ceiling(bb**(j0+jj+1)))          

                    do in=1,nnj
                       
                       if ((in*in+l*(l+1)) <= bb**(2d0*(j0+jj)+1)) then
                          factor = (phas_n(m,in,jj) + par_lm*phas_s(m,in,jj))*gln(l,in,jj)
                          alm(1:nshell-1,1, l, m) = alm(1:nshell-1,1, l, m) + factor * lam_lm *in*argpi*sinc(in*argpi*rvec(1:nshell-1))
                       endif
                       
                    enddo !end in
                 enddo !!end jj

                 if (feedback>4) print*,'summing over n and v done! ..l=',l 

!!!               if (ith.eq.2*nsmax) then
!!!                  ylmx(l,m)=lam_lm   !xxx
!!!               endif


                 lam_0 = lam_1 / a_rec
                 lam_1 = lam_2
                 a_rec = recfac(l,m)
                 lam_2 = (cth * lam_1 - lam_0) * a_rec

                 if (abs(lam_2) .gt. OVFLOW) then
                    lam_0=lam_0/OVFLOW
                    lam_1=lam_1/OVFLOW
                    lam_2 = (cth * lam_1 - lam_0) * a_rec
                    scalel=scalel+1
                    corfac=exp((scalem+scalel)*logOVFLOW)
                 else
                    if (abs(lam_2) .lt. UNFLOW) then
                       lam_0=lam_0*OVFLOW
                       lam_1=lam_1*OVFLOW
                       lam_2 = (cth * lam_1 - lam_0) * a_rec 
                       scalel=scalel-1
                       corfac=exp((scalem+scalel)*logOVFLOW)
                    endif
                 endif
              enddo
           enddo

        endif

     enddo

     !     --------------------------------------------------------------
     !     normalisation of the alm (multiplication by pixel solid angle)
     !     --------------------------------------------------------------
     alm = alm * omega_pix

!!!open(12,file='ylmx.unf',status='unknown',form='unformatted')   !xxx
!!!write(12) ylmx(0:nlmax,0:nmmax)
!!!close(12)

     !     --------------------
     !     free memory and exit
     !     --------------------
     DEALLOCATE(phas_n)
     DEALLOCATE(phas_s)
     DEALLOCATE(ring)
     !DEALLOCATE(recfac)
     DEALLOCATE(phas_x)

     RETURN
   END subroutine map2alm_ball_rall


!!===============

!=====================================
   subroutine map2alm_ball_jr(nsmax, nlmax, nmmax, nnmax, nvmaps,map,bb,jj,gln, alm, cos_theta_cut, w8ring,imin,imax)
     !=======================================================================
     !     computes the a(l,m) from a map for the HEALPIX pixelisation
     !      for the Temperature field
     !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
     !            = int dtheta lambda_lm(theta) 
     !                  * int dphi T(theta,phi) e^(-i*m*phi)
     !
     !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
     !
     !     * the recurrence of Ylm is the standard one (cf Num Rec)
     !     * the integral over phi is done by FFT
     !
     !     cos_theta_cut (>0) is the cosine of the 
     !     symmetric cut applied to the sky
     !     if it is <= 0 no cut is applied
     !
     !     NB : these a(l,m) have to be multiplied by the pixel solid angle
     !      to give the correct coefficients
     !
     !             -------------------------------
     !         precomputes the Lambda_lm recurrence factor
     !      and is therefore ~50% faster than previous versions
     !     the multiplication by omega_pix is done in the routine
     !             -------------------------------
     !
     !=======================================================================
     IMPLICIT none
     INTEGER(I4B), INTENT(IN) :: nsmax
     INTEGER(I4B), INTENT(IN) :: nlmax
     INTEGER(I4B), INTENT(IN) :: nmmax
     INTEGER(I4B), INTENT(IN) :: nnmax,nvmaps,jj
     REAL(DP), INTENT(IN) :: bb
     REAL(SP),     INTENT(IN), DIMENSION(0:12*nsmax**2-1,0:nvmaps) :: map
     REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax) :: gln
     COMPLEX(SPC), INTENT(OUT), DIMENSION(1:1,0:nlmax,0:nlmax) :: alm
     REAL(DP),     INTENT(IN) :: cos_theta_cut
     REAL(DP),     INTENT(IN), DIMENSION(1:2*nsmax,1) :: w8ring
     INTEGER(I4B), OPTIONAL, INTENT(IN) :: imin,imax

     INTEGER(I4B) :: l, m, ith, im_max, scalem, scalel       ! alm related
     INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
     INTEGER(I4B) :: i, j, irmax, irmin,iv,in,ij,nvj

     REAL(DP) :: omega_pix,vval
     REAL(DP) :: cth, sth, dth1, dth2, dst1
     REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
     REAL(DP) :: f2m, fm2, fl2, corfac,factor


     CHARACTER(LEN=7), PARAMETER :: code = 'MAP2ALM'
     COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_n
     COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_s
     COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_x
     INTEGER(I4B) :: status


     REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
     LOGICAL   :: keep_it,keep_n,keep_s

!!!      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ylmx    !xxx
     !=======================================================================


     irmax=4*nsmax
     irmin=1

     if (present(imin)) then
        irmin=imin
        irmax=imax
     endif

     nvj=ceiling(bb**(jj))

     ! --- allocates space for arrays ---
     if (allocated(recfac) .eqv. .false.)  then 
        ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)
        if (status /= 0) stop 'recfac'
     endif

     ALLOCATE(phas_n(0:nmmax,0:nvmaps),stat = status) 
     if (status /= 0) stop 'phas_n'

     ALLOCATE(phas_s(0:nmmax,0:nvmaps),stat = status) 
     if (status /= 0) stop 'phas_s'

     ALLOCATE(phas_x(0:nmmax),stat = status)
     if (status /= 0) stop 'phas_x'

     ALLOCATE(ring(0:4*nsmax-1),stat = status) 
     if (status /= 0) stop 'ring'

!!!      ALLOCATE(ylmx(0:nlmax,0:nmmax),stat=status)      !xxx
!!!      if (status /= 0) call die_alloc(code,'ylmx')

     !     ------------ initiate arrays ----------------

     if (recfac_count<2) then
        do m = 0, nmmax
           fm2 = DBLE(m) **2
           do l = m, nlmax
              fl2 = DBLE(l+1) **2
              recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
           enddo
        enddo
     endif
     recfac_count = recfac_count+1

     !     --------------------------------------------

     !     ----- set the whole alm array to zero ------
     alm = 0.0
!!!      ylmx=0.d0     !xxx
     !     -------------------------------------------

     omega_pix = pi / (3.d0 * nsmax * nsmax)
     istart_north = 0
     istart_south = 12*nsmax**2


     dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
     dth2 = 2.d0 / (3.d0*DBLE(nsmax))
     dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

     !     --------------------------------------------

     !-----------------------------------------------------------------------
     !           computes the integral in phi : phas_m(theta)
     !           for each parallele from north to south pole
     !-----------------------------------------------------------------------

     do ith = 1, 2*nsmax


        phas_x(0:nmmax) = CMPLX(0.d0, 0.d0)   ! arbitrary
        phas_n(0:nmmax,0:nvmaps) = CMPLX(0.d0, 0.d0)   ! North    m >= 0
        phas_s(0:nmmax,0:nvmaps) = CMPLX(0.d0, 0.d0)   ! South    m >= 0

        if (ith .le. nsmax-1) then      ! north polar cap
           nph = 4*ith
           kphi0 = 1 
           cth = 1.d0  - DBLE(ith)**2 * dth1
           sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
        else                            ! tropical band + equat.
           nph = 4*nsmax
           kphi0 = MOD(ith+1-nsmax,2)
           cth = DBLE(2*nsmax-ith) * dth2
           sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
        endif

        keep_it = (ABS(cth).gt.cos_theta_cut) ! part of the sky out of the symmetric cut
        keep_n = (ith.ge.irmin).and.(ith.le.irmax)
        keep_s = (((4*nsmax-ith).ge.irmin).and.((4*nsmax-ith).le.irmax))
        if (.not.(keep_n.or.keep_s)) keep_it=keep_n

        !         if (.not.keep_n) print*,'Not N'
        !         if (.not.keep_s) print*,'Not S'

        if (keep_it.and.keep_n) then

           if (feedback>4)  print*,'north:',ith
           do iv=1,nvmaps-1
              ring(0:nph-1) = map(istart_north:istart_north+nph-1,iv) * w8ring(ith,1)
              !print*,'ith, iv, sum(ring)',ith, iv,sum(ring)
              call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
              phas_n(:,iv)=phas_x
           enddo
           if (feedback>4)  print*,'ring analysis north done:'
        endif
        istart_north = istart_north + nph

        istart_south = istart_south - nph
        if (ith .lt. 2*nsmax .and. (keep_it.and.keep_s)) then

           if (feedback>4)  print*,'south:',ith
           do iv=1,nvmaps-1
              ring(0:nph-1) = map(istart_south:istart_south+nph-1,iv) * w8ring(ith,1)
              !print*,'ith, iv, sum(ring)',ith, iv,sum(ring)
              call ring_analysis_org_ball(nsmax, nlmax, nmmax, ring, nph, phas_x, kphi0)
              phas_s(:,iv)=phas_x
           enddo
           if (feedback>3)  print*,'ring analysis south done:'
        endif

        !-----------------------------------------------------------------------
        !              computes the a_lm by integrating over theta
        !                  lambda_lm(theta) * phas_m(theta)
        !                         for each m and l
        !-----------------------------------------------------------------------

        if (keep_it) then

           if (feedback>3) print*,'recursion'

           lam_mm = sq4pi_inv ! lambda_00
           scalem=1
           im_max = nmmax !16/05/1999 B. D. Wandelt -> proper delta function sampling
           do m = 0, im_max
              f2m = 2.d0 * m

              !           ---------- l = m ----------
              par_lm = 1.d0  ! = (-1)^(l+m)
              if (m .ge. 1) then ! lambda_0_0 for m>0
                 lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
              endif

              if (abs(lam_mm).lt.UNFLOW) then
                 lam_mm=lam_mm*OVFLOW
                 scalem=scalem-1
              endif
              corfac = exp(scalem*logOVFLOW)   ! Actual lambda_mm 
              lam_lm = corfac*lam_mm/OVFLOW

              if (feedback>4) print*,'summing over n and v ..m=',m 
              do in=1,nnmax                   
                 if ((in*in+l*(l+1)) < bb**(2d0*(jj))) then
                    do iv=1,nvj-1
                       vval= dble(iv)
                       factor = snrv(in, vval,jj,bb)*gln(m,in)
                       alm(1, m, m) = alm(1, m, m) + factor * lam_lm * (phas_n(m,iv) + phas_s(m,iv))
                    enddo
                 endif
              enddo
              if (feedback>4) print*,'summing over n and v done! ..m=',m 
!!!            if (ith.eq.2*nsmax) then
!!!               ylmx(m,m)=lam_lm   !xxx
!!!            endif
              !           ---------- l > m ----------
              lam_0 = 0.d0
              lam_1 = 1.d0
              scalel=0
              a_rec = recfac(m,m)
              lam_2 = cth * lam_1 * a_rec
              do l = m+1, nlmax
                 par_lm = - par_lm  ! = (-1)^(l+m)

                 lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors

              if (feedback>4) print*,'summing over n and v ..l=',l 
                 do in=1,nnmax                   
                    if ((in*in+l*(l+1)) < bb**(2d0*(jj))) then
                       do iv=1,nvj-1
                          vval= dble(iv)
                          factor = snrv(in, vval,jj,bb)*gln(m,in)
                          alm(1, l, m) = alm(1, l, m) + factor * lam_lm * (phas_n(m,iv) + par_lm*phas_s(m,iv))
                       enddo
                    endif
                 enddo
              if (feedback>4) print*,'summing over n and v done! ..l=',l 

!!!               if (ith.eq.2*nsmax) then
!!!                  ylmx(l,m)=lam_lm   !xxx
!!!               endif


                 lam_0 = lam_1 / a_rec
                 lam_1 = lam_2
                 a_rec = recfac(l,m)
                 lam_2 = (cth * lam_1 - lam_0) * a_rec

                 if (abs(lam_2) .gt. OVFLOW) then
                    lam_0=lam_0/OVFLOW
                    lam_1=lam_1/OVFLOW
                    lam_2 = (cth * lam_1 - lam_0) * a_rec
                    scalel=scalel+1
                    corfac=exp((scalem+scalel)*logOVFLOW)
                 else
                    if (abs(lam_2) .lt. UNFLOW) then
                       lam_0=lam_0*OVFLOW
                       lam_1=lam_1*OVFLOW
                       lam_2 = (cth * lam_1 - lam_0) * a_rec 
                       scalel=scalel-1
                       corfac=exp((scalem+scalel)*logOVFLOW)
                    endif
                 endif
              enddo
           enddo

        endif

     enddo

     !     --------------------------------------------------------------
     !     normalisation of the alm (multiplication by pixel solid angle)
     !     --------------------------------------------------------------
     alm = alm * omega_pix

!!!open(12,file='ylmx.unf',status='unknown',form='unformatted')   !xxx
!!!write(12) ylmx(0:nlmax,0:nmmax)
!!!close(12)

     !     --------------------
     !     free memory and exit
     !     --------------------
     DEALLOCATE(phas_n)
     DEALLOCATE(phas_s)
     DEALLOCATE(ring)
     !DEALLOCATE(recfac)
     DEALLOCATE(phas_x)

     RETURN
   END subroutine map2alm_ball_jr

!!============================
  subroutine ring_synthesis_ball(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)

    use healpix_types
    IMPLICIT none

    INTEGER(I4B), INTENT(IN) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    INTEGER(I4B), INTENT(IN) :: nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0

    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(IN) :: datain
    REAL(DP),     DIMENSION(0:nph-1), INTENT(OUT)     :: dataout
    REAL(DP),     DIMENSION(0:nph-1)     :: data
    COMPLEX(DPC),     DIMENSION(:,:), ALLOCATABLE, SAVE :: trig
    REAL(DP) :: phi0
    INTEGER(I4B) :: i,iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    REAL(DP), DIMENSION(8*8192+15), SAVE :: wsave
    INTEGER(I4B), SAVE :: npho=0, trg=-1,oldnmmax=-1
    COMPLEX(DPC) :: dat
    INTEGER(I4B) :: status

    !=======================================================================

    if (npho.ne.nph) then

       call rffti(nph,wsave)
       npho=nph

    endif

    if ((trg/=nsmax).or.(oldnmmax/=nmmax)) then
       if(allocated(trig)) deallocate(trig)

       ALLOCATE(trig(0:max(2*nsmax,nmmax),0:nsmax),stat = status) 
       if (status /= 0) stop 'trig'


       trig(:,:)=CMPLX(1.0_dp,0.0_dp)
       do i=1,nsmax
          phi0=PI/DBLE(i*4)
          do m=0,max(2*nsmax,nmmax)
             trig(m,i)= CMPLX( DCOS(m*phi0), DSIN(m*phi0), kind=DP)
          enddo
       enddo
       trg=nsmax
       oldnmmax=nmmax

    endif

    !-----------------------------------------------------------------------

    ksign = + 1
    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0, 0.0, DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, nmmax                        ! in -nmmax, nmmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo
    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    data(0)=REAL(bw(0))
    do iw = 1, nph/2-1
       m = ksign*(iw)
       if(kphi0==1) then
          dat =bw(iw) * trig(m,nph/4)
       else
          dat =bw(iw)
       endif
       data(iw*2-1) = REAL(dat)
       data(iw*2  ) = AIMAG(dat)
    enddo
    iw=nph/2
    m = ksign*(iw)
    if(kphi0==1) then
       dat =bw(iw) * trig(m,nph/4)
    else
       dat =bw(iw)
    endif
    data(iw*2-1) = REAL(dat)
    call rfftb(nph,data,wsave)

    !     ^^^^^^^^^^^^
    dataout=REAL(data)

    RETURN
  END subroutine ring_synthesis_ball

!=========================

  !**************************************************************************
  !
  !             FOURIER TRANSFORM ON RINGS
  !
  !**************************************************************************
  !=======================================================================
  subroutine ring_synthesis_complex(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
    !=======================================================================
    !     RING_SYNTHESIS
    !       called by alm2map
    !       calls     complex_fft
    !
    !     dataout(j) = sum_m datain(m) * exp(i*m*phi(j))
    !     with phi(j) = j*2pi/nph + kphi0*pi/nph and kphi0 =0 or 1
    !
    !     as the set of frequencies {m} is larger than nph,
    !     we wrap frequencies within {0..nph-1}
    !     ie  m = k*nph + m' with m' in {0..nph-1}
    !     then
    !     noting bw(m') = exp(i*m'*phi0)
    !                   * sum_k (datain(k*nph+m') exp(i*k*pi*kphi0))
    !        with bw(nph-m') = CONJ(bw(m')) (if datain(-m) = CONJ(datain(m)))
    !     dataout(j) = sum_m' [ bw(m') exp (i*j*m'*2pi/nph) ]
    !                = Fourier Transform of bw
    !        is real
    !
    !         NB nph is not necessarily a power of 2
    !
    !=======================================================================


    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(IN)  :: datain
    COMPLEX(DPC),     DIMENSION(0:nph-1), INTENT(OUT) :: dataout

    INTEGER(I4B) :: iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    type(planck_fft2_plan) :: plan
    COMPLEX(DPC) :: dat
    real(DP)     :: arg
    !=======================================================================

    ksign = + 1
    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, nmmax                        ! in -nmmax, nmmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo

    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    dataout(0)=bw(0)
    do iw = 1, nph/2-1
       m = ksign*(iw)
       if(kphi0==1) then
          arg = m * PI / dble(nph)
          dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
       else
          dat =bw(iw)
       endif
       dataout(iw*2-1) = CMPLX(REAL(dat, kind=DP),0.0_dp,kind=DP)
       dataout(iw*2  ) = CMPLX(0.0_dp,AIMAG(dat),kind=DP) 
    enddo
    iw=nph/2
    m = ksign*(iw)
    if(kphi0==1) then
       arg = m * PI / dble(nph)
       dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
    else
       dat =bw(iw)
    endif
    dataout(iw*2-1) = dat

    call make_fft2_plan(plan,length=nph,direction=fft2_backward)
    call complex_fft2 (plan, dataout)

    !call real_fft2 (plan, dataout)
    !     ^^^^^^^^^^^^
    call destroy_fft2_plan (plan)
    RETURN
  END subroutine ring_synthesis_complex



!!===========================
!=======================================================================
   subroutine ring_analysis_ball(nsmax,nmmax,datain,nph,dataout,kphi0)
     !=======================================================================
     !     ring_analysis
     !       called by map2alm
     !       calls     fft_gpd
     !
     !     integrates (data * phi-dependence-of-Ylm) over phi
     !     --> function of m can be computed by FFT
     !     with  0<= m <= npoints/2 (: Nyquist)
     !     because the data is real the negative m are the conjugate of the 
     !     positive ones
     !=======================================================================
     IMPLICIT none

     INTEGER(I4B), INTENT(IN) :: nsmax
     INTEGER(I4B), INTENT(IN) :: nmmax
     INTEGER(I4B), INTENT(IN) :: nph, kphi0

     REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
     COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(OUT) :: dataout
     COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE , SAVE :: trig
     INTEGER(I4B) :: i,m,im_max,ksign
     REAL(DP) :: phi0
     REAL(DP), DIMENSION(0:nph-1) :: data
     REAL(DP), DIMENSION(8*8192+15), SAVE :: wsave
     INTEGER(I4B), SAVE :: npho=0, trg=-1,oldnmmax=-1
     INTEGER(I4B) :: status

     !-----------------------------------------------------------------------

     if (npho.ne.nph) then

        call rffti(nph,wsave)
        npho=nph

     endif


     if ((trg/=nsmax).or.(oldnmmax/=nmmax)) then
        if(allocated(trig)) deallocate(trig)

        ALLOCATE(trig(0:nmmax,0:nsmax),stat = status) 
!        if (status /= 0) call die_alloc('ring_analysis','trig')


        trig(:,:)=CMPLX(1.0_dp,0.0_dp)
        do i=1,nsmax
           phi0=PI/DBLE(i*4)
           do m=0,nmmax
              trig(m,i)= CMPLX( DCOS(m*phi0), DSIN(m*phi0), kind=DP)
           enddo
        enddo
        trg=nsmax
        oldnmmax=nmmax

     endif

     !-----------------------------------------------------------------------


     ksign = - 1
     data=0.
     data(0:nph-1)=datain(0:nph-1)

     call rfftf(nph,data,wsave)

     im_max = MIN(nph/2,nmmax)
     dataout(0)=CMPLX(data(0),0_dp,kind=DP)

     do i = 1, im_max*2-3, 2
        dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP) 
     enddo

     if(im_max==nph/2) then
        dataout(im_max)= CMPLX( data(nph-1),0,kind=DP)
     else
        dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
     endif

     if(im_max==nmmax) goto 1000

     do i =  im_max+1,min(nph,nmmax)
        dataout(i) = conjg(dataout(2*im_max-i) )
     end do

     if(min(nph,nmmax)==nmmax) goto 1000

     do i =  2*im_max+1,nmmax
        dataout(i) = dataout(mod(i,2*im_max)) 
     end do

1000 continue

     if(kphi0==1)then
        do i =  0,nmmax
           m = ksign*i
           dataout(i)=dataout(i)* CONJG(trig(-m,nph/4))
        enddo
     end if

     RETURN
   END subroutine ring_analysis_ball

!==================================
   subroutine ring_analysis_org_ball(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
     !=======================================================================
     !     ring_analysis
     !       called by map2alm
     !       calls     fft_gpd
     !
     !     integrates (data * phi-dependence-of-Ylm) over phi
     !     --> function of m can be computed by FFT
     !     with  0<= m <= npoints/2 (: Nyquist)
     !     because the data is real the negative m are the conjugate of the 
     !     positive ones
     !=======================================================================
     IMPLICIT none

     INTEGER(I4B), INTENT(IN) :: nsmax
     INTEGER(I4B), INTENT(IN) :: nlmax
     INTEGER(I4B), INTENT(IN) :: nmmax
     INTEGER(I4B), INTENT(IN) :: nph, kphi0

     REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
     COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(OUT) :: dataout
     COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE , SAVE :: trig
     INTEGER(I4B) :: i,m,im_max,ksign
     REAL(DP) :: phi0
     REAL(DP), DIMENSION(0:nph-1) :: data
     REAL(DP), DIMENSION(8*8192+15), SAVE :: wsave
     INTEGER(I4B), SAVE :: npho=0, trg=-1,oldnmmax=-1
     INTEGER(I4B) :: status

     !-----------------------------------------------------------------------

     if (npho.ne.nph) then

        call rffti(nph,wsave)
        npho=nph

     endif


     if ((trg/=nsmax).or.(oldnmmax/=nmmax)) then
        if(allocated(trig)) deallocate(trig)

        ALLOCATE(trig(0:max(2*nsmax,nmmax),0:nsmax),stat = status) 
        if (status /= 0) stop 'ring_analysis'


        trig(:,:)=CMPLX(1.0_dp,0.0_dp)
        do i=1,nsmax
           phi0=PI/DBLE(i*4)
           do m=0,max(2*nsmax,nmmax)
              trig(m,i)= CMPLX( DCOS(m*phi0), DSIN(m*phi0), kind=DP)
           enddo
        enddo
        trg=nsmax
        oldnmmax=nmmax

     endif

     !-----------------------------------------------------------------------


     ksign = - 1
     data=0.
     data(0:nph-1)=datain(0:nph-1)

     call rfftf(nph,data,wsave)

     im_max = MIN(nph/2,nmmax)
     dataout(0)=CMPLX(data(0),0_dp,kind=DP)

     do i = 1, im_max*2-3, 2
        dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP) 
     enddo

     if(im_max==nph/2) then
        dataout(im_max)= CMPLX( data(nph-1),0,kind=DP)
     else
        dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
     endif

     if(im_max==nmmax) goto 1000

     do i =  im_max+1,min(nph,nmmax)
        dataout(i) = conjg(dataout(2*im_max-i) )
     end do

     if(min(nph,nmmax)==nmmax) goto 1000

     do i =  2*im_max+1,nmmax
        dataout(i) = dataout(mod(i,2*im_max)) 
     end do

1000 continue

     if(kphi0==1)then
        do i =  0,nmmax
           m = ksign*i
           dataout(i)=dataout(i)* CONJG(trig(-m,nph/4))
        enddo
     end if

     RETURN
   END subroutine ring_analysis_org_ball

!=================================

  subroutine alm2map_multi_ball_check(nsmax,nnmax,  nlmax, nmmax, almn,bb,j0,map, gln ,nmaps,nvmaps, der)

    use healpix_types
    IMPLICIT none
    INTEGER(I4B), INTENT(IN) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax, nnmax,j0,nvmaps
    INTEGER(I4B), INTENT(IN) :: nmmax,nmaps
    INTEGER(I4B), INTENT(IN), OPTIONAL :: der
    COMPLEX(DPC), INTENT(IN),  DIMENSION(0:nlmax,0:nmmax,0:nnmax) :: almn
    REAL(DP), INTENT(IN), DIMENSION(0:nlmax,0:nnmax, j0:j0+nmaps-1) :: gln
    REAL(DP),     INTENT(OUT), DIMENSION(0:12*nsmax**2-1,0:nvmaps, 0:nmaps-1) :: map

    INTEGER(I4B) :: l, m, ith, indl, scalem, scalel, i         ! alm related
    INTEGER(I4B) :: nph, kphi0, istart_south, istart_north  ! map related
    integer(i4b):: iv, vv,in,jj,nv,nvj

    REAL(DP) :: cth, sth, dth1, dth2, dst1, OVFLOW, UNFLOW, logOVFLOW,bb
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) :: f2m, fm2, fl2, corfac,vval
    COMPLEX(DPC) :: factor

    CHARACTER(LEN=7), PARAMETER :: code = 'ALMN2MAP_BALL'
    COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: b_north, b_south
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: b_n, b_s
    INTEGER(I4B) :: status

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: recfac
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    !=======================================================================



    !     --- allocates space for arrays ---
    ALLOCATE(recfac(0:nlmax,0:nmmax),stat = status)    
    if (status /= 0) stop 'recfac'

    ALLOCATE(b_north(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_north'

    ALLOCATE(b_n(0:nvmaps-1,0:nmaps),stat = status) 
    if (status /= 0) stop 'b_n'

    ALLOCATE(b_south(0:nmmax,0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_south'

    ALLOCATE(b_s(0:nvmaps,0:nmaps-1),stat = status) 
    if (status /= 0) stop 'b_s'

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) stop 'ring'


    if (feedback>2) print*, 'l320: total(gln): ',sum(gln)
    if (feedback>2) print*, 'l321: total(almn): ',sum(almn)
    if (feedback>1) print*,'l322: alm2map_multi_ball l321 A '

    !     ------------ initiate arrays ----------------

    do m = 0, nmmax
       fm2 = DBLE(m) **2
       do l = m, nlmax
          fl2 = DBLE(l+1) **2
          recfac(l,m) = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       enddo
    enddo
    !     --------------------------------------------

    !     ----- set the whole map to zero ------
    map = 0.0
    !     --------------------------------------

    istart_north = 0
    istart_south = 12*nsmax**2
    OVFLOW=exp(log(1.0d30))
    UNFLOW=exp(log(1.0d-30))
    logOVFLOW=log(1.0d30)

    dth1 = 1.d0 / (3.d0*DBLE(nsmax)**2)
    dth2 = 2.d0 / (3.d0*DBLE(nsmax))
    dst1 = 1.d0 / (SQRT(6.d0) * DBLE(nsmax) )

    do ith = 1,2*nsmax      ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.d0  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.d0 * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.D0-cth)*(1.D0+cth)) ! sin(theta)
       endif
       lam_mm = sq4pi_inv ! lambda_00
       scalem=1

       !print*,'alm2map_multi_ball l325 B m loop start '
       b_north=0
       b_south=0

       do m = 0, nmmax

          b_n = 0 !holds sum over n
          b_s = 0

          f2m = 2.d0 * m

          !           ---------- l = m ----------
          par_lm = 1.d0  ! = (-1)^(l+m)
          if (m .ge. 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.d0)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac=exp(scalem*logOVFLOW) !Actual lam_mm
          lam_lm = lam_mm*corfac/OVFLOW

          do in=0,nnmax
             
             factor = lam_lm * almn(m,m,in)
             if (present(der)) factor = lam_lm * almn(m,m,in)/sth**der
             
             do jj=0,nmaps-1
                nvj=ceiling(bb**(j0+jj))
                if (feedback>10) print*, 'm,jj,nvj: ',m,jj+j0,nvj
                do iv=1,nvj-1
                   vval= iv*(bb**(j0+jj)-1d0)/dble(nvj-1)
                   b_n(iv, jj) = b_n(iv,jj) +  factor*snrv(in, vval,j0+jj,bb)*gln(m,in, j0+jj)
                   b_s(iv, jj) = b_s(iv,jj) +  factor*snrv(in, vval,j0+jj,bb)*gln(m,in, j0+jj)

                   !b_n(iv, jj) = b_n(iv,jj) +  factor*snr(in,nvmaps,iv)*gln(m,in,j0+jj)**2  !snrv(in, vval,j0+jj,bb)
                   !b_s(iv, jj) = b_s(iv,jj) +  factor*snr(in,nvmaps,iv)*gln(m,in, j0+jj)**2 !snrv(in, vval,j0+jj,bb)
                   !print*, 'bb, j0, in, jj, iv, vval, lam_lm, abs(almn), snrv, gln: ',bb, j0,in, jj, iv, vval, lam_lm,abs(almn(in,m,m)),snrv(in, vval,j0+jj,bb),gln(m,in, j0+jj) 
                   
                enddo
             enddo
          enddo

          if (feedback>4) print*, 'l395: sum(b_n) ',sum(b_n)
          if (feedback>4) print*, 'l396: sum(b_s) ',sum(b_s)

          !           ---------- l > m ----------
          lam_0 = 0.d0
          lam_1 = 1.d0 
          scalel=0
          a_rec = recfac(m,m)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 

             do in=0,nnmax
                factor = lam_lm * almn(l,m,in)
                if (present(der)) factor = lam_lm * almn(l,m,in)/sth**der

                do jj=0,nmaps-1
                   nvj=ceiling(bb**(j0+jj))
                   do iv=1,nvj-1
                      vval= iv*(bb**(j0+jj)-1d0)/dble(nvj-1)
                      b_n(iv, jj) = b_n(iv,jj) +          factor*snrv(in, vval,j0+jj,bb)*gln(l,in,j0+jj)
                      b_s(iv, jj) = b_s(iv,jj) + par_lm * factor*snrv(in, vval,j0+jj,bb)*gln(l,in,j0+jj) 

                      ! b_n(iv, jj) = b_n(iv,jj) +          factor*snr(in,nvmaps,iv)*gln(l,in,j0+jj)**2 !snrv(in, vval,j0+jj,bb)
                      ! b_s(iv, jj) = b_s(iv,jj) + par_lm * factor*snr(in,nvmaps,iv)*gln(l,in,j0+jj)**2 !snrv(in, vval,j0+jj,bb)
                   enddo
                enddo !end loop over jj
             enddo !end loop over in

             if (feedback>4) print*, 'l418: sum(b_n) ',sum(b_n)
             if (feedback>4) print*, 'l419: sum(b_s) ',sum(b_s)

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = recfac(l,m)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) .gt. OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac=exp((scalel+scalem)*logOVFLOW)
             else
                if (abs(lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac=exp((scalel+scalem)*logOVFLOW)
                endif
             endif
          enddo

          do jj=0,nmaps-1
                nvj=ceiling(bb**(j0+jj))
             do iv=1,nvj-1
                b_north(m,iv,jj) = b_north(m,iv,jj)+b_n(iv,jj) !sum over n sum(almn)
                b_south(m,iv,jj) = b_south(m,iv,jj)+b_s(iv,jj)
             enddo
          enddo


       enddo !end loop over m

       if (feedback>2) print*, 'l448: sum(b_north) ',sum(b_north)
       if (feedback>2) print*, 'l449: sum(b_south) ',sum(b_south)

       !print*,'alm2map_multi_ball l325 C m loop done! '
       !        ---------------------------------------------------------------

       !        obtains the negative m of b(m,theta) (= complex conjugate)


       !        ---------------------------------------------------------------
       !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta)
       !        ---------------------------------------------------------------

       do i=0,nmaps-1
          nvj=ceiling(bb**(j0+i))
          do vv=1,nvj-1
             call ring_synthesis_ball(nsmax,nlmax,nmmax,b_north(:,vv,i),nph,ring(:),kphi0)   ! north hemisph. + equator
             map(istart_north:istart_north+nph-1,vv,i) = ring(0:nph-1)
             if (feedback>3) print*, 'l464: j, v, sum(b_north) ',j0+i, vv, sum(b_north(:,vv,i))
             if (feedback>3) print*, 'l465: j, v, sum(b_south) ',j0+i, vv, sum(b_south(:,vv,i))
          enddo
       enddo
       istart_north = istart_north + nph

       if (feedback>2) print*, 'l469: ith_ring north, sum(map(istart_north:istart_north+nph-1,:,:)) ',ith, sum(map(istart_north:istart_north+nph-1,:,:))


       !print*,'alm2map_multi_ball l325 north done, ith= ',ith

       istart_south = istart_south - nph
       if (ith .lt. 2*nsmax) then
          do i=0,nmaps-1
             nvj=ceiling(bb**(j0+i))
             do vv=1,nvj-1
                call ring_synthesis_ball(nsmax,nlmax,nmmax,b_south(:,vv,i),nph,ring(:),kphi0) ! south hemisph. w/o equat
                map(istart_south:istart_south+nph-1,vv,i) = ring(0:nph-1)
             enddo !end vv
          enddo
          if (feedback>2) print*, 'l482: ith_ring south, sum(map(istart_south:istart_south+nph-1,:,:)) ',ith,sum(map(istart_south:istart_south+nph-1,:,:))

          !print*,'alm2map_multi_ball l325 south done, ith= ',ith
       endif

    enddo    ! loop on cos(theta)

    if (feedback>1)  print*,'alm2map_multi_ball l325 E, multiball done! sum(map) = ',sum(map)
    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(recfac)
    DEALLOCATE(b_north)
    DEALLOCATE(b_south)
    DEALLOCATE(b_n)
    DEALLOCATE(b_s)
    DEALLOCATE(ring)
    return
  end subroutine alm2map_multi_ball_check


!========================

   subroutine create_almn (nshell,nnmax,nsmax,nlmax,nmmax,iseed, almn,cl,cn)
      USE ran_tools
      IMPLICIT none
      INTEGER(I4B), INTENT(IN) :: nshell
      INTEGER(I4B), INTENT(IN) :: nlmax
      INTEGER(I4B), INTENT(IN) :: nmmax
      INTEGER(I4B), INTENT(IN) :: nsmax
      INTEGER(I4B)             :: nnmax
      INTEGER(I4B),       INTENT(INOUT)      :: iseed



      !radial and angular PS
      REAL(DP), INTENT(IN), DIMENSION(0:nlmax) :: cl
      REAL(DP), INTENT(IN), DIMENSION(0:nshell-1) :: cn
      !
      COMPLEX(DPC),INTENT(OUT), DIMENSION(0:nshell-1,0:nlmax,0:nmmax) :: almn


      !working arrays
      !COMPLEX(DPC), DIMENSION(0:nshell-1,0:nlmax,0:nmmax) :: alm_r
      COMPLEX(DPC),DIMENSION(0:nshell-1)                  :: fr
      COMPLEX(DPC)                                        :: a_lm
      !REAL(DP), DIMENSION(0:12*nsmax**2-1) :: map


      !COMPLEX(DPC),DIMENSION(0:nlmax,0:nmmax) :: alm_TGC
      !COMPLEX(DPC),DIMENSION(0:nshell/2) :: an

      INTEGER(I4B) :: unit, i, j, in,ir,n,l, m, n_l
      REAL(DP) ::  hsqrt2
      REAL(DP) ::  rms_n,rms_tt, rms_g1, rms_g2, rms_cc
      REAL(DP) ::  zeta1_r, zeta1_i, zeta2_r, zeta2_i, zeta3_r, zeta3_i
      REAL(DP), DIMENSION(:),   ALLOCATABLE :: cls_tt, cls_tg, cls_gg, cls_cc

      INTEGER(I4B) :: status


      !alm_TGC = CMPLX(0.0,0.0)
      almn = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
      fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
      hsqrt2 = 1d0/ sqrt(2d0)


      ! if (nnmax>nshell/2) nnmax=nshell/2
      ! do in = 0,nnmax
      !    rms_n = 0d0
      !    if (cn(in) .gt. 1e-13) then
      !       rms_n   = sqrt(cn(in))
      !       zeta1_r = randgauss_boxmuller(iseed) * hsqrt2
      !       zeta1_i = randgauss_boxmuller(iseed) * hsqrt2
      !       an(in) = CMPLX(zeta1_r, zeta1_i) * rms_n
      !    endif
      ! enddo
      
      !print*,'sum(an)',sum(an)
      ! !     --- smoothes the initial power spectrum ---
      ! do l = 0, nlmax
      !    cls_tt(l) = cls_tt(l) * window(l,0)**2
      !    if (polar.gt.1) cls_tg(l) = cls_tg(l) * window(l,1)**2
      !    if (polar.gt.2)cls_gg(l) = cls_gg(l) * window(l,2)**2
      !    if (polar.gt.3)cls_cc(l) = cls_cc(l) * window(l,3)**2
      ! enddo

      !     --- generates randomly the alm according to their power spectrum ---
      !     alm_T = zeta1 * rms_tt
      !     alm_G = zeta1 * rms_g1 + zeta2 * rms_g2
      !     alm_C = zeta3 * rms_cc

      print*, 'almn_create: gen alm per shell '

      !do ir=1,nshell-2  !first and last shells are set 0 by hand here

         do l = 2, nlmax
            rms_tt = 0d0
            
            if (cl(l) .gt. 1e-13) then
               rms_tt   = sqrt(cl(l))
               
               !        ------ m = 0 ------
               zeta1_r = randgauss_boxmuller(iseed)  ! gaussian deviate based on ran_mwc
               zeta1_i = 0d0
               a_lm   = CMPLX(zeta1_r, zeta1_i) * rms_tt ! T
               almn(0:nshell-1,l,0) = a_lm*cn(0:nshell-1)
               
               !        ------ m > 0 ------
               do m = 1,l
                  zeta1_r = randgauss_boxmuller(iseed) * hsqrt2
                  zeta1_i = randgauss_boxmuller(iseed) * hsqrt2
                  a_lm   = CMPLX(zeta1_r, zeta1_i) * rms_tt ! T
                  almn(0:nshell-1,l,m) = a_lm*cn(0:nshell-1)
               enddo !m loop
               
            endif  !rms_tt is not zero
            
         enddo !l loop



         !call alm2map(nsmax,nlmax,nmmax,almn(20:20,:,:),map)
         !print*,'sim ballmap(0:100,20) = ',map(0:100)


         
      !enddo !end ir

      !print*,abs(almn(20,4,:))
      print*, 'almn_create: calling complex_fft '
      !print*, 'sim almnr(0:nshell-1,20,2)',almn(0:nshell-1,20,2)

      do l=0,nmmax
         do m=0,l
           fr(0:nshell-1)=almn(0:nshell-1,l,m)
           call complex_fft(fr,backward=.false.) 
           almn(0:nshell-1,l,m) = fr(0:nshell-1)/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
        enddo
     enddo

      !print*,abs(almn(20,4,:))


      !print*,'sum(alm)',sum(alm_TGC)
      
      ! do ir=2,nshell-1
      !    if (ir<=nshell/2) then
      !       almn(ir,0:nlmax,0:nmmax) = alm_TGC(0:nlmax,0:nmmax)*(an(ir))
      !    else
      !       almn(ir,0:nlmax,0:nmmax) = alm_TGC(0:nlmax,0:nmmax)*conjg(an(nshell-ir))
      !    endif
      ! enddo

      !print*,'sum(almn)',sum(almn)

      return
   end subroutine create_almn


subroutine read_ringweights(nside, nmaps, healpix_data, weights )
    implicit none

    integer(i4b),intent(in)  :: nside, nmaps
    character(len=128), intent(in) :: healpix_data
 !   real(dp), pointer, dimension(:,:) :: weights
    real(dp), dimension(1:2*nside,nmaps), intent(out) :: weights
    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

   ! allocate(weights(1:2*nside,nmaps))
   

    call int2string(nside, nside_text)
    weight_file = TRIM(healpix_data)//"/data/weight_ring_n" // nside_text // ".fits"
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
     !  call read_dbintab(weight_file, weights, 2*nside, nmaps, nullval, anynull)
       call input_map(weight_file, weights, 2*nside, nmaps)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if

  
  end subroutine read_ringweights

  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string
    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string

  function num2str(integer)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=30) :: string
    character(len=30) :: num2str
    integer(i4b)               :: temp_int, i, k

    write(string,*)integer
    num2str=trim(adjustl(string))
    return
  end function num2str
! !===================================
!   subroutine calc_smhw_ball(nside,lmax,nscales,scales,g_real_tot,loadwav,savewav,filename)

!     use pix_tools
!     use alm_tools

!     IMPLICIT NONE


!     INTEGER(I4B), INTENT(IN) :: nside,lmax,nscales
!     REAL(DP),DIMENSION(0:lmax,0:nscales-1), INTENT(OUT) :: g_real_tot
!     CHARACTER(LEN=128), INTENT(IN) :: filename
!     LOGICAL(LGT), INTENT(IN) :: loadwav,savewav
!     REAL(DP),DIMENSION(0:nscales-1), INTENT(IN) :: scales
!     COMPLEX(SPC), ALLOCATABLE, DIMENSION(:,:,:) :: almw
!     REAL(SP),ALLOCATABLE,DIMENSION(:) ::wav_pix
!     REAL(DP) :: r,nr,y,theta,phi
!     INTEGER(I4B) :: n_pix_T,i_pix,irisol,i
!     REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
!     REAL(DP), DIMENSION(1:2) :: z = 0.0d0


!     if (loadwav) then

!        open(12,file=filename,status='old',form='unformatted')
!        read(12) g_real_tot
!        close(12)

!     else


!        n_pix_T=nside**2*12

!        ALLOCATE(w8ring(1:2*nside,1:1))
!        w8ring=1d0

!        ALLOCATE(wav_pix(0:n_pix_T-1))
!        ALLOCATE(almw(1:1,0:lmax,0:lmax))

!        DO irisol=0,nscales-1

!           !!     r = (irisol+1)*r0

!           if (scales(irisol).ne.0.) then

!              r=scales(irisol)/60.d0/180.d0*pi
!              nr=r*SQRT(1.+r**2./2.+r**4./4.) 

!              DO i_pix = 0,n_pix_T-1 

!                 CALL pix2ang_ring(nside,i_pix,theta,phi)

!                 y =2.*SIN(theta/2.)/COS(theta/2.)

!                 if ((-y**2./2./r**2.).gt.(-700.d0)) then
!                    wav_pix(i_pix)=1./nr/SQRT(2.*pi)*(1.+(y/2.)**2.)**2.*(2.- (y/r)**2. )*EXP(-y**2./2./r**2.)
!                 else
!                    wav_pix(i_pix)=0.d0
!                 endif


!              END DO     !i_pix


!              CALL map2alm(nside, lmax, lmax, wav_pix, almw, z, w8ring)  

!              g_real_tot(:,irisol) = REAL(almw(1,:,0),DP)


!              DO i=0,lmax
!                 g_real_tot(i,irisol) = g_real_tot(i,irisol)/SQRT(2d0*REAL(i,DP)+1)*SQRT(4d0*pi)
!              END DO

!           else
!              g_real_tot(0:lmax,irisol)=1.
!           endif


!        END DO  !irisol

!        if (savewav) then
!           open(12,file=filename,status='unknown',form='unformatted')
!           write(12) g_real_tot
!           close(12)
!        endif

!        DEALLOCATE(wav_pix,almw,w8ring)

!     endif


!   END subroutine calc_smhw_ball



  !=================================================

end module wav_ball_mod


