  
!!!IDL extension to HEALPix [TESTVERSION] 2D needlet: Written by Frode K. Hansen, Nov. 2006
  !! Modified for 3D needlet by Y. Fantaye Sep 2013 
  
  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp


  IMPLICIT NONE

  !healpix
  INTEGER(I4B) :: nside,npix,lmax,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx,glnpow
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring,pixelw
  REAL(DP), DIMENSION(1:2) :: zbounds


  !power spectra
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cl,cn

  !ball
  INTEGER(I4B) :: nmax,nvmax,nshell,in,iv,ir
  INTEGER(I4B) :: j0,nj,first_j,last_j,first_v, last_v
  integer(i4b) :: nrpix,sn,s0

  !general 
  INTEGER(I4B) :: i,j,k,ii,jj,kk,iseed=-1234
  integer(i4b) ::oddInd, evenInd

  !needlet
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln

  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm
  REAL(DP) ::rmin,rmax,dr,rmid,y

  !timing 
  real     (SP) :: clock_time, clock_time2,time0, time1
  real     (SP) :: ptime, ptime0, ptime1

  !timing and accuracy array
  integer, parameter ::  nnsim=6
  integer::  isim,nsim
  REAL(DP), DIMENSION(nnsim,5) :: metric

  !MPI integers
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nct_pp,nmap_pp
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit,stat
  real(dp):: s2n, s2nm1

  !maps and almn
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(Dp), dimension(:,:), allocatable::ballmap,map
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map1,map2,map3
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_ball,almn_orig
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr


  CHARACTER(LEN=200) :: fname,pfile
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn
  CHARACTER(LEN=80),dimension(60)::header


  me=0
  unit=12+me

  call wall_clock_time(time0)
  call cpu_time(ptime0)

  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'

  ! if (iargc() > 0) then
  !    call getarg(1,pfile)
  ! else
  !    pfile=trim(adjustl(dir))//'params_ball2beta.unf'
  ! endif

  ! open(12,file=pfile,status='old',form='unformatted')
  ! read(12) glnpow
  ! read(12) bb
  ! read(12) j0
  ! read(12) nj
  ! read(12) nside
  ! read(12) polar
  ! read(12) nmax
  ! read(12) lmax
  ! read(12) nshell
  ! read(12) feedback
  ! close(12)

  !!necessary to be defined params
  polar=0
  wav=0
  p=0
  m1=1
  m2=0



  ! print*, '**************** ball2almn.f90 *****'
  ! print*, '---------------------------'
  ! print*, 'nside = ',nside 
  ! print*, 'polar = ',polar 
  ! print*, 'lmax = ',lmax
  ! print*, 'nshell = ',nshell
  ! print*, 'rmin, rmax = ',real(rmin),real(rmax)
  ! print*, 'bb,glnpow = ',bb,glnpow
  ! print*, 'npix, nj, j0 = ',npix,nj,j0
  ! print*, '---------------------------'



  zbounds = (/ -1.0_dp, 1.0_dp /)
  !compute needlet weight gln
  CALL calc_f2


  nsim=nnsim

  do isim=1,nsim

     feedback=1
     glnpow=1
     bb=4d0
     j0=-1

     nside=2**(isim+3)
     lmax=2**(isim+4)

     nshell = 2**(isim+5)
     nmax=nshell/2

     nj = ceiling( log( sqrt( dble(lmax*(lmax+1)+nmax**2) ) )/log(bb) + 1 )

     npix=12*nside**2


     print*,'==========================='
     print*,''
     print*, '************** isim **** ',isim
     print*, 'glnpow, j0,bb,nj: ',glnpow, j0,bb,nj
     print*, 'lmax, nshell, nside: ',lmax, nshell, nside
     print*,''
     print*,'==========================='

     !!======================== ====
     print*, 'generating almn'
     !!======================== ====

     ALLOCATE(w8ring(1:2*nside,1:1+2*polar))
     allocate(cl(0:lmax),cn(0:nshell-1))
     allocate(map(0:npix-1,1))
     ALLOCATE(beta_ballmap(0:npix-1,0:nshell-1, 0:nj-1))
     allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
     allocate(alm_TGC(1:1,0:lmax,0:lmax))
     allocate(almn_orig(0:nshell-1,0:lmax,0:lmax))
     allocate(almn_ball(0:nshell-1,0:lmax,0:lmax))
     allocate(fr(0:nshell-1))
     ALLOCATE(gln(0:lmax,0:nshell-1,0:nj-1))

     !===========

     w8ring=1d0
     call read_ringweights(nside, 1 ,healpixdir, w8ring )
     ! open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
     ! read(12) w8ring
     ! close(12)

     cl=1d0
     !qudratic scaling
     do i=0,nshell-1
        cn(i)=sin(2d0*pi*3*dble(i)/dble(nshell-1))+sin(2d0*pi*4*dble(i)/dble(nshell-1))
     enddo


     cl(0)=0d0
     cl(1)=0d0


     alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     alm_TGC_r=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     almn_orig = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)


     call create_almn (nshell,nmax,lmax,lmax,iseed, almn_orig,cl,cn)
     almn_ball=almn_orig

     fname = trim(adjustl(dir))//'almn_sim.unf'
     print*, 'saving almn sim: ',trim(fname)
     print*, 'sshape(almn) = ',shape(almn_ball)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_ball
     close(unit)

     print*, sum(almn_ball)


     !!======================== ====
     print*, 'COMPUTING almn2beta ..'
     !!======================== ====


     if (wav.eq.0) wavtyp='standard'
     CALL calc_gln(j0,nj,lmax,nshell-1,gln,bb,wavtyp,p,gln4fft=1)
     if (glnpow>1) then
        if  (me == 0) print*, 'glnpow=2 .. using gln=gln^2 ..'
        gln=gln*gln
     endif

     call alm2map_fftw_ball(nside,nshell-1,lmax,lmax, almn_ball,bb,j0,beta_ballmap, gln(:,:, 0:nj-1),nj)

     call wall_clock_time(time1)
     call cpu_time(ptime1)
     clock_time = time1 - time0
     ptime      = ptime1 - ptime0
     write(*,*) "ball2beta: Clock and CPU time [s]: ", clock_time, ptime






     !=========================================================
     !================== .......beta2aln ....==================
     print*, 'COMPUTING beta2ball ..'
     !=========================================================
     !========================================================

     call wall_clock_time(time0)
     call cpu_time(ptime0)

     !=========================================================
     print*, 'COMPUTING beta2almn ..'
     !=========================================================



     almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     do j=0,nj-1

        map=0.
        do iv=0,nshell-1
           map(0:npix-1,1) = beta_ballmap(0:npix-1,iv, j)
           alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
           if (count( map(:,1)>1e-13 )>10)  call map2alm_iterative(nside, lmax, lmax, 3, map, alm_TGC, zbounds, w8ring) !sum over k
           alm_TGC_r(iv,0:lmax,0:lmax) = alm_TGC(1,0:lmax,0:lmax)
        enddo

        do m=0,lmax        
           do l=0,lmax
              fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
              call complex_fft(fr,backward=.false.) !sum over v
              alm_TGC_r(0:nshell-1,l,m) = fr(0:nshell-1)*gln(l,0:nshell-1,j)/dble(nshell) !an for a given l,m and j 
           enddo
        enddo

        almn_ball = almn_ball + alm_TGC_r !sum over local j

     enddo !end sum over j


     fname = trim(adjustl(dir))//'almn_rec.unf'
     print*, 'saving almn rec: ',trim(fname)
     print*, 'sshape(almn) = ',shape(almn_ball)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_ball
     close(unit)

     print*, 'maxval almn_ball-almn_orig', maxval(abs(almn_ball(2:nshell/2-1,2:lmax,2:lmax)-almn_orig(2:nshell/2-1,2:lmax,2:lmax)))

     ! !=========================================================
     !   print*, 'COMPUTING almn2almr ..'
     ! !=========================================================

     !   do m=0,lmax        
     !      do l=0,lmax
     !         fr(0:nshell-1)=almn_ball(0:nshell-1,l,m)
     !         call complex_fft(fr,backward=.true.) !sum over n
     !         alm_TGC_r(0:nshell-1,l,m) = fr(0:nshell-1)
     !      enddo
     !   enddo

     ! !=========================================================
     !   print*, 'COMPUTING almr2ball ..'
     ! !=========================================================

     !   do ir=0,nshell-1
     !      alm_TGC(1,0:lmax,0:lmax) =  alm_TGC_r(ir,0:lmax,0:lmax)
     !      call alm2map(nside,lmax,lmax,alm_TGC,ballmap(:,ir)) !sum over lm
     !   enddo





     call wall_clock_time(time1)
     call cpu_time(ptime1)
     clock_time2 = time1 - time0
     ptime      = ptime1 - ptime0
     write(*,*) "beta2ball: Clock and CPU time [s]: ", clock_time2, ptime


     metric(isim,1)=dble(nshell*npix)
     metric(isim,2)=dble(lmax)
     metric(isim,3)=dble(nshell/2)
     metric(isim,4)=dble(clock_time+clock_time2)/2d0
     metric(isim,5)=maxval(abs(almn_ball(2:nshell/2-1,2:lmax,2:lmax)-almn_orig(2:nshell/2-1,2:lmax,2:lmax)))


     deALLOCATE(w8ring)
     deallocate(cl,cn)
     deallocate(map)
     deALLOCATE(beta_ballmap)
     deallocate(alm_TGC_r)
     deallocate(alm_TGC)
     deallocate(almn_orig)
     deallocate(almn_ball)
     deallocate(fr)
     deALLOCATE(gln)


  enddo

     fname = trim(adjustl(dir))//'metric_almn2beta_beta2almn.unf'
     print*, 'saving metric: ',trim(fname)
     print*, 'sshape(almn) = ',shape(metric)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) metric
     close(unit)

end program

