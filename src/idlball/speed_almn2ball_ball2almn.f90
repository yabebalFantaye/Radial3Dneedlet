  
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
  integer, parameter ::  nnsim=4
  integer::  imc, isim,nsim
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

do imc = 1,20

   print*,'******** imc = ****',imc
  k=0
  kk=0

  do isim=1,nsim


     nside=2**7  !base is nside=32
     nshell = 2**7  !base is nshell=128


     !k=k+1  !nside increases 
     if (isim<=nsim/2) then 
        k=k+1
        nside=2**(k+4)
     endif

     !kk=kk+1  !nshell increases 
     if (isim>nsim/2) then 
        kk=kk+1
        nshell = 2**(kk+6)
     endif



     lmax=2*nside
     nmax=nshell/2

     feedback=1
     glnpow=1
     bb=4d0
     j0=-1



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
     print*, 'initialize time...'
     !!======================== ====

     time0=0
     ptime0=0
     time1=0
     ptime1=0
     call wall_clock_time(time0)
     call cpu_time(ptime0)



     !!======================== ====
     Print*, 'initializing arrays...'
     !!======================== ====

     !!read pixel window
     !allocate(pixelw(0:lmax,1:1))
     !call pixel_window( pixelw, nside) 
     !print*, 'create_almn will multiply alm by pixel window'



     allocate(cl(0:lmax),cn(0:nshell-1))
     cl=1d0
     !qudratic scaling
     do i=0,nshell-1
        cn(i)=200*sin(2d0*pi*3*dble(i)/dble(nshell-1))+sin(2d0*pi*4*dble(i)/dble(nshell-1))
     enddo

     cl(0)=0d0
     cl(1)=0d0



     !!======================== ====
     print*, 'simulating almn...'
     !!======================== ====


     allocate(almn_ball(0:nshell-1,0:lmax,0:lmax))
     almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

     call create_almn (nshell,nmax,nside,lmax,lmax,iseed, almn_ball,cl,cn)

     !=========================================================
     print*, 'saving almn_ball_orig ..'
     !=========================================================

     fname = trim(adjustl(dir))//'almn_sim.unf'
     print*, 'saving almn sim: ',trim(fname)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_ball
     close(unit)

     print*, 'sum(almn_sim): ',sum(almn_ball)

     deallocate(cl,cn)

     !!======================== ====
     print*, 'COMPUTING almn to almr ..'
     !!======================== ====

     allocate(fr(0:nshell-1))
     fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     do m=0,lmax 
        do l=0,lmax
           fr(0:nshell-1)=almn_ball(0:nshell-1,l,m)
           call complex_fft(fr,backward=.true.) 
           almn_ball(0:nshell-1,l,m) = fr(0:nshell-1)
        enddo
     enddo
     deallocate(fr)

     ! print*, 'rec almnr(0:nshell-1,20,2)',almn_ball(0:nshell-1,20,2)

     !!======================== ====
     print*, 'COMPUTING almr to ball ..'
     !!======================== ====


     ALLOCATE(ballmap(0:npix-1,0:nshell-1))
     ballmap=0d0

     allocate(alm_TGC(1:1,0:lmax,0:lmax))
     alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)

     do ii=0,nshell-1
        alm_TGC(1,0:lmax,0:lmax)=almn_ball(ii,0:lmax,0:lmax)
        call alm2map(nside,lmax,lmax,alm_TGC ,ballmap(:,ii)) 
        almn_ball(ii,0:lmax,0:lmax) = alm_TGC(1,0:lmax,0:lmax) 
     enddo
     deallocate(alm_TGC)
     !print*,'rec ballmap(0:100,20) = ',ballmap(0:100,20)


     !!======================== ====
     print*, 'time it takes for almn2ball ..'
     !!======================== ====

     call wall_clock_time(time1)
     call cpu_time(ptime1)
     clock_time = time1 - time0
     ptime      = ptime1 - ptime0
     write(*,*) "almn2ball: Clock and CPU time [s]: ", clock_time, ptime



     !=========================================================
     !================== .......ball2almn ....==================
     print*, 'COMPUTING ball2almn ..'
     !=========================================================
     !========================================================

     time0=0
     ptime0=0
     time1=0
     ptime1=0
     call wall_clock_time(time0)
     call cpu_time(ptime0)

     !=========================================================
     print*, 'COMPUTING ball to almr ..'
     !=========================================================



     ALLOCATE(w8ring(1:2*nside,1:1+2*polar))
     w8ring=1d0
     call read_ringweights(nside, 1 ,healpixdir, w8ring )

     allocate(map(0:npix-1,1))
     allocate(alm_TGC(1:1,0:lmax,0:lmax))

     do iv=0,nshell-1

        alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
        map(0:npix-1,1) = ballmap(0:npix-1,iv)
        if (count( map(:,1)>1e-13 )>10)  call map2alm_iterative(nside, lmax, lmax, 4, map, alm_TGC, zbounds, w8ring) !sum over k
        almn_ball(iv,0:lmax,0:lmax) = alm_TGC(1,0:lmax,0:lmax)

     enddo

     deallocate(map)
     deallocate(alm_TGC)
     deALLOCATE(w8ring)
     deALLOCATE(ballmap)


     !print*, 'rec2 almnr(0:nshell-1,20,2)',almn_ball(0:nshell-1,20,2)

     allocate(fr(0:nshell-1))
     fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     do m=0,lmax        
        do l=0,lmax
           fr(0:nshell-1)=almn_ball(0:nshell-1,l,m)
           call complex_fft(fr,backward=.false.) !sum over v
           almn_ball(0:nshell-1,l,m) = fr(0:nshell-1)/dble(nshell) !an for a given l,m 
        enddo
     enddo
     deallocate(fr)


     !=========================================================
     print*, 'saving almn_ball_rec ..'
     !=========================================================


     fname = trim(adjustl(dir))//'almn_rec.unf'
     print*, 'saving almn rec: ',trim(fname)
     print*, 'sshape(almn) = ',shape(almn_ball)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_ball
     close(unit)



     call wall_clock_time(time1)
     call cpu_time(ptime1)
     clock_time2 = time1 - time0
     ptime      = ptime1 - ptime0
     write(*,*) "beta2ball: Clock and CPU time [s]: ", clock_time2, ptime


     !=========================================================
     print*, 'recording speed and  accuracy ..'
     !=========================================================
     

     allocate(almn_orig(0:nshell-1,0:lmax,0:lmax))
     almn_orig = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

     fname = trim(adjustl(dir))//'almn_sim.unf'
     open(unit,file=trim(fname),status='old',form='unformatted')
     read(unit) almn_orig
     close(unit)

     print*, 'sum(almn_sim): ',sum(almn_orig)


     print*, '****** maxval almn_ball-almn_orig', maxval(abs(almn_ball(2:nshell/2-1,2:lmax,2:lmax)-almn_orig(2:nshell/2-1,2:lmax,2:lmax)))


     metric(isim,1)=dble(nshell*npix)
     metric(isim,2)=dble(lmax)
     metric(isim,3)=dble(nshell/2)
     metric(isim,4)=dble(clock_time+clock_time2)/2d0
     metric(isim,5)=maxval(abs(almn_ball(2:nshell/2-1,2:lmax,2:lmax)-almn_orig(2:nshell/2-1,2:lmax,2:lmax)))


     deallocate(almn_orig)
     deallocate(almn_ball)




  enddo !end isim 

  fname = '../speed_measure/metric_almn2ball_ball2almn_r7_8_ns5_6_sim'//trim(adjustl(num2str(imc)))//'.unf'
  print*, 'saving metric: ',trim(fname)
  print*, 'sshape(almn) = ',shape(metric)
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) metric
  close(unit)


enddo !

end program

