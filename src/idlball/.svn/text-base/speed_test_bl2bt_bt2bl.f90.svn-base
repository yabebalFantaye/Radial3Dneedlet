  
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
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
  REAL(DP), DIMENSION(1:2) :: zbounds

  !ball
  INTEGER(I4B) :: nmax,nvmax,nshell,in,iv,ir
  INTEGER(I4B) :: j0,nj,first_j,last_j,first_v, last_v
  integer(i4b) :: nrpix,sn,s0

  !general 
  INTEGER(I4B) :: i,j,k,ii,jj,kk
  integer(i4b) ::oddInd, evenInd

  !needlet
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln

  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm
  REAL(DP) ::rmin,rmax,dr,rmid,y
  
  !timing 
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1

  !MPI integers
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nct_pp,nmap_pp
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit,stat
  real(dp):: s2n, s2nm1

  !maps and almn
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(Dp), dimension(:,:), allocatable::ballmap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map,map1,map2,map3
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_ball
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr


  CHARACTER(LEN=200) :: fname,pfile
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn

  me=0
  unit=12+me

  call wall_clock_time(time0)
  call cpu_time(ptime0)

  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'

  if (iargc() > 0) then
     call getarg(1,pfile)
  else
     pfile=trim(adjustl(dir))//'params_ball2beta.unf'
  endif




  open(12,file=pfile,status='old',form='unformatted')
  read(12) glnpow
  read(12) bb
  read(12) j0
  read(12) nj
  read(12) nside
  read(12) polar
  read(12) nmax
  read(12) lmax
  read(12) nshell
  read(12) feedback
  close(12)

  !!necessary to be defined params
  polar=0
  wav=0
  p=0
  m1=1
  m2=0

  npix=12*nside**2



  print*, '**************** ball2almn.f90 *****'
  print*, '---------------------------'
  print*, 'nside = ',nside 
  print*, 'polar = ',polar 
  print*, 'lmax = ',lmax
  print*, 'nshell = ',nshell
  print*, 'rmin, rmax = ',real(rmin),real(rmax)
  print*, 'bb,glnpow = ',bb,glnpow
  print*, 'npix, nj, j0 = ',npix,nj,j0
  print*, '---------------------------'
  

  !!======================== ====
  print*, 'reading ball map '
  !!======================== ====

  ALLOCATE(w8ring(1:2*nside,1:1+2*polar))
  w8ring=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)


  ALLOCATE(ballmap(0:npix-1,0:nshell-1))
  open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
  read(12) ballmap
  close(12)


  !!======================== ====
  print*, 'computing ball(:,ir) to almr ..'
  !!======================== ====

  !allocate big arrays
  allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
  alm_TGC_r=CMPLX(0.0_dp, 0.0_dp, KIND=DP)


  do i=0,nshell-1
     call map2alm(nside, lmax, lmax, ballmap(:,i), alm_TGC, zbounds, w8ring)
     alm_TGC_r(i,0:lmax,0:lmax) = alm_TGC(1,0:lmax,0:lmax)
  enddo

  !!======================== ====
  print*, 'computing almr to almn'
  !!======================== ====

  !allocate big array 
  allocate(almn_ball(0:nshell-1,0:lmax,0:lmax))
  allocate(fr(0:nshell-1))
  almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
  fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
  
  do m=0,lmax 
     do l=0,lmax
        fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
        call complex_fft(fr,backward=.false.) 
        almn_ball(0:nshell-1,l,m) = fr(0:nshell-1)/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
     enddo
  enddo

  

  !!======================== ====
  print*, 'COMPUTING almn2beta ..'
  !!======================== ====


  !compute needlet weight gln
  ALLOCATE(gln(0:lmax,0:nshell-1,0:nj-1))
  CALL calc_f2
  if (wav.eq.0) wavtyp='standard'
  CALL calc_gln(j0,nj,lmax,nshell-1,gln,bb,wavtyp,p,gln4fft=1)
  if (glnpow>1) then
     if  (me == 0) print*, 'glnpow=2 .. using gln=gln^2 ..'
     gln=gln*gln
  endif


  ALLOCATE(beta_ballmap(0:npix-1,0:nshell-1, 0:nj-1))
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

  allocate(map(0:npix-1))

  almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)
  do j=0,nj-1
     
     map=0.
     do iv=0,nshell-1
        map(0:npix-1) = beta_ballmap(0:npix-1,iv, j)
        call map2alm(nside, lmax, lmax, map, alm_TGC, zbounds, w8ring) !sum over k
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


!=========================================================
  print*, 'COMPUTING almn2almr ..'
!=========================================================
  
  do m=0,lmax        
     do l=0,lmax
        fr(0:nshell-1)=almn_ball(0:nshell-1,l,m)
        call complex_fft(fr,backward=.true.) !sum over n
        alm_TGC_r(0:nshell-1,l,m) = fr(0:nshell-1)
     enddo
  enddo

!=========================================================
  print*, 'COMPUTING almr2ball ..'
!=========================================================

  do ir=0,nshell-1
     alm_TGC(1,0:lmax,0:lmax) =  alm_TGC_r(ir,0:lmax,0:lmax)
     call alm2map(nside,lmax,lmax,alm_TGC,ballmap(:,ir)) !sum over lm
  enddo

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) "beta2ball: Clock and CPU time [s]: ", clock_time, ptime


end program

