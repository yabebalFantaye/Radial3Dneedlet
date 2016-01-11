  
!!!IDL extension to HEALPix [TESTVERSION] 2D needlet: Written by Frode K. Hansen, Nov. 2006
  !! Modified for 3D needlet by Y. Fantaye Sep 2013 
  
  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp
  use rngmod !  USE ran_tools


  IMPLICIT NONE

  !healpix
  INTEGER(I4B) :: nside,npix,lmax,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx,glnpow
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring,pixelw
  REAL(DP), DIMENSION(1:2) :: zbounds
  type(planck_rng) :: rng_handle


  !power spectra
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cn

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
  integer, parameter ::  nnsim=3
  integer::  isim,nsim
  REAL(DP), DIMENSION(nnsim,4) :: metric

  !MPI integers
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nct_pp,nmap_pp
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit,stat
  real(dp):: s2n, s2nm1

  !maps and almn
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(Dp), dimension(:,:), allocatable::ballmap

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map1,map2,map3
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC,alm1,alm2

  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_ball,almn_orig
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr


  CHARACTER(LEN=200) :: fname,pfile
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn
  CHARACTER(LEN=80),dimension(1:60)::header


  me=0
  unit=12+me

  call wall_clock_time(time0)
  call cpu_time(ptime0)

  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'


  zbounds = (/ -1.0_dp, 1.0_dp /)


  polar=0
  lmax=128
  nside=64
  npix=12*nside**2


  print*, 'allocating ..'
  allocate(pixelw(0:lmax,1:1))
  ALLOCATE(w8ring(1:2*nside,1:1+2*polar))
  allocate(cl(0:lmax,1:1))
  allocate(map(0:npix-1,1))
  allocate(alm1(1:1,0:lmax,0:lmax))
  allocate(alm2(1:1,0:lmax,0:lmax))
  print*, 'allocating done'

  !!initialize Cl
  cl=1d0
  cl(0,1)=0d0
  cl(1,1)=0d0
  !print*, 'writing cl'
  !call write_asctab (cl,lmax,1,header,60,'temp/cl.fits')
  !print*, 'cl file written to temp/cl.fits'

  !!read weight
  w8ring=1d0
  call read_ringweights(nside, 1 ,healpixdir, w8ring )

  !!read pixel window
  call pixel_window( pixelw, nside) 
  print*, 'pixel window read '

  !!generate alm
  call rand_init(rng_handle, -1) 
  call create_alm(nside, lmax, lmax, polar, '../cls/cl.fits', rng_handle, 0d0, alm1, header, trim(healpixdir)//'/data/pixel_window_n0064.fits')
  print*, 'alm created'
  print*,'sum(alm1) = ',sum(alm1)


  print*, 'lmax, nside: ',lmax,nside
  fname = trim(adjustl(dir))//'alm_sim.unf'
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) alm1
  close(unit)

  call alm2map(nside,lmax,lmax,alm1,map(:,1)) 
  call map2alm_iterative(nside, lmax, lmax, 3,map, alm2, zbounds, w8ring) 

  fname = trim(adjustl(dir))//'alm_rec.unf'
  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) alm2
  close(unit)


  print*, '-----------------'
  print*, 'max val diff(alm): ',maxval(abs(alm1-alm2))
  print*, '-----------------'

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) "beta2ball: Clock and CPU time [s]: ", clock_time, ptime
  
  
  ! fname = trim(adjustl(dir))//'metric_healpix.unf'
  ! print*, 'saving metric: ',trim(fname)
  ! print*, 'sshape(metric) = ',shape(metric)
  ! open(unit,file=trim(fname),status='unknown',form='unformatted')
  ! write(unit) metric
  ! close(unit)

  deallocate(pixelw)
  deALLOCATE(w8ring)

  deallocate(cl)
  deallocate(map)
  deallocate(alm1)
  deallocate(alm2)

  
end program

