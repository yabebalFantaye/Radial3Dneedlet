  
  
!!!IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod
  use sincm

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,j0,nj,p,wav,nside_out
  INTEGER(I4B) :: npix_out,large,m1,m2,first_j,last_j,savealm
  INTEGER(I4B) :: nv, v0, nmax,nvj,nshell,ir, in,iv,nrpix
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: map_slet
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_ball,map
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map1,map2,map3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam,rvec

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE ::gln
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almv,anlm,anlm_all
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: anlmj
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fn
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr
  REAL(DP), DIMENSION(1:2) :: zbounds

  REAL(DP),DIMENSION(:,:),allocatable :: w8ring_map

  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm,rmin,rmax,rmid,vval,blnj_sn,temp
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit

  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)


  call wall_clock_time(time0)
  call cpu_time(ptime0)

  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=trim(adjustl(dir))//'params_beta2ball.unf',status='old',form='unformatted')
  read(12) bb
  read(12) j0
  read(12) nj
  read(12) nside
  read(12) nside_out
  read(12) polar
  read(12) lmax
  read(12) nmax
  read(12) nshell
  read(12) feedback
  close(12)



  polar=0
  wav=0
  p=0
  m1=0
  m2=0

  npix=nside**2*12
  npix_out=nside_out**2*12


  if (me == 0) then
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'nside, nside_out = ',nside,nside_out
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nmax
     print*, 'nshell = ',nshell
     print*, 'nj, j0 = ',nj,j0
     print*, '---------------------------'
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'


     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) ' Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif



  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1

  if (feedback>2) print*,'proc:',me,nj_pp(me),first_j,last_j



  ALLOCATE( w8ring_map(2*nside,1+2*polar) )
  w8ring_map=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring_map
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)



  allocate(map_slet(0:npix-1,0:nshell-1,0:nj_pp(me)-1))

  do j=0,nj_pp(me)-1
     !print*,me,first_j-j0+j
     WRITE(filenum,'(I4)') first_j-j0+j !+1000
     open(12,file=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'.unf',status='old',form='unformatted')
     read(12) map_slet(:,0:nshell-1,j)
     close(12)
  enddo

  if (feedback>2) print*,'beta map read!'
  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  ALLOCATE(gln(0:lmax,0:nshell-1,j0:j0+nj-1))
  CALL calc_f2
  !print*,'beta2map1:',j0,nj,lmax,bb,wavtyp,p
  CALL calc_gln(j0,nj,lmax,nshell-1,gln,bb,wavtyp,p,gln4fft=1)
  if (feedback>2) print*,'gln computed!'



  ALLOCATE(map(0:npix-1,1))
  ALLOCATE(alm(1:1+2*polar,0:lmax,0:lmax))
  ALLOCATE(almv(0:nshell-1,0:lmax,0:lmax))
  ALLOCATE(anlm(0:nshell-1,0:lmax,0:lmax))
  ALLOCATE(fn(0:nshell-1))

  anlm=CMPLX(0.0_dp, 0.0_dp, KIND=DP)

  do j=0,nj_pp(me)-1 
     
     map=0.
     do iv=0,nshell-1
        map(0:npix-1,1) = map_slet(0:npix-1,iv, j)
        if (count( map(:,1)>1e-13 )>10) call map2alm_iterative(nside, lmax, lmax,3, map, alm, zbounds, w8ring_map) !sum over k
        almv(iv,0:lmax,0:lmax) = alm(1,0:lmax,0:lmax)
     enddo
     
     do m=0,lmax        
        do l=0,lmax
           fn(0:nshell-1)=almv(0:nshell-1,l,m)
           call complex_fft(fn,backward=.false.) !sum over v
           almv(0:nshell-1,l,m) = fn(0:nshell-1)*gln(l,0:nshell-1,first_j+j)/dble(nshell) !an for a given l,m and j 
        enddo
     enddo

     anlm = anlm + almv !sum over local j

  enddo !end local j

  deALLOCATE(w8ring_map)
  deALLOCATE(map_slet)
  deALLOCATE(gln)
  deALLOCATE(almv)
  deALLOCATE(map)


  !sum over globalj     
  allocate(anlm_all(0:nshell-1,0:lmax,0:lmax))
  cnt=(nshell)*(lmax+1)**2
  if (ntasks >1 ) then
     if (me==0) print*,'calling MPI_ALLREDUCE'
     CALL MPI_ALLREDUCE(anlm,anlm_all,cnt,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     print*,'No MPI_ALLREDUCE'
     anlm_all=anlm
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  deallocate(anlm)
     
  !save almn
  if (me==0) then
     open(12,file=trim(adjustl(dir))//'almn_from_beta.unf',status='unknown',form='unformatted')
     write(12) anlm_all
     close(12)
  endif


  !=== for now compute alm_r in all processor ===
  do m=0,lmax        
     do l=0,lmax
        fn(0:nshell-1)=anlm_all(0:nshell-1,l,m)
        call complex_fft(fn,backward=.true.) !sum over n
        anlm_all(0:nshell-1,l,m) = fn(0:nshell-1)
     enddo
  enddo

  !save alm_r
  if ( me.eq.0) then 
     open(unit,file=trim(adjustl(dir))//'almr_from_beta.unf' ,status='unknown',form='unformatted')
     write(unit) anlm_all
     close(unit)
  endif

  !=== for now compute sum over lmn in main processor ===

  if ( me.eq.0) then 

     ALLOCATE(map_ball(0:npix_out-1,0:nshell-1))

     do ir=0,nshell-1
        alm(1,0:lmax,0:lmax) =  anlm_all(ir,0:lmax,0:lmax)
        call alm2map(nside_out,lmax,lmax,alm,map_ball(:,ir)) !sum over lm
     enddo

     if (feedback>3) print*,'calling alm2map done!'

     print*,'writing file to :',trim(adjustl(dir))//'map_ball_beta.unf'

     !if (me.eq.0) then
     open(12,file=trim(adjustl(dir))//'map_ball_beta.unf',status='unknown',form='unformatted')
     write(12) npix_out,nshell,4,1
     write(12) map_ball
     close(12)

     deALLOCATE(map_ball)

  endif


  deALLOCATE(fn)
  deALLOCATE(nj_pp)
  deALLOCATE(anlm_all)



  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime


  print*, 'Calling MPI finalize'
  CALL MPI_FINALIZE(ierr)


end program

