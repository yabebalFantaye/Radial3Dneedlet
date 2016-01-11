program idl_almn2ball

  !yabebal fantaye 12 Jan 2013

  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod
  use sincm

  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,large
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_TQU_pp,map_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma
  CHARACTER(LEN=128) :: healpixdir,dir,str
  CHARACTER(LEN=200) :: fname

  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_TGC
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fn
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrpix
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1

  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me
  call wall_clock_time(time0)
  call cpu_time(ptime0)


  call getEnvironment("HEALPIX",healpixdir)
  dir='large/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  !read parameters
  open(12,file=trim(adjustl(dir))//'params_almn2ball.unf',status='old',form='unformatted')
  read(12) large
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nnmax
  read(12) nshell
  close(12)

  nnmax=nshell


  !check correctness of input parameters
  if (me == 0) then
     print*, '************** almn2ball.f90 *****'
     print*, '---------------------------'
     print*, 'large = ', large 
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, '---------------------------'

     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'


     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'Fatal error! NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif
  polar=0  


  npix = nside2npix(nside)

  !define per processor task array
  ALLOCATE(nmap_pp(0:ntasks-1), nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*npix
  sn=nmap_pp(me)

  ALLOCATE(map_TQU_pp(0:npix-1,0:sn-1))
  allocate(almn_TGC(0:nshell-1,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))


  !allocate(alm_TGC_r(0:nshell-1,1:1+2*polar,0:lmax,0:lmax))

  !load almn array
  open(unit,file=trim(adjustl(dir))//'almn.unf',status='old',form='unformatted')
  read(unit) almn_TGC
  close(unit)


  if  (me == 0) print*, 'calling alm2map in MPI distributed loop.'


  if (me==0) print*, 'starting inverse fftw integration'
  allocate(fn(0:nshell-1))
  do m=0,lmax 
     do l=0,lmax
        fn(0:nshell-1)=almn_TGC(0:nshell-1,l,m)
        call complex_fft(fn,backward=.true.) 
        almn_TGC(0:nshell-1,l,m) = fn(0:nshell-1)
     enddo
  enddo
  deallocate(fn)
  if (me==0)  print*, 'FFTW done!'


  !now alm(r) to map(r)

  map_TQU_pp=0
  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif


     !if (.not. (sum( real( almn_TGC(i,0:lmax,0:lmax) ))==0d0 )) then !set by hand the map at the 0th shell to be 0
     alm_TGC(1,0:lmax,0:lmax) = almn_TGC(i,0:lmax,0:lmax)
     call alm2map(nside,lmax,lmax,alm_TGC,map_TQU_pp(:,ii)) 
     !endif

     write(str,*)i
     fname = trim(adjustl(dir))//'almr.unf_'//trim(adjustl(str)) 
     if (me==0 .and. i==0) then 
        print*, '******* almn2ball: saving alm_r array to **'
        print*, trim(fname)
        print*, '---------------------------'
     endif

     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) alm_TGC
     close(unit)

  enddo

  !free alm memory
  deallocate(alm_TGC,almn_TGC)

  !allocate big array to collect all maps at the master node 
  if (me==0) ALLOCATE(map_TQU(0:npix-1,0:nshell-1))

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) print*, 'MPI send revcieve maps from work processors'

  if (ntasks>1) then

     dest=0
     if (me .ne. 0) then
        tag = me
        cnt=nct_pp(me)
        call  MPI_SSEND(map_TQU_pp,cnt,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
     else
        do i=0,ntasks-1

           if (i==0) then
              s0 = 0
              sn = nmap_pp(0)-1
              map_TQU(:,s0:sn)=map_TQU_pp
           else

              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1

              call MPI_RECV(map_TQU(:,s0:sn),cnt,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
        enddo
     endif

  else
     map_TQU = map_TQU_pp
  endif

  if (me==0) then


     !if (large.eq.1) 
     fname=trim(adjustl(dir))//'ballmap.unf'
     !if (large.eq.0) fname=trim(adjustl(dir))//'ballmap.unf'

     print*, '************* saving ballmap array to **'
     print*, trim(fname)
     print*, '---------------------------'

     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(12) npix,nshell,4,1
     write(12) map_TQU
     close(12)
  endif


  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime


  CALL MPI_FINALIZE(ierr)

end program idl_almn2ball

