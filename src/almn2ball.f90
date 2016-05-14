program idl_almn2ball

  !yabebal fantaye 12 Jan 2013
  use parsers
  use healpix_types
  use pix_tools
  use extension
  use alm_tools

  use mapio
  use almnutil
  use fftw_wrapper
  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_pp,map
  REAL(DP), DIMENSION(1:2) :: zbounds

  REAL(DP) :: fwhm,sigma,norm
  CHARACTER(LEN=30) :: filestr,str
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext

  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl, cln
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fn
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrpix
  real(dp) ::rmin,rmax,dr,rmid
  real(dp):: s2n, s2nm1

  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)

  if (iargc() > 0) then
     call getarg(1,parfile)
  else
     stop 'usage: box2ball box_filename.dat'
  endif


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me
  call wall_clock_time(time0)
  call cpu_time(ptime0)


  call getEnvironment("HEALPIX",healpixdir)

  !read parameter file
  call almnparser(parfile,nside, lmax, nnmax, nshell,ifirstshell,ilastshell,iwidth,norm,balltype,inputfile,outputfile,finext,tempdir)

  if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
     if (me == 0) write ( *, '(a)' ) 'Fatal error! NSHELL must be eve and > 1.'
     if (me == 0) write ( *, '(a,i8)' ) '  NTAB = ',nshell
     stop
  end if

  if ( nshell <= 2*nnmax ) then
     if (me == 0) print *,  'WARNNING! setting nnmax=NSHELL/2. Given nnmax=',nnmax
     nnmax=nshell/2
  end if



  !check correctness of input parameters
  if (me == 0) then
     print*, '************** almn2ball.f90 *****'
     print*, '---------------------------'
     print*, 'nside = ',nside 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, '---------------------------'

     if ( nshell <= 1 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'Fatal error! NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif
  polar=0  


  npix = nside2npix(nside)

  !------ simple load balancing
  !define per processor task array
  ALLOCATE(nmap_pp(0:ntasks-1), nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*npix
  sn=nmap_pp(me)

  !the big almn array in all processor
  allocate(almn(0:nshell-1,0:lmax,0:lmax))

  !=================================
  ! real almn
  !==================================
  fname=trim(adjustl(inputfile))
  call read_almx(almn, fname)

  call wall_clock_time(time1)
  if (me==0) print*, 't-t0 (sec) when file read is done: ',time1-time0

  !---- 2D power spectrum cln
  if  (me == 0) then
     allocate(cln(0:nshell-1,0:lmax))
     fname = trim(adjustl(tempdir))//'cln.unf'
     call almn2cln(almn,cln,fname)     
     deallocate(cln)
  endif

  !===========================
  !inverse FFT 
  !===========================
  if  (me == 0) print*, 'calling inverse fftw almn in openMP loop.' 
  call fftw_complex_lm(almn, .true.) !now almn is alm(r) in all nodes

  !------- mean of the power spectra on each shell
  if  (me == 0) then     
     allocate(cl(0:lmax, 2))
     fname=trim(adjustl(tempdir))//'meancl.unf'
     call almr2cl(almn, cl,fname)
     deallocate(cl)
  endif
        
  if (me==0)  print*, 'FFTW and meancl done!'

  call wall_clock_time(time1)
  if (me==0) print*, 't-t0 when fftw and almn2cln are done: ',time1-time0

  !============================
  !now alm(r) to map(r)
  !==============================

  !ALLOCATE(map_pp(0:npix-1,0:sn-1)) !if we needed to 
  allocate(map_pp(0:npix-1,1)) !as we are not interested to storeall  maps 
  allocate(alm(1:1+2*polar,0:lmax,0:lmax))
  map_pp=0

  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif


     !if (.not. (sum( real( almn(i,0:lmax,0:lmax) ))==0d0 )) then !set by hand the map at the 0th shell to be 0
     alm(1,0:lmax,0:lmax) = almn(i,0:lmax,0:lmax)
     call alm2map(nside,lmax,lmax,alm,map_pp(:,1)) 
     !endif

     
!      if (tempdir .ne. '') then
!         fname = trim(adjustl(tempdir))//'almr.unf_'//trim(adjustl(str)) 
!         open(unit,file=trim(fname),status='unknown',form='unformatted')
!         write(unit) alm
!         close(unit)
!      endif
     
     write(str,*)i
     fname = trim(adjustl(outputfile))//'_'//trim(adjustl(str)) 

     if (me==0 .and. i==0) then 
        print*, '******* almn2ball: saving mapr and alm_r array to **'
        print*, trim(fname)
        print*, '---------------------------'
     endif

     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) map_pp(:,1)
     close(unit)
     
  enddo

  !free alm memory
  deallocate(alm,almn)
  deallocate(map_pp)

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  if (me==0) print*, "Total Clock and CPU time [s] since start: ", clock_time, ptime


  CALL MPI_FINALIZE(ierr)

end program idl_almn2ball




!   !allocate big array to collect all maps at the master node 
!   if (me==0) ALLOCATE(map(0:npix-1,0:nshell-1))

!   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!   if (me==0) print*, 'MPI send revcieve maps from work processors'

!   if (ntasks>1) then

!      dest=0
!      if (me .ne. 0) then
!         tag = me
!         cnt=nct_pp(me)
!         call  MPI_SSEND(map_pp,cnt,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
!      else
!         do i=0,ntasks-1

!            if (i==0) then
!               s0 = 0
!               sn = nmap_pp(0)-1
!               map(:,s0:sn)=map_pp
!            else

!               src = i
!               tag = i
!               cnt=nct_pp(i)
!               s0 = sum(nmap_pp(0:i-1))
!               sn = sum(nmap_pp(0:i))-1

!               call MPI_RECV(map(:,s0:sn),cnt,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,stat,ierr)
!            endif
!         enddo
!      endif

!   else
!      map = map_pp
!   endif

!   if (me==0) then

!      fname=trim(adjustl(outputfile))

!      print*, '************* saving ballmap array to **'
!      print*, trim(fname)
!      print*, '---------------------------'

!      open(unit,file=trim(fname),status='unknown',form='unformatted')
!      write(12) npix,nshell,4,1
!      write(12) map
!      close(12)
!   endif
