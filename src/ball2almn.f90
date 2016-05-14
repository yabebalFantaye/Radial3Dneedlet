
program ball2almn

  !yabebal fantaye 12 Jan 2013

  use parsers
  use healpix_types
  use alm_tools
  use pix_tools
  use extension


  use fitstools
  use mapio

  use almnutil
  use fftw_wrapper
  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(LEN=30) :: filestr,str
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext
 
  INTEGER(I4B) :: nside,npix,lmax,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  INTEGER(I4B) :: i, ii, j, jj
  INTEGER(I4B) :: ifstr,ffstr
  integer(i4b) ::nshell,nnmax,oddInd, evenInd,nrsample,kk,nrpix
  !vectors
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam  
  !map
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map

  !alm
  COMPLEX(DPC) :: fr_mean
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr
  COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: anm
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm,alm1
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_rpp,alm_r
  !ps
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cln
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cl,cl_mean
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl_rpp,cl_temp
  !
  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma,norm
  !mpi
  INTEGER(I4B) :: ierr,ntasks,unit,me,MPI_DPC, MPI_p, MPI_r
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp
  !omp
  INTEGER :: nthreads, myid
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  !time
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

  !define new HEALPIX DPC type for MPI communication
!   MPI_p=precision(cmplx( 0, 0, kind = dpc))
!   MPI_r=range(cmplx( 0, 0, kind = dpc))
!   print*, 'MPI_p, MPI_r, MPI_dpc_p',MPI_p, MPI_r,MPI_DOUBLE_COMPLEX

  !call MPI_TYPE_CREATE_F90_COMPLEX (MPI_p ,MPI_r , MPI_DPC, ierr)
  MPI_DPC=MPI_DOUBLE_COMPLEX

  unit=12+me

  call wall_clock_time(time0)
  call cpu_time(ptime0)


  call getEnvironment("HEALPIX",healpixdir)

  !read parameter file
  call almnparser(parfile,nside, lmax, nnmax, nshell,ifirstshell,ilastshell,iwidth,norm,balltype,inputfile,outputfile,finext,tempdir)

  ifstr=6-iwidth+1
  ffstr=6

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
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'nside = ',nside 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, 'ifirst = ',ifirstshell
     print*, 'ilast = ',ilastshell
     print*, 'iwidth = ',iwidth
     print*, 'inputfile=',trim(adjustl(inputfile))
     print*, 'finext=',trim(adjustl(finext))
     print*, 'outputfile=',trim(adjustl(outputfile))
     print*, 'tempdir=',trim(adjustl(tempdir))
     print*, '---------------------------'

  endif

  polar=0  
  npix=nside2npix(nside) 

  !===================================================
  !----- the following is a simple load balancing -----
  !task distribution across processors: each processor recieves nmap_pp(me) shells
  !===================================================
  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks !vector
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*(lmax+1)*(lmax+1) !vector data count per each processor


  !===================================================
  !healpix pixel weights
  !===================================================
  ALLOCATE(w8ring(1:2*nside,1:1))
  w8ring=1d0
  call read_ringweights(nside, 1 ,healpixdir, w8ring )
  zbounds = (/ -1.0_dp, 1.0_dp /)

  !
  sn=nmap_pp(me)  !the number of tasks at processor me
  allocate(alm_rpp(0:sn-1,0:lmax,0:lmax))
  allocate(alm(1,0:lmax,0:lmax))
  alm=0

  if  (me == 0) then
     print*, 'nnmax, lmax, nshell = ',nnmax, lmax, nshell
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif


  ALLOCATE(map(0:npix-1,1))

  !===================================================
  !in each processor loop over radial shells, and do harmonic expansion
  !===================================================
  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     str=''
     !===================================================
     !construct filename 
     !===================================================
     if (ifstr==ffstr) then
        WRITE(str,*) i
        fname=trim(adjustl(inputfile))//trim(adjustl(str))//'.'//trim(adjustl(finext))
     else
        WRITE(str,'(I6)') i + 100000 
        fname=trim(adjustl(inputfile))//trim(adjustl(str(ifstr:ffstr)))//'.'//trim(adjustl(finext))
     endif

     if (i == 0) then !the 0th layer is an empty map
       print*, 'reading file: '//trim(fname)
     endif
     !===================================================
     !read shell map
     !===================================================
     if (trim(finext).eq.'fits') then
        call read_dble_map(fname,map,nside,1,RING) 
     elseif (trim(finext).eq.'unf' .or. trim(finext).eq.'bin0') then
        open(12,file=trim(fname),status='old',form='unformatted',ACCESS="STREAM",IOSTAT=ierr)
        if (ierr==0) then
           read(12) map
           close(12)
        else
           print*, 'open file, status: ',trim(fname),ierr
           stop
        endif
     endif
     
     !===================================================
     !Do map2alm on none zero valued shells - double check if there could be violation to the condition set
     !===================================================
     !if (count( abs(map(:,1))>1e-30 )>1)  
     call map2alm_iterative(nside, lmax, lmax, 3, map, alm, zbounds, w8ring)

     !smooth map by fwhm arcmin
     fwhm=20d0
     call alter_alm(nside,lmax,lmax,fwhm,alm)

     alm_rpp(ii,:,:) = alm(1,:,:)

  enddo

  !free some memory
  deallocate(map,w8ring)
  deallocate(alm)

  call wall_clock_time(time1)
  if (me==0) print*, 't-t0 when map2alm is done: ',time1-time0

  
  !===================================================
  !---------- collect all alm(r) from different processors into one big array -----
  !---------- This is necessary to do FFT
  !===================================================

  dest=0  !the root/master processor where all alm(r) are collected
  if (me==dest) then     
     print*, 'MPI send revciev alm_rpp from work processors'
     print*,'nshell,ntasks,me' ,nshell,ntasks,me
     print*, '----------'

     allocate(alm_r(0:nshell-1, 0:lmax,0:lmax))
     alm_r = 0
  endif

  
  if (ntasks>1) then
     if (me .ne. dest) then
        tag = me
        cnt=nct_pp(me)
        !print*, 'sending: me, dest',me,dest
        if (tag<3) print*, 'sum(alm)= sent to root with tag=',sum(alm_rpp),tag
        call  MPI_SEND(alm_rpp(0:nmap_pp(me)-1,0:lmax,0:lmax),cnt,MPI_DPC,dest,tag,MPI_COMM_WORLD,ierr)
     else           
        do i=0,ntasks-1              
           !print*, 'recv: me,dest,itask',me,dest,i
           if (i==dest) then
              s0 = 0
              sn = nmap_pp(i)-1
              alm_r(s0:sn,0:lmax,0:lmax)=alm_rpp(:,0:lmax,0:lmax)
           else                
              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1
              call MPI_RECV(alm_r(s0:sn,0:lmax,0:lmax),cnt,MPI_DPC,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
           if (tag<3) print*, 'sum(alm)= received in root with tag=',sum(alm_r(s0:sn,0:lmax,0:lmax)),tag
        enddo
     endif
  else
     alm_r(0:nshell-1,0:lmax,0:lmax) = alm_rpp(0:nshell-1,0:lmax,0:lmax)
  endif

  deallocate(alm_rpp)
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (me==0) print*, 'MPI send recv done! me, nshell,dest:',me,nshell,dest

  call wall_clock_time(time1)
  if (me==0) print*, 't-t0 when MPI send recv is done: ',time1-time0

  !===================================================
  !Once all alm(r) are collected, do forward FFT(alm(r))
  !===================================================
  if  (me == 0) then

     !print*,'after mpi_gather, alm_r:',alm_r(:,10,0)
     print*,'---------------------'

     allocate(cl_temp(0:lmax,2))
     fname = trim(adjustl(tempdir))//'meancl.unf'
     call almr2cl(alm_r,cl_temp,fname)
     deallocate(cl_temp)

     !do parallel fft forward
     call fftw_complex_lm(alm_r, .false.) !FFT forward
     !save fft output
     fname = trim(adjustl(outputfile)) 
     call write_almx(alm_r,fname)      

     !save power spectrum
     allocate(cln(0:nshell-1,0:lmax))
     fname = trim(adjustl(tempdir))//'cln.unf'
     call almn2cln(alm_r,cln,fname)

  endif
   
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  if (me==0) write(*,*) "Total Clock and CPU time [s]: ", clock_time, ptime
  
  CALL MPI_FINALIZE(ierr)
  
end program ball2almn



!   !take the sum per node and sum over all nodes of alm
!   alm(1,:,:) = sum(alm_rpp,1)   
!   call MPI_REDUCE(alm, alm1, (lmax+1)*(lmax+1), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD,ierr)  
  
!   !take the sum per node and sum all over node  of cl_rpp
!   !cl_rpp(:,0) = sum(cl_rpp,2)
!   call MPI_REDUCE(cl_rpp, cl_temp, lmax+1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)  

!   !save cl: one from mean alm, and the other from sum(cl(r))
!   if (me==0) then
!      fname = trim(adjustl(tempdir))//'cl_mean_clr.unf'
!      print*, '******* saving sum cl(r) **'
!      print*, 'cl_mean_clr(3:10)',cl_temp(3:10,1)
!      print*, trim(fname)
!      print*, '---------------------------'
!      open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
!      write(unit) cl_temp(:,1)
!      close(unit)

!      print*, 'sum(alm), max(alm_r), min(alm_r), max(alm_i),min(alm_i)',sum(alm),maxval(real(alm)),minval(real(alm)),maxval(aimag(alm)),minval(aimag(alm))
!      call alm2cl(lmax, lmax, alm1,cl_temp) 
!      fname = trim(adjustl(tempdir))//'cl_mean_almr.unf'
!      print*, '******* saving cl from mean alm **'
!      print*, 'cl_mean_alm(3:10)',cl_temp(3:10,1)
!      print*, trim(fname)
!      print*, '---------------------------'
!      open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
!      write(unit) cl_temp(:,1)
!      close(unit)
!   endif
