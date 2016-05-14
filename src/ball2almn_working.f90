
program idl_ball2almn

  !yabebal fantaye 12 Jan 2013

  use parsers
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod

  use fitstools
  use head_fits
  use mapio

  use fftw_wrapper
  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  INTEGER(I4B) :: ifstr,ffstr
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC,alm
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_TQU

  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma,norm
  CHARACTER(LEN=30) :: filestr,str
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC) :: fr_mean
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_TGC, alm_TGC_rpp,alm_TGC_r
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr
  COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: anm
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cln
  REAL(DP), DIMENSION(:), ALLOCATABLE :: cl,cl_mean
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl_rpp,cl_temp
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1

  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)

  INTEGER :: nthreads, myid
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS




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


  !task distribution across processors: each processor recieves nmap_pp(me) shells
  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks !vector
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*(lmax+1) !*(lmax+1) !vector data count per each processor



  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))
  w8ring_TQU=1d0
  call read_ringweights(nside, 1 ,healpixdir, w8ring_TQU )


  zbounds = (/ -1.0_dp, 1.0_dp /)


  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,0:lmax,0:lmax))
  allocate(cl_rpp(0:lmax,1))
  allocate(cl_temp(0:lmax,1))
  allocate(cl(0:lmax))
  allocate(alm_TGC(1:1,0:lmax,0:lmax))
  allocate(alm(1,0:lmax,0:lmax))

  cl_rpp=0d0
  cl_temp=0d0
  cl=0d0
  alm_TGC=0
  alm=0d0


  if  (me == 0) then
     print*, 'nnmax, lmax, nshell = ',nnmax, lmax, nshell
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif


  ALLOCATE(map(0:npix-1,1))
  !ALLOCATE(map_TQU(0:npix-1,0:nmap_pp(me)-1))


  !in each processor loop over radial shells, and do harmonic expansion
  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     str=''

     !print*, 'proc: ',trim(adjustl(str(ifstr:ffstr)))
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

     if (trim(finext).eq.'fits') then
        call read_dble_map(fname,map,nside,1,RING) 
     elseif (trim(finext).eq.'unf' .or. trim(finext).eq.'bin0') then
        open(12,file=trim(fname),status='old',form='unformatted',IOSTAT=ierr)
        if (ierr==0) then
           read(12) map
           close(12)
        else
           print*, 'open file, status: ',trim(fname),ierr
           stop
        endif
     endif
     
     !map_TQU(:,ii)=dble(map(:,1))
     !print*, 'nshell = , max(map) =  ', i , maxval(map(:,1))

    if (count( map(:,1)>1e-13 )>10)  call map2alm_iterative(nside, lmax, lmax, 0, map, alm_TGC, zbounds, w8ring_TQU)

    alm_TGC_rpp(ii,:,:) = alm_TGC(1,:,:)

    !compute power spectra on each shell
    call alm2cl(lmax, lmax, alm_TGC,cl_temp)
    cl_rpp(:,1) = cl_rpp(:,1)+cl_temp(:,1) !take the sum on each node

  enddo

  deallocate(map,w8ring_TQU)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) print*, 'done with map2alm on all shells'

  !take the sum per node and sum over all nodes of alm
  alm_TGC(1,:,:) = sum(alm_TGC_rpp,1)   
  call MPI_REDUCE(alm_TGC, alm, (lmax+1)*(lmax+1), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD,ierr)  
  
  !take the sum per node and sum all over node  of cl_rpp
  !cl_rpp(:,0) = sum(cl_rpp,2)
  call MPI_REDUCE(cl_rpp, cl_temp, lmax+1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)  

  !save cl: one from mean alm, and the other from sum(cl(r))
  if (me==0) then
     fname = trim(adjustl(tempdir))//'cl_mean_clr.unf'
     print*, '******* saving sum cl(r) **'
     print*, 'cl_mean_clr(3:10)',cl_temp(3:10,1)
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) cl_temp(:,1)
     close(unit)

     call alm2cl(lmax, lmax, alm,cl_temp) 
     fname = trim(adjustl(tempdir))//'cl_mean_almr.unf'
     print*, '******* saving cl from mean alm **'
     print*, 'cl_mean_alm(3:10)',cl_temp(3:10,1)
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) cl_temp(:,1)
     close(unit)
  endif

  !free some memory
  deallocate(cl, cl_temp, cl_rpp, alm, alm_TGC)

  dest=0
  if (me==dest) then     
     print*, 'MPI send revciev alm_TQU_rpp from work processors'
     print*,'nshell,ntasks,me' ,nshell,ntasks,me
     print*, '----------'
     allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
     alm_TGC_r = 0
  endif

  
  if (ntasks>1) then
     if (me .ne. dest) then
        tag = me
        cnt=nct_pp(me)
        print*, 'sending: me, dest',me,dest
        call  MPI_SSEND(alm_TGC_rpp(0:nmap_pp(me)-1,0:lmax,0:lmax),cnt,MPI_DOUBLE_COMPLEX,dest,tag,MPI_COMM_WORLD,ierr)
     else           
        do i=0,ntasks-1              
           print*, 'recv: me,dest,itask',me,dest,i
           if (i==dest) then
              s0 = dest
              sn = nmap_pp(dest)-1
              alm_TGC_r(s0:sn,0:lmax,0:lmax)=alm_TGC_rpp(:,0:lmax,0:lmax)
           else                
              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1
              call MPI_RECV(alm_TGC_r(s0:sn,0:lmax,0:lmax),cnt,MPI_DOUBLE_COMPLEX,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
        enddo
     endif
  else
     alm_TGC_r(0:nshell-1,0:lmax,0:lmax) = alm_TGC_rpp(0:nshell-1,0:lmax,0:lmax)
  endif
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  print*, 'MPI send recv done! me, nshell,dest:',me,nshell,dest


  deallocate(alm_TGC_rpp)
  
  if  (me == 0) then

     allocate(cln(0:nshell-1,0:lmax))
     allocate(fr(-3:nshell-1+3)) !pad -3:0, nshell:nshell+2
     !allocate(anm(0:nshell-1,0:lmax))
     !anm=0

     cln=0d0

     !--------- Start OpenMp -----------
     !$OMP PARALLEL DEFAULT(PRIVATE), &     
     !$OMP& SHARED(alm_TGC_r,nshell,lmax,unit,outputfile,nthreads,cln)
     myid = OMP_GET_THREAD_NUM()                        
     if (myid==0) then
        nthreads = OMP_GET_NUM_THREADS()
        print*, 'starting openmp fft. lmax, nshell,nthreads',lmax, nshell,nthreads
     endif

     !$OMP DO PRIVATE(myid,l, m,fr)
     do l=0,lmax        
        do m=0,l
           fr=0
           fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
           
           if (m==0 .and. mod(l,50)==0) print*, 'myid, nthread, l, nshell, fr before fft ',myid, nthreads, l, nshell,fr(0:5)
           !Healpix forward fft. First half positive wavenumber, second half negative wavenumber 
           call complex_fft(fr,backward=.false.) 
           if (m==0  .and. mod(l,50)==0) print*, 'myid, l, fr after fft ',myid,l, fr(0:5)           
           
           alm_TGC_r(0:nshell-1,l,m) = fr(0:nshell-1) !/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
           !DC mode n=0
           !Nyquist mode n=n_nyquist=nshell/2
           !positive frequency 0:nshell/2 -  n=0,1,2..nshell/2
           !negative frequency nshell/2+1:nshell-1 - n=[1,2,3...nshell/2-1]-nshell/2,
           cln(0:nshell-1,l) = cln(0:nshell-1,l) + abs(fr(0:nshell-1)*conjg(fr(0:nshell-1)))/(l+1d0)           
        enddo !m loop        
     enddo  !l loop
     !$OMP END DO
     !$OMP END PARALLEL     
     !-------------------------------------

     print*, 'openmp done! Saving cln.'


     fname = trim(adjustl(outputfile))  !//'_'//trim(adjustl(str))
     print*, 'save almn to file: '//trim(fname)
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) alm_TGC_r
     close(unit)
     
  
     fname = trim(adjustl(tempdir))//'cln.unf'
     print*, '******* saving 2d power spectrum cln **'
     print*, 'cln(3:10,3:10)',cln(3:10,3:10)
     print*, trim(fname)
     print*, 'size(almn), shape(almn) = ',size(cln), shape(cln)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) cln
     close(unit)
          
  endif
  
  
  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  if (me<=10) write(*,*) " Clock and CPU time [s]: ", clock_time, ptime
  
  CALL MPI_FINALIZE(ierr)
  
end program idl_ball2almn



     
!      do l = 0,lmax
!         do m=0,l
!            fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
              
!            if (m==0) print*, 'l, nshell, fr before fft ',l, nshell,fr(0:5)
!            !Healpix forward fft. First half positive wavenumber, second half negative wavenumber 
!            call complex_fft(fr,backward=.false.) 
!            if (m==0) print*, 'fr after fft ',fr(0:5)           
           
!            anm(0:nshell-1,m) = fr(0:nshell-1) !/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
!            !DC mode n=0
!            !Nyquist mode n=n_nyquist=nshell/2
!            !positive frequency 0:nshell/2 -  n=0,1,2..nshell/2
!            !negative frequency nshell/2+1:nshell-1 - n=[1,2,3...nshell/2-1]-nshell/2,
!            cln(0:nshell-1,l) = cln(0:nshell-1,l) + abs(fr(0:nshell-1)*conjg(fr(0:nshell-1)))/(l+1d0)           
!         enddo !m loop
        
!         !save anm
!         WRITE(str,*) l
!         fname = trim(adjustl(outputfile))//'_'//trim(adjustl(str))
!         open(unit,file=trim(fname),status='unknown',form='unformatted')
!         write(unit) anm
!         close(unit)
        
        
!         !info
!         if (me==dest .and. mod(l,10)==0) then
!            print *, 'mpi send rcv done. l=',l
!            print*, 'fr',fr(0:5)
!         endif
        
!      enddo !end l loop



!      !damp alm(r) to a file     
!      if (tempdir .ne. '') then
!         fname = trim(adjustl(tempdir))//'almr.unf' 
!         print*, '******* saving alm_r array to **'
!         print*, trim(fname)
!         print*, '---------------------------'
!         open(unit,file=trim(fname),status='unknown',form='unformatted')
!         write(unit) alm_TGC_r
!         close(unit)
!      endif


!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
