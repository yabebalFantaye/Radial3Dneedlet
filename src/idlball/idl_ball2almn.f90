program idl_ball2almn

  !yabebal fantaye 12 Jan 2013

  use parsers
  use healpix_types
  use alm_tools
  use pix_tools
  use extension
  use wav_ball_mod
  
  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,i,j,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  REAL(DP), DIMENSION(:), ALLOCATABLE :: beam

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map_TQU

  REAL(DP), DIMENSION(1:2) :: zbounds
  REAL(DP) :: fwhm,sigma,norm
  CHARACTER(LEN=128) :: healpixdir,tempdir,str,finext
  CHARACTER(LEN=200) :: fname,parfile,inputfilename,outputfilename


  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src,sn,s0,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp

  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd,nrsample,kk,nrpix
  real(dp) ::rmin,rmax,dr,rmid,y
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
  call ball2almn(parfile,nside, polar, nnmax, nshell,ifirstshell,ilastshell,iwidth,norm,balltype,inputfile,outputfile,finext,tempdir)

  !check correctness of input parameters
  if (me == 0) then
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nnmax
     print*, 'nshell = ',nshell
     print*, '---------------------------'
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'


     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif

  polar=0  
  npix=nside2npix(nside) 


  !task distribution across processors: each processor recieves nmap_pp(me) shells
  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*(lmax+1)*(lmax+1) !data count per each processor



  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))
  w8ring_TQU=1d0
  call read_ringweights(nside, 1 ,healpixdir, w8ring_TQU )


  zbounds = (/ -1.0_dp, 1.0_dp /)


  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  alm_TGC=0


  if  (me == 0) then
     print*, 'nnmax, lmax, nshell = ',nnmax, lmax, nshell
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif


  ALLOCATE(map_TQU(0:npix-1,0:nmap_pp(me)-1))

  !in each processor loop over radial shells, and do harmonic expansion
  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     write(str,*)i
     open(12,file=trim(adjustl(inputfile))//trim(adjustl(str))//'.'//trim(adjustl(finext)),status='old',form='unformatted')
     read(12) map_TQU
     close(12)

     !print*, 'nshell = , max(map) =  ', i, maxval(map_TQU(:,1,i))
     !if (i == 0) then !the 0th layer is an empty map
     !   alm_TGC=0.
     !else
    if (count( map_TQU(:,i)>1e-13 )>10)  call map2alm_iterative(nside, lmax, lmax, 4, map_TQU(:,ii:ii), alm_TGC, zbounds, w8ring_TQU)
     !endif
     alm_TGC_rpp(ii,:,:) = alm_TGC(1,:,:)

     !print*, 'done with nshell =  from jj:kk= ',i,(nrsample+1)*i,(nrsample+1)*(i+1)-1
  enddo

  deallocate(map_TQU)


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) then
     !allocate big arrays
     allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
     alm_TGC_r = 0
     
     print*, 'MPI send revciev alm_TQU_rpp from work processors'
     print*, ntasks,me
  endif
  

  if (ntasks>1) then

     dest=0
     if (me .ne. 0) then
        tag = me
        cnt=nct_pp(me)
        print*, 'sending proc=',me
        call  MPI_SSEND(alm_TGC_rpp,cnt,MPI_DOUBLE_COMPLEX,dest,tag,MPI_COMM_WORLD,ierr)
     else

        do i=0,ntasks-1

           if (i==0) then
              s0 = 0
              sn = nmap_pp(0)-1
              alm_TGC_r(s0:sn,:,:)=alm_TGC_rpp
           else

              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1
              print*, 'recieving from proc=',i
              call MPI_RECV(alm_TGC_r(s0:sn,:,:),cnt,MPI_DOUBLE_COMPLEX,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
        enddo
        print*, 'MPI send revciev alm_TQU_rpp done!'        
     endif

  else
     alm_TGC_r = alm_TGC_rpp
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !print*,'deallocating alm_TGC_rpp ',me
  deallocate(alm_TGC_rpp)

  if  (me == 0) then

     !damp alm(r) to a file
     
     if (tempdir ~= '') then
        fname = trim(adjustl(tempdir))//'almr.unf' 
        print*, '******* saving alm_r array to **'
        print*, trim(fname)
        print*, '---------------------------'
        open(unit,file=trim(fname),status='unknown',form='unformatted')
        write(unit) alm_TGC_r
        close(unit)
     endif

     print*, 'lmax, nshell = ',lmax,nshell
     !print*, 'size(almn), shape(almn) = ',size(alm_TGC), shape(alm_TGC)

     !allocate big array 
     allocate(almn_TGC(0:nshell-1,0:lmax,0:lmax))
     almn_TGC = 0


     print*, 'starting fftw integration'
     allocate(fr(0:nshell-1))
     do m=0,lmax 
        do l=0,lmax
           fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
           call complex_fft(fr,backward=.false.) 
           almn_TGC(0:nshell-1,l,m) = fr(0:nshell-1)/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
        enddo
     enddo
     deallocate(fr)



     fname = trim(adjustl(outputfile))
     print*, '******* saving almn array to **'
     print*, trim(fname)
     print*, 'size(almn), shape(almn) = ',size(almn_TGC), shape(almn_TGC)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_TGC
     close(unit)

  endif


  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime

  CALL MPI_FINALIZE(ierr)

end program idl_ball2almn



!central r value for the ith radial pixel
!i = ii-1
!rmid = (rvec(i+1) + rvec(i))/2d0
!almn_TGC(jj,:,:,:) =  almn_TGC(jj,:,:,:) + alm_TGC_r(i,:,:,:)*sin(PI*jj*rmid)*rmid*dr
