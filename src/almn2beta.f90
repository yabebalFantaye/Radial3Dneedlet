program almn2beta

  !yabebal fantaye 12 Jan 2013
  use parsers
  use healpix_types
  use pix_tools
  use extension
  use alm_tools

  use needlet_kernel
  use almnutil
  use almn_tools
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

  IMPLICIT NONE

  INCLUDE 'mpif.h'


  INTEGER(I4B) :: nside,npix,nside_jmap, npix_jmap 
  INTEGER(I4B) :: lmax,i,j,k,l,m,polar,balltype
  INTEGER(I4B) :: ifirstshell,ilastshell,iwidth
  integer(i4b) ::nshell,ii,nnmax,jj,oddInd, evenInd

  !needlet
  CHARACTER(LEN=30) :: wavtyp
  INTEGER(I4B) :: j0,nj,glnpow,wav, first_j,last_j,first_v, last_v
  REAL(DP) :: fwhm,sigma,norm, bb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln

  !alms, maps, and cls
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fn
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_TGC
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: cl, cln
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: beta_ballmap

  !strings
  CHARACTER(LEN=*),parameter :: code='ALMN2BETA_PROGRAM'
  CHARACTER(LEN=30) :: filestr,str,strn,filenum
  CHARACTER(LEN=FILELEN) :: fname,parfile,inputfile,outputfile,healpixdir,tempdir,finext

  !MPI params
  INTEGER(I4B) :: ierr,ntasks,unit,me
  INTEGER(I4B) :: cnt,dest,tag,src, s0, sn,stat
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nmap_pp,nct_pp,nj_pp

  !timing
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1
  !type(planck_fft2_plan) :: plan_ball2beta
  !call make_fft2_plan (plan_ball2beta,length,fft2_forward)

  if (iargc() > 0) then
     call getarg(1,parfile)
  else
     stop 'usage: '//code//' params_filename.dat'
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
  call betaparser(parfile, bb, j0, nj, glnpow, wav,nside_jmap)
  if (nside_jmap==0) nside_jmap=nside

  if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
     if (me == 0) write ( *, '(a)' ) 'Fatal error! NSHELL must be eve and > 1.'
     if (me == 0) write ( *, '(a,i8)' ) '  NTAB = ',nshell
     stop
  end if

  if ( 2*nnmax>nshell ) then
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
  npix_jmap = nside2npix(nside_jmap)

  allocate(almn_TGC(0:nshell-1,0:lmax,0:lmax))


  fname=trim(adjustl(inputfile))
  call read_almx(almn_TGC, fname)

  call wall_clock_time(time1)
  if (me==0) print*, 't-t0 (sec) when file read is done: ',time1-time0

  if  (me == 0) then
     allocate(cln(0:nshell-1,0:lmax))
     fname = trim(adjustl(tempdir))//'cln.unf'
     call almn2cln(almn_TGC,cln,fname)     
     deallocate(cln)
  endif


  !-------------- now to computing beta_ballmap --------

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1


  
  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  if  (me == 0) print*, 'calling .. calc_f2  ..'
  CALL calc_f2

  if  (me == 0) print*, 'calling .. calc_gln  ..'

  ALLOCATE(gln(0:lmax,0:nshell-1,j0:j0+nj-1))  

  CALL calc_gln(j0,nj,lmax,nnmax,gln,bb,wavtyp,1,gln4fft=1)

  if (me==0) then 
     fname = trim(adjustl(tempdir))//'gln.unf'
     print*, 'saving gln array to: '//trim(fname)
     print*, 'sshape(gln) = ',shape(gln)
     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) gln
     close(unit)
  endif

  if (glnpow>1) then
     if  (me == 0) print*, 'glnpow=2 .. using gln=gln^2 ..'
     gln=gln*gln

     fname = trim(adjustl(tempdir))//'gln2.unf'
     open(unit,file=trim(fname),status='unknown',form='unformatted', ACCESS="STREAM")
     write(unit) gln
     close(unit)
  endif


  ALLOCATE(beta_ballmap(0:npix_jmap-1,0:nshell-1, 0:nj_pp(me)-1))

  !if  (me == 0) 
  print*, 'calling almn_to_beta .. inode=', me
  call almn_to_beta(nside,nshell,lmax,lmax, almn_TGC,beta_ballmap, gln(:,:, first_j:last_j),nj_pp(me))

  print*, 'me, beta_ball.sum()',sum(beta_ballmap)

  if (balltype == slice_balltype) then
     do j=0,nj_pp(me)-1
        do k=0,nshell-1
           WRITE(filenum,*) first_j-j0+j !+1000
           WRITE(strn,*) k
           fname=trim(adjustl(outputfile))//'_j'//trim(adjustl(filenum))//'_r'//trim(adjustl(strn))//'.unf'
           if (glnpow>1) fname=trim(adjustl(outputfile))//'_gln2_j'//trim(adjustl(filenum))//'_r'//trim(adjustl(strn))//'.unf'
           
           open(12,file=fname,status='unknown',form='unformatted', ACCESS="STREAM")
           write(12) beta_ballmap(0:npix-1,k,j)
           close(12)
        enddo
     enddo
  else !full ball map
     do j=0,nj_pp(me)-1
        WRITE(filenum,*) first_j-j0+j !+1000
        fname=trim(adjustl(outputfile))//'_j'//trim(adjustl(filenum))//'.unf'
        if (glnpow>1) fname=trim(adjustl(outputfile))//'_gln2_j'//trim(adjustl(filenum))//'.unf'
        
        open(12,file=fname,status='unknown',form='unformatted', ACCESS="STREAM")
        write(12) beta_ballmap(0:npix-1,0:nshell-1,j)
        close(12)
     enddo
  endif

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  if (me==0) print*, "Total Clock and CPU time [s] since start: ", clock_time, ptime


  CALL MPI_FINALIZE(ierr)

end program almn2beta


