  
!!!IDL extension to HEALPix [TESTVERSION] 2D needlet: Written by Frode K. Hansen, Nov. 2006
  !! Modified for 3D needlet by Y. Fantaye Sep 2013 
  
  use healpix_types
  use pix_tools
  use extension
  use alm_tools

  use mapio
  use healpix_fft
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !healpix
  INTEGER(I4B) :: nside,npix,lmax,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx,glnpow
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds

  !ball
  INTEGER(I4B) :: nmax,nvmax,nshell,in
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
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: almn_ball
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: fr


  CHARACTER(LEN=200) :: fname,pfile
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn




  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me

  call wall_clock_time(time0)
  call cpu_time(ptime0)

  if (iargc() > 0) then
     call getarg(1,pfile)
  else
     pfile=trim(adjustl(dir))//'params_ball2beta.unf'
  endif


  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

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


  !if (nmax>nshell/2) nmax=nshell/2

  !!necessary to be defined params
  polar=0
  wav=0
  p=0
  m1=1
  m2=0

  npix=12*nside**2


  if (me == 0) then
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
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'



     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) ' NSHELL must be even and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  !!======= ball to almn =====




  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))
  w8ring_TQU=1d0
  call read_ringweights(nside, 1 ,healpixdir, w8ring )

  zbounds = (/ -1.0_dp, 1.0_dp /)

  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nshell/ntasks
  if(mod(nshell,ntasks).ne.0) nmap_pp(1:mod(nshell,ntasks)-1)=nmap_pp(1:mod(nshell,ntasks)-1)+1
  nct_pp = nmap_pp*(lmax+1)*(lmax+1) !*(1+2*polar) !data count per each processor



  if  (me == 0) then
     print*, 'nmax, lmax, nshell = ',nmax, lmax, nshell
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif


  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))

  alm_TGC=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
  alm_TGC_rpp=CMPLX(0.0_dp, 0.0_dp, KIND=DP)

  ALLOCATE(ballmap(0:npix-1,0:nmap_pp(me)-1))

  !in each processor loop over radial shells, and do harmonic expansion
  do ii=0,nmap_pp(me)-1

     open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
     read(12) ballmap
     close(12)


     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     !print*, 'nshell = , max(map) =  ', i, maxval(map_TQU(:,1,i))
     !if (.not. (sum(ballmap(:,i))==0d0)) then !the 0th layer is an empty map
     if (count( ballmap(:,1)>1e-13 )>10) call map2alm_iterative(nside, lmax, lmax,3, ballmap(:,i:i), alm_TGC, zbounds, w8ring_TQU)
     alm_TGC_rpp(ii,0:lmax,0:lmax) = alm_TGC(1,0:lmax,0:lmax)
     !endif

  enddo

  deallocate(ballmap)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) then
     !allocate big arrays
     allocate(alm_TGC_r(0:nshell-1,0:lmax,0:lmax))
     alm_TGC_r = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

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

     fname = trim(adjustl(dir))//'almr_from_ball.unf' 
     print*, '******* saving alm_r array to **'
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) alm_TGC_r
     close(unit)

     print*, 'lmax, nshell, nshell = ',lmax,nshell,nshell
     !print*, 'size(almn), shape(almn) = ',size(alm_TGC), shape(alm_TGC)
     print*, 'sum(almn) = ',sum(alm_TGC_r)


     !allocate big array 
     allocate(almn_ball(0:nshell-1,0:lmax,0:lmax))
     almn_ball = CMPLX(0.0_dp, 0.0_dp, KIND=DP)


     print*, 'starting fftw integration'
     allocate(fr(0:nshell-1))
     fr=CMPLX(0.0_dp, 0.0_dp, KIND=DP)
     do m=0,lmax 
        do l=0,lmax
           fr(0:nshell-1)=alm_TGC_r(0:nshell-1,l,m)
           call complex_fft(fr,backward=.false.) 
           almn_ball(0:nshell-1,l,m) = fr(0:nshell-1)/dble(nshell) !normalized FFTW -> FFTW(IFFTW(f))=f
        enddo
     enddo
     deallocate(fr)

    ! do iv=0,nmax
    !    print*, 'main: iv, sum(almn)',sum(almn(iv,:,:))
    ! enddo


     fname = trim(adjustl(dir))//'almn_from_ball.unf'
     !if (large.eq.0)  fname = trim(adjustl(dir))//'almn.unf' 
     print*, '******* saving almn array to **'
     print*, trim(fname)
     print*, 'size(almn) = ',size(almn_ball)
     print*, 'sshape(almn) = ',shape(almn_ball)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) almn_ball
     close(unit)

     deallocate(almn_ball)
     deallocate(alm_TGC_r)
  endif

  !!========= almn to beta ====

  !allocate big array 
  allocate(almn_ball(0:nshell-1,0:lmax,0:lmax))
  almn_ball = 0

  fname = trim(adjustl(dir))//'almn_from_ball.unf'
  open(unit,file=trim(fname),status='old',form='unformatted')
  read(unit) almn_ball
  close(unit)



  if  (me == 0) print*, 'almn_ball write and read from file ..'

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1



  !print*,'ntasks:',ntasks,nj_pp(me),first_j,last_j

  ALLOCATE(beta_ballmap(0:npix-1,0:nshell-1, 0:nj_pp(me)-1))
  ALLOCATE(gln(0:lmax,0:nshell-1,j0:j0+nj-1))



  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  if  (me == 0) print*, 'calling .. calc_f2  ..'
  CALL calc_f2

  if  (me == 0) print*, 'calling .. calc_gln  ..'
  !print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p
  CALL calc_gln(j0,nj,lmax,nshell-1,gln,bb,wavtyp,p,gln4fft=1)

  if (glnpow>1) then
     if  (me == 0) print*, 'glnpow=2 .. using gln=gln^2 ..'
     gln=gln*gln
  endif

  if (me==0) then 
     fname = trim(adjustl(dir))//'gln.unf'
     if (glnpow>1) fname = trim(adjustl(dir))//'gln2.unf'
     !if (large.eq.0)  fname = trim(adjustl(dir))//'almn.unf' 
     print*, '******* saving gln array to **'
     print*, trim(fname)
     print*, 'size(gln) = ',size(gln)
     print*, 'sshape(gln) = ',shape(gln)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) gln
     close(unit)
  endif


  if  (me == 0) print*, 'calling .. synmap_multi_ball ..'
  call almn2beta_dp(nside,nshell,lmax,lmax, almn_ball,beta_ballmap, gln(:,:, first_j:last_j),nj_pp(me))


  WRITE(strn,*) glnpow

  do j=0,nj_pp(me)-1
     WRITE(filenum,*) first_j-j0+j !+1000

     fname=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'.unf'
     if (glnpow>1) fname=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'_glnpow'//trim(adjustl(strn))//'.unf'
     
     open(12,file=fname,status='unknown',form='unformatted')
     write(12) beta_ballmap(0:npix-1,0:nshell-1,j)
     close(12)
     !enddo
  enddo

  ! do j=0,nshell
  !    print*, 'j, sum(beta_ballmap)',sum(beta_ballmap(:,j,:))
  ! enddo

  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime

  !print*,'done',me

  CALL MPI_FINALIZE(ierr)

end program

