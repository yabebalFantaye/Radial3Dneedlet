  
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

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,nmax,i,j,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx
  INTEGER(I4B) :: j0,nj,first_j,last_j, v0,nv,nvj,first_v, last_v, in
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almn_ball
  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum
  INTEGER(I4B) :: ierr,ntasks
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nj_pp
  INTEGER(I4B) :: cnt,dest,tag,src,me,unit,stat

  integer(i4b) :: nrpix,sn,s0
  REAL(Dp), dimension(:,:,:), allocatable::ballmap
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU
  REAL(DP), DIMENSION(1:2) :: zbounds
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nct_pp,nmap_pp
  COMPLEX(DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm_TGC
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: alm_TGC_r,almn_TGC, alm_TGC_rpp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili
  integer(i4b) ::ii,jj,oddInd, evenInd,nshell,kk
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1
  CHARACTER(LEN=200) :: fname,pfile

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

  if (iargc() > 0) then
     call getarg(1,pfile)
  else
     pfile=trim(adjustl(dir))//'params_ball2beta.unf'
  endif


  call getEnvironment("HEALPIX",healpixdir)
  dir='temp/'!trim(adjustl(healpixdir))//'/src/idl/almtools/'

  open(12,file=pfile,status='old',form='unformatted')
  read(12) bb
  read(12) j0
  read(12) nj
  read(12) nv
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nmax
  read(12) nrpix
  read(12) rmin,rmax
  read(12) feedback
  close(12)


  !nv = ceiling(bb**(j0+nj))


  nshell=nrpix-1
  if (nv < nshell) nv=nshell

  nmax=nshell


  !!necessary to be defined params
  polar=0
  wav=0
  p=0
  m1=1
  m2=0

  npix=12*nside**2
  !nv=int(bb**(j0+nj-1))

  if (me == 0) then
     print*, '**************** ball2almn.f90 *****'
     print*, '---------------------------'
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, 'bb = ',bb
     print*, 'npix, nv, nj, j0 = ',npix, nv,nj,j0
     print*, '---------------------------'
     if (polar>0) print*, 'WARNING: using polar=0, idl_ball2almn can not handle polarization yet.'

     if ( rmin == rmax ) then
        stop 'idl_ball2beta : rmin=rmax'
     end if

     if ( nshell <= 2 .or. mod ( nshell, 2 ) /= 0 ) then
        write ( *, '(a)' ) 'FILON_SIN - Fatal error! NSHELL must be eve and > 1.'
        write ( *, '(a,i8)' ) '  NTAB = ',nshell
        stop
     end if
  endif


  !!======= ball to almn =====




  ALLOCATE(w8ring_TQU(1:2*nside,1:1+2*polar))

  w8ring_TQU=1d0
  open(12,file=trim(adjustl(dir))//'w8.unf',status='old',form='unformatted')
  read(12) w8ring_TQU
  close(12)
  zbounds = (/ -1.0_dp, 1.0_dp /)

  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nrpix/ntasks
  if(mod(nrpix,ntasks).ne.0) nmap_pp(1:mod(nrpix,ntasks)-1)=nmap_pp(1:mod(nrpix,ntasks)-1)+1
  nct_pp = nmap_pp*(1+2*polar)*(lmax+1)*(lmax+1) !data count per each processor

  ALLOCATE(ballmap(0:npix-1,1:1+2*polar,0:nshell))

  open(12,file=trim(adjustl(dir))//'ballmap.unf',status='old',form='unformatted')
  read(12) ballmap
  close(12)

  sn=nmap_pp(me)
  allocate(alm_TGC_rpp(0:sn-1,1:1+2*polar,0:lmax,0:lmax))
  allocate(alm_TGC(1:1+2*polar,0:lmax,0:lmax))
  alm_TGC=0

  if  (me == 0) then
     print*, 'nmax, lmax, nrpix = ',nmax, lmax, nrpix
     print*, 'ntasks, nshell, nmap_pp = ',ntasks, nshell, nmap_pp
     print*, 'calling map2alm in MPI distributed loop'
  endif

  !in each processor loop over radial shells, and do harmonic expansion
  do ii=0,nmap_pp(me)-1

     if (me.ne.0) then
        i=ii+sum(nmap_pp(0:me-1))
     else
        i=ii
     endif

     !print*, 'nshell = , max(map) =  ', i, maxval(map_TQU(:,1,i))
     if (i == 0) then !the 0th layer is an empty map
        alm_TGC=0.
     else
        if (polar.eq.0) call map2alm(nside, lmax, lmax, ballmap(:,1,i), alm_TGC, zbounds, w8ring_TQU)
        if (polar.eq.1) call map2alm(nside, lmax, lmax, ballmap(:,:,i), alm_TGC, zbounds, w8ring_TQU)
     endif
     alm_TGC_rpp(ii,:,:,:) = alm_TGC(:,:,:)

     !print*, 'done with nshell =  from jj:kk= ',i,(nrsample+1)*i,(nrsample+1)*(i+1)-1
  enddo

  deallocate(ballmap)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (me==0) then
     !allocate big arrays
     allocate(alm_TGC_r(0:nshell,1:1+2*polar,0:lmax,0:lmax))
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
              alm_TGC_r(s0:sn,:,:,:)=alm_TGC_rpp
           else

              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nmap_pp(0:i-1))
              sn = sum(nmap_pp(0:i))-1
              print*, 'recieving from proc=',i
              call MPI_RECV(alm_TGC_r(s0:sn,:,:,:),cnt,MPI_DOUBLE_COMPLEX,src,tag,MPI_COMM_WORLD,stat,ierr)
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

     ! nmap_pp=nrpix/ntasks
     ! if(mod(nrpix,ntasks).ne.0) nmap_pp(1:mod(nrpix,ntasks)-1)=nmap_pp(1:mod(nrpix,ntasks)-1)+1
     ! nct_pp = nmap_pp*(1+2*polar)*(lmax+1)*(lmax+1) !data count per each processor
     
     ! !in each processor loop over single, (m,l)  arrays, and do 1d FFTW
     ! do ii=0,nmap_pp(me)-1
     
     !damp alm(r) to a file

     fname = trim(adjustl(dir))//'almr_from_ball.unf' 
     print*, '******* saving alm_r array to **'
     print*, trim(fname)
     print*, '---------------------------'
     open(unit,file=trim(fname),status='unknown',form='unformatted')
     write(unit) alm_TGC_r
     close(unit)

     print*, 'lmax, nrpix, nshell = ',lmax,nrpix,nshell
     !print*, 'size(almn), shape(almn) = ',size(alm_TGC), shape(alm_TGC)


     !allocate Filon-simpson quadrature nodes and weights
     allocate(rfili(0:nshell),wfili(0:nshell))

     !allocate big array 
     allocate(almn_ball(1:1+2*polar,0:lmax,0:lmax,0:nshell/2)) !fftw will assume nmax=nshell/2+1 --> 0:nshell/2
     almn_ball = 0

     print*, 'starting fftw integration'

     do ii=0,lmax 
        do jj=0,lmax
           wfili=alm_TGC_r(0:nshell,1,ii,jj)
           call complex_fft(wfili,backward=.false.) 
           almn_ball(1,ii,jj,0:nshell) = wfili(0:nshell)
        enddo
     enddo

     !almn_ball = 2d0*almn_ball  !2d0 normalization 

     if (me==0) print*, 'r integration: ',rfili

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
  allocate(almn_ball(1:1+2*polar,0:lmax,0:lmax,0:nmax))
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

     ALLOCATE(beta_ballmap(0:npix-1,1:1+polar*2,0:nv, 0:nj_pp(me)-1))
     ALLOCATE(gln(0:lmax,0:nmax,j0:j0+nj-1))


     !sigma=fwhm/60./180.*pi/sqrt(8.*log(2.))
     !allocate(beam(0:lmax,1:1+polar*2))
     ! open(12,file=trim(adjustl(dir))//'beam_beta.unf',status='old',form='unformatted')
     ! read(12) beam
     ! close(12)


     ! do l=0,lmax
     !    alm_TGC(1,l,0:l,:)=alm_TGC(1,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,1)
     !    if (polar.eq.1) then
     !       alm_TGC(2,l,0:l,:)=alm_TGC(2,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,2)
     !       alm_TGC(3,l,0:l,:)=alm_TGC(3,l,0:l,:)*exp(-sigma**2*l*(l+1)/2d0)*beam(l,3)
     !    endif
     ! enddo

     if (wav.eq.0) wavtyp='standard'
     if (wav.eq.1) wavtyp='mexican'

  if  (me == 0) print*, 'calling .. calc_f2  ..'
     CALL calc_f2

  if  (me == 0) print*, 'calling .. calc_gln  ..'
     !print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p
     CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p)


     if (me==0) then 
        fname = trim(adjustl(dir))//'gln.unf'
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

! do ii=0,lmax
! do kk=0,nmax
! do jj=j0,j0+nv-1
!    print*,'l, n, j, glnj = ',ii,kk,jj,gln(ii,kk,jj)
! enddo
! enddo
! enddo


  if  (me == 0) print*, 'calling .. synmap_multi_ball ..'
     CALL synmap_multi_ball(almn_ball,beta_ballmap,nside,nmax, lmax, bb,first_j,gln(:,:, first_j:last_j),nj_pp(me),nv)

     do j=0,nj_pp(me)-1
        WRITE(filenum,*) first_j-j0+j !+1000

        print*, 'writing file:'
        print*, trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'.unf'

        open(12,file=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'.unf',status='unknown',form='unformatted')
        write(12) beta_ballmap(0:npix-1,1,0:nv,j)
        close(12)
     enddo


  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " Clock and CPU time [s]: ", clock_time, ptime

     !print*,'done',me

     CALL MPI_FINALIZE(ierr)

   end program

