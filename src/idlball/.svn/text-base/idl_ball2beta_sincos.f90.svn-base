  
!!!IDL extension to HEALPix [TESTVERSION] 2D needlet: Written by Frode K. Hansen, Nov. 2006
  !! Modified for 3D needlet by Y. Fantaye Sep 2013 
  
  use healpix_types
  use pix_tools
  use extension
  use alm_tools
  use wav_ball_mod

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: nside,npix,lmax,nmax,i,j,l,m,polar,p,wav,large,m1,m2,multialm,lmaxx
  INTEGER(I4B) :: j0,nj,first_j,last_j, v0,nv,nvj,first_v, last_v, in,glnpow,itype
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almn_ball
  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm,anorm
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn,stype
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
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rvec,alpha,beta, gamma,rfili,wfili,w2fili
  integer(i4b) ::ii,jj,oddInd, evenInd,nshell,kk
  real(dp) ::rmin,rmax,dr,rmid,y
  real(dp):: s2n, s2nm1
  CHARACTER(LEN=200) :: fname,pfile

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me

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
  read(12) nv
  read(12) nside
  read(12) polar
  read(12) lmax
  read(12) nmax
  read(12) nrpix
  read(12) rmin,rmax
  read(12) feedback
  close(12)


  nshell=nrpix-1
  if (nv < nshell) nv=nshell

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
     print*, 'itype (1-sinc, 2-sin, 3-cos, 4-sin+cos)= ',4 
     print*, 'nside = ',nside 
     print*, 'polar = ',polar 
     print*, 'lmax = ',lmax
     print*, 'nmax = ',nmax
     print*, 'nshell = ',nshell
     print*, 'rmin, rmax = ',real(rmin),real(rmax)
     print*, 'glnpow,bb = ',glnpow,bb
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
     allocate(rfili(0:nshell),wfili(0:nshell),w2fili(0:nshell))

     !allocate big array 
     allocate(almn_ball(1:1+2*polar,0:lmax,0:lmax,0:nmax))


     print*, 'starting filon_sin integration'

     do itype=2,3

        intflag=itype
        write(stype,*) itype
        if (me==0) print*, 'computing almn_ball for itype=',itype

        almn_ball = 0

        do jj = 0,nmax

           anorm=1d0
           if (jj==0) anorm=2d0

           y = argpi*jj

           call fil_simp_points(rmin,rmax,y,nshell,1,rfili,wfili,itype=itype) !sin or cos
           !alm_TGC = 0
           do ii=0,nshell-1
              almn_ball(1,0:lmax,0:lmax,jj) = almn_ball(1,0:lmax,0:lmax,jj) + wfili(ii)*alm_TGC_r(ii,1,0:lmax,0:lmax)/anorm
           enddo

        enddo

        almn_ball = 2d0*almn_ball  !2d0 normalization 

        !if (me==0) print*, 'r integration: ',rfili

        fname = trim(adjustl(dir))//'almn_from_ball_stype'//trim(adjustl(stype))//'.unf'
        !if (large.eq.0)  fname = trim(adjustl(dir))//'almn.unf' 
        print*, '******* saving almn array to **'
        print*, trim(fname)
        !print*, 'size(almn) = ',size(almn_ball)
        !print*, 'sshape(almn) = ',shape(almn_ball)
        print*, '---------------------------'
        open(unit,file=trim(fname),status='unknown',form='unformatted')
        write(unit) almn_ball
        close(unit)

     enddo !end itype loope

     deallocate(almn_ball)
     deallocate(alm_TGC_r)
  endif

  !!========= almn to beta ====


  if  (me == 0) print*, 'almn_ball write and read from file ..'

  ALLOCATE(nj_pp(0:ntasks-1))
  nj_pp=nj/ntasks
  if(mod(nj,ntasks).ne.0) nj_pp(0:mod(nj,ntasks)-1)=nj_pp(0:mod(nj,ntasks)-1)+1
  first_j=j0
  if (me.ne.0) first_j=j0+sum(nj_pp(0:me-1))
  last_j=first_j+nj_pp(me)-1


  ALLOCATE(beta_ballmap(0:npix-1,1:1+polar*2,0:nv, 0:nj_pp(me)-1))
  ALLOCATE(gln(0:lmax,0:nmax,j0:j0+nj-1))


  if (wav.eq.0) wavtyp='standard'
  if (wav.eq.1) wavtyp='mexican'

  if  (me == 0) print*, 'calling .. calc_f2  ..'
  CALL calc_f2

  if  (me == 0) print*, 'calling .. calc_gln  ..'
  !print*,'alm2beta1:',j0,nj,lmax,bb,wavtyp,p
  CALL calc_gln(j0,nj,lmax,nmax,gln,bb,wavtyp,p)

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



  !allocate big array 
  allocate(almn_ball(1:1+2*polar,0:lmax,0:lmax,0:nmax))
  almn_ball = 0

  do itype=2,3 

     intflag=itype
     write(stype,*) itype
     if (me==0) print*, 'computing beta_ball for itype=',itype     


     !read the corresponding almn
     fname = trim(adjustl(dir))//'almn_from_ball_stype'//trim(adjustl(stype))//'.unf'
     open(unit,file=trim(fname),status='old',form='unformatted')
     read(unit) almn_ball
     close(unit)


     if  (me == 0) print*, 'calling .. alm2map_sincos_ball ..'
     call alm2map_sincos_ball(itype,nside,nmax,lmax,lmax, almn_ball,bb,first_j,beta_ballmap, gln(:,:, first_j:last_j),nj_pp(me),nv)

     WRITE(strn,*) glnpow

     do j=0,nj_pp(me)-1

        WRITE(filenum,*) first_j-j0+j !+1000
        fname=trim(adjustl(dir))//'beta_ball_stype'//trim(adjustl(stype))//'_j'//trim(adjustl(filenum))//'.unf'
        if (glnpow>1) fname=trim(adjustl(dir))//'beta_ball_stype'//trim(adjustl(stype))//'_j'//trim(adjustl(filenum))//'_glnpow'//trim(adjustl(strn))//'.unf'

        !print*, 'writing file:'

        !print*,fname 

        open(12,file=fname,status='unknown',form='unformatted')
        write(12) beta_ballmap(0:npix-1,1,0:nv,j)
        close(12)
     enddo

  enddo !end itype loop


  !print*,'done',me

  CALL MPI_FINALIZE(ierr)

end program

