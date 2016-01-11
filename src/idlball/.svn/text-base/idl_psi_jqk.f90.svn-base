  
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
  INTEGER(I4B) :: j0,nj,first_j,last_j, v0,nv,nvj,first_v, last_v, in,glnpow
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: beta_ballmap
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beam,gl
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gln
  COMPLEX(DPC), DIMENSION(:,:,:,:), ALLOCATABLE :: almn_ball
  REAL(DP) :: fwhm,sigma
  REAL(DP) :: bb,norm
  CHARACTER(LEN=128) :: healpixdir,dir,wavtyp,filenum,strn
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

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me

  if (iargc() > 0) then
     call getarg(1,pfile)
  else
     pfile=trim(adjustl(dir))//'params_psi_jqk.unf'
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


  !nv = ceiling(bb**(j0+nj))


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


  endif



  zbounds = (/ -1.0_dp, 1.0_dp /)

  ALLOCATE(nmap_pp(0:ntasks-1),nct_pp(0:ntasks-1))
  nmap_pp=nrpix/ntasks
  if(mod(nrpix,ntasks).ne.0) nmap_pp(1:mod(nrpix,ntasks)-1)=nmap_pp(1:mod(nrpix,ntasks)-1)+1
  nct_pp = nmap_pp*(1+2*polar)*(lmax+1)*(lmax+1) !data count per each processor

  !!========= almn to beta ====

  !allocate big array 
  allocate(almn_ball(1:1+2*polar,0:lmax,0:lmax,0:nmax))
  almn_ball = 0

  almn_ball=1


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


  if  (me == 0) print*, 'calling .. synmap_multi_ball ..'
     CALL synmap_multi_ball(almn_ball,beta_ballmap,nside,nmax, lmax, bb,first_j,gln(:,:, first_j:last_j),nj_pp(me),nv)
     
   WRITE(strn,*) glnpow

     do j=0,nj_pp(me)-1
        WRITE(filenum,*) first_j-j0+j !+1000

        fname=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'.unf'
        if (glnpow>1) fname=trim(adjustl(dir))//'beta_ball_j'//trim(adjustl(filenum))//'_glnpow'//trim(adjustl(strn))//'.unf'

        !print*, 'writing file:'

        !print*,fname 

        open(12,file=fname,status='unknown',form='unformatted')
        write(12) beta_ballmap(0:npix-1,1,0:nv,j)
        close(12)
     enddo

     !print*,'done',me

     CALL MPI_FINALIZE(ierr)

   end program

